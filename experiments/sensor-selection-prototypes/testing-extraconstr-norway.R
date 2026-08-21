#pak::pak("trafikkdata/xdtkit")
#install.packages("jsonlite")
#install.packages("sf")
# install.packages(
#   "INLA",
#   repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"),
#   dep = TRUE)
# inla.binary.install(os = c("Ubuntu-22.04"))
library(INLA)
library(xdtkit)
library(jsonlite)
source("R/scripts_inla_sensor_placement.R")

year <- 2024

norway_directed_traffic_links <- jsonlite::read_json(
  path = "data/directed-traffic-links-2024.json",
  simplifyVector = TRUE,
  simplifyDataFrame = TRUE
)
norway_nodes_raw <- sf::st_read("data/traffic-nodes-2024.geojson")


preprocessed_traffic_links_norway <- preprocess_traffic_links(norway_directed_traffic_links, year = year)

stops_on_traffic_links_norway <- read.csv(paste0("data/trafikklenker_med_holdeplasser_", year, ".csv"))
bus_counts_norway <- read.csv(paste0("data/holdeplasspasseringer_entur_", year, ".csv"))

bus_aadt_norway <- calculate_bus_aadt(stops_on_traffic_links_norway, bus_counts_norway, year = year)

prepared_traffic_links_norway <- fill_missing_values(
  df = preprocessed_traffic_links_norway,
  unknown_impute_columns = c("functionClass", "highestSpeedLimit", "lowestSpeedLimit","maxLanes", "minLanes"),
  mode_impute_columns = c("hasOnlyPublicTransportLanes"),
  median_impute_columns = c("lastYearAadt_aadt", "lastYearAadt_heavyRatio", 
                            "lastYearAadt_heavyAadt")) |>
  remove_negative_aadt() |> 
  add_logLastYear() |>
  join_bus_to_traffic(bus_aadt_norway)

adjacency_matrix_norway <- build_adjacency_matrix(
  prepared_traffic_links_norway,
  exclude_public_transport = TRUE)

clusters <- strategic_network_clustering(
  data = prepared_traffic_links_norway,
  year = year,
  boundary_links = c("Trafikkdata_continuous")
)

nodes_norway <- identify_unbalanceable_nodes(norway_nodes_raw, prepared_traffic_links_norway)
covariates <- ~ functionalRoadClass:maxLanes +
  functionalRoadClass:roadCategory +
  minLanes:roadCategory + functionalRoadClass +
  maxLanes + roadCategory +
  hasOnlyPublicTransportLanes +
  functionalRoadClass*isRamp +
  lastYearAadt_logAadt

bal_node_ids <- nodes_norway$id[!nodes_norway$unbalanceable_node]

# Direct match: prepared_traffic_links_norway$id and norway_directed_traffic_links$id
# share the same format (e.g. "0.0-1.0@384223-AGAINST"), so no regex parsing is needed.
link_direction <- data.frame(
  id         = as.character(norway_directed_traffic_links$id),
  start_node = as.character(norway_directed_traffic_links$startTrafficNodeId),
  end_node   = as.character(norway_directed_traffic_links$endTrafficNodeId),
  stringsAsFactors = FALSE
)

link_idx <- match(prepared_traffic_links_norway$id, link_direction$id)

# Each directed link leaves start_node (outgoing) and enters end_node (incoming).
n_links <- nrow(prepared_traffic_links_norway)
incidence <- data.frame(
  row_idx   = c(seq_len(n_links), seq_len(n_links)),
  node_id   = c(link_direction$start_node[link_idx], link_direction$end_node[link_idx]),
  direction = c(rep("out", n_links), rep("in", n_links)),
  stringsAsFactors = FALSE
)

incidence <- incidence[!is.na(incidence$node_id) &
                         incidence$node_id %in% as.character(bal_node_ids), ]
message("Incidence rows: ", nrow(incidence),
        " | unique balanceable nodes: ", length(unique(incidence$node_id)))

build_flow_conservation_extraconstr <- function(data, incidence,
                                                weight_col = "lastYearAadt_aadt") {
  n         <- nrow(data)
  bal_nodes <- unique(incidence$node_id)
  m         <- length(bal_nodes)
  message("Flow conservation constraints: ", m, " balanceable nodes")
  
  # Build sparse matrix via (i, j, x) triplets — avoids allocating a dense m×n matrix
  node_row <- setNames(seq_along(bal_nodes), bal_nodes)
  rows  <- node_row[incidence$node_id]
  cols  <- incidence$row_idx
  ws    <- data[[weight_col]][cols]
  ws[is.na(ws)] <- 1
  signs <- ifelse(incidence$direction == "in", 1L, -1L)
  
  A <- Matrix::sparseMatrix(
    i    = rows,
    j    = cols,
    x    = as.numeric(ws * signs),
    dims = c(m, n)
  )
  list(A = A, e = rep(0, m))
}

flow_constr <- build_flow_conservation_extraconstr(
  data      = prepared_traffic_links_norway,
  incidence = incidence
)

# ── Shared model settings ─────────────────────────────────────────────────────
ordinal_levels_road <- list(
  functionClass     = c("unknown", "E", "D", "C", "B", "A"),
  highestSpeedLimit = c("unknown", "20", "30", "40", "50", "60", "70", "80", "90", "100", "110")
)
similarity_covariates_norway <- c("minLanes", "highestSpeedLimit", "functionClass", "lastYearAadt_logAadt")

# ── 80/20 holdout split ───────────────────────────────────────────────────────
set.seed(123)
measured_idx <- which(!is.na(prepared_traffic_links_norway$aadt))
test_idx     <- sample(measured_idx, size = round(0.2 * length(measured_idx)))

data_test <- prepared_traffic_links_norway
data_test$aadt[test_idx]      <- NA
data_test$heavyAadt[test_idx] <- NA

true_aadt <- prepared_traffic_links_norway$aadt[test_idx]

# ── Model 1: Weighted without flow conservation ──────────────────────────
message("Fitting baseline weighted model...")
fit_baseline <- fit_weighted_inla_model(
  data                  = data_test,
  adjacency_matrix      = adjacency_matrix_norway,
  spatial_term          = "besagproper_rbf",
  fixed_effects         = covariates,
  iid_effects           = "roadSystem",
  ordinal_levels        = ordinal_levels_road,
  similarity_covariates = similarity_covariates_norway,
  family                = "nbinomial"
)

# ── Model 2: Weighted with flow conservation extraconstr ─────────────────
message("Fitting flow-constrained weighted model...")
fit_constr <- fit_weighted_inla_model(
  data                  = data_test,
  adjacency_matrix      = adjacency_matrix_norway,
  spatial_term          = "besagproper_rbf",
  fixed_effects         = covariates,
  iid_effects           = "roadSystem",
  ordinal_levels        = ordinal_levels_road,
  similarity_covariates = similarity_covariates_norway,
  family                = "nbinomial",
  extraconstr           = flow_constr
)

# ── Metrics ───────────────────────────────────────────────────────────────────
pred_baseline <- fit_baseline$predictions$inla_pred[test_idx]
pred_constr   <- fit_constr$predictions$inla_pred[test_idx]

data.frame(
  model = c("Weighted (baseline)", "Weighted + flow conservation"),
  MALE  = c(.male(pred_baseline, true_aadt), .male(pred_constr, true_aadt)),
  MAPE  = c(.mape(pred_baseline, true_aadt), .mape(pred_constr, true_aadt)),
  RMSE  = c(.rmse(pred_baseline, true_aadt), .rmse(pred_constr, true_aadt)),
  MAE   = c(.mae( pred_baseline, true_aadt), .mae( pred_constr, true_aadt))
) |> print(digits = 4, row.names = FALSE)

# ── Plots ─────────────────────────────────────────────────────────────────────
plot_inla_model(fit_baseline, prepared_traffic_links_norway$aadt, test_idx,
                type = "pred_vs_obs", title = "Weighted (baseline)")
plot_inla_model(fit_constr,   prepared_traffic_links_norway$aadt, test_idx,
                type = "pred_vs_obs", title = "Weighted + flow conservation")




