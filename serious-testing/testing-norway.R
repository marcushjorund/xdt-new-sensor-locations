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

bus_aadt_norway <- calculate_bus_aadt(stops_on_traffic_links_norway, bus_counts, year = year)

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
#covariates <- ~ maxLanes + hasOnlyPublicTransportLanes + functionalRoadClass + lastYearAadt_logAadt

inla_model_norway <- fit_inla_model(
  data = prepared_traffic_links_norway,
  adjacency_matrix_norway,
  fixed_effects = covariates,
  iid_effects = "roadSystem",
  family = "nbinomial")

predictions_total_norway <- dplyr::full_join(prepared_traffic_links_norway, inla_model_norway$predictions)

plot(predictions_total_norway$lastYearAadt_aadt, predictions_total_norway$inla_pred)

set.seed(123)
measured_idx_norway <- which(!is.na(prepared_traffic_links_norway$aadt))
test_idx_norway     <- sample(measured_idx_norway, size = round(0.2 * length(measured_idx_norway)))

data_cv_norway <- prepared_traffic_links_norway
data_cv_norway$aadt[test_idx_norway]      <- NA
data_cv_norway$heavyAadt[test_idx_norway] <- NA

true_aadt_norway <- prepared_traffic_links_norway$aadt[test_idx_norway]

inla_model_norway_test <- fit_inla_model(
  data = data_cv_norway,
  adjacency_matrix_norway,
  fixed_effects = covariates,
  iid_effects = "roadSystem",
  family = "nbinomial")

source("serious-testing/scripts_inla_sensor_placement.R")

ordinal_levels_road <- list(
  functionalRoadClass = as.character(c(7, 6, 5, 4, 3, 2, 1, 0)),
  roadCategory        = c("KOMMUNAL_VEG", "FYLKESVEG", "RIKSVEG", "EUROPAVEG")
)

similarity_covariates_norway <- c(
  "maxLanes", "minLanes",
  "highestSpeedLimit", "lowestSpeedLimit",
  "hasOnlyPublicTransportLanes", "isRamp",
  "roadCategory", "functionalRoadClass", "roadSystem",
  "lastYearAadt_logAadt"
)

inla_rbf_norway <- fit_inla_rbf_model(
  data                  = data_cv_norway,
  adjacency_matrix      = adjacency_matrix_norway,
  spatial_term          = "besagproper_rbf",
  fixed_effects         = ~ functionalRoadClass + roadCategory + functionalRoadClass*isRamp + hasOnlyPublicTransportLanes + maxLanes + minLanes + lastYearAadt_logAadt,
  iid_effects           = "roadSystem",
  ordinal_levels        = ordinal_levels_road,
  similarity_covariates = c("lastYearAadt_logAadt", "functionalRoadClass", "roadCategory"),
  family                = "nbinomial")

# ── Diagnostic plots ──────────────────────────────────────────────────────────
plot_inla_model(inla_model_norway_test, prepared_traffic_links_norway$aadt, test_idx_norway,
                type = "pred_vs_obs",  title = "Besagproper")
plot_inla_model(inla_rbf_norway,       prepared_traffic_links_norway$aadt, test_idx_norway,
                type = "pred_vs_obs",  title = "RBF-weighted Besagproper")

plot_inla_model(inla_model_norway_test, prepared_traffic_links_norway$aadt, test_idx_norway,
                type = "qq",           title = "Besagproper")
plot_inla_model(inla_rbf_norway,        prepared_traffic_links_norway$aadt, test_idx_norway,
                type = "qq",           title = "RBF-weighted Besagproper")

# ── Metrics comparison ────────────────────────────────────────────────────────
pred_besag <- inla_model_norway_test$predictions$inla_pred[test_idx_norway]
pred_rbf   <- inla_rbf_norway$predictions$inla_pred[test_idx_norway]

# Exclude zero-AADT observations from MAPE (division by zero → Inf)
nonzero <- true_aadt_norway > 0

metrics_norway <- data.frame(
  model = c("Besagproper", "RBF-weighted Besagproper"),
  RMSE  = c(.rmse(pred_besag,           true_aadt_norway),         .rmse(pred_rbf,           true_aadt_norway)),
  MAE   = c(.mae( pred_besag,           true_aadt_norway),         .mae( pred_rbf,           true_aadt_norway)),
  MAPE  = c(.mape(pred_besag[nonzero],  true_aadt_norway[nonzero]),.mape(pred_rbf[nonzero],  true_aadt_norway[nonzero]))
)
print(metrics_norway, digits = 4, row.names = FALSE)
