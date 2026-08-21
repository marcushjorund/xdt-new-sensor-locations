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

# ── Load raw data ─────────────────────────────────────────────────────────────
norway_directed_traffic_links <- jsonlite::read_json(
  path = "data/directed-traffic-links-2024.json",
  simplifyVector = TRUE,
  simplifyDataFrame = TRUE
)
norway_nodes_raw <- sf::st_read("data/traffic-nodes-2024.geojson")
# ── Pre-simplified geometry file (create once, reuse across sessions) ─────────────
# Reduces coordinate precision to ~1 m (UTM 33N) before writing, so
# add_geometries() loads a much smaller file on every subsequent run.
simplified_geojson_path <- "data/directed_traffic_links_2024_simplified.geojson"
if (!file.exists(simplified_geojson_path)) {
  message("Creating simplified geometry file (one-time) ...")
  geo_full <- sf::st_read("data/directed_traffic_links_2024.geojson", quiet = TRUE)
  geo_simplified <- geo_full[, "id"] |>
    sf::st_transform(25833) |>                              # project to UTM 33N (metres)
    sf::st_simplify(dTolerance = 5, preserveTopology = TRUE) |>  # simplify to 5 m
    sf::st_transform(4326)                                  # back to WGS84
  # Round coordinates to 5 decimal places (~1 m at Norwegian latitudes)
  geo_simplified <- sf::st_set_precision(geo_simplified, precision = 1e5)
  sf::st_write(geo_simplified, simplified_geojson_path, delete_dsn = TRUE, quiet = TRUE)
  message("Saved: ", simplified_geojson_path,
          " (", round(file.size(simplified_geojson_path) / 1024^2, 1), " MB)")
  rm(geo_full, geo_simplified)
}
# ── Preprocess ────────────────────────────────────────────────────────────────
preprocessed_traffic_links_norway <- preprocess_traffic_links(norway_directed_traffic_links, year = year)

stops_on_traffic_links_norway <- read.csv(paste0("data/trafikklenker_med_holdeplasser_", year, ".csv"))
bus_counts_norway <- read.csv(paste0("data/holdeplasspasseringer_entur_", year, ".csv"))

bus_aadt_norway <- calculate_bus_aadt(stops_on_traffic_links_norway, bus_counts_norway, year = year)

prepared_traffic_links_norway <- fill_missing_values(
  df = preprocessed_traffic_links_norway,
  unknown_impute_columns = c("functionClass", "highestSpeedLimit", "lowestSpeedLimit", "maxLanes", "minLanes"),
  mode_impute_columns = c("hasOnlyPublicTransportLanes"),
  median_impute_columns = c("lastYearAadt_aadt", "lastYearAadt_heavyRatio",
                            "lastYearAadt_heavyAadt")) |>
  remove_negative_aadt() |>
  add_logLastYear() |>
  join_bus_to_traffic(bus_aadt_norway)

adjacency_matrix_norway <- build_adjacency_matrix(
  prepared_traffic_links_norway,
  exclude_public_transport = TRUE)

nodes_norway <- identify_unbalanceable_nodes(norway_nodes_raw, prepared_traffic_links_norway)

# ── INLA model ────────────────────────────────────────────────────────────────
covariates <- ~ functionalRoadClass:maxLanes +
  functionalRoadClass:roadCategory +
  minLanes:roadCategory + functionalRoadClass +
  maxLanes + roadCategory +
  hasOnlyPublicTransportLanes +
  functionalRoadClass*isRamp +
  lastYearAadt_logAadt

ordinal_levels_road <- list(
  functionClass       = c("unknown", "E", "D", "C", "B", "A"),
  highestSpeedLimit   = c("unknown", "20", "30", "40", "50", "60", "70", "80", "90", "100", "110")
)

similarity_covariates_norway <- c(
  "minLanes", "highestSpeedLimit", "functionClass",
  "lastYearAadt_logAadt")

inla_weighted_norway_full <- fit_weighted_inla_model(
  data                  = prepared_traffic_links_norway,
  adjacency_matrix      = adjacency_matrix_norway,
  spatial_term          = "besagproper_rbf",
  fixed_effects         = covariates,
  iid_effects           = "roadSystem",
  ordinal_levels        = ordinal_levels_road,
  similarity_covariates = similarity_covariates_norway,
  family                = "nbinomial")


spatial_hyperpar <- inla_weighted_norway_full$spatial_hyperpar[2:4, "mean"]
tau   <- spatial_hyperpar[1]
d     <- spatial_hyperpar[2]
sigma <- spatial_hyperpar[3]

distances <- inla_weighted_norway_full$distances

# ── Identify derivable nodes (flow conservation) ──────────────────────────────
# Enriches the data frame with enables_derivable_links, derivable, n_derivable,
# enables_derivable, and n_enables_derivable columns.
incidence_matrix_norway <- build_incidence_matrix(nodes_norway, prepared_traffic_links_norway, nodes_to_balance = "complete nodes")

prepared_traffic_links_norway_deriv <- identify_derivable_nodes(
  incidence_matrix = incidence_matrix_norway,
  traffic_links    = prepared_traffic_links_norway,
  nodes            = nodes_norway)

# ── Join model predictions onto enriched data ─────────────────────────────────
prepared_traffic_links_norway_deriv <- dplyr::full_join(
  prepared_traffic_links_norway_deriv,
  inla_weighted_norway_full$predictions
)

# ── Sensor 1: Unweighted, bundle scoring (raw, no added information) ────────────
norway_sensors_unweighted_bundle <- greedy_mi_sensor_selection_norway(
  data                = prepared_traffic_links_norway_deriv,
  adjacency_matrix    = adjacency_matrix_norway,
  distances           = distances,
  tau                 = tau,
  d                   = d,
  sigma               = sigma,
  hops                = 3,
  k                   = 100,
  r                   = 100,
  weighting_vars      = NULL,
  neighbourhood_decay = c(1, 0.5, 0.25, 0.125),
  neighbour_hops      = 2,
  bundle_scoring      = TRUE
)

saveRDS(norway_sensors_unweighted_bundle,
        file.path("results", "unweighted_bundle_100percounty_100overall.rds"))
# Fixed algorithm (dynamic-sum neighbourhood weighting) — for Shiny comparison.
saveRDS(norway_sensors_unweighted_bundle,
        file.path("results", "unweighted_bundle_distsum_100percounty_100overall.rds"))

plot_sensor_selection_map(norway_sensors_unweighted_bundle)

# ── Sensor 2: FRC (60%) + log-prop AADT (30%) + heavy ratio (10%), bundle scoring ──
norway_sensors_frc_aadt_heavy <- greedy_mi_sensor_selection_norway(
  data                = prepared_traffic_links_norway_deriv,
  adjacency_matrix    = adjacency_matrix_norway,
  distances           = distances,
  tau                 = tau,
  d                   = d,
  sigma               = sigma,
  hops                = 3,
  k                   = 100,
  r                   = 100,
  weighting_vars      = list(
    functionalRoadClass = setNames(sqrt(seq(8, 1, by = -1)), as.character(0:7)),
    inla_pred           = "log_proportional",
    lastYearAadt_heavyRatio = "proportional"
  ),
  weighting_vars_alpha = c(0.6, 0.3, 0.1),
  neighbourhood_decay  = c(1, 0.5, 0.25, 0.125),
  neighbour_hops       = 2,
  bundle_scoring       = TRUE
)

saveRDS(norway_sensors_frc_aadt_heavy,
        file.path("results", "frc_0.6_logaadt_0.3_heavyratio_0.1_100percounty_100overall.rds"))

plot_sensor_selection_map(norway_sensors_frc_aadt_heavy)
