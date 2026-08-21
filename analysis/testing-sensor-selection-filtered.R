# testing-sensor-selection-filtered.R
#
# Tests the filtering argument in greedy_mi_sensor_selection_norway().
# Use case: optimal sensor placement on Fylkesvei roads in Nordland county.
#
# Weighting: functionalRoadClass (sqrt 0.9→0, 80%)
#            + inla_pred log_proportional (15%)
#            + lastYearAadt_heavyRatio identity (5%)
#
# Two tests:
#   Test 1 — Norway-wide data + filtering = list(county, roadCategory)
#   Test 2 — Nordland-only data + filtering = list(roadCategory)
#   Both should select the same sensor locations.

library(INLA)
library(xdtkit)
library(jsonlite)
source("R/scripts_inla_sensor_placement.R")

year <- 2024

# ── Load / build prepared traffic links ──────────────────────────────────────
# On a cold start, rebuild the full pipeline (mirrors testing-sensor-selection-norway.R).
# After the first run the INLA cache file makes subsequent runs fast.

norway_directed_traffic_links <- jsonlite::read_json(
  path = "data/directed-traffic-links-2024.json",
  simplifyVector    = TRUE,
  simplifyDataFrame = TRUE
)
norway_nodes_raw <- sf::st_read("data/traffic-nodes-2024.geojson")

preprocessed_traffic_links_norway <- preprocess_traffic_links(
  norway_directed_traffic_links,
  year = year
)

stops_on_traffic_links_norway <- read.csv(
  paste0("data/trafikklenker_med_holdeplasser_", year, ".csv")
)
bus_counts_norway <- read.csv(
  paste0("data/holdeplasspasseringer_entur_", year, ".csv")
)

bus_aadt_norway <- calculate_bus_aadt(
  stops_on_traffic_links_norway,
  bus_counts_norway,
  year = year
)

prepared_traffic_links_norway <- fill_missing_values(
  df = preprocessed_traffic_links_norway,
  unknown_impute_columns = c("functionClass", "highestSpeedLimit", "lowestSpeedLimit",
                             "maxLanes", "minLanes"),
  mode_impute_columns    = c("hasOnlyPublicTransportLanes"),
  median_impute_columns  = c("lastYearAadt_aadt", "lastYearAadt_heavyRatio",
                             "lastYearAadt_heavyAadt")
) |>
  remove_negative_aadt() |>
  add_logLastYear() |>
  join_bus_to_traffic(bus_aadt_norway)

adjacency_matrix_norway <- build_adjacency_matrix(
  prepared_traffic_links_norway,
  exclude_public_transport = TRUE
)

nodes_norway <- identify_unbalanceable_nodes(
  norway_nodes_raw,
  prepared_traffic_links_norway
)

# ── INLA model (cache to disk after first run) ────────────────────────────────
inla_cache_path <- "data/inla_norway_hyperpar_cache.rds"

if (file.exists(inla_cache_path)) {
  message("Loading cached INLA hyperparameters from ", inla_cache_path)
  inla_cache <- readRDS(inla_cache_path)
  tau       <- inla_cache$tau
  d         <- inla_cache$d
  sigma     <- inla_cache$sigma
  distances <- inla_cache$distances
  prepared_traffic_links_norway <- dplyr::full_join(
    prepared_traffic_links_norway,
    inla_cache$predictions
  )
} else {
  message("Fitting INLA model (this may take 20–60 minutes)...")
  covariates <- ~ functionalRoadClass:maxLanes +
    functionalRoadClass:roadCategory +
    minLanes:roadCategory + functionalRoadClass +
    maxLanes + roadCategory +
    hasOnlyPublicTransportLanes +
    functionalRoadClass * isRamp +
    lastYearAadt_logAadt

  ordinal_levels_road <- list(
    functionClass     = c("unknown", "E", "D", "C", "B", "A"),
    highestSpeedLimit = c("unknown", "20", "30", "40", "50", "60", "70", "80", "90", "100", "110")
  )

  similarity_covariates_norway <- c(
    "minLanes", "highestSpeedLimit", "functionClass", "lastYearAadt_logAadt"
  )

  inla_weighted_norway_full <- fit_weighted_inla_model(
    data                  = prepared_traffic_links_norway,
    adjacency_matrix      = adjacency_matrix_norway,
    spatial_term          = "besagproper_rbf",
    fixed_effects         = covariates,
    iid_effects           = "roadSystem",
    ordinal_levels        = ordinal_levels_road,
    similarity_covariates = similarity_covariates_norway,
    family                = "nbinomial"
  )

  spatial_hyperpar <- inla_weighted_norway_full$spatial_hyperpar[2:4, "mean"]
  tau       <- spatial_hyperpar[1]
  d         <- spatial_hyperpar[2]
  sigma     <- spatial_hyperpar[3]
  distances <- inla_weighted_norway_full$distances

  prepared_traffic_links_norway <- dplyr::full_join(
    prepared_traffic_links_norway,
    inla_weighted_norway_full$predictions
  )

  saveRDS(
    list(
      tau         = tau,
      d           = d,
      sigma       = sigma,
      distances   = distances,
      predictions = inla_weighted_norway_full$predictions
    ),
    inla_cache_path
  )
  message("INLA cache saved to ", inla_cache_path)
}

# ── Identify derivable links via flow conservation ────────────────────────────
incidence_matrix_norway <- build_incidence_matrix(
  nodes_norway,
  prepared_traffic_links_norway,
  nodes_to_balance = "complete nodes"
)
prepared_traffic_links_norway <- identify_derivable_nodes(
  incidence_matrix = incidence_matrix_norway,
  traffic_links    = prepared_traffic_links_norway,
  nodes            = nodes_norway
)

# ── Shared weighting specification ───────────────────────────────────────────
# 80% functionalRoadClass, 15% log-proportional predicted AADT, 5% heavy ratio
weighting_filtered <- list(
  functionalRoadClass     = setNames(sqrt(seq(0.9,0, by = -0.1)), as.character(0:9)),
  inla_pred               = "proportional",
  lastYearAadt_heavyRatio = "identity"
)
weighting_alpha_filtered <- c(0.6, 0.30, 0.10)

# ══════════════════════════════════════════════════════════════════════════════
# Test 1: Norway-wide data + filtering = list(county, roadCategory)
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Test 1: Norway-wide data, filter to Nordland FYLKESVEG ──────────────")

norway_sensors_nordland_fylkes_1 <- greedy_mi_sensor_selection_norway(
  data                 = prepared_traffic_links_norway,
  adjacency_matrix     = adjacency_matrix_norway,
  distances            = distances,
  tau                  = tau,
  d                    = d,
  sigma                = sigma,
  hops                 = 3,
  k                    = 15,
  r                    = 15,
  filtering            = list(
    county       = "Nordland",
    roadCategory = "FYLKESVEG"
  ),
  weighting_vars       = weighting_filtered,
  weighting_vars_alpha = weighting_alpha_filtered,
  neighbourhood_decay  = c(1, 0.5, 0.25, 0.125),
  neighbour_hops       = 2
)
library(dplyr)
saveRDS(
  prepared_traffic_links_norway,
  "data/prepared_traffic_links_norway.rds"
)
saveRDS(
  norway_sensors_nordland_fylkes_1,
  "results/analysis/norway_sensors_15_nordland_fylkes_1.rds"
)

saveRDS(
  prepared_traffic_links_norway %>% filter(county == "Nordland"),
  "data/prepared_traffic_links_nordland.rds"
)

plot_sensor_selection_map(
  norway_sensors_nordland_fylkes_1,
  all_data = prepared_traffic_links_norway %>%filter(county == "Nordland")
)

norway_sensors_nordland_fylkes_1$summary
norway_sensors_nordland_fylkes_1$summary$counts_table
norway_sensors_nordland_fylkes_1$selected_data_entries

saveRDS(
  norway_sensors_nordland_fylkes_1,
  "results/analysis/nordland_fylkesveg_frc_logaadt_heavyratio_100k100r_norway_data.rds"
)
