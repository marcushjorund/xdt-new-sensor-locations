# ═══════════════════════════════════════════════════════════════════════════════
# PART 1: Data Preparation & Saving
# Run this section once to build and persist all intermediate objects needed for
# the sensor selection experiments in Part 2.
# ═══════════════════════════════════════════════════════════════════════════════

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

nodes_norway <- identify_unbalanceable_nodes(norway_nodes_raw, prepared_traffic_links_norway)
covariates <- ~ functionalRoadClass:maxLanes +
  functionalRoadClass:roadCategory +
  minLanes:roadCategory + functionalRoadClass +
  maxLanes + roadCategory +
  hasOnlyPublicTransportLanes +
  functionalRoadClass*isRamp +
  lastYearAadt_logAadt

#Defining the ordinal levels and their ordering
ordinal_levels_road <- list(
  functionClass       = c("unknown", "E", "D", "C", "B", "A"),
  highestSpeedLimit   = c("unknown","20", "30", "40", "50", "60", "70", "80", "90", "100", "110")
)
#Best similarity covariates
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


# ── Join model predictions onto traffic links for AADT-based weighting ────────
prepared_traffic_links_norway <- dplyr::full_join(
  prepared_traffic_links_norway,
  inla_weighted_norway_full$predictions
)

# ── Identify derivable links via flow conservation ────────────────────────────
# Enriches prepared_traffic_links_norway with:
#   enables_derivable, enables_derivable_links, n_enables_derivable  (sensor placement value)
#   derivable, derivable_flow_nodes, n_derivable                     (already derivable now)
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

#Saving data so that we don't have to rerun the preprocessing multiple times
saveRDS(prepared_traffic_links_norway, "analysis/model_cache/traffic_links_norway_with_derivable.rds")
saveRDS(adjacency_matrix_norway, "analysis/model_cache/adjacency_matrix_norway.rds")
saveRDS(inla_weighted_norway_full, "analysis/model_cache/inla_rbf_norway_full.rds")


# ── Save minimal measured links for Shiny app overlay ────────────────────────
saveRDS(
  prepared_traffic_links_norway[
    !is.na(prepared_traffic_links_norway$aadt),
    intersect(c("id", "aadt", "traffic_volume_source"),
              names(prepared_traffic_links_norway)),
    drop = FALSE
  ],
  "data/measured_traffic_links_minimal.rds"
)

# ═══════════════════════════════════════════════════════════════════════════════
# PART 2: Sensor Selection Experiments
# Loads the objects saved in Part 1 and runs greedy MI sensor selection under
# various weighting schemes. This section can be run independently of Part 1
# once its output files exist.
# ═══════════════════════════════════════════════════════════════════════════════

library(INLA)
library(xdtkit)
source("R/scripts_inla_sensor_placement.R")

prepared_traffic_links_norway <- readRDS("analysis/model_cache/traffic_links_norway_with_derivable.rds")
adjacency_matrix_norway       <- readRDS("analysis/model_cache/adjacency_matrix_norway.rds")
inla_weighted_norway_full          <- readRDS("analysis/model_cache/inla_rbf_norway_full.rds")

spatial_hyperpar <- inla_weighted_norway_full$spatial_hyperpar[2:4, "mean"]
tau   <- spatial_hyperpar[1]
d     <- spatial_hyperpar[2]
sigma <- spatial_hyperpar[3]

distances <- inla_weighted_norway_full$distances

# ── Select top-k sensors per county, top-r across all of Norway ──────────────
norway_sensors_1 <- greedy_mi_sensor_selection_norway(
  data                = prepared_traffic_links_norway,
  adjacency_matrix    = adjacency_matrix_norway,
  distances           = distances,
  tau                 = tau,
  d                   = d,
  sigma               = sigma,
  hops                = 3,
  k                   = 70,           # sensors selected per county
  r                   = 100,           # set e.g. r = 500 for Norway-wide top-r
  weighting_vars      = list(
    functionalRoadClass = setNames(sqrt(seq(0.9, 0, by = -0.1)), as.character(0:9))
  ),
  neighbourhood_decay = c(1, 0.5,0.25),
  neighbour_hops      = 2
)

selected_sensor_placements_1 <- norway_sensors_1$selected_data_entries

# ── Plot sensor selection map ─────────────────────────────────────────────────
plot_sensor_selection_map(norway_sensors_1, all_data = prepared_traffic_links_norway)

# Statistics of selected points
norway_sensors_1$summary
norway_sensors_1$summary$counts_table

# Saving results
saveRDS(norway_sensors_1, "results/production/frc_sqrt(0.9_to_0)_100percounty_100overall.rds")

# ── Select top-k sensors per county, top-r across all of Norway ──────────────
norway_sensors_2 <- greedy_mi_sensor_selection_norway(
  data                = prepared_traffic_links_norway,
  adjacency_matrix    = adjacency_matrix_norway,
  distances           = distances,
  tau                 = tau,
  d                   = d,
  sigma               = sigma,
  hops                = 3,
  k                   = 100,           # sensors selected per county
  r                   = 100,           # set e.g. r = 500 for Norway-wide top-r
  weighting_vars      = list(
    functionalRoadClass = setNames(1/seq(1, 10, by = 1), as.character(0:9))
  ),
  neighbourhood_decay = c(1, 0.5, 0.25),
  neighbour_hops      = 2
)

selected_sensor_placements_2 <- norway_sensors_2$selected_data_entries


# ── Plot sensor selection map ─────────────────────────────────────────────────
plot_sensor_selection_map(norway_sensors_2, all_data = prepared_traffic_links_norway)

# Statistics of selected points
norway_sensors_2$summary

# Saving results
saveRDS(norway_sensors_2, "results/production/frc_inverse_100percounty_100overall.rds")

# The inverse frc weighting selected exclusively points of high frc, while the
# sqrt(0.9_to_0) weighting selected a more balanced mix of points across frc
# levels, with a clear preference for higher frc but not exclusively so.
# This suggests that the inverse weighting may be too extreme in prioritizing
# high-frc links, while the sqrt(0.9_to_0) weighting provides a more nuanced
# balance that still favors higher-frc links without completely excluding
# lower-frc ones.


# ── Select top-k sensors per county, top-r across all of Norway ──────────────
norway_sensors_3 <- greedy_mi_sensor_selection_norway(
  data                 = prepared_traffic_links_norway,
  adjacency_matrix     = adjacency_matrix_norway,
  distances            = distances,
  tau                  = tau,
  d                    = d,
  sigma                = sigma,
  hops                 = 3,
  k                    = 100,           # sensors selected per county
  r                    = 100,           # set e.g. r = 500 for Norway-wide top-r
  weighting_vars       = list(
    functionalRoadClass = setNames(sqrt(seq(0.9, 0, by = -0.1)), as.character(0:9)),
    inla_pred           = "log_proportional"
  ),
  weighting_vars_alpha = c(0.9, 0.1),   # 90% FRC, 10% predicted AADT
  neighbourhood_decay  = c(1, 0.5, 0.25),
  neighbour_hops       = 2
)

selected_sensor_placements_3 <- norway_sensors_3$selected_data_entries


# ── Plot sensor selection map ─────────────────────────────────────────────────
plot_sensor_selection_map(norway_sensors_3, all_data = prepared_traffic_links_norway)

# Statistics of selected points
norway_sensors_3$summary

# Saving results
saveRDS(norway_sensors_3, "results/production/frc_sqrt(0.9_to_0)_and_logAADT_0.9_to_0.1_weight_100percounty_100overall.rds")

# ── Select top-k sensors per county, top-r across all of Norway ──────────────
norway_sensors_4 <- greedy_mi_sensor_selection_norway(
  data                = prepared_traffic_links_norway,
  adjacency_matrix    = adjacency_matrix_norway,
  distances           = distances,
  tau                 = tau,
  d                   = d,
  sigma               = sigma,
  hops                = 3,
  k                   = 10,           # sensors selected per county
  r                   = 150,           # set e.g. r = 500 for Norway-wide top-r
  weighting_vars      = list(
    functionalRoadClass = setNames(sqrt(seq(0.9, 0, by = -0.1)), as.character(0:9))
  ),
  neighbourhood_decay = c(1, 0.5, 0.25 ),
  neighbour_hops      = 2
)

selected_sensor_placements_4 <- norway_sensors_4$selected_data_entries


# ── Plot sensor selection map ─────────────────────────────────────────────────
plot_sensor_selection_map(norway_sensors_4, all_data = prepared_traffic_links_norway)

# Statistics of selected points
norway_sensors_4$summary

# Saving results
saveRDS(norway_sensors_4, "results/production/frc_sqrt(0.9_to_0)_10percounty_150overall.rds")

# ── Select top-k sensors per county, top-r across all of Norway ──────────────
norway_sensors_5 <- greedy_mi_sensor_selection_norway(
  data                = prepared_traffic_links_norway,
  adjacency_matrix    = adjacency_matrix_norway,
  distances           = distances,
  tau                 = tau,
  d                   = d,
  sigma               = sigma,
  hops                = 3,
  k                   = 100,           # sensors selected per county
  r                   = 100,           # set e.g. r = 500 for Norway-wide top-r
  weighting_vars      = NULL,
  neighbourhood_decay = c(1, 0.5, 0.25),
  neighbour_hops      = 2
)

selected_sensor_placements_5 <- norway_sensors_5$selected_data_entries


# ── Plot sensor selection map ─────────────────────────────────────────────────
plot_sensor_selection_map(norway_sensors_5, all_data = prepared_traffic_links_norway)

# Statistics from selected sensors
norway_sensors_5$summary

# Saving results
saveRDS(norway_sensors_5, "results/production/unweighted_100percounty100overall.rds")

# ── Select top-k sensors per county, top-r across all of Norway ──────────────
norway_sensors_6 <- greedy_mi_sensor_selection_norway(
  data                = prepared_traffic_links_norway,
  adjacency_matrix    = adjacency_matrix_norway,
  distances           = distances,
  tau                 = tau,
  d                   = d,
  sigma               = sigma,
  hops                = 3,
  k                   = 10,           # sensors selected per county
  r                   = 150,           # set e.g. r = 500 for Norway-wide top-r
  weighting_vars      = list(
    functionalRoadClass = setNames(sqrt(seq(0.9, 0, by = -0.1)), as.character(0:9)),
    inla_pred         = "log_proportional"
  ),
  weighting_vars_alpha = c(0.9, 0.1),   # 90% FRC, 10% predicted AADT
  neighbourhood_decay = c(1, 0.5, 0.25),
  neighbour_hops      = 2
)

selected_sensor_placements_6 <- norway_sensors_6$selected_data_entries


# ── Plot sensor selection map ─────────────────────────────────────────────────
plot_sensor_selection_map(norway_sensors_6, all_data = prepared_traffic_links_norway)

# Statistics of selected points
norway_sensors_6$summary

# Saving results
saveRDS(norway_sensors_6, "results/production/frc_sqrt(0.9_to_0)_logAADT_0.9_0.1_weight_10percounty_150overall.rds")

# ── Select top-k sensors per county, top-r across all of Norway ──────────────
norway_sensors_7 <- greedy_mi_sensor_selection_norway(
  data                = prepared_traffic_links_norway,
  adjacency_matrix    = adjacency_matrix_norway,
  distances           = distances,
  tau                 = tau,
  d                   = d,
  sigma               = sigma,
  hops                = 3,
  k                   = 100,           # sensors selected per county
  r                   = 100,           # set e.g. r = 500 for Norway-wide top-r
  weighting_vars      = list(
    inla_pred         = "log_proportional"
  ),
  neighbourhood_decay = c(1, 0.5, 0.25),
  neighbour_hops      = 2
)

selected_sensor_placements_7 <- norway_sensors_7$selected_data_entries


# ── Plot sensor selection map ─────────────────────────────────────────────────
plot_sensor_selection_map(norway_sensors_7, all_data = prepared_traffic_links_norway)

# Statistics from selected points
norway_sensors_7$summary

# Saving results
saveRDS(norway_sensors_7, "results/production/AADT_logprop_100percounty_100overall.rds")

# ── Select top-k sensors per county, top-r across all of Norway ──────────────
norway_sensors_8 <- greedy_mi_sensor_selection_norway(
  data                = prepared_traffic_links_norway,
  adjacency_matrix    = adjacency_matrix_norway,
  distances           = distances,
  tau                 = tau,
  d                   = d,
  sigma               = sigma,
  hops                = 3,
  k                   = 100,           # sensors selected per county
  r                   = 100,           # set e.g. r = 500 for Norway-wide top-r
  weighting_vars      = list(
    inla_pred         = "identity"
  ),
  neighbourhood_decay = c(1, 0.5, 0.25),
  neighbour_hops      = 2
)

selected_sensor_placements_8 <- norway_sensors_8$selected_data_entries


# ── Plot sensor selection map ─────────────────────────────────────────────────
plot_sensor_selection_map(norway_sensors_8, all_data = prepared_traffic_links_norway)

# Statistics from selected points
norway_sensors_8$summary

# Saving results
saveRDS(norway_sensors_8, "results/production/AADT_identity_100percounty_100overall.rds")


# ── Select top-k sensors per county, top-r across all of Norway ──────────────
norway_sensors_9 <- greedy_mi_sensor_selection_norway(
  data                = prepared_traffic_links_norway,
  adjacency_matrix    = adjacency_matrix_norway,
  distances           = distances,
  tau                 = tau,
  d                   = d,
  sigma               = sigma,
  hops                = 3,
  k                   = 100,           # sensors selected per county
  r                   = 100,           # set e.g. r = 500 for Norway-wide top-r
  weighting_vars      = list(
    inla_pred         = "proportional"
  ),
  weighting_vars_alpha = NULL,   # scale down aadt influence
  neighbourhood_decay = c(1, 0.5, 0.25),
  neighbour_hops      = 2
)

selected_sensor_placements_9 <- norway_sensors_9$selected_data_entries


# ── Plot sensor selection map ─────────────────────────────────────────────────
plot_sensor_selection_map(norway_sensors_9, all_data = prepared_traffic_links_norway)

# Statistics from selected points
norway_sensors_9$summary

# Saving results
saveRDS(norway_sensors_9, "results/production/AADT_prop_100percounty_100overall.rds")

# ── Select top-k sensors per county, top-r across all of Norway ──────────────
norway_sensors_10 <- greedy_mi_sensor_selection_norway(
  data                = prepared_traffic_links_norway,
  adjacency_matrix    = adjacency_matrix_norway,
  distances           = distances,
  tau                 = tau,
  d                   = d,
  sigma               = sigma,
  hops                = 3,
  k                   = 100,           # sensors selected per county
  r                   = 100,           # set e.g. r = 500 for Norway-wide top-r
  weighting_vars      = list(
    functionClass      = setNames(sqrt(seq(0.6, 0, length.out = 6)), c("A", "B", "C", "D", "E", "unknown"))
  ),
  weighting_vars_alpha = NULL,
  neighbourhood_decay = c(1, 0.5, 0.25),
  neighbour_hops      = 2
)
selected_sensor_placements_10 <- norway_sensors_10$selected_data_entries

# ── Plot sensor selection map ─────────────────────────────────────────────────
plot_sensor_selection_map(norway_sensors_10, all_data = prepared_traffic_links_norway)

# Statistics from selected points
norway_sensors_10$summary

# Saving results
saveRDS(norway_sensors_10, "results/production/functionClass_100percounty_100overall.rds")

# ── Select top-k sensors per county, top-r across all of Norway ──────────────
norway_sensors_11 <- greedy_mi_sensor_selection_norway(
  data                = prepared_traffic_links_norway,
  adjacency_matrix    = adjacency_matrix_norway,
  distances           = distances,
  tau                 = tau,
  d                   = d,
  sigma               = sigma,
  hops                = 3,
  k                   = 100,           # sensors selected per county
  r                   = 100,           # set e.g. r = 500 for Norway-wide top-r
  weighting_vars      = list(
    functionClass           = setNames(sqrt(seq(0.6, 0, by = -0.1)), c("A", "B", "C", "D", "E", "unknown")),
    inla_pred               = "log_proportional",
    lastYearAadt_heavyRatio = "identity"
  ),
  weighting_vars_alpha = c(0.8,0.15,0.05),   # scale down aadt influence
  neighbourhood_decay = c(1, 0.5, 0.25 ),
  neighbour_hops      = 2
)
selected_sensor_placements_11 <- norway_sensors_11$selected_data_entries
# ── Plot sensor selection map ─────────────────────────────────────────────────
plot_sensor_selection_map(norway_sensors_11, all_data = prepared_traffic_links_norway)

# Statistics from selected points
norway_sensors_11$summary

# Saving results
saveRDS(norway_sensors_11, "results/production/functionClass_sqrt(0.6_to_0)logpropaadt_heavy_08015005_100percounty_100overall.rds")

# ── Select top-k sensors per county, top-r across all of Norway ──────────────
norway_sensors_12 <- greedy_mi_sensor_selection_norway(
  data                = prepared_traffic_links_norway,
  adjacency_matrix    = adjacency_matrix_norway,
  distances           = distances,
  tau                 = tau,
  d                   = d,
  sigma               = sigma,
  hops                = 3,
  k                   = 100,           # sensors selected per county
  r                   = 100,           # set e.g. r = 500 for Norway-wide top-r
  weighting_vars      = list(
    functionalRoadClass           = setNames(sqrt(seq(0.9, 0, by = -0.1)), as.character(0:9)),
    inla_pred               = "log_proportional",
    lastYearAadt_heavyRatio = "identity"
  ),
  weighting_vars_alpha = c(0.8,0.15,0.05),   # scale down aadt influence
  neighbourhood_decay = c(1, 0.5, 0.25 ),
  neighbour_hops      = 2
)
selected_sensor_placements_12 <- norway_sensors_12$selected_data_entries

# ── Plot sensor selection map ─────────────────────────────────────────────────
plot_sensor_selection_map(norway_sensors_12, all_data = prepared_traffic_links_norway)

# Statistics from selected points
norway_sensors_12$summary

# Saving results
saveRDS(norway_sensors_12, "results/production/frc_sqrt(0.9_to_0)logpropaadt_heavy_08015005_100percounty_100overall.rds")

# ── Experiment: treat non-continuous sources as unmeasured ───────────────────
continuous_sources <- c("Trafikkdata_continuous", "AutoPASS")

prepared_traffic_links_experiment <- prepared_traffic_links_norway
non_continuous_mask <- !is.na(prepared_traffic_links_experiment$aadt) &
  !(prepared_traffic_links_experiment$traffic_volume_source %in% continuous_sources)
prepared_traffic_links_experiment$aadt[non_continuous_mask] <- NA
# traffic_volume_source is intentionally kept — visible in selected_data_entries

norway_sensors_experiment <- greedy_mi_sensor_selection_norway(
  data             = prepared_traffic_links_experiment,
  adjacency_matrix = adjacency_matrix_norway,
  distances        = distances,
  tau = tau, d = d, sigma = sigma,
  hops = 3, k = 100, r = 100,
  weighting_vars      = list(
    functionalRoadClass           = setNames(seq(0.9, 0, by = -0.1), as.character(0:9)),
    inla_pred               = "log_proportional",
    lastYearAadt_heavyRatio = "identity"
  ),
  weighting_vars_alpha = c(0.8,0.15,0.05),
  neighbourhood_decay = c(1, 0.5, 0.25 ),
  neighbour_hops = 2
)

# Inspect which selected locations already have periodic registrations:
norway_sensors_experiment$selected_data_entries |>
  dplyr::filter(selected, !is.na(traffic_volume_source)) |>
  dplyr::select(id, traffic_volume_source, mi_score)

plot_sensor_selection_map(norway_sensors_experiment, all_data = prepared_traffic_links_norway)

# Saving results
saveRDS(norway_sensors_experiment, "results/production/frc_0.8_logaadt_0.15_heavyratio_0.05_noncontinuous_100percounty_100overall.rds")
