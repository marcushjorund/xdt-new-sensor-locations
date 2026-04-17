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
source("serious-testing/scripts_inla_sensor_placement.R")

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

#Defining the ordinal levels and their ordering
ordinal_levels_road <- list(
  functionClass       = c("unknown", "E", "D", "C", "B", "A"),
  highestSpeedLimit   = c("unknown","20", "30", "40", "50", "60", "70", "80", "90", "100", "110")
)
#Best similarity covariates
similarity_covariates_norway <- c(
  "minLanes", "highestSpeedLimit", "functionClass",
  "lastYearAadt_logAadt")

inla_rbf_norway_full <- fit_inla_rbf_model(
  data                  = prepared_traffic_links_norway,
  adjacency_matrix      = adjacency_matrix_norway,
  spatial_term          = "besagproper_rbf",
  fixed_effects         = covariates,
  iid_effects           = "roadSystem",
  ordinal_levels        = ordinal_levels_road,
  similarity_covariates = similarity_covariates_norway,
  family                = "nbinomial")


spatial_hyperpar <- inla_rbf_norway_full$spatial_hyperpar[2:4, "mean"]
tau   <- spatial_hyperpar[1]
d     <- spatial_hyperpar[2]
sigma <- spatial_hyperpar[3]

distances <- inla_rbf_norway_full$distances

# ── Select top-k sensors per county, top-r across all of Norway ──────────────
norway_sensors <- greedy_mi_sensor_selection_norway(
  data               = prepared_traffic_links_norway,
  adjacency_matrix   = adjacency_matrix_norway,
  distances          = distances,
  tau                = tau,
  d                  = d,
  sigma              = sigma,
  hops               = 3,
  k                  = 100,          # sensors selected per county
  r                  = 200,         # set e.g. r = 500 for Norway-wide top-r
  weighting_bias     = 1 / seq_len(8),
  include_neighbours = TRUE          # append adjacency context rows for plotting
)

selected_sensor_placements <- norway_sensors$selected_data_entries



# ── Plot sensor selection map ─────────────────────────────────────────────────
plot_sensor_selection_map(norway_sensors)
