# pak::pak("trafikkdata/xdtkit")
# install.packages("jsonlite")
# install.packages("sf")
# install.packages(
#   "INLA",
#   repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"),
#   dep = TRUE)
library(INLA)
# inla.binary.install(os = c("Ubuntu-22.04"))
library(xdtkit)
library(jsonlite)
library(sf)
source("R/scripts_inla_sensor_placement.R")

year <- 2024

norway_directed_traffic_links <- jsonlite::fromJSON("data/directed-traffic-links-2024.json")
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

ordinal_levels_road <- list(
  functionClass       = c("unknown", "E", "D", "C", "B", "A"),
  highestSpeedLimit   = c("unknown","20", "30", "40", "50", "60", "70", "80", "90", "100", "110")
)
#Best similarity covariates
similarity_covariates_norway <- c(
  "minLanes", "highestSpeedLimit", "functionClass",
  "lastYearAadt_logAadt"
)
inla_weighted_norway <- fit_weighted_inla_model(
  data                  = prepared_traffic_links_norway,
  adjacency_matrix      = adjacency_matrix_norway,
  spatial_term          = "besagproper_rbf",
  fixed_effects         = covariates,
  iid_effects           = "roadSystem",
  ordinal_levels        = ordinal_levels_road,
  similarity_covariates = similarity_covariates_norway,
  family                = "nbinomial")

predictions_total_norway <- dplyr::left_join(prepared_traffic_links_norway, inla_weighted_norway$predictions)
balanced_results <- balance_predictions(
  data = predictions_total_norway,
  nodes = nodes_norway,
  balancing_grouping_variable = clusters,
  nodes_to_balance = "complete nodes",
  year = year)
predictions_total_norway <- dplyr::full_join(predictions_total_norway, balanced_results$balanced_res)


library(leaflet)

map_sf <- add_geometries(predictions_total_norway) |>
  sf::st_transform(4326)
map_sf <- map_sf[!sf::st_is_empty(map_sf) & !is.na(map_sf$balanced_pred), ]

# Log-scale palette so high-AADT links don't drown out the rest
pal <- colorNumeric("viridis", domain = log1p(map_sf$balanced_pred))

leaflet(map_sf, options = leafletOptions(preferCanvas = TRUE)) |>
  addTiles() |>
  addPolylines(
    color   = ~pal(log1p(balanced_pred)),
    weight  = 2,
    opacity = 0.9
  ) |>
  addLegend(
    "bottomright",
    pal      = pal,
    values   = ~log1p(balanced_pred),
    title    = "Balanced predictions",
    labFormat = labelFormat(transform = function(x) round(expm1(x)))
  )

map_sf_measured <- add_geometries(predictions_total_norway) |>
  sf::st_transform(4326)
map_sf_measured <- map_sf_measured[
  !sf::st_is_empty(map_sf_measured) & !is.na(map_sf_measured$aadt), 
]

pal_aadt <- colorNumeric("viridis", domain = log1p(map_sf_measured$aadt))

leaflet(map_sf_measured, options = leafletOptions(preferCanvas = TRUE)) |>
  addTiles() |>
  addPolylines(
    color   = ~pal_aadt(log1p(aadt)),
    weight  = 2,
    opacity = 0.9
  ) |>
  addLegend(
    "bottomright",
    pal       = pal_aadt,
    values    = ~log1p(aadt),
    title     = "Measured AADT",
    labFormat = labelFormat(transform = function(x) round(expm1(x)))
  )

map_sf_unc <- add_geometries(predictions_total_norway) |>
  sf::st_transform(4326)

map_sf_unc <- map_sf_unc[
  !sf::st_is_empty(map_sf_unc) &
    !is.na(map_sf_unc$balanced_pred) &
    !is.na(map_sf_unc$balanced_sd) &
    map_sf_unc$balanced_pred > 0,
]

map_sf_unc$relative_uncertainty <- map_sf_unc$balanced_sd / map_sf_unc$balanced_pred
map_sf_unc <- map_sf_unc[map_sf_unc$relative_uncertainty < 0.5, ]

pal_unc <- colorNumeric(
  palette = c("green", "yellow", "red"),
  domain  = log1p(map_sf_unc$relative_uncertainty),
  reverse = FALSE
)
leaflet(map_sf_unc, options = leafletOptions(preferCanvas = TRUE)) |>
  addTiles() |>
  addPolylines(
    color   = ~pal_unc(log1p(relative_uncertainty)),
    weight  = 2,
    opacity = 0.9
  ) |>
  addLegend(
    "bottomright",
    pal       = pal_unc,
    values    = ~log1p(relative_uncertainty),
    title     = "Relative uncertainty",
    labFormat = labelFormat(transform = function(x) round(expm1(x), 3))
  )