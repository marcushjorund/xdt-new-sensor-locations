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
incidence_matrix_norway <- build_incidence_matrix(nodes_norway, prepared_traffic_links_norway, nodes_to_balance = "complete nodes", sparse = TRUE)
derivable_nodes_data <- identify_derivable_nodes(incidence_matrix_norway, prepared_traffic_links_norway, nodes = nodes_norway)


# Links that are currently derivable from existing measurements (no sensor needed)
map_traffic_links(derivable_nodes_data, color_by = "derivable",
                  title = "Currently derivable")

# Links where placing a sensor unlocks at least one other link via exact flow conservation
map_traffic_links(derivable_nodes_data, color_by = "enables_derivable",
                  title = "Enables derivable if measured")

# How many links would be unlocked per sensor placement
map_traffic_links(derivable_nodes_data, color_by = "n_enables_derivable",
                  title = "N links unlocked by sensor")

# All traffic links in black — baseline to check if any links are missing from the plots above
all_links_sf <- add_geometries(prepared_traffic_links_norway) |>
  sf::st_transform(crs = 4326)
all_links_sf <- all_links_sf[!sf::st_is_empty(all_links_sf) & sf::st_is_valid(all_links_sf), ]

leaflet::leaflet(all_links_sf, options = leaflet::leafletOptions(preferCanvas = TRUE)) |>
  leaflet::addTiles() |>
  leaflet::addPolylines(
    color   = "#000000",
    weight  = 1,
    opacity = 0.6,
    popup   = ~id
  )

