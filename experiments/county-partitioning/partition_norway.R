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

# County partitioning 

expand_border <- function(seed_idx, adjacency_matrix, hops = 2) {
  current <- seed_idx
  for (i in seq_len(hops)) {
    adj_sub <- adjacency_matrix[current, , drop = FALSE]
    new_idx <- which(Matrix::colSums(adj_sub != 0) > 0)
    current <- union(current, new_idx)
  }
  setdiff(current, seed_idx)
}

#' Partition a nationwide traffic-link dataset into per-county sub-datasets,
#' each widened by \code{hops} adjacency-matrix steps across county borders.
#'
#' @param data             Data frame with a \code{county} factor column.
#' @param adjacency_matrix Sparse square adjacency matrix matching \code{nrow(data)}.
#' @param distances        Numeric vector of edge distances (upper-triangle of adjacency_matrix),
#'   as returned by \code{fit_weighted_inla_model()$distances}.
#' @param hops             Number of hops to expand beyond the county boundary (default 2).
#' @return Named list (one element per county level), each containing:
#'   \code{$data}, \code{$adjacency_matrix}, and \code{$distances} (all subsetted consistently).
#'   Core rows are identified by \code{$data$county == county_name}.
partition_by_county <- function(data, adjacency_matrix, distances, hops = 2) {
  # Pre-compute upper-triangle edge indices from the full adjacency matrix once.
  # distances[k] corresponds to edge (ui[k], uj[k]); same logic as
  # create_covariance_and_precision_matrix().
  adj_trip   <- Matrix::summary(as(adjacency_matrix, "dgCMatrix"))
  mask_upper <- adj_trip$i < adj_trip$j
  ui         <- adj_trip$i[mask_upper]
  uj         <- adj_trip$j[mask_upper]

  counties   <- levels(data$county)
  partitions <- lapply(counties, function(county) {
    county_idx <- which(data$county == county)
    border_idx <- expand_border(county_idx, adjacency_matrix, hops = hops)
    subset_idx <- sort(union(county_idx, border_idx))
    # Keep only edges whose both endpoints are in subset_idx.
    # subset_idx is sorted, so original_i < original_j  ↔  i_sub < j_sub,
    # preserving the column-major order that Matrix::summary uses.
    edge_mask  <- ui %in% subset_idx & uj %in% subset_idx
    list(
      data             = data[subset_idx, ],
      adjacency_matrix = adjacency_matrix[subset_idx, subset_idx],
      distances        = distances[edge_mask]
    )
  })
  setNames(partitions, counties)
}

county_partitions <- partition_by_county(
  data             = prepared_traffic_links_norway,
  adjacency_matrix = adjacency_matrix_norway,
  distances        = inla_weighted_norway_full$distances,
  hops             = 2
)

county_partitions
