library(INLA)
library(xdtkit)
library(igraph)
library(Matrix)
library(dplyr)

# =============================================================================
# Data preparation
# =============================================================================

year <- 2025

preprocessed_traffic_links <- preprocess_traffic_links(buskerud_directed_traffic_links, year = year)
bus_aadt <- calculate_bus_aadt(stops_on_traffic_links, bus_counts, year = year)

prepared_traffic_links <- fill_missing_values(
  df = preprocessed_traffic_links,
  unknown_impute_columns = c("functionClass", "highestSpeedLimit", "lowestSpeedLimit", "maxLanes", "minLanes"),
  mode_impute_columns = c("hasOnlyPublicTransportLanes"),
  median_impute_columns = c("lastYearAadt_aadt", "lastYearAadt_heavyRatio", "lastYearAadt_heavyAadt")
) |>
  remove_negative_aadt() |>
  add_logLastYear() |>
  join_bus_to_traffic(bus_aadt)

prepared_traffic_links_emma <- readRDS("prepared_traffic_links_emma.rds")
prepared_traffic_links <- subset(prepared_traffic_links_emma, select = -c(spatial.idx, inla_pred))

adjacency_matrix <- build_adjacency_matrix(
  prepared_traffic_links,
  exclude_public_transport = TRUE)

clusters <- strategic_network_clustering(
  data = prepared_traffic_links,
  year = year,
  boundary_links = c("Trafikkdata_continuous")
)

nodes <- identify_unbalanceable_nodes(buskerud_nodes, prepared_traffic_links)

# =============================================================================
# INLA model
# =============================================================================

covariates <- ~ maxLanes + hasOnlyPublicTransportLanes + roadCategory +
  minLanes:roadCategory + lastYearAadt_logAadt

inla_model_total <- fit_inla_model(
  data = prepared_traffic_links,
  adjacency_matrix,
  fixed_effects = covariates,
  iid_effects = "roadSystem",
  family = "poisson")

predictions_total <- dplyr::full_join(prepared_traffic_links, inla_model_total$predictions)

balanced_model_total <- balance_predictions(
  data = predictions_total,
  nodes = nodes,
  balancing_grouping_variable = clusters,
  nodes_to_balance = "complete nodes",
  year = year)

predictions_total <- dplyr::full_join(predictions_total, balanced_model_total$balanced_res)
predictions_total$rel_cv <- predictions_total$balanced_sd / predictions_total$balanced_pred
predictions_total$logrel_cv <- log(predictions_total$rel_cv)

# =============================================================================
# Precision/covariance matrices from Besag proper hyperparameters
# =============================================================================

build_precision_matrix <- function(adjacency_matrix, tau, d) {
  n <- dim(adjacency_matrix)[1]
  precision_matrix <- matrix(0, ncol = n, nrow = n)
  for (i in seq_len(n)) {
    diag(precision_matrix)[i] <- tau * (d + sum(adjacency_matrix[i, ]))
    non_zero_idx <- which(adjacency_matrix[i, ] != 0)
    precision_matrix[i, non_zero_idx] <- -tau
  }
  precision_matrix
}

d   <- inla_model_total$model_summary$hyperpar$mean[2]
tau <- inla_model_total$model_summary$hyperpar$mean[1]

precision_matrix   <- build_precision_matrix(adjacency_matrix, tau, d)
covariance_matrix  <- solve(precision_matrix)
correlation_matrix <- cov2cor(covariance_matrix)

# =============================================================================
# Hop-distance matrix
# =============================================================================

build_hop_distance_matrix <- function(adjacency_matrix) {
  g <- igraph::graph_from_adjacency_matrix(
    adjacency_matrix,
    mode     = "undirected",
    weighted = NULL,
    diag     = FALSE)
  distance_matrix <- igraph::distances(g, algorithm = "unweighted")
  rownames(distance_matrix) <- rownames(adjacency_matrix)
  colnames(distance_matrix) <- colnames(adjacency_matrix)
  distance_matrix
}

# Summarise covariance (or correlation) by hop distance.
# Returns a data frame with columns: hop, mean, median, sd, n.
# Set use_correlation = TRUE to summarise the correlation matrix instead.
covariance_by_hops <- function(distance_matrix, covariance_matrix, use_correlation = FALSE) {
  mat <- if (use_correlation) cov2cor(covariance_matrix) else covariance_matrix
  hops <- distance_matrix[upper.tri(distance_matrix) & is.finite(distance_matrix)]
  vals <- mat[upper.tri(distance_matrix) & is.finite(distance_matrix)]
  hop_levels <- sort(unique(hops[hops > 0]))
  result <- do.call(rbind, lapply(hop_levels, function(h) {
    v <- vals[hops == h]
    data.frame(hop = h, mean = mean(v), median = median(v), sd = sd(v), n = length(v))
  }))
  result
}

distance_matrix <- build_hop_distance_matrix(adjacency_matrix)

# =============================================================================
# Sensor selection: naive MI algorithm (individual links)
# =============================================================================

greedy_sensor_selection <- function(covariance_matrix, measured_idxs, unmeasured_idxs, k) {
  new_sensor_location_idxs <- rep(0, k)
  for (i in seq_len(k)) {
    delta <- rep(0, length(unmeasured_idxs))
    for (j in seq_along(unmeasured_idxs)) {
      y <- unmeasured_idxs[j]
      temp_unmeasured_idxs <- unmeasured_idxs[unmeasured_idxs != y]
      delta[j] <- (
        covariance_matrix[y, y] -
          t(covariance_matrix[y, measured_idxs]) %*%
          solve(covariance_matrix[measured_idxs, measured_idxs]) %*%
          covariance_matrix[measured_idxs, y]
      ) / (
        covariance_matrix[y, y] -
          t(covariance_matrix[y, temp_unmeasured_idxs]) %*%
          solve(covariance_matrix[temp_unmeasured_idxs, temp_unmeasured_idxs]) %*%
          covariance_matrix[temp_unmeasured_idxs, y]
      )
    }
    best <- which(delta == max(delta))[1]
    new_sensor_location_idxs[i] <- unmeasured_idxs[best]
    measured_idxs   <- c(measured_idxs, unmeasured_idxs[best])
    unmeasured_idxs <- unmeasured_idxs[unmeasured_idxs != unmeasured_idxs[best]]
  }
  new_sensor_location_idxs
}

unmeasured_idxs <- which(is.na(predictions_total$aadt))
measured_idxs   <- which(!is.na(predictions_total$aadt))

new_sensor_location_idxs <- greedy_sensor_selection(covariance_matrix, measured_idxs, unmeasured_idxs, k = 10)
new_sensor_ids <- predictions_total$id[new_sensor_location_idxs]

# =============================================================================
# WITH + AGAINST linear combination
# =============================================================================

build_link_combination_matrix <- function(prepared_traffic_links) {
  with_idx_int    <- which(grepl("-WITH$",    prepared_traffic_links$id))
  against_idx_int <- which(grepl("-AGAINST$", prepared_traffic_links$id))

  base_with    <- sub("-WITH$",    "", prepared_traffic_links$id[with_idx_int])
  base_against <- sub("-AGAINST$", "", prepared_traffic_links$id[against_idx_int])

  paired_bases   <- intersect(base_with, base_against)
  unpaired_bases <- setdiff(base_with, base_against)

  paired_with_idx    <- match(paste0(paired_bases,   "-WITH"),    prepared_traffic_links$id)
  paired_against_idx <- match(paste0(paired_bases,   "-AGAINST"), prepared_traffic_links$id)
  unpaired_with_idx  <- match(paste0(unpaired_bases, "-WITH"),    prepared_traffic_links$id)

  n <- nrow(prepared_traffic_links)
  m <- length(paired_bases) + length(unpaired_bases)

  row_idx <- c(
    seq_along(paired_bases),
    seq_along(paired_bases),
    length(paired_bases) + seq_along(unpaired_bases)
  )
  col_idx <- c(paired_with_idx, paired_against_idx, unpaired_with_idx)
  values  <- rep(1, length(row_idx))

  C <- sparseMatrix(i = row_idx, j = col_idx, x = values, dims = c(m, n))

  list(C = C, paired_bases = paired_bases, unpaired_bases = unpaired_bases,
       paired_with_idx = paired_with_idx, unpaired_with_idx = unpaired_with_idx)
}

link_combo           <- build_link_combination_matrix(prepared_traffic_links)
C                    <- link_combo$C
paired_bases         <- link_combo$paired_bases
unpaired_bases       <- link_combo$unpaired_bases
covariance_matrix_sum <- C %*% covariance_matrix %*% t(C)

# Measured/unmeasured indices in the combined covariance matrix
measured_paired_idxs     <- which(!is.na(prepared_traffic_links[link_combo$paired_with_idx, ]$aadt))
unmeasured_paired_idxs   <- which( is.na(prepared_traffic_links[link_combo$paired_with_idx, ]$aadt))
measured_unpaired_idxs   <- which(!is.na(prepared_traffic_links[link_combo$unpaired_with_idx, ]$aadt))
unmeasured_unpaired_idxs <- which( is.na(prepared_traffic_links[link_combo$unpaired_with_idx, ]$aadt))

measured_idxs_sum   <- c(measured_paired_idxs,   measured_unpaired_idxs   + length(paired_bases))
unmeasured_idxs_sum <- c(unmeasured_paired_idxs, unmeasured_unpaired_idxs + length(paired_bases))

# =============================================================================
# Sensor selection on combined WITH+AGAINST covariance
# =============================================================================

new_sensor_location_idxs_sum <- greedy_sensor_selection(
  covariance_matrix_sum, measured_idxs_sum, unmeasured_idxs_sum, k = 30)

paired_new_sensor_locations   <- new_sensor_location_idxs_sum[new_sensor_location_idxs_sum <= length(paired_bases)]
unpaired_new_sensor_locations <- new_sensor_location_idxs_sum[new_sensor_location_idxs_sum >  length(paired_bases)]

paired_new_sensor_ids   <- paired_bases[paired_new_sensor_locations]
unpaired_new_sensor_ids <- unpaired_bases[unpaired_new_sensor_locations - length(paired_bases)]
