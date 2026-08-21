#County partitioning 

#' Expand a seed set of link indices outward across an adjacency matrix
#'
#' Repeatedly grows `seed_idx` by `hops` adjacency steps and returns only the newly
#' reached border indices (the original seed set is excluded from the result). Used by
#' `partition_by_county()` to widen each county partition across its borders.
#'
#' @param seed_idx         Integer vector of starting row/column indices.
#' @param adjacency_matrix Sparse square adjacency matrix.
#' @param hops             Number of adjacency hops to expand (default 2).
#' @return Integer vector of border indices reached within `hops` steps, excluding `seed_idx`.
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
#' @param hops             Number of hops to expand beyond the county boundary (default 2).
#' @return Named list (one element per county level), each containing:
#'   \code{$data} (subset of \code{data}) and
#'   \code{$adjacency_matrix} (subsetted and reindexed).
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

# ════════════════════════════════════════════════════════════════════════════════
# categorise_aadt() ----
# ════════════════════════════════════════════════════════════════════════════════

#' Categorise raw AADT values into ordered traffic-volume bands.
#'
#' @param x       Numeric vector of AADT values (NAs passed through as NA level).
#' @param breaks  Numeric break points passed to \code{cut()}.
#'   Default: \code{c(0, 1000, 5000, 15000, 35000, Inf)}.
#' @param labels  Character labels for each interval (length = length(breaks)-1).
#'   Default: \code{c("Very Low","Low","Medium","High","Very High")}.
#' @return An ordered factor of the same length as \code{x}.
categorise_aadt <- function(
    x,
    breaks = c(0, 1000, 5000, 15000, 35000, Inf),
    labels = c("Very Low", "Low", "Medium", "High", "Very High")) {
  cut(x, breaks = breaks, labels = labels,
      include.lowest = TRUE, right = FALSE, ordered_result = TRUE)
}

# ════════════════════════════════════════════════════════════════════════════════
# summarise_sensor_selection() ----
# ════════════════════════════════════════════════════════════════════════════════

#' Compute summary statistics for a completed sensor selection.
#'
#' Produces a single tidy long data frame (\code{counts_table}) with one row
#' per covariate-level combination for all selected (non-derived, non-neighbour)
#' sensors, plus a handful of scalar fields.
#'
#' @param selected_data  The \code{selected_data_entries} data frame returned by
#'   a Norway-level sensor-selection wrapper.
#' @param params         Named list of call parameters to record verbatim.
#' @param aadt_breaks    Forwarded to \code{categorise_aadt()}.
#' @param aadt_labels    Forwarded to \code{categorise_aadt()}.
#' @return A list with elements:
#'   \describe{
#'     \item{counts_table}{Long data frame: \code{covariate}, \code{level},
#'       \code{n}, \code{pct}.}
#'     \item{n_selected}{Total number of selected sensor rows.}
#'     \item{parameters}{The \code{params} list, verbatim.}
#'     \item{enables_derivable_count}{Number of selected sensors with
#'       \code{enables_derivable == TRUE}, or \code{NULL} if column absent.}
#'     \item{total_derivable_links_enabled}{Sum of \code{n_derivable} across
#'       selected sensors, or \code{NULL} if column absent.}
#'   }
summarise_sensor_selection <- function(
    selected_data,
    params       = list(),
    aadt_breaks  = c(0, 1000, 5000, 15000, 35000, Inf),
    aadt_labels  = c("Very Low", "Low", "Medium", "High", "Very High")) {

  # Filter to actual placed sensors (not flow-derived rows or neighbour context).
  is_sel     <- isTRUE(selected_data$selected) | identical(selected_data$selected, TRUE)
  is_sel     <- selected_data$selected %in% TRUE
  not_deriv  <- if ("is_derived" %in% names(selected_data))
    !selected_data$is_derived %in% TRUE else rep(TRUE, nrow(selected_data))
  not_neighbour <- if ("neighbour_hop" %in% names(selected_data))
    is.na(selected_data$neighbour_hop) else rep(TRUE, nrow(selected_data))

  df         <- selected_data[is_sel & not_deriv & not_neighbour, ]
  # Deduplicate by physical sensor location: directed link IDs end with
  # -WITH or -AGAINST; strip that suffix so both directions count as one sensor.
  df         <- df[!duplicated(sub("-WITH$|-AGAINST$", "", df$id)), ]
  n_selected <- nrow(df)

  # Helper: build one block of the long table from a column or derived vector.
  make_block <- function(covariate_name, values) {
    tbl <- table(as.character(values), useNA = "ifany")
    block <- data.frame(
      covariate = covariate_name,
      level     = names(tbl),
      n         = as.integer(tbl),
      stringsAsFactors = FALSE
    )
    block$pct <- if (n_selected > 0) round(100 * block$n / n_selected, 1) else NA_real_
    block
  }

  blocks <- list()

  if ("functionalRoadClass" %in% names(df))
    blocks[["functionalRoadClass"]] <- make_block("functionalRoadClass", df$functionalRoadClass)

  if ("functionClass" %in% names(df))
    blocks[["functionClass"]] <- make_block("functionClass", df$functionClass)

  if ("roadCategory" %in% names(df))
    blocks[["roadCategory"]] <- make_block("roadCategory", df$roadCategory)

  if ("county" %in% names(df))
    blocks[["county"]] <- make_block("county", df$county)

  if ("lastYearAadt_aadt" %in% names(df))
    blocks[["aadt_category"]] <- make_block(
      "aadt_category",
      categorise_aadt(df$lastYearAadt_aadt, breaks = aadt_breaks, labels = aadt_labels)
    )

  counts_table <- if (length(blocks) > 0L)
    do.call(rbind, blocks) else
    data.frame(covariate = character(), level = character(),
               n = integer(), pct = numeric(), stringsAsFactors = FALSE)
  rownames(counts_table) <- NULL

  # Derivable scalars (only when enriched by identify_derivable_nodes()).
  enables_derivable_count <- if ("enables_derivable" %in% names(df))
    sum(df$enables_derivable %in% TRUE, na.rm = TRUE) else NULL

  total_derivable_links_enabled <- if ("n_derivable" %in% names(df))
    sum(df$n_derivable, na.rm = TRUE) else NULL

  list(
    counts_table                  = counts_table,
    n_selected                    = n_selected,
    parameters                    = params,
    enables_derivable_count       = enables_derivable_count,
    total_derivable_links_enabled = total_derivable_links_enabled
  )
}

# ════════════════════════════════════════════════════════════════════════════════
# add_geometries() ----
# ════════════════════════════════════════════════════════════════════════════════

# Cache environment so the GeoJSON is only read once per session.
.geom_cache <- new.env(parent = emptyenv())

#' Attach LineString geometries to a traffic-link data frame
#'
#' Reads \code{data/directed_traffic_links_2024.geojson} (once, then cached) and
#' left-joins its geometries onto \code{traffic_data} by the link ID column.
#' This overrides \code{xdtkit::add_geometries}, which uses a bundled geometry
#' dataset with an older ID scheme incompatible with the 2024 data.
#'
#' @param traffic_data  A data frame containing at least the column named by
#'   \code{id_col}.
#' @param id_col        Name of the ID column in \code{traffic_data}
#'   (default \code{"id"}).
#' @param geojson_path  Path to the GeoJSON file.  Defaults to the canonical
#'   2024 file relative to the project root.
#'
#' @return An \code{sf} data frame with LineString geometries attached.
add_geometries <- function(
    traffic_data,
    id_col      = "id",
    geojson_path = if (file.exists("data/directed_traffic_links_2024_simplified.geojson"))
                     "data/directed_traffic_links_2024_simplified.geojson"
                   else
                     "data/directed_traffic_links_2024.geojson") {

  if (!id_col %in% names(traffic_data))
    stop("Column '", id_col, "' not found in traffic_data.")

  # ── Read (and cache) geometry lookup ──────────────────────────────────────────
  cache_key <- normalizePath(geojson_path, mustWork = FALSE)
  if (!exists(cache_key, envir = .geom_cache, inherits = FALSE)) {
    if (!file.exists(geojson_path))
      stop("GeoJSON file not found: ", geojson_path,
           "\nSet geojson_path to the correct path.")
    geo <- sf::st_read(geojson_path, quiet = TRUE)
    geo <- sf::st_transform(geo[, "id"], 25833) |>   # project to UTM 33N (metres)
      sf::st_simplify(dTolerance = 1, preserveTopology = TRUE) |>
      sf::st_transform(4326)                         # back to WGS84
    assign(cache_key, geo, envir = .geom_cache)
  }
  geo <- get(cache_key, envir = .geom_cache, inherits = FALSE)

  # ── Join and return sf ────────────────────────────────────────────────────────
  sf::st_as_sf(
    dplyr::left_join(traffic_data, geo,
                     by = stats::setNames("id", id_col))
  )
}


