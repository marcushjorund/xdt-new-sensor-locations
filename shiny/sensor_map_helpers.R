# =============================================================================
# sensor_map_helpers.R
#
# Self-contained mapping helpers for the sensor placement Shiny app.
# Extracted from scripts_inla_sensor_placement.R (add_geometries,
# plot_sensor_selection_map, map_traffic_links) and xdtkit (nvdb_objects).
#
# No dependency on INLA, xdtkit, or any other heavy package.
# Required packages: sf, dplyr, leaflet
#
# The GeoJSON path is read from options("sensor_app.geojson_path").
# Set it at app startup:
#   options(sensor_app.geojson_path = "/absolute/path/to/file.geojson")
# =============================================================================

library(sf)
library(dplyr)
library(leaflet)


# ════════════════════════════════════════════════════════════════════════════════
# nvdb_objects() ----
# Extracted from xdtkit — returns the NVDB tile URL, attribution, and CRS.
# ════════════════════════════════════════════════════════════════════════════════

nvdb_objects <- function() {
  nvdb_map_url <- paste0(
    "https://nvdbcache.geodataonline.no/arcgis/rest/services/",
    "Trafikkportalen/GeocacheTrafikkJPG/MapServer/tile/{z}/{y}/{x}"
  )
  nvdb_map_attribution <- paste0(
    "NVDB, Geovekst, kommunene og Open Street Map contributors (utenfor Norge)"
  )
  nvdb_crs <- leaflet::leafletCRS(
    crsClass   = "L.Proj.CRS",
    code       = "EPSG:25833",
    proj4def   = "+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs",
    resolutions = c(
      21674.7100160867, 10837.3550080434,  5418.67750402168,
       2709.33875201084,  1354.66937600542,   677.334688002709,
        338.667344001355,   169.333672000677,    84.6668360003387,
         42.3334180001693,    21.1667090000847,    10.5833545000423,
          5.29167725002117,     2.64583862501058,     1.32291931250529,
          0.661459656252646,    0.330729828126323,    0.165364914063161
    ),
    origin = c(-2500000, 9045984)
  )
  list(
    nvdb_url         = nvdb_map_url,
    nvdb_attribution = nvdb_map_attribution,
    nvdb_crs         = nvdb_crs
  )
}


# ════════════════════════════════════════════════════════════════════════════════
# .geom_cache ----
# Persistent geometry cache — the GeoJSON is read once and reused across calls.
# ════════════════════════════════════════════════════════════════════════════════

.geom_cache <- new.env(parent = emptyenv())


# ════════════════════════════════════════════════════════════════════════════════
# add_geometries() ----
# ════════════════════════════════════════════════════════════════════════════════

#' Attach LineString geometries to a traffic-link data frame
#'
#' Reads the GeoJSON specified by \code{options("sensor_app.geojson_path")}
#' (once per session, then cached) and left-joins its geometries onto
#' \code{traffic_data} by the link ID column.
#'
#' @param traffic_data  A data frame containing at least the column named by
#'   \code{id_col}.
#' @param id_col        Name of the ID column in \code{traffic_data}
#'   (default \code{"id"}).
#' @param geojson_path  Path to the GeoJSON file.  Defaults to
#'   \code{getOption("sensor_app.geojson_path")} which must be set at app
#'   startup to an absolute path.
#'
#' @return An \code{sf} data frame with LineString geometries attached.
add_geometries <- function(
    traffic_data,
    id_col       = "id",
    geojson_path = getOption("sensor_app.geojson_path")) {

  if (is.null(geojson_path))
    stop("geojson_path is NULL. ",
         "Set options(sensor_app.geojson_path = '/path/to/file.geojson') ",
         "at app startup.")

  if (!id_col %in% names(traffic_data))
    stop("Column '", id_col, "' not found in traffic_data.")

  # ── Read (and cache) geometry lookup ──────────────────────────────────────────
  cache_key <- normalizePath(geojson_path, mustWork = FALSE)
  if (!exists(cache_key, envir = .geom_cache, inherits = FALSE)) {
    if (!file.exists(geojson_path))
      stop("GeoJSON file not found: ", geojson_path,
           "\nSet options(sensor_app.geojson_path = '/correct/path').")
    geo <- sf::st_read(geojson_path, quiet = TRUE)
    geo <- sf::st_transform(geo[, "id"], 25833) |>
      sf::st_simplify(dTolerance = 1, preserveTopology = TRUE) |>
      sf::st_transform(4326)
    assign(cache_key, geo, envir = .geom_cache)
  }
  geo <- get(cache_key, envir = .geom_cache, inherits = FALSE)

  # ── Join and return sf ────────────────────────────────────────────────────────
  sf::st_as_sf(
    dplyr::left_join(traffic_data, geo,
                     by = stats::setNames("id", id_col))
  )
}


# ════════════════════════════════════════════════════════════════════════════════
# plot_sensor_selection_map() ----
# ════════════════════════════════════════════════════════════════════════════════

plot_sensor_selection_map <- function(
    sensor_results,
    color_selected    = "#E63946",
    color_neighbour   = "#457B9D",
    color_derived     = "#F4A261",
    weight_selected   = 6,
    weight_neighbour  = 3,
    opacity_selected  = 0.95,
    opacity_neighbour = 0.7,
    all_data          = NULL,
    color_measured    = "#FF00FF") {

  # ── 1. Resolve input ──────────────────────────────────────────────────────────
  if (is.list(sensor_results) && !is.data.frame(sensor_results)) {
    if (!"selected_data_entries" %in% names(sensor_results))
      stop("sensor_results must contain '$selected_data_entries'.")
    sel_df <- sensor_results$selected_data_entries
  } else if (is.data.frame(sensor_results)) {
    sel_df <- sensor_results
  } else {
    stop("sensor_results must be a list from greedy_mi_sensor_selection_norway() or a data frame.")
  }
  if (!"selected" %in% names(sel_df))
    stop("sel_df must contain a 'selected' column (TRUE = selected, FALSE = neighbour).")

  sel_df <- sel_df[!duplicated(sel_df$id), , drop = FALSE]

  # ── 2. Geometry ───────────────────────────────────────────────────────────────
  sel_sf <- add_geometries(sel_df)
  sel_sf <- sf::st_transform(sel_sf, crs = 4326)

  valid  <- !sf::st_is_empty(sel_sf) & sf::st_is_valid(sel_sf)
  sel_sf <- sel_sf[valid, ]

  # ── 3. Rank + popup ───────────────────────────────────────────────────────────
  sel_sf$base_id <- sub("-WITH$|-AGAINST$", "", sel_sf$id)
  sel_sf$rank    <- NA_integer_
  true_idx       <- which(sel_sf$selected %in% TRUE)

  has_derived_col <- "is_derived"      %in% names(sel_sf)
  has_sel_order   <- "selection_order" %in% names(sel_sf)
  is_derived_vec  <- if (has_derived_col) sel_sf$is_derived %in% TRUE else rep(FALSE, nrow(sel_sf))
  use_bundle_rank <- has_derived_col && has_sel_order &&
                     length(true_idx) > 0 && any(is_derived_vec[true_idx])

  if (length(true_idx) > 0) {
    true_sf <- sel_sf[true_idx, ]
    if (use_bundle_rank) {
      bundle_scores <- tapply(true_sf$mi_score, true_sf$selection_order,
                              max, na.rm = TRUE)
      bundle_order  <- as.integer(names(sort(bundle_scores, decreasing = TRUE)))
      rank_lookup   <- stats::setNames(seq_along(bundle_order),
                                       as.character(bundle_order))
      sel_sf$rank[true_idx] <- rank_lookup[as.character(sel_sf$selection_order[true_idx])]
    } else if ("mi_score" %in% names(true_sf)) {
      base_order  <- unique(true_sf$base_id[order(-true_sf$mi_score, na.last = TRUE)])
      rank_lookup <- stats::setNames(seq_along(base_order), base_order)
      sel_sf$rank[true_idx] <- rank_lookup[sel_sf$base_id[true_idx]]
    } else {
      base_order  <- unique(true_sf$base_id)
      rank_lookup <- stats::setNames(seq_along(base_order), base_order)
      sel_sf$rank[true_idx] <- rank_lookup[sel_sf$base_id[true_idx]]
    }
  }

  has_mi  <- "mi_score"                %in% names(sel_sf)
  has_ly  <- "lastYearAadt_aadt"       %in% names(sel_sf)
  has_co  <- "county"                  %in% names(sel_sf)
  has_src <- "traffic_volume_source"   %in% names(sel_sf)
  has_frc <- "functionalRoadClass"     %in% names(sel_sf)
  has_fc  <- "functionClass"           %in% names(sel_sf)
  has_hr  <- "lastYearAadt_heavyRatio" %in% names(sel_sf)

  has_hop_col <- "neighbour_hop" %in% names(sel_sf)
  hop_labels  <- ifelse(sel_sf$neighbour_hop == 1L, "Nabo",
                 ifelse(sel_sf$neighbour_hop == 2L, "Nabo til nabo",
                        paste0("Nabo (hop ", sel_sf$neighbour_hop, ")")))
  type_label  <- ifelse(
    sel_sf$selected %in% TRUE & !is_derived_vec, "Valgt trafikkregistreringspunkt",
    ifelse(sel_sf$selected %in% TRUE & is_derived_vec, "Utledet (bevaring av trafikkmengde)",
    ifelse(has_hop_col & !is.na(sel_sf$neighbour_hop),
           hop_labels,
           "Nabo")))
  pop <- paste0("<strong>ID:</strong> ", sel_sf$id,
                "<br><strong>Type:</strong> ", type_label)
  if (any(!is.na(sel_sf$rank)))
    pop <- paste0(pop, "<br><strong>Sensor rank:</strong> ",
                  ifelse(is.na(sel_sf$rank), "\u2014", as.character(sel_sf$rank)))
  if (has_mi)
    pop <- paste0(pop, "<br><strong>MI score:</strong> ",
                  ifelse(is.na(sel_sf$mi_score), "\u2014",
                         sprintf("%.4f", sel_sf$mi_score)))
  if (has_co)
    pop <- paste0(pop, "<br><strong>County:</strong> ",
                  ifelse(is.na(sel_sf$county), "\u2014", as.character(sel_sf$county)))
  if (has_src)
    pop <- paste0(pop, "<br><strong>Source:</strong> ",
                  ifelse(is.na(sel_sf$traffic_volume_source), "Unmeasured",
                         sel_sf$traffic_volume_source))
  if (has_ly)
    pop <- paste0(pop, "<br><strong>Last year AADT:</strong> ",
                  format(round(sel_sf$lastYearAadt_aadt),
                         big.mark = "\u00a0", scientific = FALSE))
  if (has_frc)
    pop <- paste0(pop, "<br><strong>Functional road class:</strong> ",
                  ifelse(is.na(sel_sf$functionalRoadClass), "\u2014",
                         as.character(sel_sf$functionalRoadClass)))
  if (has_fc)
    pop <- paste0(pop, "<br><strong>Function class:</strong> ",
                  ifelse(is.na(sel_sf$functionClass), "\u2014",
                         as.character(sel_sf$functionClass)))
  if (has_hr)
    pop <- paste0(pop, "<br><strong>Heavy ratio:</strong> ",
                  ifelse(is.na(sel_sf$lastYearAadt_heavyRatio), "\u2014",
                         sprintf("%.2f", sel_sf$lastYearAadt_heavyRatio)))
  sel_sf$popup_text <- pop

  # ── 4. Split ──────────────────────────────────────────────────────────────────
  sel_sensor     <- sel_sf[sel_sf$selected %in% TRUE & !is_derived_vec,  ]
  sel_derived_sf <- sel_sf[sel_sf$selected %in% TRUE &  is_derived_vec,  ]
  sel_neighbour  <- sel_sf[sel_sf$selected %in% FALSE, ]

  # ── 5. Midpoint markers ───────────────────────────────────────────────────────
  markers_sf <- NULL
  if (nrow(sel_sensor) > 0) {
    one_per_sensor <- sel_sensor[!duplicated(sel_sensor$base_id), ]
    suppressWarnings({
      mid <- sf::st_line_sample(sf::st_transform(one_per_sensor, 25833), sample = 0.5)
      mid <- sf::st_cast(mid, "POINT")
      mid <- sf::st_transform(mid, 4326)
    })
    markers_sf <- sf::st_set_geometry(one_per_sensor, mid)
    coords     <- sf::st_coordinates(markers_sf)
    markers_sf <- markers_sf[!is.na(coords[, 1]) & !is.na(coords[, 2]), ]
    if (nrow(markers_sf) == 0) markers_sf <- NULL
  }

  derived_markers_sf <- NULL
  if (nrow(sel_derived_sf) > 0) {
    one_per_derived <- sel_derived_sf[!duplicated(sel_derived_sf$base_id), ]
    suppressWarnings({
      mid_d <- sf::st_line_sample(sf::st_transform(one_per_derived, 25833), sample = 0.5)
      mid_d <- sf::st_cast(mid_d, "POINT")
      mid_d <- sf::st_transform(mid_d, 4326)
    })
    derived_markers_sf <- sf::st_set_geometry(one_per_derived, mid_d)
    coords_d           <- sf::st_coordinates(derived_markers_sf)
    derived_markers_sf <- derived_markers_sf[!is.na(coords_d[, 1]) & !is.na(coords_d[, 2]), ]
    if (nrow(derived_markers_sf) == 0) derived_markers_sf <- NULL
  }

  # ── 6. Build leaflet map ──────────────────────────────────────────────────────
  nvdb <- nvdb_objects()

  m <- leaflet::leaflet(
    sel_sf,
    options = leaflet::leafletOptions(crs = nvdb$nvdb_crs, zoomControl = TRUE)
  ) |>
    leaflet::addTiles(urlTemplate = nvdb$nvdb_url,
                      attribution = nvdb$nvdb_attribution)

  has_hop_data <- "neighbour_hop" %in% names(sel_neighbour)
  hops_in_data <- if (has_hop_data && nrow(sel_neighbour) > 0)
    sort(unique(stats::na.omit(as.integer(sel_neighbour$neighbour_hop))))
  else if (nrow(sel_neighbour) > 0) 1L else integer(0)
  hops_to_show <- hops_in_data
  neigh_groups <- ifelse(hops_to_show == 1L, "Naboer",
                  ifelse(hops_to_show == 2L, "Naboer til naboer",
                         paste0("Naboer (hop ", hops_to_show, ")")))

  for (h in hops_to_show) {
    hop_sf <- if (has_hop_data)
      sel_neighbour[!is.na(sel_neighbour$neighbour_hop) & sel_neighbour$neighbour_hop == h, ]
    else sel_neighbour
    if (nrow(hop_sf) == 0L) next
    hop_op  <- opacity_neighbour
    hop_wt  <- weight_neighbour
    hop_grp <- if (h == 1L) "Naboer" else if (h == 2L) "Naboer til naboer" else paste0("Naboer (hop ", h, ")")
    m <- m |> leaflet::addPolylines(
      data             = hop_sf,
      color            = color_neighbour,
      weight           = hop_wt,
      opacity          = hop_op,
      popup            = hop_sf$popup_text,
      group            = hop_grp,
      highlightOptions = leaflet::highlightOptions(
        weight = hop_wt + 2, color = "white", bringToFront = FALSE))
  }

  if (nrow(sel_derived_sf) > 0)
    m <- m |> leaflet::addPolylines(
      data             = sel_derived_sf,
      color            = color_derived,
      weight           = weight_selected,
      opacity          = opacity_selected,
      popup            = sel_derived_sf$popup_text,
      group            = "Utledede trafikklenker",
      highlightOptions = leaflet::highlightOptions(
        weight = weight_selected + 2, color = "white", bringToFront = TRUE))

  if (nrow(sel_sensor) > 0)
    m <- m |> leaflet::addPolylines(
      data             = sel_sensor,
      color            = color_selected,
      weight           = weight_selected,
      opacity          = opacity_selected,
      popup            = sel_sensor$popup_text,
      group            = "Valgte trafikkregistreringspunkter",
      highlightOptions = leaflet::highlightOptions(
        weight = weight_selected + 2, color = "white", bringToFront = TRUE))

  if (!is.null(markers_sf))
    m <- m |> leaflet::addCircleMarkers(
      data         = markers_sf,
      radius       = 8,
      color        = "#000000",
      weight       = 1.5,
      fillColor    = color_selected,
      fillOpacity  = 0.95,
      popup        = markers_sf$popup_text,
      label        = as.character(markers_sf$rank),
      labelOptions = leaflet::labelOptions(
        permanent = TRUE,
        textOnly  = TRUE,
        style     = list("font-weight" = "bold",
                         "font-size"   = "10px",
                         "color"       = "white",
                         "text-shadow" = "0 0 3px #000")))

  if (!is.null(derived_markers_sf))
    m <- m |> leaflet::addCircleMarkers(
      data         = derived_markers_sf,
      radius       = 8,
      color        = "#000000",
      weight       = 1.5,
      fillColor    = color_derived,
      fillOpacity  = 0.95,
      popup        = derived_markers_sf$popup_text,
      label        = as.character(derived_markers_sf$rank),
      labelOptions = leaflet::labelOptions(
        permanent = TRUE,
        textOnly  = TRUE,
        style     = list("font-weight" = "bold",
                         "font-size"   = "10px",
                         "color"       = "white",
                         "text-shadow" = "0 0 3px #000")))

  # ── 7. Measured traffic links (dashed overlay) ──────────────────────────────
  measured_sf <- NULL
  if (!is.null(all_data)) {
    measured_df <- all_data[!is.na(all_data$aadt), , drop = FALSE]
    measured_df <- measured_df[!duplicated(measured_df$id), ]
    if (nrow(measured_df) > 0) {
      measured_sf <- add_geometries(measured_df)
      measured_sf <- sf::st_transform(measured_sf, 4326)
      measured_sf <- measured_sf[
        !sf::st_is_empty(measured_sf) & sf::st_is_valid(measured_sf), ]
      if (nrow(measured_sf) == 0) measured_sf <- NULL
    }
  }

  legend_colors <- c(color_selected, color_neighbour)
  legend_labels <- c("Valgt trafikkregistreringspunkt", "Nabolenker")
  ctrl_groups   <- c(rev(neigh_groups), "Valgte trafikkregistreringspunkter")
  if (nrow(sel_derived_sf) > 0) {
    legend_colors <- c(color_selected, color_derived, color_neighbour)
    legend_labels <- c("Valgt trafikkregistreringspunkt",
                       "Utledet (bevaring av trafikkmengde)", "Nabolenker")
    ctrl_groups   <- c(rev(neigh_groups), "Utledede trafikklenker",
                       "Valgte trafikkregistreringspunkter")
  }
  if (!is.null(measured_sf)) {
    legend_colors <- c(legend_colors, color_measured)
    legend_labels <- c(legend_labels, "M\u00e5lt trafikklenke")
    ctrl_groups   <- c(ctrl_groups, "M\u00e5lte trafikklenker")
  }
  m <- m |>
    leaflet::addLayersControl(
      overlayGroups = ctrl_groups,
      options       = leaflet::layersControlOptions(collapsed = FALSE)) |>
    leaflet::addLegend(
      position = "bottomright",
      colors   = legend_colors,
      labels   = legend_labels,
      title    = "Nye trafikkregistreringspunkter",
      opacity  = 0.9)

  if (!is.null(measured_sf)) {
    has_src_m <- "traffic_volume_source" %in% names(measured_sf)
    meas_pop  <- paste0("<strong>ID:</strong> ", measured_sf$id,
                        "<br><strong>AADT:</strong> ",
                        format(round(measured_sf$aadt), big.mark = "\u00a0",
                               scientific = FALSE))
    if (has_src_m)
      meas_pop <- paste0(meas_pop,
                         "<br><strong>Kilde:</strong> ",
                         ifelse(is.na(measured_sf$traffic_volume_source), "\u2014",
                                measured_sf$traffic_volume_source))
    m <- m |> leaflet::addPolylines(
      data             = measured_sf,
      color            = color_measured,
      weight           = 2,
      opacity          = 0.75,
      dashArray        = "5, 10",
      popup            = meas_pop,
      group            = "M\u00e5lte trafikklenker",
      highlightOptions = leaflet::highlightOptions(
        weight = 4, color = "white", bringToFront = TRUE))
  }

  m
}


# ════════════════════════════════════════════════════════════════════════════════
# map_traffic_links() ----
# ════════════════════════════════════════════════════════════════════════════════

map_traffic_links <- function(
    data,
    color_by,
    title   = color_by,
    weight  = 2,
    opacity = 0.9) {

  if (!color_by %in% names(data))
    stop("Column '", color_by, "' not found in data.")

  map_sf <- add_geometries(data) |> sf::st_transform(crs = 4326)
  map_sf <- map_sf[!sf::st_is_empty(map_sf) & sf::st_is_valid(map_sf), ]
  if (nrow(map_sf) == 0L)
    stop("No valid geometries found after add_geometries().")

  col_values <- map_sf[[color_by]]
  if (is.logical(col_values) || is.character(col_values) || is.factor(col_values)) {
    map_sf[[color_by]] <- as.character(col_values)
    pal <- leaflet::colorFactor(palette = "Set1", domain = map_sf[[color_by]],
                                na.color = "#888888")
  } else {
    pal <- leaflet::colorNumeric(palette = "viridis", domain = col_values,
                                 na.color = "#888888")
  }

  has <- function(col) col %in% names(map_sf)
  pred_col <- if (has("balanced_pred")) "balanced_pred" else if (has("inla_pred")) "inla_pred" else NULL
  sd_col   <- if (has("balanced_sd"))   "balanced_sd"   else if (has("inla_sd"))   "inla_sd"   else NULL

  if (!is.null(pred_col) && !is.null(sd_col)) {
    map_sf$.rel_unc <- ifelse(
      !is.na(map_sf[[pred_col]]) & map_sf[[pred_col]] > 0,
      map_sf[[sd_col]] / map_sf[[pred_col]], NA_real_)
  }

  fmt <- function(x) format(round(x), big.mark = "\u00a0", scientific = FALSE)
  pop <- paste0("<strong>Link ID:</strong> ", map_sf$id)
  if (!is.null(pred_col))
    pop <- paste0(pop, "<br><strong>Predicted AADT:</strong> ",
                  ifelse(is.na(map_sf[[pred_col]]), "\u2014", fmt(map_sf[[pred_col]])))
  if (!is.null(sd_col))
    pop <- paste0(pop, "<br><strong>Standard deviation:</strong> ",
                  ifelse(is.na(map_sf[[sd_col]]), "\u2014", fmt(map_sf[[sd_col]])))
  if (!is.null(pred_col) && !is.null(sd_col))
    pop <- paste0(pop, "<br><strong>Relative uncertainty:</strong> ",
                  ifelse(is.na(map_sf$.rel_unc), "\u2014",
                         sprintf("%.3f", map_sf$.rel_unc)))
  if (has("aadt")) {
    src <- if (has("traffic_volume_source"))
      paste0(" (", ifelse(is.na(map_sf$traffic_volume_source), "unknown",
                          map_sf$traffic_volume_source), ")")
    else ""
    pop <- paste0(pop, "<br><strong>Measured AADT:</strong> ",
                  ifelse(is.na(map_sf$aadt), "unmeasured", fmt(map_sf$aadt)), src)
  }
  if (has("lastYearAadt_aadt"))
    pop <- paste0(pop, "<br><strong>Last year AADT:</strong> ",
                  ifelse(is.na(map_sf$lastYearAadt_aadt), "\u2014",
                         fmt(map_sf$lastYearAadt_aadt)))
  pop <- paste0(pop, "<br><strong>", color_by, ":</strong> ",
                as.character(map_sf[[color_by]]))
  map_sf$popup_text <- pop

  leaflet::leaflet(map_sf, options = leaflet::leafletOptions(preferCanvas = TRUE)) |>
    leaflet::addTiles() |>
    leaflet::addPolylines(
      color            = ~pal(map_sf[[color_by]]),
      weight           = weight,
      opacity          = opacity,
      popup            = ~popup_text,
      highlightOptions = leaflet::highlightOptions(
        weight = weight + 2, color = "white", bringToFront = TRUE)) |>
    leaflet::addLegend(
      position = "bottomright",
      pal      = pal,
      values   = ~map_sf[[color_by]],
      title    = title,
      opacity  = 0.9)
}
