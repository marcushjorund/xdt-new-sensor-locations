# ════════════════════════════════════════════════════════════════════════════════
# plot_inla_model() ----
# ════════════════════════════════════════════════════════════════════════════════

#' Diagnostic plot for a fitted weighted_inla_traffic_model
#'
#' @param model   An \code{weighted_inla_traffic_model} object.
#' @param observed Numeric vector of true AADT values (full-length, same row order
#'   as the data passed to \code{fit_weighted_inla_model}).
#' @param test_idx Optional integer index vector.  When supplied, only these
#'   positions are used (e.g. a held-out test set).
#' @param type    One of \code{"pred_vs_obs"}, \code{"residuals_vs_fitted"},
#'   \code{"qq"}, \code{"residuals_hist"}.  Produces a single ggplot.
#' @param title   Optional title string.  Auto-generated if \code{NULL}.
#' @param col     Hex colour for points.  Defaults to a blue.
#' @return A \code{ggplot} object (printed automatically when not assigned).
#'   Also returns \code{invisible(list(fitted, observed, residuals))}.
plot_inla_model <- function(model,
                            observed,
                            test_idx = NULL,
                            type     = c("pred_vs_obs", "residuals_vs_fitted",
                                         "qq", "residuals_hist"),
                            title    = NULL,
                            col      = "#0066CC",
                            x_label  = NULL,
                            y_label  = NULL,
                            ylim     = NULL) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required for plot_inla_model(). Install it with install.packages('ggplot2').")
  
  type <- match.arg(type)
  
  pred_col <- if (!is.null(model$heavy_vehicle) && isTRUE(model$heavy_vehicle))
    "inla_pred_heavy" else "inla_pred"
  
  fitted <- model$predictions[[pred_col]]
  
  if (!is.null(test_idx)) {
    fitted   <- fitted[test_idx]
    observed <- observed[test_idx]
  }
  
  resid <- fitted - observed
  
  if (is.null(title)) {
    family_lab <- switch(model$family,
                         gaussian  = "Gaussian",
                         nbinomial = "NB",
                         poisson   = "Poisson",
                         model$family)
    spatial_lab <- switch(model$spatial_term,
                          besagproper     = "Besag",
                          besagproper_rbf = "Weighted Besag",
                          model$spatial_term)
    title <- paste(family_lab, "\u2014", spatial_lab)
  }
  
  pt_col   <- col
  pt_alpha <- 0.45
  pt_size  <- 1.4
  
  p <- switch(type,
              
              pred_vs_obs = {
                lims <- if (!is.null(ylim)) ylim else range(c(observed, fitted), na.rm = TRUE)
                df   <- data.frame(obs = observed, pred = fitted)
                ggplot2::ggplot(df, ggplot2::aes(x = obs, y = pred)) +
                  ggplot2::geom_abline(slope = 1, intercept = 0,
                                       colour = "firebrick", linetype = "dashed", linewidth = 0.8) +
                  ggplot2::geom_point(colour = pt_col, alpha = pt_alpha, size = pt_size) +
                  ggplot2::coord_fixed(xlim = lims, ylim = lims) +
                  ggplot2::annotate("label", x = lims[1], y = lims[2],
                                    label  = .metric_label(fitted, observed),
                                    hjust  = 0, vjust = 1, size = 3,
                                    family = "mono") +
                  ggplot2::labs(x = "Observed AADT", y = "Predicted AADT", title = title) +
                  ggplot2::theme_bw(base_size = 11)
              },
              
              residuals_vs_fitted = {
                df <- data.frame(fitted = fitted, resid = resid)
                p_rv <- ggplot2::ggplot(df, ggplot2::aes(x = fitted, y = resid)) +
                  ggplot2::geom_hline(yintercept = 0,
                                      colour = "firebrick", linetype = "dashed", linewidth = 0.8) +
                  ggplot2::geom_point(colour = pt_col, alpha = pt_alpha, size = pt_size) +
                  ggplot2::geom_smooth(method = "loess", formula = y ~ x,
                                       se = FALSE, colour = "black", linewidth = 0.9) +
                  ggplot2::labs(x = "Fitted AADT", y = "Residual (Fitted \u2212 Observed)",
                                title = title) +
                  ggplot2::theme_bw(base_size = 11)
                if (!is.null(ylim)) p_rv <- p_rv + ggplot2::coord_cartesian(ylim = ylim)
                p_rv
              },
              
              qq = {
                df <- data.frame(resid = resid[!is.na(resid)])
                p_qq <- ggplot2::ggplot(df, ggplot2::aes(sample = resid)) +
                  ggplot2::stat_qq(colour = pt_col, alpha = pt_alpha, size = pt_size) +
                  ggplot2::stat_qq_line(colour = "firebrick", linetype = "dashed", linewidth = 0.8) +
                  ggplot2::labs(x = "Theoretical quantiles", y = "Sample quantiles",
                                title = paste(title, "\u2014 Q-Q plot")) +
                  ggplot2::theme_bw(base_size = 11)
                if (!is.null(ylim)) p_qq <- p_qq + ggplot2::coord_cartesian(ylim = ylim)
                p_qq
              },
              
              residuals_hist = {
                df <- data.frame(resid = resid[!is.na(resid)])
                ggplot2::ggplot(df, ggplot2::aes(x = resid)) +
                  ggplot2::geom_histogram(fill = pt_col, colour = "white", alpha = 0.8,
                                          bins = 30) +
                  ggplot2::geom_rug(colour = pt_col, alpha = 0.3) +
                  ggplot2::geom_vline(xintercept = 0,
                                      colour = "firebrick", linetype = "dashed", linewidth = 0.8) +
                  ggplot2::labs(x = "Residual (Fitted \u2212 Observed)", y = "Count",
                                title = paste(title, "\u2014 Residual distribution")) +
                  ggplot2::theme_bw(base_size = 11)
              }
  )
  
  if (!is.null(x_label) || !is.null(y_label))
    p <- p + ggplot2::labs(x = x_label, y = y_label)

  print(p)
  invisible(p)
}

# ════════════════════════════════════════════════════════════════════════════════
# plot_kfold_cv() ----
# ════════════════════════════════════════════════════════════════════════════════

#' Diagnostic plots for k-fold CV results
#'
#' @param cv_result An \code{inla_kfold_cv} object from \code{kfold_cv_inla()}.
#' @param type      One of \code{"pred_vs_obs"}, \code{"residuals_vs_fitted"},
#'   \code{"qq"}, \code{"residuals_hist"}, \code{"metrics_by_fold"}.
#' @return A \code{ggplot} object (printed automatically).
plot_kfold_cv <- function(cv_result,
                          type = c("pred_vs_obs", "residuals_vs_fitted",
                                   "qq", "residuals_hist", "metrics_by_fold")) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required. Install with install.packages('ggplot2').")
  
  type <- match.arg(type)
  
  oof   <- cv_result$oof_predictions
  oof$fold  <- factor(oof$fold)
  oof$resid <- oof$pred_aadt - oof$true_aadt
  
  family_lab <- switch(cv_result$family,
                       gaussian  = "Gaussian",
                       nbinomial = "NB",
                       poisson   = "Poisson",
                       cv_result$family)
  spatial_lab <- switch(cv_result$spatial_term,
                        besagproper     = "Besag",
                        besagproper_rbf = "Weighted Besag",
                        cv_result$spatial_term)
  base_title <- paste0(family_lab, " \u2014 ", spatial_lab,
                       "  (", cv_result$k, "-fold CV)")
  
  p <- switch(type,
              
              pred_vs_obs = {
                lims <- range(c(oof$true_aadt, oof$pred_aadt), na.rm = TRUE)
                
                # Per-fold metrics for legend labels
                fold_labs <- vapply(levels(oof$fold), function(f) {
                  sub <- oof[oof$fold == f, ]
                  sprintf("Fold %s  RMSE=%.0f  MAE=%.0f  MAPE=%.1f%%  MALE=%.3f",
                          f,
                          .rmse(sub$pred_aadt, sub$true_aadt),
                          .mae( sub$pred_aadt, sub$true_aadt),
                          .mape(sub$pred_aadt, sub$true_aadt),
                          .male(sub$pred_aadt, sub$true_aadt))
                }, character(1))
                
                # Mean ± SD across folds for annotation
                fm    <- cv_result$fold_metrics
                annot <- sprintf("Mean \u00b1 SD\nRMSE = %.0f \u00b1 %.0f\nMAE  = %.0f \u00b1 %.0f\nMAPE = %.1f \u00b1 %.1f%%\nMALE = %.3f \u00b1 %.3f",
                                 mean(fm$RMSE), sd(fm$RMSE),
                                 mean(fm$MAE),  sd(fm$MAE),
                                 mean(fm$MAPE), sd(fm$MAPE),
                                 mean(fm$MALE), sd(fm$MALE))
                
                ggplot2::ggplot(oof, ggplot2::aes(x = true_aadt, y = pred_aadt,
                                                  colour = fold)) +
                  ggplot2::geom_abline(slope = 1, intercept = 0,
                                       colour = "firebrick", linetype = "dashed", linewidth = 0.8) +
                  ggplot2::geom_point(alpha = 0.45, size = 1.4) +
                  ggplot2::coord_fixed(xlim = lims, ylim = lims) +
                  ggplot2::annotate("label", x = lims[1], y = lims[2],
                                    label = annot, hjust = 0, vjust = 1, size = 3,
                                    family = "mono") +
                  ggplot2::scale_colour_discrete(labels = fold_labs) +
                  ggplot2::labs(x = "Observed AADT", y = "Predicted AADT",
                                colour = NULL, title = base_title) +
                  ggplot2::theme_bw(base_size = 11) +
                  ggplot2::theme(legend.text = ggplot2::element_text(family = "mono", size = 8))
              },
              
              residuals_vs_fitted = {
                ggplot2::ggplot(oof, ggplot2::aes(x = pred_aadt, y = resid,
                                                  colour = fold)) +
                  ggplot2::geom_hline(yintercept = 0,
                                      colour = "firebrick", linetype = "dashed", linewidth = 0.8) +
                  ggplot2::geom_point(alpha = 0.45, size = 1.4) +
                  ggplot2::geom_smooth(ggplot2::aes(group = 1),
                                       method = "loess", formula = y ~ x,
                                       se = FALSE, colour = "black", linewidth = 0.9) +
                  ggplot2::labs(x = "Fitted AADT", y = "Residual (Fitted \u2212 Observed)",
                                colour = "Fold", title = base_title) +
                  ggplot2::theme_bw(base_size = 11)
              },
              
              qq = {
                oof_clean <- oof[!is.na(oof$resid), ]
                ggplot2::ggplot(oof_clean, ggplot2::aes(sample = resid, colour = fold)) +
                  ggplot2::stat_qq(alpha = 0.45, size = 1.4) +
                  ggplot2::stat_qq_line(ggplot2::aes(group = 1),
                                        colour = "firebrick", linetype = "dashed",
                                        linewidth = 0.8) +
                  ggplot2::labs(x = "Theoretical quantiles", y = "Sample quantiles",
                                colour = "Fold",
                                title = paste(base_title, "\u2014 Q-Q plot")) +
                  ggplot2::theme_bw(base_size = 11)
              },
              
              residuals_hist = {
                oof_clean <- oof[!is.na(oof$resid), ]
                ggplot2::ggplot(oof_clean, ggplot2::aes(x = resid, fill = fold)) +
                  ggplot2::geom_histogram(colour = "white", alpha = 0.7,
                                          bins = 30, position = "identity") +
                  ggplot2::geom_vline(xintercept = 0,
                                      colour = "firebrick", linetype = "dashed",
                                      linewidth = 0.8) +
                  ggplot2::labs(x = "Residual (Fitted \u2212 Observed)", y = "Count",
                                fill = "Fold",
                                title = paste(base_title, "\u2014 Residual distribution")) +
                  ggplot2::theme_bw(base_size = 11)
              },
              
              metrics_by_fold = {
                fm      <- cv_result$fold_metrics
                fm$fold <- factor(fm$fold)
                # Pivot to long form without tidyr dependency
                long <- rbind(
                  data.frame(fold = fm$fold, metric = "RMSE", value = fm$RMSE),
                  data.frame(fold = fm$fold, metric = "MAE",  value = fm$MAE),
                  data.frame(fold = fm$fold, metric = "MAPE", value = fm$MAPE),
                  data.frame(fold = fm$fold, metric = "MALE", value = fm$MALE)
                )
                long$metric <- factor(long$metric, levels = c("RMSE", "MAE", "MAPE", "MALE"))
                
                means <- aggregate(value ~ metric, data = long, FUN = mean)
                
                ggplot2::ggplot(long, ggplot2::aes(x = fold, y = value, fill = fold)) +
                  ggplot2::geom_col(width = 0.65, show.legend = FALSE) +
                  ggplot2::geom_hline(data = means,
                                      ggplot2::aes(yintercept = value),
                                      colour = "firebrick", linetype = "dashed",
                                      linewidth = 0.8) +
                  ggplot2::facet_wrap(~ metric, scales = "free_y") +
                  ggplot2::labs(x = "Fold", y = NULL,
                                title = paste(base_title, "\u2014 Metrics by fold")) +
                  ggplot2::theme_bw(base_size = 11) +
                  ggplot2::theme(strip.background = ggplot2::element_blank(),
                                 strip.text = ggplot2::element_text(face = "bold"))
              }
  )
  
  print(p)
  invisible(p)
}

# ════════════════════════════════════════════════════════════════════════════════
# plot_sensor_selection_map() ----
# ════════════════════════════════════════════════════════════════════════════════

#' Interactive leaflet map of sensor selection results
#'
#' Minimal version: plots only the selected sensors (red) and their neighbour
#' context links (steel blue) from \code{selected_data_entries}.
#' Permanently visible numbered circle markers are placed at the midpoint of
#' each selected link, ranked by \code{mi_score} descending.
#'
#' @param sensor_results  List from \code{greedy_mi_sensor_selection_norway()}
#'   (must contain \code{$selected_data_entries}), or the data frame directly.
#' @param color_selected  Hex colour for selected links. Default \code{"#E63946"}.
#' @param color_neighbour Hex colour for neighbour links. Default \code{"#457B9D"}.
#' @param weight_selected   Line weight for selected links (default 6).
#' @param weight_neighbour  Line weight for neighbour links (default 3).
#' @param opacity_selected  Opacity for selected links (default 0.95).
#' @param opacity_neighbour Base opacity for neighbour links (default 0.7).
#' @param all_data          Optional data frame of all traffic links. When provided, links with
#'   a measured \code{aadt} value are overlaid as dashed lines (layer group
#'   "M\u00e5lte trafikklenker"), including those outside the selected neighbourhood.
#' @param color_measured    Hex colour for the measured-links dash overlay (default \code{"#333333"}).
#'
#' @return A \code{leaflet} map object.
#' @export
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
  
  # ── 2. Geometry: add_geometries() directly on sel_df ─────────────────────────
  # sel_df from greedy_mi_sensor_selection already carries 'selected' and
  # 'mi_score' as plain columns — call add_geometries() once, no joins on sf.
  sel_sf <- add_geometries(sel_df)
  sel_sf <- sf::st_transform(sel_sf, crs = 4326)
  
  # Drop rows with empty or invalid geometry (prevents lat/lon warnings)
  valid <- !sf::st_is_empty(sel_sf) & sf::st_is_valid(sel_sf)
  sel_sf <- sel_sf[valid, ]
  
  # ── 3. Rank + popup (before split so subsets inherit both columns) ────────────
  # Rank is per SENSOR (base ID = strip -WITH/-AGAINST), so both directions of
  # the same sensor share the same rank number.  With r=100, ranks run 1..100.
  sel_sf$base_id <- sub("-WITH$|-AGAINST$", "", sel_sf$id)
  sel_sf$rank <- NA_integer_
  true_idx <- which(sel_sf$selected %in% TRUE)

  # Detect derivable-algorithm output: use selection_order as bundle key when
  # is_derived and selection_order columns are both present with any derived rows.
  has_derived_col <- "is_derived"      %in% names(sel_sf)
  has_sel_order   <- "selection_order" %in% names(sel_sf)
  is_derived_vec  <- if (has_derived_col) sel_sf$is_derived %in% TRUE else rep(FALSE, nrow(sel_sf))
  use_bundle_rank <- has_derived_col && has_sel_order &&
                     length(true_idx) > 0 && any(is_derived_vec[true_idx])

  if (length(true_idx) > 0) {
    true_sf <- sel_sf[true_idx, ]
    if (use_bundle_rank) {
      # Rank by bundle (selection_order); derived links share order with their sensor.
      bundle_scores <- tapply(true_sf$mi_score, true_sf$selection_order,
                              max, na.rm = TRUE)
      bundle_order  <- as.integer(names(sort(bundle_scores, decreasing = TRUE)))
      rank_lookup   <- stats::setNames(seq_along(bundle_order),
                                       as.character(bundle_order))
      sel_sf$rank[true_idx] <- rank_lookup[as.character(sel_sf$selection_order[true_idx])]
    } else if ("mi_score" %in% names(true_sf)) {
      # Unique base IDs ordered by mi_score descending (same score for WITH+AGAINST)
      base_order  <- unique(true_sf$base_id[order(-true_sf$mi_score, na.last = TRUE)])
      rank_lookup <- stats::setNames(seq_along(base_order), base_order)
      sel_sf$rank[true_idx] <- rank_lookup[sel_sf$base_id[true_idx]]
    } else {
      base_order  <- unique(true_sf$base_id)
      rank_lookup <- stats::setNames(seq_along(base_order), base_order)
      sel_sf$rank[true_idx] <- rank_lookup[sel_sf$base_id[true_idx]]
    }
  }

  has_mi  <- "mi_score"                  %in% names(sel_sf)
  has_ly  <- "lastYearAadt_aadt"         %in% names(sel_sf)
  has_co  <- "county"                    %in% names(sel_sf)
  has_src <- "traffic_volume_source"     %in% names(sel_sf)
  has_frc <- "functionalRoadClass"       %in% names(sel_sf)
  has_fc  <- "functionClass"             %in% names(sel_sf)
  has_hr  <- "lastYearAadt_heavyRatio"   %in% names(sel_sf)

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
  
  # ── 4. Split (inherits popup_text, rank, base_id) ────────────────────────────
  sel_sensor     <- sel_sf[sel_sf$selected %in% TRUE & !is_derived_vec,  ]
  sel_derived_sf <- sel_sf[sel_sf$selected %in% TRUE & is_derived_vec,   ]
  sel_neighbour  <- sel_sf[sel_sf$selected %in% FALSE, ]

  # ── 5. Midpoint markers: one per sensor base_id (red) + one per derived base_id (orange) ──
  markers_sf <- NULL
  if (nrow(sel_sensor) > 0) {
    # Keep one representative directed link per base_id (first occurrence = -WITH)
    one_per_sensor <- sel_sensor[!duplicated(sel_sensor$base_id), ]
    # Sample at 50% arc length (true midpoint) rather than st_point_on_surface,
    # which tends to land near one vertex/end of the LineString.
    suppressWarnings({
      mid <- sf::st_line_sample(sf::st_transform(one_per_sensor, 25833), sample = 0.5)
      mid <- sf::st_cast(mid, "POINT")
      mid <- sf::st_transform(mid, 4326)
    })
    markers_sf <- sf::st_set_geometry(one_per_sensor, mid)
    coords <- sf::st_coordinates(markers_sf)
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
    coords_d <- sf::st_coordinates(derived_markers_sf)
    derived_markers_sf <- derived_markers_sf[!is.na(coords_d[, 1]) & !is.na(coords_d[, 2]), ]
    if (nrow(derived_markers_sf) == 0) derived_markers_sf <- NULL
  }
  
  # ── 6. Build leaflet map ──────────────────────────────────────────────────────
  nvdb <- nvdb_objects()
  
  m <- leaflet::leaflet(
    sel_sf,
    options = leaflet::leafletOptions(crs = nvdb$nvdb_crs, zoomControl = TRUE)
  ) |>
    leaflet::addTiles(
  urlTemplate = "https://services.geodataonline.no/arcgis/rest/services/Geocache_UTM33_EUREF89/GeocacheGraatone/MapServer/tile/{z}/{y}/{x}",
  attribution = "\u00a9 Geodata AS, Kartverket, Geovekst og kommunene, OpenStreetMap"
)
  
  # Neighbour links — one addPolylines per hop, uniform opacity and weight
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
  
  # Derived links — orange, same weight as selected
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

  # Selected links — red, thick
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

  # ── . Measured traffic links (dashed overlay) ──────────────────────────────
  # Computed and added here — before circle markers — so it renders above all
  # polyline layers (neighbours, derived, selected) but below the rank markers.
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

  # Layer controls + legend
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

  # Derived link markers — smaller orange circles, grouped so they hide with
  # "Utledede trafikklenker" layer toggle. Added before red markers so red
  # renders on top.
  if (!is.null(derived_markers_sf))
    m <- m |> leaflet::addCircleMarkers(
      data         = derived_markers_sf,
      radius       = 6,
      color        = "#000000",
      weight       = 1.5,
      fillColor    = color_derived,
      fillOpacity  = 0.95,
      popup        = derived_markers_sf$popup_text,
      group        = "Utledede trafikklenker",
      label        = as.character(derived_markers_sf$rank),
      labelOptions = leaflet::labelOptions(
        permanent = TRUE,
        textOnly  = TRUE,
        style     = list("font-weight" = "bold",
                         "font-size"   = "10px",
                         "color"       = "white",
                         "text-shadow" = "0 0 3px #000")))

  # Numbered circle markers — no group → always visible; added last so they
  # render on top of everything including derived (yellow) markers.
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
  m
}


# ════════════════════════════════════════════════════════════════════════════════
# map_traffic_links() ----
# ════════════════════════════════════════════════════════════════════════════════

#' Interactive leaflet map of traffic links
#'
#' General-purpose map that mirrors the popup content of
#' \code{xdtkit::plot_traffic_links_map} (link ID, predicted AADT, SD, measured
#' AADT + source, relative uncertainty, last-year AADT) while accepting any
#' \code{color_by} column for flexible colouring.  Fields that are absent from
#' \code{data} are silently omitted from the popup.
#'
#' @param data     A data frame with at least an \code{id} column compatible
#'   with \code{add_geometries()}.
#' @param color_by Name of the column used for colouring links.
#' @param title    Legend title. Defaults to \code{color_by}.
#' @param weight   Polyline weight (default 2).
#' @param opacity  Polyline opacity (default 0.9).
#'
#' @return A \code{leaflet} map object.
map_traffic_links <- function(
    data,
    color_by,
    title   = color_by,
    weight  = 2,
    opacity = 0.9) {

  if (!color_by %in% names(data))
    stop("Column '", color_by, "' not found in data.")

  # ── 1. Add geometries and reproject to WGS84 ─────────────────────────────────
  map_sf <- add_geometries(data) |>
    sf::st_transform(crs = 4326)

  map_sf <- map_sf[!sf::st_is_empty(map_sf) & sf::st_is_valid(map_sf), ]

  if (nrow(map_sf) == 0L)
    stop("No valid geometries found after add_geometries().")

  # ── 2. Build colour palette (factor/logical/char → colorFactor, numeric → viridis) ──
  col_values <- map_sf[[color_by]]
  if (is.logical(col_values) || is.character(col_values) || is.factor(col_values)) {
    map_sf[[color_by]] <- as.character(col_values)
    pal <- leaflet::colorFactor(
      palette  = "Set1",
      domain   = map_sf[[color_by]],
      na.color = "#888888"
    )
  } else {
    pal <- leaflet::colorNumeric(
      palette  = "viridis",
      domain   = col_values,
      na.color = "#888888"
    )
  }

  # ── 3. Build popup — same fields as plot_traffic_links_map, present columns only ──
  has <- function(col) col %in% names(map_sf)

  # Resolve prediction / SD column names (prefer balanced over raw INLA)
  pred_col <- if (has("balanced_pred")) "balanced_pred" else if (has("inla_pred")) "inla_pred" else NULL
  sd_col   <- if (has("balanced_sd"))   "balanced_sd"   else if (has("inla_sd"))   "inla_sd"   else NULL

  # Relative uncertainty
  if (!is.null(pred_col) && !is.null(sd_col)) {
    map_sf$.rel_unc <- ifelse(
      !is.na(map_sf[[pred_col]]) & map_sf[[pred_col]] > 0,
      map_sf[[sd_col]] / map_sf[[pred_col]],
      NA_real_
    )
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

  # ── 4. Build leaflet map ──────────────────────────────────────────────────────
  m <- leaflet::leaflet(map_sf, options = leaflet::leafletOptions(preferCanvas = TRUE)) |>
    leaflet::addTiles() |>
    leaflet::addPolylines(
      color            = ~pal(map_sf[[color_by]]),
      weight           = weight,
      opacity          = opacity,
      popup            = ~popup_text,
      highlightOptions = leaflet::highlightOptions(
        weight       = weight + 2,
        color        = "white",
        bringToFront = TRUE
      )
    ) |>
    leaflet::addLegend(
      position = "bottomright",
      pal      = pal,
      values   = ~map_sf[[color_by]],
      title    = title,
      opacity  = 0.9
    )

  m
}


