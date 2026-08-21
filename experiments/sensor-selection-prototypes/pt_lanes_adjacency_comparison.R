library(INLA)
library(xdtkit)
library(jsonlite)
library(ggplot2)

source("R/scripts_inla_sensor_placement.R")

# ════════════════════════════════════════════════════════════════════════════════
# Data loading and preprocessing  ----
# (identical pipeline to testing-norway.R)
# ════════════════════════════════════════════════════════════════════════════════

year <- 2024

norway_directed_traffic_links <- jsonlite::read_json(
  path = "data/directed-traffic-links-2024.json",
  simplifyVector   = TRUE,
  simplifyDataFrame = TRUE
)

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

data <- fill_missing_values(
  df = preprocessed_traffic_links_norway,
  unknown_impute_columns  = c("functionClass", "highestSpeedLimit",
                               "lowestSpeedLimit", "maxLanes", "minLanes"),
  mode_impute_columns     = c("hasOnlyPublicTransportLanes"),
  median_impute_columns   = c("lastYearAadt_aadt", "lastYearAadt_heavyRatio",
                               "lastYearAadt_heavyAadt")
) |>
  remove_negative_aadt() |>
  add_logLastYear()       |>
  join_bus_to_traffic(bus_aadt_norway)

# ════════════════════════════════════════════════════════════════════════════════
# Shared model inputs  ----
# ════════════════════════════════════════════════════════════════════════════════

covariates <- ~ functionalRoadClass:maxLanes +
  functionalRoadClass:roadCategory +
  minLanes:roadCategory +
  functionalRoadClass +
  maxLanes +
  roadCategory +
  hasOnlyPublicTransportLanes +
  functionalRoadClass * isRamp +
  lastYearAadt_logAadt

similarity_covariates_weighted <- c(
  "lastYearAadt_logAadt",
  "functionalRoadClass",
  "roadCategory",
  "hasOnlyPublicTransportLanes"
)

ordinal_levels_weighted <- list(
  functionalRoadClass = as.character(c(7, 6, 5, 4, 3, 2, 1, 0)),
  roadCategory        = c("KOMMUNAL_VEG", "FYLKESVEG", "RIKSVEG", "EUROPAVEG")
)

# ────────────────────────────────────────────────────────────────────────────────
# Two adjacency matrices
# ────────────────────────────────────────────────────────────────────────────────

adj_excl_pt <- build_adjacency_matrix(data, exclude_public_transport = TRUE)
adj_incl_pt <- build_adjacency_matrix(data, exclude_public_transport = FALSE)

message("adj_excl_pt non-zeros: ", Matrix::nnzero(adj_excl_pt))
message("adj_incl_pt non-zeros: ", Matrix::nnzero(adj_incl_pt))

# ────────────────────────────────────────────────────────────────────────────────
# PT-link indices
# ────────────────────────────────────────────────────────────────────────────────

pt_idx <- which(data$hasOnlyPublicTransportLanes == TRUE)

# PT-link neighbours: non-PT links adjacent to at least one PT link in adj_incl_pt.
# These are at risk of deflated AADT estimates when PT links are included in the
# adjacency matrix because spatial smoothing toward low-traffic PT links can pull
# their estimates down.
adj_incl_trip    <- Matrix::summary(as(adj_incl_pt, "dgCMatrix"))
pt_neighbour_idx <- unique(c(
  adj_incl_trip$i[adj_incl_trip$j %in% pt_idx],
  adj_incl_trip$j[adj_incl_trip$i %in% pt_idx]
))
pt_neighbour_idx <- setdiff(pt_neighbour_idx, pt_idx)   # exclude PT links themselves

message(
  "PT-only links     : ", length(pt_idx),
  " (", sum(!is.na(data$aadt[pt_idx])), " with measured aadt)"
)
message(
  "PT-link neighbours: ", length(pt_neighbour_idx),
  " (", sum(!is.na(data$aadt[pt_neighbour_idx])), " with measured aadt)"
)

# Best-available truth: measured aadt when present, else lastYearAadt_aadt.
# Ensures ground-truth sensor data is always preferred over the imputed proxy.
truth_aadt <- ifelse(!is.na(data$aadt), data$aadt, data$lastYearAadt_aadt)

# ════════════════════════════════════════════════════════════════════════════════
# TEST 1 — Full-data fit; evaluate PT-link predictions vs lastYearAadt_aadt  ----
# ════════════════════════════════════════════════════════════════════════════════

message("\n", strrep("═", 60))
message("  TEST 1 — Full-data fit, PT-link and PT-neighbour evaluation")
message(strrep("═", 60), "\n")

# Fit on ALL data with each adjacency matrix
fit_excl <- fit_weighted_inla_model(
  data                  = data,
  adjacency_matrix      = adj_excl_pt,
  spatial_term          = "besagproper_rbf",
  fixed_effects         = covariates,
  iid_effects           = "roadSystem",
  similarity_covariates = similarity_covariates_weighted,
  ordinal_levels        = ordinal_levels_weighted,
  family                = "nbinomial",
  weight_fn             = "laplacian",
  distance_type         = "gower"
)

fit_incl <- fit_weighted_inla_model(
  data                  = data,
  adjacency_matrix      = adj_incl_pt,
  spatial_term          = "besagproper_rbf",
  fixed_effects         = covariates,
  iid_effects           = "roadSystem",
  similarity_covariates = similarity_covariates_weighted,
  ordinal_levels        = ordinal_levels_weighted,
  family                = "nbinomial",
  weight_fn             = "laplacian",
  distance_type         = "gower"
)

# ── Build per-group comparison tables ────────────────────────────────────────

# truth_aadt = measured aadt when available, else lastYearAadt_aadt (see above).
pt_compare <- data.frame(
  id          = data$id[pt_idx],
  group       = "PT-only",
  truth_aadt  = truth_aadt[pt_idx],
  meas_source = ifelse(!is.na(data$aadt[pt_idx]), "measured", "lastYearAadt"),
  pred_excl   = fit_excl$predictions$inla_pred[pt_idx],
  pred_incl   = fit_incl$predictions$inla_pred[pt_idx]
)

nb_compare <- data.frame(
  id          = data$id[pt_neighbour_idx],
  group       = "PT-neighbour",
  truth_aadt  = truth_aadt[pt_neighbour_idx],
  meas_source = ifelse(!is.na(data$aadt[pt_neighbour_idx]), "measured", "lastYearAadt"),
  pred_excl   = fit_excl$predictions$inla_pred[pt_neighbour_idx],
  pred_incl   = fit_incl$predictions$inla_pred[pt_neighbour_idx]
)

cat("\n── PT-link predictions ───────────────────────────────────────────────\n")
print(pt_compare[, c("id", "truth_aadt", "meas_source", "pred_excl", "pred_incl")],
      row.names = FALSE)

cat("\n── PT-neighbour predictions (first 20 rows) ──────────────────────────\n")
print(head(nb_compare[, c("id", "truth_aadt", "meas_source", "pred_excl", "pred_incl")], 20),
      row.names = FALSE)

# ── Metrics helper ────────────────────────────────────────────────────────────

eval_metrics <- function(df, pred_col, config_label) {
  ok   <- !is.na(df$truth_aadt) & df$truth_aadt > 0
  pred <- df[[pred_col]][ok]
  obs  <- df$truth_aadt[ok]
  data.frame(
    group  = df$group[1],
    config = config_label,
    n      = sum(ok),
    RMSE   = .rmse(pred, obs),
    MAE    = .mae(pred,  obs),
    MALE   = .male(pred, obs)
  )
}

metrics_t1 <- do.call(rbind, list(
  eval_metrics(pt_compare, "pred_excl", "excl_PT (adj)"),
  eval_metrics(pt_compare, "pred_incl", "incl_PT (adj)"),
  eval_metrics(nb_compare, "pred_excl", "excl_PT (adj)"),
  eval_metrics(nb_compare, "pred_incl", "incl_PT (adj)")
))

cat("\n── Test 1: metrics by group and adjacency config ─────────────────────\n")
print(metrics_t1, digits = 4, row.names = FALSE)

# ── Scatter plots ─────────────────────────────────────────────────────────────

plot_group_scatter <- function(df, group_label) {
  long <- rbind(
    data.frame(truth = df$truth_aadt, pred = df$pred_excl,
               config = "excl_PT (adj)", source = df$meas_source),
    data.frame(truth = df$truth_aadt, pred = df$pred_incl,
               config = "incl_PT (adj)", source = df$meas_source)
  )
  long <- long[!is.na(long$truth) & long$truth > 0, ]
  lims <- range(c(long$truth, long$pred), na.rm = TRUE)

  ggplot(long, aes(x = truth, y = pred, colour = config, shape = source)) +
    geom_abline(slope = 1, intercept = 0,
                colour = "firebrick", linetype = "dashed", linewidth = 0.8) +
    geom_point(alpha = 0.75, size = 2.5) +
    coord_fixed(xlim = lims, ylim = lims) +
    scale_shape_manual(values = c(measured = 16, lastYearAadt = 1)) +
    labs(
      x      = "Truth AADT (measured if available, else lastYearAadt_aadt)",
      y      = "Predicted AADT",
      colour = NULL,
      shape  = "Truth source",
      title  = paste0("Test 1 — ", group_label, ": predicted vs truth")
    ) +
    theme_bw(base_size = 11)
}

print(plot_group_scatter(pt_compare, "PT-only links"))
print(plot_group_scatter(nb_compare, "PT-link neighbours"))

# ════════════════════════════════════════════════════════════════════════════════
# TEST 2 — K-fold CV; overall MALE comparison  ----
# ════════════════════════════════════════════════════════════════════════════════

message("\n", strrep("═", 60))
message("  TEST 2 — 5-fold CV, overall MALE (adj_excl vs adj_incl)")
message(strrep("═", 60), "\n")

# k = 5: PT links are few; higher k risks folds with no measurable links.
# Same seed for both to ensure identical fold assignments on measured links.

cv_excl <- kfold_cv_inla(
  data                  = data,
  adjacency_matrix      = adj_excl_pt,
  k                     = 5,
  seed                  = 42,
  fixed_effects         = covariates,
  spatial_term          = "besagproper_rbf",
  similarity_covariates = similarity_covariates_weighted,
  ordinal_levels        = ordinal_levels_weighted,
  iid_effects           = "roadSystem",
  family                = "nbinomial",
  weight_fn             = "laplacian",
  distance_type         = "gower"
)

cv_incl <- kfold_cv_inla(
  data                  = data,
  adjacency_matrix      = adj_incl_pt,
  k                     = 5,
  seed                  = 42,
  fixed_effects         = covariates,
  spatial_term          = "besagproper_rbf",
  similarity_covariates = similarity_covariates_weighted,
  ordinal_levels        = ordinal_levels_weighted,
  iid_effects           = "roadSystem",
  family                = "nbinomial",
  weight_fn             = "laplacian",
  distance_type         = "gower"
)

# ── Side-by-side fold metrics ─────────────────────────────────────────────────

fold_comparison <- merge(
  cv_excl$fold_metrics,
  cv_incl$fold_metrics,
  by   = "fold",
  suffixes = c("_excl", "_incl")
)

cat("\n── Test 2: fold-level MALE ───────────────────────────────────────────\n")
print(fold_comparison[, c("fold", "MALE_excl", "MALE_incl")],
      digits = 4, row.names = FALSE)

summary_comparison <- data.frame(
  config   = c("excl_PT (adj)", "incl_PT (adj)"),
  mean_MALE = c(mean(cv_excl$fold_metrics$MALE), mean(cv_incl$fold_metrics$MALE)),
  sd_MALE   = c(sd(cv_excl$fold_metrics$MALE),   sd(cv_incl$fold_metrics$MALE)),
  mean_RMSE = c(mean(cv_excl$fold_metrics$RMSE),  mean(cv_incl$fold_metrics$RMSE)),
  mean_MAE  = c(mean(cv_excl$fold_metrics$MAE),   mean(cv_incl$fold_metrics$MAE)),
  mean_MAPE = c(mean(cv_excl$fold_metrics$MAPE),  mean(cv_incl$fold_metrics$MAPE))
)

cat("\n── Test 2: summary metrics ───────────────────────────────────────────\n")
print(summary_comparison, digits = 4, row.names = FALSE)

# ── Pred vs obs plots ─────────────────────────────────────────────────────────

plot_kfold_cv(cv_excl, type = "pred_vs_obs")
plot_kfold_cv(cv_incl, type = "pred_vs_obs")

# ── MALE per fold — side-by-side bar chart ────────────────────────────────────

male_long <- rbind(
  data.frame(fold = factor(cv_excl$fold_metrics$fold),
             MALE = cv_excl$fold_metrics$MALE,
             config = "excl_PT (adj)"),
  data.frame(fold = factor(cv_incl$fold_metrics$fold),
             MALE = cv_incl$fold_metrics$MALE,
             config = "incl_PT (adj)")
)

p_male <- ggplot(male_long, aes(x = fold, y = MALE, fill = config)) +
  geom_col(position = "dodge", width = 0.65) +
  geom_hline(
    data = summary_comparison,
    aes(yintercept = mean_MALE, colour = config),
    linetype = "dashed", linewidth = 0.8, show.legend = FALSE
  ) +
  labs(
    x     = "Fold",
    y     = "MALE",
    fill  = NULL,
    title = "Test 2 — MALE per fold: excl vs incl PT links in adjacency"
  ) +
  theme_bw(base_size = 11)

print(p_male)

# ════════════════════════════════════════════════════════════════════════════════
# Save results  ----
# ════════════════════════════════════════════════════════════════════════════════

saveRDS(
  list(
    fit_excl           = fit_excl,
    fit_incl           = fit_incl,
    pt_compare         = pt_compare,
    nb_compare         = nb_compare,
    metrics_t1         = metrics_t1,
    cv_excl            = cv_excl,
    cv_incl            = cv_incl,
    summary_comparison = summary_comparison
  ),
  file = "data/pt_adjacency_comparison.rds"
)

message("\nResults saved to data/pt_adjacency_comparison.rds")
