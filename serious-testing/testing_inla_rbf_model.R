#install.packages("pak")
#install.packages("INLA", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)
#inla.binary.install(os = c("Ubuntu-22.04"))
#pak::pak("trafikkdata/xdtkit")
source("serious-testing/scripts_inla_sensor_placement.R")
library(xdtkit)
library(ggplot2)

# =============================================================================
# Data preparation ----
# =============================================================================

year <- 2025

preprocessed_traffic_links <- preprocess_traffic_links(buskerud_directed_traffic_links, year = year)
bus_aadt <- calculate_bus_aadt(stops_on_traffic_links, bus_counts, year = year)

prepared_traffic_links <- fill_missing_values(
  df                     = preprocessed_traffic_links,
  unknown_impute_columns = c("functionClass", "highestSpeedLimit", "lowestSpeedLimit", "maxLanes", "minLanes"),
  mode_impute_columns    = c("hasOnlyPublicTransportLanes"),
  median_impute_columns  = c("lastYearAadt_aadt", "lastYearAadt_heavyRatio", "lastYearAadt_heavyAadt")
) |>
  remove_negative_aadt() |>
  add_logLastYear() |>
  join_bus_to_traffic(bus_aadt)

adjacency_matrix <- build_adjacency_matrix(
  prepared_traffic_links,
  exclude_public_transport = TRUE
)

# clusters <- strategic_network_clustering(
#   data           = prepared_traffic_links,
#   year           = year,
#   boundary_links = c("Trafikkdata_continuous")
# )
# nodes <- identify_unbalanceable_nodes(buskerud_nodes, prepared_traffic_links)

# =============================================================================
# Shared configuration ----
# =============================================================================

ordinal_levels_road <- list(
  functionalRoadClass = as.character(c(7, 6, 5, 4, 3, 2, 1, 0)),
  roadCategory        = c("KOMMUNAL_VEG", "FYLKESVEG", "RIKSVEG", "EUROPAVEG")
)

# Two experiment specs: with and without last-year AADT as a covariate
experiments <- list(

  with_last_year = list(
    label = "With last-year AADT",
    fixed_effects = ~ functionalRoadClass:maxLanes +
      functionalRoadClass:roadCategory +
      minLanes:roadCategory + functionalRoadClass +
      maxLanes + roadCategory +
      hasOnlyPublicTransportLanes +
      functionalRoadClass * isRamp +
      lastYearAadt_logAadt,
    similarity_covariates = c(
      "maxLanes", "minLanes",
      "highestSpeedLimit", "lowestSpeedLimit",
      "hasOnlyPublicTransportLanes", "isRamp",
      "roadCategory", "functionClass", "functionalRoadClass", "roadSystem",
      "lastYearAadt_logAadt"
    )
  ),

  no_last_year = list(
    label = "Without last-year AADT",
    fixed_effects = ~ functionalRoadClass:maxLanes +
      functionalRoadClass:roadCategory +
      minLanes:roadCategory + functionalRoadClass +
      maxLanes + roadCategory +
      hasOnlyPublicTransportLanes +
      functionalRoadClass * isRamp,
    similarity_covariates = c(
      "maxLanes", "minLanes",
      "highestSpeedLimit", "lowestSpeedLimit",
      "hasOnlyPublicTransportLanes", "isRamp",
      "roadCategory", "functionClass", "functionalRoadClass", "roadSystem"
    )
  )
)

# =============================================================================
# Train / test split ----
# =============================================================================

set.seed(420)
measured_idx <- which(!is.na(prepared_traffic_links$aadt))
test_idx     <- sample(measured_idx, size = round(0.2 * length(measured_idx)))

data_cv <- prepared_traffic_links
data_cv$aadt[test_idx]      <- NA
data_cv$heavyAadt[test_idx] <- NA

true_aadt <- prepared_traffic_links$aadt[test_idx]

# =============================================================================
# Experiment loop: Besagproper vs RBF-Besag x Gaussian vs NB ----
# =============================================================================

# Fit all 4 model variants for one experiment spec
fit_experiment <- function(exp, data_cv, adjacency_matrix, ordinal_levels) {
  base_args <- list(
    data             = data_cv,
    adjacency_matrix = adjacency_matrix,
    fixed_effects    = exp$fixed_effects,
    iid_effects      = "roadSystem",
    ordinal_levels   = ordinal_levels
  )
  list(
    gauss_besag = do.call(fit_inla_rbf_model, c(base_args,
      list(spatial_term = "besagproper",     family = "gaussian"))),
    gauss_rbf   = do.call(fit_inla_rbf_model, c(base_args,
      list(spatial_term = "besagproper_rbf", family = "gaussian",
           similarity_covariates = exp$similarity_covariates))),
    nb_besag    = do.call(fit_inla_rbf_model, c(base_args,
      list(spatial_term = "besagproper",     family = "nbinomial"))),
    nb_rbf      = do.call(fit_inla_rbf_model, c(base_args,
      list(spatial_term = "besagproper_rbf", family = "nbinomial",
           similarity_covariates = exp$similarity_covariates)))
  )
}

# Print metrics table and RBF lift summary for a named list of fitted models
print_metrics <- function(fits, test_idx, true_aadt, label) {
  preds <- lapply(fits, function(m) m$predictions$inla_pred[test_idx])
  metrics <- data.frame(
    Model = names(fits),
    RMSE  = vapply(preds, .rmse, numeric(1), o = true_aadt),
    MAE   = vapply(preds, .mae,  numeric(1), o = true_aadt),
    MAPE  = vapply(preds, .mape, numeric(1), o = true_aadt)
  )
  cat("\n", label, "\n", strrep("=", nchar(label)), "\n", sep = "")
  print(metrics, digits = 4, row.names = FALSE)
  cat(sprintf("\nRBF lift (Besag - RBF, positive = RBF better):\n"))
  cat(sprintf("  Gaussian  RMSE: %+.0f  MAE: %+.0f  MAPE: %+.1f%%\n",
              metrics$RMSE[1] - metrics$RMSE[2],
              metrics$MAE[1]  - metrics$MAE[2],
              metrics$MAPE[1] - metrics$MAPE[2]))
  cat(sprintf("  NB        RMSE: %+.0f  MAE: %+.0f  MAPE: %+.1f%%\n",
              metrics$RMSE[3] - metrics$RMSE[4],
              metrics$MAE[3]  - metrics$MAE[4],
              metrics$MAPE[3] - metrics$MAPE[4]))
  invisible(metrics)
}

# ── Run experiments ────────────────────────────────────────────────────────────

for (exp_name in names(experiments)) {
  exp  <- experiments[[exp_name]]
  fits <- fit_experiment(exp, data_cv, adjacency_matrix, ordinal_levels_road)

  print_metrics(fits, test_idx, true_aadt, exp$label)

  for (model_name in names(fits)) {
    m <- fits[[model_name]]
    for (plt_type in c("pred_vs_obs", "residuals_vs_fitted", "qq", "residuals_hist")) {
      plot_inla_model(
        model    = m,
        observed = prepared_traffic_links$aadt,
        test_idx = test_idx,
        type     = plt_type,
        title    = paste(exp$label, "|", model_name)
      )
    }
  }
}

# =============================================================================
# K-fold cross-validation example ----
# =============================================================================

# Adjust k, seed, or experiment spec as needed
cv_result <- kfold_cv_inla(
  data                  = prepared_traffic_links,
  adjacency_matrix      = adjacency_matrix,
  k                     = 5,
  seed                  = 42,
  fixed_effects         = experiments$with_last_year$fixed_effects,
  spatial_term          = "besagproper_rbf",
  similarity_covariates = experiments$with_last_year$similarity_covariates,
  ordinal_levels        = ordinal_levels_road,
  iid_effects           = "roadSystem",
  family                = "gaussian"
)

# Diagnostic plots for CV result
for (plt_type in c("pred_vs_obs", "residuals_vs_fitted", "qq",
                   "residuals_hist", "metrics_by_fold")) {
  plot_kfold_cv(cv_result, type = plt_type)
}

