source("R/scripts_inla_sensor_placement.R")
library(xdtkit)

# ── Data loading and preprocessing ───────────────────────────────────────────
year <- 2025
preprocessed_traffic_links <- preprocess_traffic_links(buskerud_directed_traffic_links, year = year)
bus_aadt <- calculate_bus_aadt(stops_on_traffic_links, bus_counts, year = year)

prepared_traffic_links <- fill_missing_values(
  df = preprocessed_traffic_links,
  unknown_impute_columns  = c("functionClass", "highestSpeedLimit", "lowestSpeedLimit", "maxLanes", "minLanes"),
  mode_impute_columns     = c("hasOnlyPublicTransportLanes"),
  median_impute_columns   = c("lastYearAadt_aadt", "lastYearAadt_heavyRatio", "lastYearAadt_heavyAadt")
) |>
  remove_negative_aadt() |>
  add_logLastYear() |>
  join_bus_to_traffic(bus_aadt)

adjacency_matrix <- build_adjacency_matrix(prepared_traffic_links, exclude_public_transport = TRUE)

# ── Model specs ───────────────────────────────────────────────────────────────
fixed_effects_no_last <- ~ functionalRoadClass:maxLanes +
  functionalRoadClass:roadCategory +
  minLanes:roadCategory + functionalRoadClass +
  maxLanes + roadCategory +
  hasOnlyPublicTransportLanes +
  functionalRoadClass * isRamp

fixed_effects_with_last <- ~ functionalRoadClass:maxLanes +
  functionalRoadClass:roadCategory +
  minLanes:roadCategory + functionalRoadClass +
  maxLanes + roadCategory +
  hasOnlyPublicTransportLanes +
  functionalRoadClass * isRamp + lastYearAadt_logAadt

similarity_covariates_no_last <- c("maxLanes", "minLanes",
                                   "highestSpeedLimit", "lowestSpeedLimit",
                                   "hasOnlyPublicTransportLanes", "isRamp",
                                   "roadCategory", "functionClass", "functionalRoadClass", "roadSystem")

similarity_covariates_with_last <- c("maxLanes", "minLanes",
                                     "highestSpeedLimit", "lowestSpeedLimit",
                                     "hasOnlyPublicTransportLanes", "isRamp",
                                     "roadCategory", "functionClass", "functionalRoadClass", "roadSystem",
                                     "lastYearAadt_logAadt")

similarity_interaction_pairs <- list(
  c("functionalRoadClass", "roadCategory"),
  c("maxLanes",            "roadCategory")
)

# ── CV split ──────────────────────────────────────────────────────────────────
set.seed(420)
measured_idx <- which(!is.na(prepared_traffic_links$aadt))
test_idx     <- sample(measured_idx, size = round(0.2 * length(measured_idx)))

data_cv                    <- prepared_traffic_links
data_cv$aadt[test_idx]     <- NA
data_cv$heavyAadt[test_idx]<- NA
true_aadt                  <- prepared_traffic_links$aadt[test_idx]

rmse <- function(p, o) sqrt(mean((p - o)^2, na.rm = TRUE))
mae  <- function(p, o) mean(abs(p - o), na.rm = TRUE)
mape <- function(p, o) mean(abs((p - o) / o), na.rm = TRUE) * 100

# ── 8 model fits ─────────────────────────────────────────────────────────────
message("\n=== [1/8] Gaussian Weighted — WITH lastYear ===")
cv_g_weighted_with <- fit_weighted_inla_model(
  data                         = data_cv,
  adjacency_matrix             = adjacency_matrix,
  spatial_term                 = "besagproper_rbf",
  fixed_effects                = fixed_effects_with_last,
  iid_effects                  = "roadSystem",
  similarity_covariates        = similarity_covariates_with_last,
  similarity_interaction_pairs = similarity_interaction_pairs,
  family                       = "gaussian"
)

message("\n=== [2/8] Gaussian Besagproper — WITH lastYear ===")
cv_g_besag_with <- fit_weighted_inla_model(
  data             = data_cv,
  adjacency_matrix = adjacency_matrix,
  spatial_term     = "besagproper",
  fixed_effects    = fixed_effects_with_last,
  iid_effects      = "roadSystem",
  family           = "gaussian"
)

message("\n=== [3/8] NB Weighted — WITH lastYear ===")
cv_nb_weighted_with <- fit_weighted_inla_model(
  data                         = data_cv,
  adjacency_matrix             = adjacency_matrix,
  spatial_term                 = "besagproper_rbf",
  fixed_effects                = fixed_effects_with_last,
  iid_effects                  = "roadSystem",
  similarity_covariates        = similarity_covariates_with_last,
  similarity_interaction_pairs = similarity_interaction_pairs,
  family                       = "nbinomial"
)

message("\n=== [4/8] NB Besagproper — WITH lastYear ===")
cv_nb_besag_with <- fit_weighted_inla_model(
  data             = data_cv,
  adjacency_matrix = adjacency_matrix,
  spatial_term     = "besagproper",
  fixed_effects    = fixed_effects_with_last,
  iid_effects      = "roadSystem",
  family           = "nbinomial"
)

message("\n=== [5/8] Gaussian Weighted — NO lastYear ===")
cv_g_weighted_no <- fit_weighted_inla_model(
  data                         = data_cv,
  adjacency_matrix             = adjacency_matrix,
  spatial_term                 = "besagproper_rbf",
  fixed_effects                = fixed_effects_no_last,
  iid_effects                  = "roadSystem",
  similarity_covariates        = similarity_covariates_no_last,
  similarity_interaction_pairs = similarity_interaction_pairs,
  family                       = "gaussian"
)

message("\n=== [6/8] Gaussian Besagproper — NO lastYear ===")
cv_g_besag_no <- fit_weighted_inla_model(
  data             = data_cv,
  adjacency_matrix = adjacency_matrix,
  spatial_term     = "besagproper",
  fixed_effects    = fixed_effects_no_last,
  iid_effects      = "roadSystem",
  family           = "gaussian"
)

message("\n=== [7/8] NB Weighted — NO lastYear ===")
cv_nb_weighted_no <- fit_weighted_inla_model(
  data                         = data_cv,
  adjacency_matrix             = adjacency_matrix,
  spatial_term                 = "besagproper_rbf",
  fixed_effects                = fixed_effects_no_last,
  iid_effects                  = "roadSystem",
  similarity_covariates        = similarity_covariates_no_last,
  similarity_interaction_pairs = similarity_interaction_pairs,
  family                       = "nbinomial"
)

message("\n=== [8/8] NB Besagproper — NO lastYear ===")
cv_nb_besag_no <- fit_weighted_inla_model(
  data             = data_cv,
  adjacency_matrix = adjacency_matrix,
  spatial_term     = "besagproper",
  fixed_effects    = fixed_effects_no_last,
  iid_effects      = "roadSystem",
  family           = "nbinomial"
)

# ── Metrics ───────────────────────────────────────────────────────────────────
extract <- function(model) model$predictions$inla_pred[test_idx]

pred_g_weighted_with   <- extract(cv_g_weighted_with)
pred_g_besag_with <- extract(cv_g_besag_with)
pred_nb_weighted_with  <- extract(cv_nb_weighted_with)
pred_nb_besag_with<- extract(cv_nb_besag_with)
pred_g_weighted_no     <- extract(cv_g_weighted_no)
pred_g_besag_no   <- extract(cv_g_besag_no)
pred_nb_weighted_no    <- extract(cv_nb_weighted_no)
pred_nb_besag_no  <- extract(cv_nb_besag_no)

metrics <- data.frame(
  lastYear     = rep(c("with", "without"), each = 4),
  Spatial      = rep(c("Weighted", "Besagproper", "Weighted", "Besagproper"), 2),
  Family       = rep(c("Gaussian", "Gaussian", "NB", "NB"), 2),
  RMSE         = c(rmse(pred_g_weighted_with,    true_aadt), rmse(pred_g_besag_with, true_aadt),
                   rmse(pred_nb_weighted_with,   true_aadt), rmse(pred_nb_besag_with,true_aadt),
                   rmse(pred_g_weighted_no,      true_aadt), rmse(pred_g_besag_no,   true_aadt),
                   rmse(pred_nb_weighted_no,     true_aadt), rmse(pred_nb_besag_no,  true_aadt)),
  MAE          = c(mae(pred_g_weighted_with,     true_aadt), mae(pred_g_besag_with,  true_aadt),
                   mae(pred_nb_weighted_with,    true_aadt), mae(pred_nb_besag_with, true_aadt),
                   mae(pred_g_weighted_no,       true_aadt), mae(pred_g_besag_no,    true_aadt),
                   mae(pred_nb_weighted_no,      true_aadt), mae(pred_nb_besag_no,   true_aadt)),
  MAPE         = c(mape(pred_g_weighted_with,    true_aadt), mape(pred_g_besag_with, true_aadt),
                   mape(pred_nb_weighted_with,   true_aadt), mape(pred_nb_besag_with,true_aadt),
                   mape(pred_g_weighted_no,      true_aadt), mape(pred_g_besag_no,   true_aadt),
                   mape(pred_nb_weighted_no,     true_aadt), mape(pred_nb_besag_no,  true_aadt))
)

cat("\n===============================================================\n")
cat("  CV results — WITH vs WITHOUT lastYearAadt_logAadt\n")
cat("===============================================================\n")
print(metrics, digits = 4, row.names = FALSE)
cat("===============================================================\n")
