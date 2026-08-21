#install.packages("pak")
#install.packages("INLA", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)
#inla.binary.install(os = c("Ubuntu-22.04"))
#pak::pak("trafikkdata/xdtkit")
source("R/scripts_inla_sensor_placement.R")
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

covariates <- ~ functionalRoadClass:maxLanes +
  functionalRoadClass:roadCategory +
  minLanes:roadCategory + functionalRoadClass +
  maxLanes + roadCategory +
  hasOnlyPublicTransportLanes +
  functionalRoadClass * isRamp +
  lastYearAadt_logAadt


similarity_covariates = c(
      "roadCategory", "functionalRoadClass",
      "lastYearAadt_logAadt"
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

inla_model_weighted_1 <- fit_weighted_inla_model(
  data                         = data_cv,
  adjacency_matrix             = adjacency_matrix,
  fixed_effects                = covariates,
  similarity_covariates        = similarity_covariates,
  iid_effects                  = "roadSystem",
  family                       = "nbinomial",
  distance_type = "gower",
  weight_fn = "laplacian",
  ordinal_levels = ordinal_levels_road
)
results_dataframe <- as.data.frame(dplyr::full_join(data_cv, inla_model_weighted_1$predictions))
results_dataframe[test_idx, "aadt"] <- true_aadt

plot_inla_model(inla_model_weighted_1, results_dataframe$aadt, test_idx = test_idx, title = "INLA weighted model - true vs predicted AADT",
                type = c("pred_vs_obs"))
# plot_inla_model(inla_model_weighted_1, results_dataframe$aadt, test_idx, title = "QQ-plot weighted model",
#                 type = c("qq"))
# plot_inla_model(inla_model_weighted_1, results_dataframe$aadt, title = "Residuals versus fitted",
#                 type = c("residuals_vs_fitted"))

#Gower distance metric outperforms "Gower squared" distance metric. Laplacian and Gaussian kernels er almost equal, but both are better than "linear"
laplacian_kfold <- kfold_cv_inla(data = prepared_traffic_links, adjacency_matrix = adjacency_matrix,
                                 k = 10, seed = 420,fixed_effects                = covariates,
                                 similarity_covariates        = similarity_covariates,
                                 iid_effects                  = "roadSystem",
                                 family                       = "nbinomial",
                                 distance_type = "gower",
                                 weight_fn = "laplacian",
                                 ordinal_levels = ordinal_levels_road
                                 )

gaussian_kfold <- kfold_cv_inla(
  data                         = prepared_traffic_links,
  adjacency_matrix             = adjacency_matrix,
  fixed_effects                = covariates,
  k = 10,
  seed = 420,
  similarity_covariates        = similarity_covariates,
  iid_effects                  = "roadSystem",
  family                       = "nbinomial",
  distance_type = "gower",
  weight_fn = "gaussian",
  ordinal_levels = ordinal_levels_road
)

linear_kfold <- kfold_cv_inla(
  data                         = prepared_traffic_links,
  adjacency_matrix             = adjacency_matrix,
  fixed_effects                = covariates,
  k = 10,
  seed = 420,
  similarity_covariates        = similarity_covariates,
  iid_effects                  = "roadSystem",
  family                       = "nbinomial",
  distance_type = "gower",
  weight_fn = "linear",
  ordinal_levels = ordinal_levels_road
)
plot_kfold_cv(linear_kfold, type = c("pred_vs_obs"))
#after running tests, both Gaussian and Laplacian kernel work similarly, and are equally good.

#Testing for the optimal subset of similarity covariates

gaussian_kfold_similarity_aadt <- kfold_cv_inla(
  data                         = prepared_traffic_links,
  adjacency_matrix             = adjacency_matrix,
  fixed_effects                = covariates,
  k = 10,
  seed = 420,
  similarity_covariates        = c("minLanes", "isRamp"),
  iid_effects                  = "roadSystem",
  family                       = "nbinomial",
  distance_type = "gower",
  weight_fn = "gaussian",
  ordinal_levels = NULL
)
plot_kfold_cv(gaussian_kfold_similarity_aadt, type = c("pred_vs_obs"))

laplacian_kfold_similarity_minlanes_isramp <- kfold_cv_inla(
  data                         = prepared_traffic_links,
  adjacency_matrix             = adjacency_matrix,
  fixed_effects                = covariates,
  k = 10,
  seed = 420,
  similarity_covariates        = c("minLanes", "isRamp"),
  iid_effects                  = "roadSystem",
  family                       = "nbinomial",
  distance_type = "gower",
  weight_fn = "laplacian",
  ordinal_levels = NULL
)
plot_kfold_cv(laplacian_kfold_similarity_minlanes_isramp, type = c("pred_vs_obs"))


#Using function selecting best subset of similarity_covariates
sel <- select_similarity_covariates(data = prepared_traffic_links, adjacency_matrix = adjacency_matrix,
                                    weight_fn = "gaussian")
sel
