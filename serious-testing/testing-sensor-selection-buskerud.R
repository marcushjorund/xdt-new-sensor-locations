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
source("serious-testing/scripts_inla_sensor_placement.R")

year <- 2025

preprocessed_traffic_links_norway <- preprocess_traffic_links(buskerud_directed_traffic_links, year = year)

bus_aadt<- calculate_bus_aadt(stops_on_traffic_links, bus_counts, year = year)

prepared_traffic_links <- fill_missing_values(
  df = preprocessed_traffic_links_norway,
  unknown_impute_columns = c("functionClass", "highestSpeedLimit", "lowestSpeedLimit","maxLanes", "minLanes"),
  mode_impute_columns = c("hasOnlyPublicTransportLanes"),
  median_impute_columns = c("lastYearAadt_aadt", "lastYearAadt_heavyRatio", 
                            "lastYearAadt_heavyAadt")) |>
  remove_negative_aadt() |> 
  add_logLastYear() |>
  join_bus_to_traffic(bus_aadt)

adjacency_matrix <- build_adjacency_matrix(
  prepared_traffic_links,
  exclude_public_transport = TRUE)

clusters <- strategic_network_clustering(
  data = prepared_traffic_links,
  year = year,
  boundary_links = c("Trafikkdata_continuous")
)

nodes <- identify_unbalanceable_nodes(buskerud_nodes, prepared_traffic_links)
covariates <- ~ functionalRoadClass:maxLanes +
  functionalRoadClass:roadCategory +
  minLanes:roadCategory + functionalRoadClass +
  maxLanes + roadCategory +
  hasOnlyPublicTransportLanes +
  functionalRoadClass*isRamp +
  lastYearAadt_logAadt

#Defining the ordinal levels and their ordering
ordinal_levels_road <- list(
  functionClass       = c("unknown", "E", "D", "C", "B", "A"),
  highestSpeedLimit   = c("unknown","20", "30", "40", "50", "60", "70", "80", "90", "100", "110")
)
#Best similarity covariates
similarity_covariates <- c(
  "minLanes", "highestSpeedLimit", "functionClass",
  "lastYearAadt_logAadt")

inla_rbf_buskerud <- fit_inla_rbf_model(
  data                  = prepared_traffic_links,
  adjacency_matrix      = adjacency_matrix,
  spatial_term          = "besagproper_rbf",
  fixed_effects         = covariates,
  iid_effects           = "roadSystem",
  ordinal_levels        = ordinal_levels_road,
  similarity_covariates = similarity_covariates,
  family                = "nbinomial")

predictions_total_buskerud <- dplyr::full_join(prepared_traffic_links, inla_rbf_buskerud$predictions)

spatial_hyperpar <- inla_rbf_buskerud$spatial_hyperpar[2:4,"mean"]
tau <- spatial_hyperpar[1]
d   <- spatial_hyperpar[2]
sigma <- spatial_hyperpar[3]

distances <- inla_rbf_buskerud$distances

prec_and_cov <- create_covariance_and_precision_matrix(adjacency_matrix, tau = tau, d = d, sigma = sigma, 
                                                       distances = distances, data = predictions_total_buskerud,
                                                       with_and_against = TRUE)

#cov_id_sum <- create_covariance_and_precision_matrix_both_directions(data = prepared_traffic_links_norway, covariance_matrix = prec_and_cov$covariance_matrix)

new_sensors_buskerud <- greedy_mi_sensor_selection(data = predictions_total_buskerud, covariance_matrix = prec_and_cov$covariance_matrix_sum, ids = prec_and_cov$ids,
                                                   k = 10, adjacency_matrix = adjacency_matrix, weighting_bias = 1/seq_len(8))
plot_traffic_links_simple_map(new_sensors_buskerud$selected_data_entries, color_by = "selected")


balanced_model_total_buskerud <- balance_predictions(data = predictions_total_buskerud,
                                                    nodes = nodes,
                                                    balancing_grouping_variable = clusters,
                                                    nodes_to_balance = "complete nodes", 
                                                    year = year)
predictions_total_buskerud <- dplyr::full_join(predictions_total_buskerud, balanced_model_total_buskerud$balanced_res)
predictions_total_buskerud$relative_uncertainty <- predictions_total_buskerud$balanced_sd/predictions_total_buskerud$balanced_pred
predictions_total_buskerud$log_relative_uncertainty <- log(predictions_total_buskerud$relative_uncertainty)
plot_traffic_links_map(predictions_total_buskerud[predictions_total_buskerud$relative_uncertainty<1 & predictions_total_buskerud$relative_uncertainty> 0.001,], color_by = "log_relative_uncertainty")

#Plotting the "real" estimates using the old model
inla_model_total <- fit_inla_model(
  data = prepared_traffic_links,
  adjacency_matrix,
  fixed_effects = covariates,
  iid_effects = "roadSystem",
  family = "nbinomial")
#Join predictions and uncertainty with traffic link data
predictions_total <- dplyr::full_join(prepared_traffic_links, inla_model_total$predictions)
balanced_model_total <- balance_predictions(data = predictions_total,
                                            nodes = nodes,
                                            balancing_grouping_variable = clusters,
                                            nodes_to_balance = "complete nodes", 
                                            year = year)
predictions_total <- dplyr::full_join(predictions_total, balanced_model_total$balanced_res)
predictions_total$relative_uncertainty <- predictions_total$balanced_sd/predictions_total$balanced_pred
predictions_total$log_relative_uncertainty <- log(predictions_total$relative_uncertainty)
plot_traffic_links_map(predictions_total[predictions_total$relative_uncertainty<1 & predictions_total$relative_uncertainty >0.001,], color_by = "log_relative_uncertainty")
