# install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("graph", "Rgraphviz"), dep=TRUE)
library(INLA)
INLA:::inla.binary.install(os = c("Ubuntu-22.04"))

#install.packages("pak")
pak::pak("trafikkdata/xdtkit")
library(xdtkit)
library(INLA)
set.seed(2)
year <- 2025
#Load Buskerud traffic links
preprocessed_traffic_links <- preprocess_traffic_links(buskerud_directed_traffic_links, year = year)
#Join bus-count data, for traffic links with only bus traffic
bus_aadt <- calculate_bus_aadt(stops_on_traffic_links, bus_counts, year = year)

#Notice that many of the columns in the directed traffic link data set are missing values
missing_counts <- colSums(is.na(preprocessed_traffic_links))
missing_counts[missing_counts > 0]
# functionClass                    maxLanes                    minLanes hasOnlyPublicTransportLanes 
# 18                           1                           1                           1 
# lastYearAadt_aadt     lastYearAadt_heavyRatio      lastYearAadt_heavyAadt                        aadt 
# 1                           1                           1                        1129 
# coverage                  heavyRatio                   heavyAadt       traffic_volume_source 
# 1129                        1460                        1460                        1129 
# traffic_volume_year 
# 1129 

#INLA does not like missing values, we deal with this in multiple ways:

prepared_traffic_links <- fill_missing_values(
  df = preprocessed_traffic_links,
  unknown_impute_columns = c("functionClass", "highestSpeedLimit", "lowestSpeedLimit", "maxLanes", "minLanes"),
  mode_impute_columns = c("hasOnlyPublicTransportLanes"),
  median_impute_columns = c("lastYearAadt_aadt", "lastYearAadt_heavyRatio",
                            "lastYearAadt_heavyAadt")) |>
  remove_negative_aadt() |>
  add_logLastYear() |>
  join_bus_to_traffic(bus_aadt) 
prepared_traffic_links_emma <- readRDS("prepared_traffic_links_emma.rds")
# prepared_traffic_links$municipality <- prepared_traffic_links_emma$municipality
# prepared_traffic_links$lastYearAadt_logAadt <- prepared_traffic_links_emma$lastYearAadt_logAadt
prepared_traffic_links <- subset(prepared_traffic_links_emma, select = -c(spatial.idx,inla_pred))
# prepared_traffic_links <- prepared_traffic_links_emma
#Check to see if this worked
missing_counts <- colSums(is.na(prepared_traffic_links))
missing_counts[missing_counts > 0]
# aadt              coverage            heavyRatio             heavyAadt 
# 1129                  1129                  1460                  1460 
# traffic_volume_source   traffic_volume_year          stopPointRef         stopCertainty 
# 1129                  1129                  1774                  1774 

#Now only non-covariate columns are missing values

#Now we build the adjacency matrix
adjacency_matrix <- build_adjacency_matrix(
  prepared_traffic_links,
  exclude_public_transport = TRUE)
#Output

# Building adjacency matrix for 1774 links...
# Finding adjacent links...
# Building sparse matrix from 15018 adjacency pairs...
# Excluding 1 public transport links...
# Adjacency matrix complete: 13326 non-zero entries

#Clustering - sub-groups of traffic links that are separated by traffic links that have registrations points

clusters <- strategic_network_clustering(
  data = prepared_traffic_links,
  year = year, 
  boundary_links = c("Trafikkdata_continuous")
)
#Output
# Joining with `by = join_by(parentTrafficLinkId)`
# Identifying mainland and island components...
# 
# Creating base clusters on mainland...
# 
# Merging small clusters...
# 
# Assigning barrier links to neighboring mainland clusters...
# 
# Assigning island components...
# 
# 
# === Clustering Summary ===
#   
#   Network Overview:
#   Total links: 955
# Mainland: 954 links
# Islands: 1 components, 1 links
# 
# Clustering Results:
#   Initial mainland clusters: 1
# After merging: 1 mainland clusters
# Total final clusters: 2 ( 1 mainland + 1 islands + 0 singletons)
# 
# Boundary Handling:
#   Duplicate assignments (boundaries): 0 
# 
# Cluster Size Distribution:
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.0   239.2   477.5   477.5   715.8   954.0

#Now we load the traffic nodes for Buskerud, and identify which nodes are balanceable. 
nodes <- identify_unbalanceable_nodes(buskerud_nodes, prepared_traffic_links)
#Joining with `by = join_by(id, roadSystems)`

#INLA model
# NOTE: The full interaction formula requires all-Norway data where every
# functionalRoadClass × roadCategory / maxLanes / minLanes cell is populated.
# On the Buskerud subset, 56 of 86 design-matrix columns are zero/collinear,
# causing a degenerate INLA fit (Hessian saddle point, VB corrections aborted).
# Use main effects only when running on a regional subset.
covariates <- ~ functionalRoadClass:maxLanes +
  functionalRoadClass:roadCategory +
  minLanes:roadCategory + functionalRoadClass +
  maxLanes + roadCategory +
  hasOnlyPublicTransportLanes +
  functionalRoadClass*isRamp +
  lastYearAadt_logAadt
# covariates <- ~ maxLanes + hasOnlyPublicTransportLanes + roadCategory +minLanes:roadCategory +
#     lastYearAadt_logAadt

inla_model_total <- fit_inla_model(
  data = prepared_traffic_links,
  adjacency_matrix,
  fixed_effects = covariates,
  iid_effects = "roadSystem",
  family = "poisson")
#Join predictions and uncertainty with traffic link data
predictions_total <- dplyr::full_join(prepared_traffic_links, inla_model_total$predictions)

#Diagnostic: check INLA predictions correlate with measured AADT before balancing
measured <- predictions_total[!is.na(predictions_total$aadt), c("aadt", "inla_pred", "inla_sd", "lastYearAadt_aadt")]
cat("INLA vs measured AADT correlation:", cor(measured$aadt, measured$inla_pred), "\n")
plot(measured$lastYearAadt_aadt, measured$inla_pred,
     xlab = "Measured AADT", ylab = "INLA predicted AADT",
     main = "INLA predictions vs last year (should lie near diagonal)")
abline(0, 1, col = "red")


#Balancing model
balanced_model_total <- balance_predictions(data = predictions_total,
                                            nodes = nodes,
                                            balancing_grouping_variable = clusters,
                                            nodes_to_balance = "complete nodes", 
                                            year = year)
predictions_total <- dplyr::full_join(predictions_total, balanced_model_total$balanced_res)
#Adding coefficient of variation to the data frame for plotting
predictions_total$rel_cv <- predictions_total$balanced_sd/predictions_total$balanced_pred
predictions_total$logrel_cv <- log(predictions_total$rel_cv)
plot_traffic_links_map(predictions_total, color_by = "inla_pred")
library(dplyr)
predictions_total[,c("id", "inla_pred", "inla_sd", "balanced_pred", "balanced_sd", "rel_cv")] %>%
  arrange(rel_cv) %>%
  slice(1:20)
dim(predictions_total)
inla_model_total$model_summary

print(predictions_total[!is.na(predictions_total$aadt),c("aadt","inla_pred","balanced_pred", "inla_sd", "balanced_sd")], n = 30)

print(head(predictions_total[order(predictions_total$aadt, decreasing = TRUE), c("id", "aadt", "inla_pred", "inla_sd")], 20))
# marcus_rds <- readRDS("prepared_traffic_links.rds")
# emma_rds <- readRDS("prepared_traffic_links_emma.rds")
# 
# rds_comparison_dataframe <- data.frame("marcus_ids" = marcus_rds$id, "marcus_muni" = as.integer(marcus_rds$municipality), "emma_ids" = emma_rds$id, "emma_muni" = as.integer(emma_rds$municipality))
# rds_comparison_dataframe[1:10,]
