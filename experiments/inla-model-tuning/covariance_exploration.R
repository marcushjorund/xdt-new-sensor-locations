# install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("graph", "Rgraphviz"), dep=TRUE)
# INLA:::inla.binary.install(os = c("Ubuntu-22.04"))
# 
# #install.packages("pak")
# pak::pak("trafikkdata/xdtkit")
library(xdtkit)
library(INLA)

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
# covariates <- ~ functionalRoadClass:maxLanes +
#   functionalRoadClass:roadCategory +
#   minLanes:roadCategory + functionalRoadClass +
#   maxLanes + roadCategory +
#   hasOnlyPublicTransportLanes +
#   functionalRoadClass*isRamp +
#   lastYearAadt_logAadt
covariates <- ~ maxLanes + hasOnlyPublicTransportLanes + roadCategory +minLanes:roadCategory +
    lastYearAadt_logAadt

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
     main = "INLA predictions vs last year (should lie near diag)")
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
max(predictions_total$rel_cv, na.rm = TRUE)
plot_traffic_links_map(predictions_total[predictions_total$rel_cv < 0.35,], color_by = "rel_cv")
library(dplyr)
predictions_total[,c("id", "inla_pred", "inla_sd", "balanced_pred", "balanced_sd", "rel_cv")] %>%
  arrange(rel_cv) %>%
  slice(1:20)
dim(predictions_total)
inla_model_total$model_summary


#Creating precision matrix as defined in besag proper
#extracting diag and precision parameter from inla_model_total
# Besag proper model has two hyperparameters: the precision parameter (tau) and the diag parameter (d). The precision matrix is defined as Q = tau * (D - A), where D is a diag matrix with D[i,i] = d + sum(A[i,]) and A is the adjacency matrix.
# In INLA, the hyperparameters are typically ordered with the precision parameter first, followed by the diag parameter. Therefore, we can extract them as follows:

d <- inla_model_total$model_summary$hyperpar$mean[2]
tau <-  inla_model_total$model_summary$hyperpar$mean[1]
# 
precision_matrix <- matrix(0, ncol = dim(adjacency_matrix)[1], nrow = dim(adjacency_matrix)[1])
for (i in 1:dim(adjacency_matrix)[1]) {
  diag(precision_matrix)[i] <- tau*(d + sum(adjacency_matrix[i,]))
  non_zero_idx <- which(adjacency_matrix[i,] != 0)
  precision_matrix[i, non_zero_idx] <- -tau
}
covariance_matrix <- solve(precision_matrix)
# correlation_matrix <- cov2cor(covariance_matrix)
# correlation_matrix[1:5, 1:5]
# correlation_matrix[1,]
# # Build hop-distance matrix via BFS from every node.
# # igraph::distances() with algorithm = "unweighted" runs BFS in C — fast for 1774 nodes.
# # distance_matrix[i,j] = number of links to traverse from i to j.
# # diag = 0. Disconnected pairs = Inf.
# library(igraph)
# g <- igraph::graph_from_adjacency_matrix(
#   adjacency_matrix,
#   mode     = "undirected",
#   weighted = NULL,
#   diag     = FALSE)
# 
# distance_matrix <- igraph::distances(g, algorithm = "unweighted")
# distance_matrix
# rownames(distance_matrix) <- rownames(adjacency_matrix)
# colnames(distance_matrix) <- colnames(adjacency_matrix)
# 
# cat("Max finite distance (diameter):", max(distance_matrix[is.finite(distance_matrix)]), "\n")
# cat("Disconnected pairs (Inf):      ", sum(!is.finite(distance_matrix)), "\n")
# cat("\nHop-distance distribution (off-diag pairs):\n")
# print(table(distance_matrix[distance_matrix > 0 & is.finite(distance_matrix)]))
# 
# #Make histogram of hop distances
# hist(distance_matrix[distance_matrix > 0 & is.finite(distance_matrix)],
#      breaks = seq(0.5, max(distance_matrix[is.finite(distance_matrix)]) + 0.5, by = 1),
#      main = "Histogram of hop distances",
#      xlab = "Hop distance", ylab = "Frequency")
# 
# # Scatter: correlation vs hop distance (sample up to 50k pairs to keep it readable)
# set.seed(42)
# all_pairs <- which(distance_matrix > 0 & is.finite(distance_matrix), arr.ind = TRUE)
# if (nrow(all_pairs) > 50000) all_pairs <- all_pairs[sample(nrow(all_pairs), 50000), ]
# hop_dist <- distance_matrix[all_pairs]
# corr_val  <- correlation_matrix[all_pairs]
# plot(jitter(hop_dist, amount = 0.2), corr_val,
#      xlab = "Hop distance", ylab = "Besag correlation",
#      main = "Correlation vs graph distance",
#      pch = ".", col = rgb(0, 0, 0, 0.15))
# abline(h = 0, lty = 2, col = "red")
# 
# # Mean correlation at each hop distance
# mean_corr_by_dist <- tapply(corr_val, hop_dist, mean)
# cat("\nMean Besag correlation by hop distance:\n")
# print(round(mean_corr_by_dist, 6))
# 
# 
# 1/diag(precision_matrix)[1:10]
# diag(covariance_matrix)[1:10]
# for (i in 1:10){
#   print(sum(adjacency_matrix[i,]))
# }
# 
# covariance_matrix[1,1] - t(covariance_matrix[1,2:dim(covariance_matrix)[1]])%*%
#   solve(covariance_matrix[2:dim(covariance_matrix)[1], 2:dim(covariance_matrix)[1]])%*%covariance_matrix[2:dim(covariance_matrix)[1],1]
# 1/precision_matrix[1,1]
# 
# #Trying the algorithm
# unmeasured_idxs <- which(is.na(predictions_total$aadt))
# measured_idxs <- which(!is.na(predictions_total$aadt))
# org_unmeasured_idxs <- unmeasured_idxs
# org_measured_idxs <- measured_idxs
# 
# k <- 10
# new_sensor_location_idxs <- rep(0, k)
# for (i in 1:k){
#   delta <- rep(0, length(unmeasured_idxs)) 
#   j <- 1
#   for (y in unmeasured_idxs){
#     temp_unmeasured_idxs <- unmeasured_idxs[which(unmeasured_idxs != y)]
#     delta[j] <-(covariance_matrix[y,y]- t(covariance_matrix[y,measured_idxs])%*%solve(covariance_matrix[measured_idxs, measured_idxs])%*%covariance_matrix[measured_idxs,y])/
#       (covariance_matrix[y,y] - t(covariance_matrix[y,temp_unmeasured_idxs])%*%solve(covariance_matrix[temp_unmeasured_idxs, temp_unmeasured_idxs])%*%covariance_matrix[temp_unmeasured_idxs,y])
#     j <- j + 1
#     }
#   best <- which(delta == max(delta))[1]
#   print(unmeasured_idxs[best])
#   measured_idxs <- c(measured_idxs,unmeasured_idxs[best])
#   print(length(measured_idxs))
#   unmeasured_idxs <- unmeasured_idxs[which(unmeasured_idxs != unmeasured_idxs[best])]
#   print(length(unmeasured_idxs))
#   cat("Iteration", i, ": added index", unmeasured_idxs[best], "with delta =", delta[best], "\n")
#   new_sensor_location_idxs[i] <- unmeasured_idxs[best]
# }
# new_sensor_ids <- predictions_total$id[new_sensor_location_idxs]
# new_sensor_ids
# 
# #In our data we have traffic links with id's "WITH" and "AGAINST"
# #For a given sensor, we measure the aadt on both the with and against traffic links. Therefore we need to account for this in our sensor selection model
# #We propose the linear combination traffic_link_i = traffic_link_with_i + traffic_link_against_i, where traffic_link_with_i and traffic_link_against_i are the two links corresponding to the with and against directions of the same road segment. We can then modify our sensor selection algorithm to select pairs of links together, ensuring that we capture the total traffic on that road segment. 
# #Now the distribution of this linear combination, is log gaussian, and we derive that the distribution has a mean of mu_with + mu_against and covariance matrix
# prepared_traffic_links$id[1:10]
# dim(prepared_traffic_links)
# #traffic link ids are strings, with the form "0.0@3134734-1.0@3233685-WITH" and "0.0@3134734-1.0@3233685-AGAINST"
# #we therefore need to split the id's to get the corresponding with and against links and sort them such that we add the correct with and against links
with_idx_int    <- which(grepl("-WITH$",    prepared_traffic_links$id))
against_idx_int <- which(grepl("-AGAINST$", prepared_traffic_links$id))
length(with_idx_int)
length(against_idx_int)

# Strip the suffix to get base IDs
base_with    <- sub("-WITH$",    "", prepared_traffic_links$id[with_idx_int])
base_against <- sub("-AGAINST$", "", prepared_traffic_links$id[against_idx_int])

# Find paired base IDs
paired_bases <- intersect(base_with, base_against)
unpaired_bases <- setdiff(base_with, base_against)

cat("Paired links:  ", length(paired_bases), "\n")   # expect 820
cat("Unpaired WITH: ", length(unpaired_bases), "\n")  # expect 134

# Get the original row indices in prepared_traffic_links for each group
paired_with_idx    <- match(paste0(paired_bases,    "-WITH"),    prepared_traffic_links$id)
paired_against_idx <- match(paste0(paired_bases,    "-AGAINST"), prepared_traffic_links$id)
unpaired_with_idx  <- match(paste0(unpaired_bases,  "-WITH"),    prepared_traffic_links$id)

library(Matrix)
n <- nrow(prepared_traffic_links)
m <- length(paired_bases) + length(unpaired_bases)  # = 954

# Build sparse C matrix
row_idx <- c(
  seq_along(paired_bases),                          # rows 1:820 for WITH
  seq_along(paired_bases),                          # rows 1:820 for AGAINST
  length(paired_bases) + seq_along(unpaired_bases)  # rows 821:954 for unpaired
)
col_idx <- c(
  paired_with_idx,
  paired_against_idx,
  unpaired_with_idx
)
values <- rep(1, length(row_idx))

C <- sparseMatrix(i = row_idx, j = col_idx, x = values,
                  dims = c(m, n))
dim(C)
#adding with + against by linear transformation of covariance_matrix by C
covariance_matrix_sum <- C %*% covariance_matrix %*% t(C)
#extracting corresponding indices for measured/unmeasured traffic link in the new covariance matrix

#Trying the algorithm with WITH + AGAINST
measured_paired_idxs <- which(!is.na(prepared_traffic_links[paired_with_idx,]$aadt))
unmeasured_paired_idxs <- which(is.na(prepared_traffic_links[paired_with_idx,]$aadt))
measured_unpaired_idxs <- which(!is.na(prepared_traffic_links[unpaired_with_idx,]$aadt))
unmeasured_unpaired_idxs <- which(is.na(prepared_traffic_links[unpaired_with_idx,]$aadt))

measured_idxs_sum <- c(measured_paired_idxs, measured_unpaired_idxs+length(paired_bases))
unmeasured_idxs_sum <- c(unmeasured_paired_idxs, unmeasured_unpaired_idxs+length(paired_bases))

org_unmeasured_idxs_sum <- unmeasured_idxs_sum
org_measured_idxs_sum <- measured_idxs_sum

k <- 10
new_sensor_location_idxs <- rep(0, k)
for (i in 1:k){
  delta <- rep(0, length(unmeasured_idxs_sum))
  j <- 1
  for (y in unmeasured_idxs_sum){
    temp_unmeasured_idxs_sum <- unmeasured_idxs_sum[which(unmeasured_idxs_sum != y)]
    delta[j] <-(covariance_matrix_sum[y,y]- t(covariance_matrix_sum[y,measured_idxs_sum])%*%solve(covariance_matrix_sum[measured_idxs_sum, measured_idxs_sum])%*%covariance_matrix_sum[measured_idxs_sum,y])/
      (covariance_matrix_sum[y,y] - t(covariance_matrix_sum[y,temp_unmeasured_idxs_sum])%*%solve(covariance_matrix_sum[temp_unmeasured_idxs_sum, temp_unmeasured_idxs_sum])%*%covariance_matrix_sum[temp_unmeasured_idxs_sum,y])
    j <- j + 1
  }
  best <- which(delta == max(delta))
  cat("Iteration", i, ": added index", unmeasured_idxs_sum[best], "with delta =", delta[best], "\n")
  new_sensor_location_idxs[i] <- unmeasured_idxs_sum[best]
  measured_idxs_sum <- c(measured_idxs_sum,unmeasured_idxs_sum[best])
  print(length(measured_idxs_sum))
  unmeasured_idxs_sum <- unmeasured_idxs_sum[which(unmeasured_idxs_sum != unmeasured_idxs_sum[best])]
  print(length(unmeasured_idxs_sum))
}
new_sensor_location_idxs
paired_new_sensor_locations <- new_sensor_location_idxs[new_sensor_location_idxs <= length(paired_bases)]
unpaired_new_sensor_locations <-new_sensor_location_idxs[new_sensor_location_idxs > length(paired_bases)]
length(paired_bases)
paired_new_sensor_ids <- paired_bases[paired_new_sensor_locations]
unpaired_new_sensor_ids <- unpaired_bases[unpaired_new_sensor_locations]
paired_new_sensor_ids
unpaired_new_sensor_ids

# ── Build normalised feature matrix (one row per traffic link) ────────────────
# Columns: standardised numerics + binaries + one-hot categoricals + interactions.
# Each row is L2-normalised so that row_i · row_j = cosine similarity, and
# ||row_i - row_j||² = squared Euclidean distance in the normalised feature space.
build_normalised_feature_matrix <- function(
  df,
  numeric_features     = c("maxLanes", "minLanes",
                            "highestSpeedLimit", "lowestSpeedLimit",
                            "lastYearAadt_logAadt"),
  binary_features      = c("hasOnlyPublicTransportLanes", "isRamp"),
  categorical_features = c("roadCategory", "functionClass", "roadSystem"),
  interaction_pairs    = list(
    c("functionClass",        "roadCategory"),
    c("functionClass",        "isRamp"),
    c("maxLanes",             "roadCategory"),
    c("lastYearAadt_logAadt", "roadCategory")
  )) {

  one_hot <- function(x, name, levels) {
    mat <- outer(as.character(x), levels, `==`) * 1L
    colnames(mat) <- paste0(name, "_", levels)
    mat
  }

  cat_levels <- lapply(
    intersect(categorical_features, names(df)),
    function(col) sort(unique(as.character(df[[col]])))
  )
  names(cat_levels) <- intersect(categorical_features, names(df))

  parts <- list()

  for (col in intersect(numeric_features, names(df)))
    parts[[col]] <- as.numeric(df[[col]])

  for (col in intersect(binary_features, names(df)))
    parts[[col]] <- as.numeric(df[[col]])

  for (col in intersect(categorical_features, names(df))) {
    oh <- one_hot(df[[col]], col, cat_levels[[col]])
    for (k in seq_len(ncol(oh))) parts[[colnames(oh)[k]]] <- oh[, k]
  }

  for (pair in interaction_pairs) {
    a <- pair[1]; b <- pair[2]
    if (!a %in% names(df) || !b %in% names(df)) next
    a_is_cat <- a %in% categorical_features
    b_is_cat <- b %in% categorical_features
    a_is_num <- a %in% c(numeric_features, binary_features)
    b_is_num <- b %in% c(numeric_features, binary_features)

    if (a_is_cat && b_is_cat) {
      cross <- paste0(as.character(df[[a]]), ":", as.character(df[[b]]))
      lvls  <- sort(unique(cross))
      oh    <- one_hot(cross, paste0(a, "x", b), lvls)
      for (k in seq_len(ncol(oh))) parts[[colnames(oh)[k]]] <- oh[, k]
    } else if (a_is_cat && b_is_num) {
      oh       <- one_hot(df[[a]], a, cat_levels[[a]])
      num_vals <- as.numeric(df[[b]])
      for (k in seq_len(ncol(oh)))
        parts[[paste0(colnames(oh)[k], "x", b)]] <- oh[, k] * num_vals
    } else if (a_is_num && b_is_cat) {
      oh       <- one_hot(df[[b]], b, cat_levels[[b]])
      num_vals <- as.numeric(df[[a]])
      for (k in seq_len(ncol(oh)))
        parts[[paste0(colnames(oh)[k], "x", a)]] <- oh[, k] * num_vals
    } else if (a_is_num && b_is_num) {
      parts[[paste0(a, "x", b)]] <- as.numeric(df[[a]]) * as.numeric(df[[b]])
    }
  }

  mat <- do.call(cbind, parts)

  # Column-standardise; leave constant columns at 0
  col_means <- colMeans(mat, na.rm = TRUE)
  col_sds   <- apply(mat, 2, sd, na.rm = TRUE)
  non_const <- col_sds > 0
  mat[, non_const]  <- sweep(mat[, non_const, drop = FALSE], 2, col_means[non_const], `-`)
  mat[, non_const]  <- sweep(mat[, non_const, drop = FALSE], 2, col_sds[non_const],   `/`)
  mat[, !non_const] <- 0
  mat[is.na(mat)]   <- 0

  # L2-normalise each row
  row_norms <- sqrt(rowSums(mat^2))
  row_norms[row_norms == 0] <- 1
  mat / row_norms
}

cosine_similarity_traffic_links <- function(df, i, j,
  numeric_features     = c("maxLanes", "minLanes",
                           "highestSpeedLimit", "lowestSpeedLimit",
                           "lastYearAadt_logAadt"),
  binary_features      = c("hasOnlyPublicTransportLanes", "isRamp"),
  categorical_features = c("roadCategory", "functionClass", "roadSystem"),
  interaction_pairs    = list(
    c("functionClass",        "roadCategory"),
    c("functionClass",        "isRamp"),
    c("maxLanes",             "roadCategory"),
    c("lastYearAadt_logAadt", "roadCategory")
  )) {

  one_hot <- function(x, name, levels) {
    mat <- outer(as.character(x), levels, `==`) * 1L
    colnames(mat) <- paste0(name, "_", levels)
    mat
  }

  # Compute levels and scaling stats from the FULL data frame (all rows)
  cat_levels <- lapply(
    intersect(categorical_features, names(df)),
    function(col) sort(unique(as.character(df[[col]])))
  )
  names(cat_levels) <- intersect(categorical_features, names(df))

  make_feature_matrix <- function(rows, scaling_stats = NULL) {
    parts <- list()

    # Numeric: standardise with full-data mean/sd
    for (col in intersect(numeric_features, names(rows))) {
      parts[[col]] <- as.numeric(rows[[col]])
    }

    # Binary: keep as 0/1, standardise with full-data mean/sd
    for (col in intersect(binary_features, names(rows))) {
      parts[[col]] <- as.numeric(rows[[col]])
    }

    # Categorical -> one-hot (levels fixed from full data)
    for (col in intersect(categorical_features, names(rows))) {
      oh <- one_hot(rows[[col]], col, cat_levels[[col]])
      for (k in seq_len(ncol(oh))) parts[[colnames(oh)[k]]] <- oh[, k]
    }

    # Interactions
    for (pair in interaction_pairs) {
      a <- pair[1]; b <- pair[2]
      if (!a %in% names(rows) || !b %in% names(rows)) next

      a_is_cat <- a %in% categorical_features
      b_is_cat <- b %in% categorical_features
      a_is_num <- a %in% c(numeric_features, binary_features)
      b_is_num <- b %in% c(numeric_features, binary_features)

      if (a_is_cat && b_is_cat) {
        cross  <- paste0(as.character(rows[[a]]), ":", as.character(rows[[b]]))
        lvls   <- sort(unique(paste0(
          as.character(df[[a]]), ":", as.character(df[[b]]))))
        oh <- one_hot(cross, paste0(a, "x", b), lvls)
        for (k in seq_len(ncol(oh))) parts[[colnames(oh)[k]]] <- oh[, k]

      } else if (a_is_cat && b_is_num) {
        oh <- one_hot(rows[[a]], a, cat_levels[[a]])
        num_vals <- as.numeric(rows[[b]])
        for (k in seq_len(ncol(oh)))
          parts[[paste0(colnames(oh)[k], "x", b)]] <- oh[, k] * num_vals

      } else if (a_is_num && b_is_cat) {
        oh <- one_hot(rows[[b]], b, cat_levels[[b]])
        num_vals <- as.numeric(rows[[a]])
        for (k in seq_len(ncol(oh)))
          parts[[paste0(colnames(oh)[k], "x", a)]] <- oh[, k] * num_vals

      } else if (a_is_num && b_is_num) {
        parts[[paste0(a, "x", b)]] <-
          as.numeric(rows[[a]]) * as.numeric(rows[[b]])
      }
    }

    do.call(cbind, parts)
  }

  # Build feature matrix for ALL rows to compute scaling stats
  full_feat <- make_feature_matrix(df)
  col_means <- colMeans(full_feat, na.rm = TRUE)
  col_sds   <- apply(full_feat, 2, sd, na.rm = TRUE)

  # Build feature matrix for only the two rows of interest
  pair_feat <- make_feature_matrix(df[c(i, j), ])

  # Drop columns where either observation is NA
  complete_cols <- apply(pair_feat, 2, function(x) all(!is.na(x)))
  pair_feat  <- pair_feat[,  complete_cols, drop = FALSE]
  col_means  <- col_means[complete_cols]
  col_sds    <- col_sds[complete_cols]

  if (ncol(pair_feat) == 0) return(NA_real_)

  # Standardise using full-data stats; skip constant columns (sd == 0)
  non_const <- col_sds > 0
  pair_feat[, non_const] <-
    sweep(pair_feat[, non_const, drop = FALSE], 2, col_means[non_const], `-`)
  pair_feat[, non_const] <-
    sweep(pair_feat[, non_const, drop = FALSE], 2, col_sds[non_const],   `/`)
  # Constant columns (e.g. a one-hot level not present in full data) -> 0
  pair_feat[, !non_const] <- 0

  v1 <- as.numeric(pair_feat[1, ])
  v2 <- as.numeric(pair_feat[2, ])

  sum(v1 * v2) / (sqrt(sum(v1^2)) * sqrt(sum(v2^2)))
}

#Doing cosine similarity between nodes 1 and 2
library(Matrix)

# ── Fix: sparse-safe upper triangle extraction ───────────────────────────────
adj_trip <- Matrix::summary(as(adjacency_matrix, "dgCMatrix"))
mask_upper <- adj_trip$i < adj_trip$j
ui <- adj_trip$i[mask_upper]
uj <- adj_trip$j[mask_upper]

feat_norm <- build_normalised_feature_matrix(prepared_traffic_links)

sim_vals <- rowSums(feat_norm[ui, ] * feat_norm[uj, ])
sim_vals  <- pmax(sim_vals, 0)   # clamp: dissimilar neighbours get weight 0

n <- nrow(prepared_traffic_links)
W <- sparseMatrix(
  i    = c(ui, uj),
  j    = c(uj, ui),
  x    = rep(sim_vals, 2),
  dims = c(n, n)
)

# Weighted Laplacian Q0 = D_W - W  (positive semi-definite)
Q0 <- diag(x = rowSums(W)) - W


# ── rgeneric: weighted Besag proper  Q(tau, d) = tau * (d*I + Q0) ────────────
# theta[1] = log(tau),  theta[2] = log(d)   (both estimated by INLA)

inla.rgeneric.weighted.besag <- function(cmd, theta) {
  # n and Q0 are injected by inla.rgeneric.define()

  graph <- function() {
    # Sparsity structure of Q — values irrelevant, only 0/nonzero matters
    Q0 + diag(n)
  }

  Q <- function() {
    tau <- exp(theta[1])
    d   <- exp(theta[2])
    tau * (d * diag(n) + Q0)
  }

  mu <- function() numeric(n)

  log.norm.const <- function() numeric(0)   # let INLA compute via Cholesky

  log.prior <- function() {
    # Log-Gamma(shape=1, rate=5e-5) prior on tau and d —
    # same vague prior INLA uses internally for besagproper.
    # log p(theta_k) = a * theta_k - b * exp(theta_k)
    # (Jacobian for theta = log(param) is already absorbed: a*theta not (a-1)*theta)
    lp <- function(theta_k, a = 1, b = 5e-5) a * theta_k - b * exp(theta_k)
    lp(theta[1]) + lp(theta[2])
  }

  initial <- function() c(0, 0)   # tau = 1, d = 1

  quit <- function() invisible()

  switch(cmd,
    graph          = graph(),
    Q              = Q(),
    mu             = mu(),
    log.norm.const = log.norm.const(),
    log.prior      = log.prior(),
    initial        = initial(),
    quit           = quit()
  )
}

weighted_besag_model <- inla.rgeneric.define(
  inla.rgeneric.weighted.besag,
  n  = n,
  Q0 = Q0
)

# ── Fit ───────────────────────────────────────────────────────────────────────
prepared_traffic_links$spatial.idx <- seq_len(n)

formula_weighted <- aadt ~
  maxLanes + hasOnlyPublicTransportLanes + roadCategory +
  minLanes:roadCategory + lastYearAadt_logAadt +
  f(spatial.idx,
    model  = weighted_besag_model,
    constr = FALSE)             # proper model is full-rank; no sum-to-zero needed

inla_model_weighted <- inla(
  formula_weighted,
  data   = prepared_traffic_links,
  family = "poisson",
  control.predictor = list(compute = TRUE, link = 1),
  control.compute   = list(dic = TRUE, waic = TRUE)
)

# Hyperparameters are stored on log scale ("Theta1", "Theta2")
hp <- inla_model_weighted$summary.hyperpar
cat("tau:", exp(hp[1, "mean"]),
    " 95% CI [", exp(hp[1, "0.025quant"]), ",", exp(hp[1, "0.975quant"]), "]\n")
cat("d:  ", exp(hp[2, "mean"]),
    " 95% CI [", exp(hp[2, "0.025quant"]), ",", exp(hp[2, "0.975quant"]), "]\n")

# Reconstruct the fitted precision matrix at posterior mean estimates
tau_hat <- exp(hp[1, "mean"])
d_hat   <- exp(hp[2, "mean"])
Q_fitted <- tau_hat * (d_hat * diag(n) + Q0)

# Invert to covariance and apply WITH+AGAINST linear combination
covariance_matrix_weighted     <- solve(Q_fitted)
covariance_matrix_weighted_sum <- C %*% covariance_matrix_weighted %*% t(C)


# ── RBF-kernel weighted Besag proper ─────────────────────────────────────────

# Precompute squared covariate distances for all neighbour pairs (fixed)
dist_sq  <- rowSums((feat_norm[ui, ] - feat_norm[uj, ])^2)
ell_init <- median(sqrt(dist_sq))   # median heuristic starting value

# Helper: given a bandwidth ell, return the per-neighbour-pair RBF similarity
# values (one entry per upper-triangle pair, matching ui/uj ordering) and the
# resulting sparse weight matrix W.
rbf_similarity_weights <- function(ell, dist_sq, ui, uj, n) {
  sim_vals <- pmax(exp(-dist_sq / (2 * ell^2)), 0)
  W <- Matrix::sparseMatrix(
    i    = c(ui, uj),
    j    = c(uj, ui),
    x    = rep(sim_vals, 2),
    dims = c(n, n)
  )
  list(sim_vals = sim_vals, W = W)
}

# rgeneric model: Q(tau, d, ell) = tau * (d*I + D_W - W)
# theta[1] = log(tau), theta[2] = log(d), theta[3] = log(ell)
inla.rgeneric.weighted.besag.RBF <- function(cmd, theta) {
  # Injected via inla.rgeneric.define(): n, ui, uj, dist_sq, ell_init

  graph <- function() {
    Matrix::sparseMatrix(
      i    = c(ui, uj, seq_len(n)),
      j    = c(uj, ui, seq_len(n)),
      x    = 1,
      dims = c(n, n)
    )
  }

  Q <- function() {
    tau      <- exp(theta[1])
    d        <- exp(theta[2])
    ell      <- exp(theta[3])
    sim_vals <- pmax(exp(-dist_sq / (2 * ell^2)), 0)
    W        <- Matrix::sparseMatrix(i = c(ui, uj), j = c(uj, ui),
                                     x = rep(sim_vals, 2), dims = c(n, n))
    Q0_local <- Matrix::Diagonal(x = Matrix::rowSums(W)) - W
    tau * (d * Matrix::Diagonal(n) + Q0_local)
  }

  mu             <- function() numeric(n)
  log.norm.const <- function() numeric(0)
  quit           <- function() invisible()

  log.prior <- function() {
    lp_gamma <- function(theta_k, a = 1, b = 5e-5) a * theta_k - b * exp(theta_k)
    lp_ell <- dnorm(theta[3], mean = log(ell_init), sd = 0.3, log = TRUE)
    lp_gamma(theta[1]) + lp_gamma(theta[2]) + lp_ell
  }

  initial <- function() c(0, 0, log(ell_init))

  switch(cmd,
    graph          = graph(),
    Q              = Q(),
    mu             = mu(),
    log.norm.const = log.norm.const(),
    log.prior      = log.prior(),
    initial        = initial(),
    quit           = quit()
  )
}

weighted_besag_model_RBF <- inla.rgeneric.define(
  inla.rgeneric.weighted.besag.RBF,
  n        = n,
  ui       = ui,
  uj       = uj,
  dist_sq  = dist_sq,
  ell_init = ell_init
)

# ── Fit ───────────────────────────────────────────────────────────────────────
prepared_traffic_links$spatial.idx <- seq_len(n)

formula_weighted_RBF <- aadt ~
  maxLanes + hasOnlyPublicTransportLanes + roadCategory +
  minLanes:roadCategory + lastYearAadt_logAadt +
  f(spatial.idx,
    model  = weighted_besag_model_RBF,
    n      = n,
    constr = FALSE)

inla_model_weighted_RBF <- inla(
  formula_weighted_RBF,
  data              = prepared_traffic_links,
  family            = "poisson",
  control.predictor = list(compute = TRUE, link = 1),
  control.compute   = list(dic = TRUE, waic = TRUE)
)

# ── Hyperparameter estimates ───────────────────────────────────────────────────
hp_RBF <- inla_model_weighted_RBF$summary.hyperpar
cat("tau:", exp(hp_RBF[1, "mean"]),
    " 95% CI [", exp(hp_RBF[1, "0.025quant"]), ",", exp(hp_RBF[1, "0.975quant"]), "]\n")
cat("d:  ", exp(hp_RBF[2, "mean"]),
    " 95% CI [", exp(hp_RBF[2, "0.025quant"]), ",", exp(hp_RBF[2, "0.975quant"]), "]\n")
cat("ell:", exp(hp_RBF[3, "mean"]),
    " 95% CI [", exp(hp_RBF[3, "0.025quant"]), ",", exp(hp_RBF[3, "0.975quant"]), "]\n")

# ── Reconstruct precision matrix at posterior mean estimates ──────────────────
tau_hat <- exp(hp_RBF[1, "mean"])
d_hat   <- exp(hp_RBF[2, "mean"])
ell_hat <- exp(hp_RBF[3, "mean"])

rw_hat       <- rbf_similarity_weights(ell_hat, dist_sq, ui, uj, n)
Q0_hat       <- Matrix::Diagonal(x = Matrix::rowSums(rw_hat$W)) - rw_hat$W
Q_fitted_RBF <- tau_hat * (d_hat * Matrix::Diagonal(n) + Q0_hat)

# Invert to covariance and apply WITH+AGAINST linear combination
covariance_matrix_RBF     <- as.matrix(Matrix::solve(Q_fitted_RBF))
covariance_matrix_RBF_sum <- C %*% covariance_matrix_RBF %*% t(C)

inla_model_weighted_RBF$summary.fitted.values

fv <- inla_model_weighted_RBF$summary.fitted.values
pathological <- which(fv$mean / fv$`0.5quant` > 100)
cat("Links with exploding mean/median ratio:", length(pathological), "\n")
cat("Their median predictions:\n")
print(fv[pathological, c("mean", "sd", "0.025quant", "0.5quant", "0.975quant")])

# Identify row indices of pathological links in the full data
path_idx <- as.integer(sub("fitted.Predictor.", "", rownames(fv)[pathological]))

# Convert path_idx (1:1774 full-link space) to summed index space (1:954)
# Rows 1:length(paired_bases)  → paired segments (keyed by WITH link position)
# Rows (length(paired_bases)+1):954 → unpaired WITH segments
path_sum_idx <- vapply(path_idx, function(idx) {
  p <- match(idx, paired_with_idx)
  if (!is.na(p)) return(p)
  p <- match(idx, paired_against_idx)   # AGAINST of a pair → same summed row
  if (!is.na(p)) return(p)
  p <- match(idx, unpaired_with_idx)
  if (!is.na(p)) return(length(paired_bases) + p)
  NA_integer_
}, integer(1))

cat("Pathological links in summed-index space:\n")
print(data.frame(
  full_idx  = path_idx,
  sum_idx   = path_sum_idx,
  road_id   = prepared_traffic_links$id[path_idx]
))

cat("\nWhich pathological links appear in new_sensor_location_idxs:\n")
hit <- match(path_sum_idx, new_sensor_location_idxs)
print(data.frame(
  sum_idx         = path_sum_idx,
  sensor_rank     = hit,
  road_id         = prepared_traffic_links$id[path_idx]
))

# ── Greedy sensor selection using RBF covariance matrix ──────────────────────

measured_paired_idxs    <- which(!is.na(prepared_traffic_links[paired_with_idx, ]$aadt))
unmeasured_paired_idxs  <- which( is.na(prepared_traffic_links[paired_with_idx, ]$aadt))
measured_unpaired_idxs  <- which(!is.na(prepared_traffic_links[unpaired_with_idx, ]$aadt))
unmeasured_unpaired_idxs <- which( is.na(prepared_traffic_links[unpaired_with_idx, ]$aadt))

measured_idxs_sum_rbf   <- c(measured_paired_idxs,
                              measured_unpaired_idxs + length(paired_bases))
unmeasured_idxs_sum_rbf <- c(unmeasured_paired_idxs,
                              unmeasured_unpaired_idxs + length(paired_bases))

k <- 10
new_sensor_location_idxs_RBF <- rep(0L, k)

for (i in seq_len(k)) {
  delta <- numeric(length(unmeasured_idxs_sum_rbf))
  for (j in seq_along(unmeasured_idxs_sum_rbf)) {
    y    <- unmeasured_idxs_sum_rbf[j]
    temp <- unmeasured_idxs_sum_rbf[-j]
    delta[j] <-
      (covariance_matrix_RBF_sum[y, y] -
         t(covariance_matrix_RBF_sum[y, measured_idxs_sum_rbf]) %*%
         solve(covariance_matrix_RBF_sum[measured_idxs_sum_rbf, measured_idxs_sum_rbf]) %*%
         covariance_matrix_RBF_sum[measured_idxs_sum_rbf, y]) /
      (covariance_matrix_RBF_sum[y, y] -
         t(covariance_matrix_RBF_sum[y, temp]) %*%
         solve(covariance_matrix_RBF_sum[temp, temp]) %*%
         covariance_matrix_RBF_sum[temp, y])
  }
  best <- which.max(delta)
  cat("Iteration", i, ": added sum_idx", unmeasured_idxs_sum_rbf[best],
      "  delta =", delta[best], "\n")
  new_sensor_location_idxs_RBF[i] <- unmeasured_idxs_sum_rbf[best]
  measured_idxs_sum_rbf   <- c(measured_idxs_sum_rbf,   unmeasured_idxs_sum_rbf[best])
  unmeasured_idxs_sum_rbf <- unmeasured_idxs_sum_rbf[-best]
}

# Decode sum-space indices back to road IDs
paired_new_RBF   <- new_sensor_location_idxs_RBF[new_sensor_location_idxs_RBF <= length(paired_bases)]
unpaired_new_RBF <- new_sensor_location_idxs_RBF[new_sensor_location_idxs_RBF >  length(paired_bases)]

cat("\nSelected paired road segment base IDs (RBF model, k=10):\n")
print(paired_bases[paired_new_RBF])
cat("\nSelected unpaired WITH road IDs (RBF model, k=10):\n")
print(unpaired_bases[unpaired_new_RBF - length(paired_bases)])

# Check whether previously missed pathological links now appear
hit_RBF <- match(path_sum_idx, new_sensor_location_idxs_RBF)
cat("\nPathological links in RBF sensor selection:\n")
print(data.frame(
  sum_idx         = path_sum_idx,
  sensor_rank_RBF = hit_RBF,
  road_id         = prepared_traffic_links$id[path_idx]
))

# ── Train/test comparison: ordinary besagproper vs RBF-weighted Besag ─────────
# Hold out 20 % of the measured links as a test set; mask their AADT values,
# refit both models, and compare out-of-sample predictive performance.

set.seed(42)
measured_all <- which(!is.na(prepared_traffic_links$aadt))
test_idx_cv  <- sample(measured_all, size = round(0.2 * length(measured_all)))
cat("Measured links:", length(measured_all),
    "  Train:", length(measured_all) - length(test_idx_cv),
    "  Test:", length(test_idx_cv), "\n")

# Mask test responses
data_cv            <- prepared_traffic_links
data_cv$aadt[test_idx_cv] <- NA
data_cv$spatial.idx       <- seq_len(n)

# ── Ordinary besagproper ──────────────────────────────────────────────────────
formula_besag_cv <- aadt ~
  maxLanes + hasOnlyPublicTransportLanes + roadCategory +
  minLanes:roadCategory + lastYearAadt_logAadt +
  f(spatial.idx,
    model  = "besagproper",
    graph  = adjacency_matrix,
    constr = FALSE) +f(roadSystem, model = "iid")

inla_cv_besag <- inla(
  formula_besag_cv,
  data              = data_cv,
  family            = "poisson",
  control.predictor = list(compute = TRUE, link = 1),
  control.compute   = list(dic = TRUE, waic = TRUE)
)

# ── RBF-weighted Besag ────────────────────────────────────────────────────────
# weighted_besag_model_RBF depends only on the graph structure (ui, uj, dist_sq)
# which is unchanged; only the response column differs.
formula_rbf_cv <- aadt ~
  maxLanes + hasOnlyPublicTransportLanes + roadCategory +
  minLanes:roadCategory + lastYearAadt_logAadt +
  f(spatial.idx,
    model  = weighted_besag_model_RBF,
    n      = n,
    constr = FALSE) + 
  f(roadSystem,model = "iid")

inla_cv_RBF <- inla(
  formula_rbf_cv,
  data              = data_cv,
  family            = "poisson",
  control.predictor = list(compute = TRUE, link = 1),
  control.compute   = list(dic = TRUE, waic = TRUE)
)

# ── Extract test-set predictions ──────────────────────────────────────────────
true_aadt  <- prepared_traffic_links$aadt[test_idx_cv]
fv_cv_besag <- inla_cv_besag$summary.fitted.values[test_idx_cv, ]
fv_cv_RBF   <- inla_cv_RBF$summary.fitted.values[test_idx_cv, ]

pred_besag <- fv_cv_besag$mean
pred_RBF   <- fv_cv_RBF$mean

# ── Metrics ───────────────────────────────────────────────────────────────────
rmse     <- function(p, o) sqrt(mean((p - o)^2))
mae      <- function(p, o) mean(abs(p - o))
coverage <- function(lo, hi, o) mean(o >= lo & o <= hi)

cat("\n=== Test-set performance (n =", length(test_idx_cv), "held-out measured links) ===\n")
cat(sprintf("%-28s %12s %12s\n", "Metric", "Besagproper", "RBF Besag"))
cat(strrep("-", 54), "\n")
cat(sprintf("%-28s %12.1f %12.1f\n", "RMSE",
    rmse(pred_besag, true_aadt), rmse(pred_RBF, true_aadt)))
cat(sprintf("%-28s %12.1f %12.1f\n", "MAE",
    mae(pred_besag,  true_aadt), mae(pred_RBF,  true_aadt)))
cat(sprintf("%-28s %12.3f %12.3f\n", "Corr (pred vs obs)",
    cor(pred_besag, true_aadt), cor(pred_RBF, true_aadt)))
cat(sprintf("%-28s %12.3f %12.3f\n", "95 %% CI coverage",
    coverage(fv_cv_besag$`0.025quant`, fv_cv_besag$`0.975quant`, true_aadt),
    coverage(fv_cv_RBF$`0.025quant`,   fv_cv_RBF$`0.975quant`,   true_aadt)))
cat(sprintf("%-28s %12.1f %12.1f\n", "Mean PI width (95%%)",
    mean(fv_cv_besag$`0.975quant` - fv_cv_besag$`0.025quant`),
    mean(fv_cv_RBF$`0.975quant`   - fv_cv_RBF$`0.025quant`)))
cat(sprintf("%-28s %12.1f %12.1f\n", "DIC",
    inla_cv_besag$dic$dic,  inla_cv_RBF$dic$dic))
cat(sprintf("%-28s %12.1f %12.1f\n", "WAIC",
    inla_cv_besag$waic$waic, inla_cv_RBF$waic$waic))

# ── Scatter plot: predicted vs observed ───────────────────────────────────────
lims <- range(c(true_aadt, pred_besag, pred_RBF), na.rm = TRUE)
par(mfrow = c(1, 2))
plot(true_aadt, pred_besag, pch = 16, col = rgb(0, 0, 1, 0.4),
     xlim = lims, ylim = lims,
     xlab = "Observed AADT", ylab = "Predicted AADT",
     main = "Ordinary besagproper")
abline(0, 1, col = "red")
plot(true_aadt, pred_RBF, pch = 16, col = rgb(0, 0.6, 0, 0.4),
     xlim = lims, ylim = lims,
     xlab = "Observed AADT", ylab = "Predicted AADT",
     main = "RBF-weighted Besag")
abline(0, 1, col = "red")
par(mfrow = c(1, 1))
