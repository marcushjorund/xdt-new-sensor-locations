library(xdtkit)
library(INLA)
library(ggplot2)
year <- 2025

# Traffic links: Load and preprocess
preprocessed_traffic_links <- preprocess_traffic_links(buskerud_directed_traffic_links, year = year)

# Bus data: Load and preprocess
bus_aadt <- calculate_bus_aadt(stops_on_traffic_links, bus_counts, year = year)

# Fill missing values and add bus data
prepared_traffic_links <- fill_missing_values(
  df = preprocessed_traffic_links,
  unknown_impute_columns = c("functionClass", "highestSpeedLimit", "lowestSpeedLimit","maxLanes", "minLanes"),
  mode_impute_columns = c("hasOnlyPublicTransportLanes"),
  median_impute_columns = c("lastYearAadt_aadt", "lastYearAadt_heavyRatio",
                            "lastYearAadt_heavyAadt")) |>
  remove_negative_aadt() |>
  add_logLastYear() |>
  join_bus_to_traffic(bus_aadt)
#prepared_traffic_links <- readRDS("prepared_traffic_links_emma.rds")
#q: how do I save prepared_traffic_links as an rds file
#a: 


# Adjacency matrix (may take several minutes to run)
adjacency_matrix <- build_adjacency_matrix(
  prepared_traffic_links,
  exclude_public_transport = TRUE)
# Balancing clusters
clusters <- strategic_network_clustering(
  data = prepared_traffic_links,
  year = year,
  boundary_links = c("Trafikkdata_continuous"))

# Nodes: Load and preprocess (may take a minute to run)
nodes <- identify_unbalanceable_nodes(buskerud_nodes, prepared_traffic_links)

prepared_traffic_links$spatial.idx <- 1:nrow(prepared_traffic_links)

formula <- aadt ~ f(spatial.idx, model = "besagproper", graph = adjacency_matrix,
                    adjust.for.con.comp = FALSE, constr = TRUE) + 
  f(roadSystem, model = "iid") +
  lastYearAadt_logAadt +isRamp + hasOnlyPublicTransportLanes + maxLanes + minLanes + functionalRoadClass
  
# formula <- aadt ~ lastYearAadt_logAadt + functionalRoadClass

# formula <- aadt ~ f(spatial.idx, model = "besagproper", graph = adjacency_matrix,
#                     adjust.for.con.comp = FALSE, constr = TRUE) + f(roadSystem,
#                                                                     model = "iid") + functionalRoadClass + maxLanes + roadCategory +
#   hasOnlyPublicTransportLanes + isRamp + lastYearAadt_logAadt +
#   functionalRoadClass:maxLanes + functionalRoadClass:roadCategory +
#   roadCategory:minLanes + functionalRoadClass:isRamp

mod_test <- inla(formula, data = prepared_traffic_links, family = "poisson")
?inla
summary(mod_test)
mod_test$
plot(prepared_traffic_links$lastYearAadt_aadt, mod_test$summary.fitted.values$mean)


prepared_traffic_links$inla_pred <- mod_test$summary.fitted.values$mean

ggplot(prepared_traffic_links, aes(x = aadt, y = inla_pred)) +
  geom_point()
names(prepared_traffic_links)

factor_variables <- c("functionalRoadClass", "roadCategory", "isRamp", "hasOnlyPublicTransportLanes", "minLanes", "maxLanes")
frequency_table <- as.data.frame(table(prepared_traffic_links[,factor_variables]))
sum(frequency_table[which(frequency_table$Freq != 0),"Freq"])
dim(frequency_table)
sum(frequency_table$Freq == 0)
sum(frequency_table$Freq)
frequency_table$Freq[fre]
#Creating a 5x5 precision matrix, with 5 nodes. Node 1 is neighbour with node 2 and 4, node 2 is neighbour with node 1 and 3, node 3 is neighbour with node 2,
#node 4 is neiighbour with node 1 and 5 and node 5 is neighbour with node 4
d <- 1
tau <- 3
adjacency_matrix_test <- matrix(0, nrow = 5, ncol = 5)
adjacency_matrix_test[1,] <- c(2+d,-1,0,-1,0)
adjacency_matrix_test[2,] <- c(-1,2+d,-1,0,0)
adjacency_matrix_test[3,] <- c(0,-1,1+d,0,0)
adjacency_matrix_test[4,] <- c(-1,0,0,2+d,-1)
adjacency_matrix_test[5,] <- c(0,0,0,-1,1+d)
adjacency_matrix_test <- adjacency_matrix_test
adjacency_matrix_test
solve(adjacency_matrix_test, sparse = FALSE)
mod_test$model.matrix
inla.graph2matrix(adjacency_matrix)
inla.spy(adjacency_matrix)
adjacency_matrix[2,]
dim(adjacency_matrix)
prepared_traffic_links$aadt
