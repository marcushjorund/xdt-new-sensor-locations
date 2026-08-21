marcus_rds <- readRDS("prepared_traffic_links.rds")
emma_rds <- readRDS("prepared_traffic_links_emma.rds")
install.packages("comparedf")
install.packages("arsenal")
library(arsenal)
comparedf(marcus_rds, emma_rds)
dim(marcus_rds)
dim(emma_rds)
colnames(marcus_rds)
colnames(emma_rds)
install.packages("waldo")
library(waldo)
comparison <- compare(marcus_rds, emma_rds)
print(comparison, n = 50)
marcus_rds$municipality
#Levels: Ål Drammen Flå Flesberg Gol Hemsedal Hol Hole Kongsberg Krødsherad Lier Modum Nesbyen Nore og Uvdal Øvre Eiker Ringerike Rollag Sigdal
emma_rds$municipality
#Levels: Drammen Flesberg Flå Gol Hemsedal Hol Hole Kongsberg Krødsherad Lier Modum Nesbyen Nore og Uvdal Ringerike Rollag Sigdal Øvre Eiker Ål
#I want to check if the levels of the municipality factor variable are in the same order in both datasets, as this could potentially cause issues when fitting the INLA model. The levels are the same but in a different order, which could be a source of discrepancies between the two datasets.
#q: how do i do this?
levels(marcus_rds$municipality)
levels(emma_rds$municipality)

