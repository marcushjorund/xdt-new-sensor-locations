# install.packages("shiny")
# install.packages("bslib")
# install.packages("leaflet")
# install.packages("DT")
# install.packages("rsconnect")
# getwd()
# For local development, run from project root:
# shiny::runApp("shiny/")

library(rsconnect)
rsconnect::setAccountInfo(
  name   = "nyeregistreringspunkter",
  token  = Sys.getenv("SHINYAPPS_TOKEN"),
  secret = Sys.getenv("SHINYAPPS_SECRET")
)

source("deploy.R")

