# =============================================================================
# deploy.R — Deploy the Shiny app to shinyapps.io
#
# Run this from the project root instead of rsconnect::deployApp("shiny/")
# directly. It copies the required data files into shiny/ so they are bundled,
# then cleans them up afterwards.
#
# Usage (from project root):
#   source("deploy.R")
# =============================================================================

library(rsconnect)

# ── Data files to bundle ──────────────────────────────────────────────────────
rds_files <- c(
  "frc_sqrt(0.9_to_0)_100percounty_100overall.rds",
  "functionClass_sqrt(0.6_to_0)_100percounty_100overall.rds",
  "unweighted_100percounty100overall.rds",
  "unweighted_bundle_100percounty_100overall.rds",
  "AADT_prop_100percounty_100overall.rds",
  "AADT_logprop_100percounty_100overall.rds",
  "AADT_identity_100percounty_100overall.rds",
  "functionClass_sqrt(0.6_to_0)logpropaadt_heavy_08015005_100percounty_100overall.rds",
  "frc_sqrt(0.9_to_0)logpropaadt_heavy_08015005_100percounty_100overall.rds",
  "frc_0.6_logaadt_0.3_heavyratio_0.1_100percounty_100overall.rds",
  # Fixed-algorithm comparison files (dynamic-sum neighbourhood weighting)
  "unweighted_distsum_100percounty_100overall.rds",
  "frc_sqrt(0.9_to_0)_distsum_100percounty_100overall.rds"
)
geojson_file <- "directed_traffic_links_2024_simplified.geojson"

# ── Stage data files into shiny/ ──────────────────────────────────────────────
message("Staging data files into shiny/ ...")
dir.create(file.path("shiny", "results", "production"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("shiny", "data"), showWarnings = FALSE)

for (f in rds_files) {
  file.copy(file.path("results", "production", f),
            file.path("shiny", "results", "production", f),
            overwrite = TRUE)
}
file.copy(file.path("data", geojson_file),
          file.path("shiny", "data", geojson_file),
          overwrite = TRUE)
file.copy(file.path("data", "measured_traffic_links_minimal.rds"),
          file.path("shiny", "data", "measured_traffic_links_minimal.rds"),
          overwrite = TRUE)
file.copy(file.path("data", "prepared_traffic_links_norway.rds"),
          file.path("shiny", "data", "prepared_traffic_links_norway.rds"),
          overwrite = TRUE)

# ── Deploy, then always clean up staged files ─────────────────────────────────
tryCatch(
  rsconnect::deployApp("shiny/"),
  finally = {
    message("Removing staged data files from shiny/ ...")
    unlink(file.path("shiny", "results"), recursive = TRUE)
    unlink(file.path("shiny", "data"),    recursive = TRUE)
  }
)
