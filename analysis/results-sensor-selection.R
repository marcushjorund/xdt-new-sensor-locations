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
library(dplyr)
source("R/scripts_inla_sensor_placement.R")

# ── Load saved outputs from the Norway-wide preprocessing pipeline ────────────
prepared_traffic_links_norway <- readRDS("analysis/model_cache/traffic_links_norway_with_derivable.rds")
adjacency_matrix_norway       <- readRDS("analysis/model_cache/adjacency_matrix_norway.rds")
inla_weighted_norway_full     <- readRDS("analysis/model_cache/inla_rbf_norway_full.rds")

# ── Extract spatial hyperparameters and edge distances ────────────────────────
spatial_hyperpar <- inla_weighted_norway_full$spatial_hyperpar[2:4, "mean"]
tau       <- spatial_hyperpar[1]
d         <- spatial_hyperpar[2]
sigma     <- spatial_hyperpar[3]
distances <- inla_weighted_norway_full$distances

# ═══════════════════════════════════════════════════════════════════════════════
# PART 2: Sensor Selection Experiments
# Unweighted selection for all of Norway and Trøndelag.
# Results are saved to results/report/.
# ═══════════════════════════════════════════════════════════════════════════════

# ── Experiment 1: Norway — unweighted, k = 100 per county, r = 100 overall ───
norway_unweighted <- greedy_mi_sensor_selection_norway(
  data                = prepared_traffic_links_norway,
  adjacency_matrix    = adjacency_matrix_norway,
  distances           = distances,
  tau                 = tau,
  d                   = d,
  sigma               = sigma,
  hops                = 3,
  k                   = 100,
  r                   = 100,
  weighting_vars      = NULL,
  neighbourhood_decay = c(1, 0.5, 0.25),
  neighbour_hops      = 2
)

plot_sensor_selection_map(norway_unweighted, all_data = prepared_traffic_links_norway)
norway_unweighted$summary
norway_unweighted$summary$counts_table

saveRDS(norway_unweighted, "results/report/unweighted_100percounty_100overall.rds")

# ── Experiment 2: Trøndelag — unweighted, k = 20 per county, r = 20 overall ──
trondelag_unweighted <- greedy_mi_sensor_selection_norway(
  data                = prepared_traffic_links_norway,
  adjacency_matrix    = adjacency_matrix_norway,
  distances           = distances,
  tau                 = tau,
  d                   = d,
  sigma               = sigma,
  hops                = 3,
  k                   = 20,
  r                   = 20,
  weighting_vars      = NULL,
  neighbourhood_decay = c(1, 0.5, 0.25),
  neighbour_hops      = 2,
  filtering           = list(county = "Trøndelag")
)

plot_sensor_selection_map(trondelag_unweighted, all_data = prepared_traffic_links_norway[prepared_traffic_links_norway$county == "Trøndelag", ])
trondelag_unweighted$summary
trondelag_unweighted$summary$counts_table

saveRDS(trondelag_unweighted, "results/report/trondelag_unweighted_20percounty_20overall.rds")

# ═══════════════════════════════════════════════════════════════════════════════
# PART 3: Derivable MI Sensor Selection Experiments
# Uses bundle_scoring = TRUE, which runs greedy_mi_sensor_selection_derivable()
# per county partition. Each bundle scores the placed sensor + its WITH/AGAINST
# directions + any flow-conservation derived links jointly.
# ═══════════════════════════════════════════════════════════════════════════════

# ── Experiment 3: Norway — unweighted derivable, k = 100 per county, r = 100 overall ──
norway_unweighted_derivable <- greedy_mi_sensor_selection_norway(
  data                = prepared_traffic_links_norway,
  adjacency_matrix    = adjacency_matrix_norway,
  distances           = distances,
  tau                 = tau,
  d                   = d,
  sigma               = sigma,
  hops                = 3,
  k                   = 100,
  r                   = 100,
  weighting_vars      = NULL,
  neighbourhood_decay = c(1, 0.5, 0.25),
  neighbour_hops      = 2,
  bundle_scoring      = TRUE
)

plot_sensor_selection_map(norway_unweighted_derivable, all_data = prepared_traffic_links_norway)
norway_unweighted_derivable$summary
norway_unweighted_derivable$summary$counts_table

saveRDS(norway_unweighted_derivable, "results/report/unweighted_derivable_100percounty_100overall.rds")

# ── Experiment 4: Trøndelag — unweighted derivable, k = 20 per county, r = 20 overall ──
trondelag_unweighted_derivable <- greedy_mi_sensor_selection_norway(
  data                = prepared_traffic_links_norway,
  adjacency_matrix    = adjacency_matrix_norway,
  distances           = distances,
  tau                 = tau,
  d                   = d,
  sigma               = sigma,
  hops                = 3,
  k                   = 20,
  r                   = 20,
  weighting_vars      = NULL,
  neighbourhood_decay = c(1, 0.5, 0.25),
  neighbour_hops      = 2,
  bundle_scoring      = TRUE,
  filtering           = list(county = "Trøndelag")
)

plot_sensor_selection_map(trondelag_unweighted_derivable, all_data = prepared_traffic_links_norway[prepared_traffic_links_norway$county == "Trøndelag", ])
trondelag_unweighted_derivable$summary
trondelag_unweighted_derivable$summary$counts_table

saveRDS(trondelag_unweighted_derivable, "results/report/trondelag_unweighted_derivable_20percounty_20overall.rds")

# ═══════════════════════════════════════════════════════════════════════════════
# PART 4: AADT-Weighted Sensor Selection Experiments
# Uses predicted AADT from the INLA model (inla_pred) as the weighting variable.
# Two scaling modes: identity (raw predicted AADT) and log_proportional.
# Results are saved to results/report/.
# ═══════════════════════════════════════════════════════════════════════════════

# ── Experiment 5: Norway — AADT identity, k = 100 per county, r = 100 overall ─
norway_aadt_identity <- greedy_mi_sensor_selection_norway(
  data                = prepared_traffic_links_norway,
  adjacency_matrix    = adjacency_matrix_norway,
  distances           = distances,
  tau                 = tau,
  d                   = d,
  sigma               = sigma,
  hops                = 3,
  k                   = 100,
  r                   = 100,
  weighting_vars      = list(inla_pred = "identity"),
  neighbourhood_decay  = c(1, 0.5,0.25,0.125),
  neighbour_hops       = 2,
  neighbourhood_scoring = "log_proportional"
)

plot_sensor_selection_map(norway_aadt_identity, all_data = prepared_traffic_links_norway)
norway_aadt_identity$summary
norway_aadt_identity$summary$counts_table

saveRDS(norway_aadt_identity, "results/report/norway_aadt_identity_100percounty_100overall.rds")

# ── Experiment 6: Trøndelag — AADT identity, k = 20 per county, r = 20 overall ─
trondelag_aadt_identity <- greedy_mi_sensor_selection_norway(
  data                = prepared_traffic_links_norway,
  adjacency_matrix    = adjacency_matrix_norway,
  distances           = distances,
  tau                 = tau,
  d                   = d,
  sigma               = sigma,
  hops                = 3,
  k                   = 20,
  r                   = 20,
  weighting_vars      = list(inla_pred = "identity"),
  neighbourhood_decay  = c(1, 0.5,0.25,0.125),
  neighbour_hops       = 2,
  neighbourhood_scoring = "log_proportional",
  filtering           = list(county = "Trøndelag")
)

plot_sensor_selection_map(trondelag_aadt_identity, all_data = prepared_traffic_links_norway[prepared_traffic_links_norway$county == "Trøndelag", ])
trondelag_aadt_identity$summary
trondelag_aadt_identity$summary$counts_table

saveRDS(trondelag_aadt_identity, "results/report/trondelag_aadt_identity_20percounty_20overall.rds")

# ── Experiment 7: Norway — AADT log-proportional, k = 100 per county, r = 100 overall ─
norway_aadt_logprop <- greedy_mi_sensor_selection_norway(
  data                = prepared_traffic_links_norway,
  adjacency_matrix    = adjacency_matrix_norway,
  distances           = distances,
  tau                 = tau,
  d                   = d,
  sigma               = sigma,
  hops                = 3,
  k                   = 100,
  r                   = 100,
  weighting_vars      = list(inla_pred = "log_proportional"),
  neighbourhood_decay  = c(1, 0.5,0.25,0.125),
  neighbour_hops       = 2,
  neighbourhood_scoring = "log_proportional",
)

plot_sensor_selection_map(norway_aadt_logprop, all_data = prepared_traffic_links_norway)
norway_aadt_logprop$summary
norway_aadt_logprop$summary$counts_table

saveRDS(norway_aadt_logprop, "results/report/norway_aadt_logprop_100percounty_100overall.rds")

# ── Experiment 8: Trøndelag — AADT log-proportional, k = 20 per county, r = 20 overall ─
trondelag_aadt_logprop <- greedy_mi_sensor_selection_norway(
  data                = prepared_traffic_links_norway,
  adjacency_matrix    = adjacency_matrix_norway,
  distances           = distances,
  tau                 = tau,
  d                   = d,
  sigma               = sigma,
  hops                = 3,
  k                   = 20,
  r                   = 20,
  weighting_vars      = list(inla_pred = "log_proportional"),
  neighbourhood_decay  = c(1, 0.5,0.25,0.125),
  neighbour_hops       = 2,
  neighbourhood_scoring = "log_proportional",
  filtering           = list(county = "Trøndelag")
)

plot_sensor_selection_map(trondelag_aadt_logprop, all_data = prepared_traffic_links_norway[prepared_traffic_links_norway$county == "Trøndelag", ])
trondelag_aadt_logprop$summary
trondelag_aadt_logprop$summary$counts_table

saveRDS(trondelag_aadt_logprop, "results/report/trondelag_aadt_logprop_20percounty_20overall.rds")

# ═══════════════════════════════════════════════════════════════════════════════
# PART 5: functionClass sqrt(5-to-1) Weighted Sensor Selection Experiments
# Uses functionClass with sqrt(seq(5, 1, by = -1)) weights mapped to classes
# A–E; every measured class retains a non-zero weight (floor = sqrt(1) = 1).
# The "unknown" class receives weight 0.
# Results are saved to results/report/.
# ═══════════════════════════════════════════════════════════════════════════════

# ── Experiment 9: Norway — functionClass sqrt(5_to_1), k = 100 per county, r = 100 overall ──
norway_fc_sqrt5to1 <- greedy_mi_sensor_selection_norway(
  data                = prepared_traffic_links_norway,
  adjacency_matrix    = adjacency_matrix_norway,
  distances           = distances,
  tau                 = tau,
  d                   = d,
  sigma               = sigma,
  hops                = 3,
  k                   = 100,
  r                   = 100,
  weighting_vars      = list(
    functionClass = setNames(c(sqrt(seq(1.5, 1, by = -0.1)), 0), c("A", "B", "C", "D", "E", "unknown"))
  ),
  neighbourhood_decay  = c(1, 0.5,0.25,0.125),
  neighbour_hops       = 2,
  neighbourhood_scoring = "log_proportional",
)

plot_sensor_selection_map(norway_fc_sqrt5to1, all_data = prepared_traffic_links_norway)
norway_fc_sqrt5to1$summary
norway_fc_sqrt5to1$summary$counts_table

saveRDS(norway_fc_sqrt5to1, "results/report/functionClass_sqrt(5_to_1)_100percounty_100overall.rds")

# ── Experiment 10: Trøndelag — functionClass sqrt(5_to_1), k = 20 per county, r = 20 overall ──
trondelag_fc_sqrt5to1 <- greedy_mi_sensor_selection_norway(
  data                = prepared_traffic_links_norway,
  adjacency_matrix    = adjacency_matrix_norway,
  distances           = distances,
  tau                 = tau,
  d                   = d,
  sigma               = sigma,
  hops                = 3,
  k                   = 20,
  r                   = 20,
  weighting_vars      = list(
    functionClass = setNames(c(sqrt(seq(1.5, 1, by = -0.1)), 0), c("A", "B", "C", "D", "E", "unknown"))
  ),
  neighbourhood_decay  = c(1, 0.5,0.25,0.125),
  neighbour_hops       = 2,
  neighbourhood_scoring = "log_proportional",
  filtering           = list(county = "Trøndelag")
)

plot_sensor_selection_map(trondelag_fc_sqrt5to1, all_data = prepared_traffic_links_norway[prepared_traffic_links_norway$county == "Trøndelag", ])
trondelag_fc_sqrt5to1$summary
trondelag_fc_sqrt5to1$summary$counts_table

saveRDS(trondelag_fc_sqrt5to1, "results/report/trondelag_functionClass_sqrt(5_to_1)_20percounty_20overall.rds")

# ═══════════════════════════════════════════════════════════════════════════════
# PART 6: functionClass sqrt(5-to-1) + AADT log-proportional + Heavy Ratio
# Combines road-class priority (floor weight = 1) with predicted traffic volume
# and heavy-vehicle ratio. Alpha split: 80% functionClass / 15% AADT / 5% heavy ratio.
# Results are saved to results/report/.
# ═══════════════════════════════════════════════════════════════════════════════

# ── Experiment 11: Norway — functionClass sqrt(5_to_1) + logprop AADT + heavy ratio, k = 100 per county, r = 100 overall ──
norway_fc_sqrt5to1_aadt_heavy <- greedy_mi_sensor_selection_norway(
  data                 = prepared_traffic_links_norway,
  adjacency_matrix     = adjacency_matrix_norway,
  distances            = distances,
  tau                  = tau,
  d                    = d,
  sigma                = sigma,
  hops                 = 3,
  k                    = 100,
  r                    = 100,
  weighting_vars       = list(
    functionClass           = setNames(c(sqrt(seq(1.5, 1, by = -0.1)), 0), c("A", "B", "C", "D", "E", "unknown")),
    inla_pred               = "log_proportional",
    lastYearAadt_heavyRatio = "identity"
  ),
  weighting_vars_alpha = c(0.7, 0.25, 0.05),
  neighbourhood_decay  = c(1,0.5,0.25,0.125),
  neighbourhood_scoring = "log_proportional",
  neighbour_hops       = 2
)

plot_sensor_selection_map(norway_fc_sqrt5to1_aadt_heavy, all_data = prepared_traffic_links_norway)
norway_fc_sqrt5to1_aadt_heavy$summary
norway_fc_sqrt5to1_aadt_heavy$summary$counts_table

saveRDS(norway_fc_sqrt5to1_aadt_heavy, "results/report/functionClass_sqrt(5_to_1)_logpropaadt_heavyratio_100percounty_100overall.rds")

# ── Experiment 12: Trøndelag — functionClass sqrt(5_to_1) + logprop AADT + heavy ratio, k = 20 per county, r = 20 overall ──
trondelag_fc_sqrt5to1_aadt_heavy <- greedy_mi_sensor_selection_norway(
  data                 = prepared_traffic_links_norway,
  adjacency_matrix     = adjacency_matrix_norway,
  distances            = distances,
  tau                  = tau,
  d                    = d,
  sigma                = sigma,
  hops                 = 3,
  k                    = 20,
  r                    = 20,
  weighting_vars       = list(
    functionClass           = setNames(c(sqrt(seq(1.5, 1, by = -0.1)), 0), c("A", "B", "C", "D", "E", "unknown")),
    inla_pred               = "log_proportional",
    lastYearAadt_heavyRatio = "identity"
  ),
  weighting_vars_alpha = c(0.7, 0.25, 0.05),
  neighbourhood_decay  = c(1, 0.5,0.25,0.125),
  neighbour_hops       = 2,
  neighbourhood_scoring = "log_proportional",
  filtering            = list(county = "Trøndelag")
)

plot_sensor_selection_map(trondelag_fc_sqrt5to1_aadt_heavy, all_data = prepared_traffic_links_norway[prepared_traffic_links_norway$county == "Trøndelag", ])
trondelag_fc_sqrt5to1_aadt_heavy$summary
trondelag_fc_sqrt5to1_aadt_heavy$summary$counts_table

saveRDS(trondelag_fc_sqrt5to1_aadt_heavy, "results/report/trondelag_functionClass_sqrt(5_to_1)_logpropaadt_heavyratio_20percounty_20overall.rds")

# ═══════════════════════════════════════════════════════════════════════════════
# PART 7: Verifisering av submodularitet – avtakende MI-bidrag
# Marginalt MI-bidrag per valgt trafikkregistreringspunkt skal være ikke-stigende
# innenfor hvert fylke, noe som bekrefter at submodularitetsbetingelsen holder.
# Plottet vises for både uvektet og vektet utvalg.
# ═══════════════════════════════════════════════════════════════════════════════

library(ggplot2)

# ── Last inn resultater ───────────────────────────────────────────────────────
res_uvektet <- readRDS("results/report/unweighted_100percounty_100overall.rds")
res_vektet  <- readRDS("results/report/functionClass_sqrt(5_to_1)_logpropaadt_heavyratio_100percounty_100overall.rds")

farger <- c("Uvektet"                               = "#2166ac",
            "Vektet (vegklasse + AADT + tungandel)" = "#d6604d")

# ── Hjelpefunksjon: bygg datasett fra mi_scores ───────────────────────────────
lag_mi_df <- function(res, etikett) {
  skarer <- res$mi_scores
  data.frame(
    nr          = seq_along(skarer),
    mi_marginal = skarer,
    modell      = etikett
  )
}

mi_df <- rbind(
  lag_mi_df(res_uvektet, "Uvektet"),
  lag_mi_df(res_vektet,  "Vektet (vegklasse + AADT + tungandel)")
)
mi_df$modell <- factor(mi_df$modell,
                       levels = c("Uvektet", "Vektet (vegklasse + AADT + tungandel)"))

# ── Figur: marginalt MI-bidrag (submodularitetsplot) ─────────────────────────
figur_submodular <- ggplot(mi_df, aes(x = nr, y = mi_marginal, colour = modell)) +
  geom_line(linewidth = 0.8, alpha = 0.9) +
  scale_colour_manual(values = farger) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title   = "Avtakende MI-bidrag ved grådig plassering av trafikkregistreringspunkt",
    subtitle = paste0(
      "Ikke-stigende marginalgevinst innenfor hvert fylke bekrefter at ",
      "submodularitetsbetingelsen holder,\n",
      "også når vektingsskjemaet er aktivt."
    ),
    x       = "Rekkefølge for valgte trafikkregistreringspunkt",
    y       = "\u0394 Gjensidig informasjon (marginalt bidrag)",
    colour  = NULL,
    caption = paste0(
      "Grådig plassering av trafikkregistreringspunkt basert på gjensidig informasjon (Krause m.fl., 2008). ",
      "k\u00a0=\u00a0100 trafikkregistreringspunkt per fylke, r\u00a0=\u00a0100 totalt. ",
      "Vektet modell bruker vegklasse (functionClass), predikert AADT og tungbilandel ",
      "med vektfordeling 80\u00a0%\u00a0/\u00a015\u00a0%\u00a0/\u00a05\u00a0%."
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "bottom",
    legend.text      = element_text(size = 10),
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(colour = "grey40", size = 10),
    plot.caption     = element_text(colour = "grey45", size = 9, hjust = 0),
    panel.grid.minor = element_blank()
  )

print(figur_submodular)

