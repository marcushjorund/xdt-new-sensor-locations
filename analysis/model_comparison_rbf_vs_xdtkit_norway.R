# =============================================================================
# Sammenligning av modeller — Norge, ÅDT (negativ binomial)
#
#   Ny modell    — RBF-vektet Besagproper (Gower-avstand, laplaciansk kjerne)
#                  likhetsvariable: minLanes, highestSpeedLimit, functionClass,
#                                   lastYearAadt_logAadt
#   Gammel modell — xdtkit Besagproper (uten RBF-vekting)
#
# Del 1: Enkelt tren/test-split (80/20) → diagnostiske plott
#   • Modellert vs. observert ÅDT
#   • Residualer mot tilpassede verdier
#   • KK-plott
#
# Del 2: 10-fold kryssvalidering → metrikksammenligning
#   • MALE (Median absolutt logaritmisk feil), RMSE, MAE, MAPE per fold
#   • Boksplott-sammenligning + trykt oppsummeringstabell
#
# Plott lagres til results-figures/
# =============================================================================

library(INLA)
library(xdtkit)
library(ggplot2)
library(patchwork)
source("R/scripts_inla_sensor_placement.R")

# ── Data ─────────────────────────────────────────────────────────────────────
# Laster fra samme kilde som adjacency_matrix_norway ble bygget fra, slik at
# dimensjonene matcher. Fjerner kolonner lagt til etter basisforbehandlingen
# (modellprediksjoner og utledbare-noder-kolonner).
prepared_traffic_links_norway <- readRDS(
  "analysis/model_cache/traffic_links_norway_with_derivable.rds"
)

cols_to_remove <- intersect(
  c("inla_pred", "inla_sd",
    "enables_derivable", "enables_derivable_links", "n_enables_derivable",
    "derivable", "derivable_flow_nodes", "n_derivable"),
  names(prepared_traffic_links_norway)
)
if (length(cols_to_remove) > 0)
  prepared_traffic_links_norway <- prepared_traffic_links_norway[
    , setdiff(names(prepared_traffic_links_norway), cols_to_remove), drop = FALSE
  ]

stopifnot(
  "lastYearAadt_logAadt" %in% names(prepared_traffic_links_norway),
  "aadt"                 %in% names(prepared_traffic_links_norway)
)

adjacency_matrix_norway <- readRDS(
  "analysis/model_cache/adjacency_matrix_norway.rds"
)

if (!dir.exists("results-figures"))
  dir.create("results-figures")

# ── Felles konfigurasjon ──────────────────────────────────────────────────────
covariates <- ~ functionalRoadClass:maxLanes +
  functionalRoadClass:roadCategory +
  minLanes:roadCategory + functionalRoadClass +
  maxLanes + roadCategory +
  hasOnlyPublicTransportLanes +
  functionalRoadClass * isRamp +
  lastYearAadt_logAadt

ordinal_levels_road <- list(
  functionClass     = c("unknown", "E", "D", "C", "B", "A"),
  highestSpeedLimit = c("unknown", "20", "30", "40", "50",
                        "60", "70", "80", "90", "100", "110")
)

similarity_covariates_norway <- c(
  "minLanes", "highestSpeedLimit", "functionClass",
  "lastYearAadt_logAadt"
)

# =============================================================================
# DEL 1: Enkelt tren/test-split — diagnostiske plott
# =============================================================================

set.seed(123)
measured_idx_split <- which(!is.na(prepared_traffic_links_norway$aadt))
test_idx           <- sample(measured_idx_split,
                             size = round(0.2 * length(measured_idx_split)))

data_cv                    <- prepared_traffic_links_norway
data_cv$aadt[test_idx]     <- NA
data_cv$heavyAadt[test_idx] <- NA

# ── Tilpasse modeller på treningssettet ───────────────────────────────────────
message("\n=== [1/2] Tilpasser gammel modell (xdtkit Besagproper) ===")
modell_gammel <- fit_inla_model(
  data          = data_cv,
  adjacency_matrix_norway,
  fixed_effects = covariates,
  iid_effects   = "roadSystem",
  family        = "nbinomial"
)

message("\n=== [2/2] Tilpasser ny modell (RBF-vektet Besagproper) ===")
modell_ny <- fit_weighted_inla_model(
  data                  = data_cv,
  adjacency_matrix      = adjacency_matrix_norway,
  spatial_term          = "besagproper_rbf",
  fixed_effects         = covariates,
  iid_effects           = "roadSystem",
  ordinal_levels        = ordinal_levels_road,
  similarity_covariates = similarity_covariates_norway,
  weight_fn             = "laplacian",
  distance_type         = "gower",
  family                = "nbinomial"
)

# ── Plottfunksjon: side om side ───────────────────────────────────────────────
# plot_inla_model() returnerer plott-objektet usynlig, og akseetiketter og
# y-akse-omfang kan sendes inn direkte via x_label / y_label / ylim.
.komponer_plott <- function(type, tittel_no, x_no, y_no) {
  # Hent prediksjoner og beregn felles y-omfang direkte fra modelldata
  pred_g <- modell_gammel$predictions$inla_pred[test_idx]
  pred_n <- modell_ny$predictions$inla_pred[test_idx]
  obs    <- prepared_traffic_links_norway$aadt[test_idx]
  resid  <- c(pred_g - obs, pred_n - obs)

  y_lim <- switch(type,
    pred_vs_obs         = range(c(pred_g, pred_n, obs), na.rm = TRUE),
    residuals_vs_fitted = range(resid,                   na.rm = TRUE),
    qq                  = range(resid,                   na.rm = TRUE),
    NULL
  )

  p_gammel <- plot_inla_model(
    modell_gammel, prepared_traffic_links_norway$aadt, test_idx,
    type    = type,
    title   = "Ikke-intrinsisk Besag modell (gammel)",
    x_label = x_no,
    y_label = y_no,
    ylim    = y_lim
  )

  p_ny <- plot_inla_model(
    modell_ny, prepared_traffic_links_norway$aadt, test_idx,
    type    = type,
    title   = "Laplace-kernel vektet Besag modell (ny)",
    x_label = x_no,
    y_label = y_no,
    ylim    = y_lim
  )

  (p_gammel | p_ny) +
    patchwork::plot_annotation(title = tittel_no)
}

# ── Modellert vs. observert ───────────────────────────────────────────────────
p_pred_vs_obs <- .komponer_plott(
  type      = "pred_vs_obs",
  tittel_no = "Modellert vs. observert \u00c5DT \u2014 enkelt tren/test-split (20\u00a0% testset)",
  x_no      = "Observert \u00c5DT",
  y_no      = "Modellert \u00c5DT"
)
ggsave("results-figures/sammenligning_modellert_vs_observert.png",
       p_pred_vs_obs, width = 14, height = 6, dpi = 150)
message("Lagret: results-figures/sammenligning_modellert_vs_observert.png")

# ── Residualer mot tilpassede verdier ─────────────────────────────────────────
p_residualer <- .komponer_plott(
  type      = "residuals_vs_fitted",
  tittel_no = "Residualer mot tilpassede verdier \u2014 enkelt tren/test-split",
  x_no      = "Tilpasset \u00c5DT",
  y_no      = "Residual (Tilpasset \u2212 Observert)"
)
ggsave("results-figures/sammenligning_residualer.png",
       p_residualer, width = 14, height = 6, dpi = 150)
message("Lagret: results-figures/sammenligning_residualer.png")

# ── KK-plott ──────────────────────────────────────────────────────────────────
p_qq <- .komponer_plott(
  type      = "qq",
  tittel_no = "KK-plott \u2014 residualer fra testsettet",
  x_no      = "Teoretiske kvantiler",
  y_no      = "Observerte kvantiler"
)
ggsave("results-figures/sammenligning_qq.png",
       p_qq, width = 14, height = 6, dpi = 150)
message("Lagret: results-figures/sammenligning_qq.png")

# =============================================================================
# DEL 2: 10-fold kryssvalidering — metrikksammenligning
# =============================================================================

k_cv    <- 10
seed_cv <- 420

# ── Ny modell: kfold_cv_inla() ────────────────────────────────────────────────
message("\n=== 10-fold KV [1/2]: Ny modell (RBF-vektet Besagproper) ===")
kv_ny <- kfold_cv_inla(
  data                  = prepared_traffic_links_norway,
  adjacency_matrix      = adjacency_matrix_norway,
  k                     = k_cv,
  seed                  = seed_cv,
  fixed_effects         = covariates,
  spatial_term          = "besagproper_rbf",
  similarity_covariates = similarity_covariates_norway,
  ordinal_levels        = ordinal_levels_road,
  iid_effects           = "roadSystem",
  family                = "nbinomial",
  weight_fn             = "laplacian",
  distance_type         = "gower"
)

# ── Gammel modell: manuell fold-løkke ────────────────────────────────────────
# Replikerer fold-tildelingen fra kfold_cv_inla():
#   set.seed(seed_cv) → sample(rep_len(1:k, length(measured_idx)))
message("\n=== 10-fold KV [2/2]: Gammel modell (xdtkit Besagproper) ===")

measured_idx_kv <- which(!is.na(prepared_traffic_links_norway$aadt))
set.seed(seed_cv)
fold_id <- sample(rep_len(seq_len(k_cv), length(measured_idx_kv)))

fold_statistikk_gammel <- vector("list", k_cv)

for (fold in seq_len(k_cv)) {
  message(sprintf("\n\u2500\u2500 Fold %d / %d \u2500\u2500", fold, k_cv))

  test_pos  <- measured_idx_kv[fold_id == fold]
  data_fold <- prepared_traffic_links_norway
  data_fold$aadt[test_pos]      <- NA
  data_fold$heavyAadt[test_pos] <- NA

  tilpasning <- fit_inla_model(
    data          = data_fold,
    adjacency_matrix_norway,
    fixed_effects = covariates,
    iid_effects   = "roadSystem",
    family        = "nbinomial"
  )

  pred_vals <- tilpasning$predictions$inla_pred[test_pos]
  true_vals <- prepared_traffic_links_norway$aadt[test_pos]

  fold_statistikk_gammel[[fold]] <- data.frame(
    fold = fold,
    RMSE = .rmse(pred_vals, true_vals),
    MAE  = .mae( pred_vals, true_vals),
    MAPE = .mape(pred_vals, true_vals),
    MALE = .male(pred_vals, true_vals)
  )

  message(sprintf("  Fold %d  RMSE = %.0f  MAE = %.0f  MAPE = %.1f%%  MALE = %.3f",
                  fold,
                  fold_statistikk_gammel[[fold]]$RMSE,
                  fold_statistikk_gammel[[fold]]$MAE,
                  fold_statistikk_gammel[[fold]]$MAPE,
                  fold_statistikk_gammel[[fold]]$MALE))
}

fold_metrikker_gammel <- do.call(rbind, fold_statistikk_gammel)

# ── Metrikkplott: boksplott per fold ─────────────────────────────────────────
fold_data <- rbind(
  transform(kv_ny$fold_metrics,    Modell = "Laplace-kernel vektet Besag modell (ny)"),
  transform(fold_metrikker_gammel, Modell = "Ikke-intrinsisk Besag modell (gammel)")
)

farger <- c("Laplace-kernel vektet Besag modell (ny)" = "#0066CC",
            "Ikke-intrinsisk Besag modell (gammel)" = "#CC3300")

.metrikk_boks <- function(df, kolonne, y_etikett) {
  ggplot2::ggplot(df, ggplot2::aes(x = Modell, y = .data[[kolonne]], fill = Modell)) +
    ggplot2::geom_boxplot(alpha = 0.65, outlier.shape = NA, width = 0.5) +
    ggplot2::geom_jitter(ggplot2::aes(colour = Modell),
                         width = 0.12, size = 2.2, alpha = 0.85) +
    ggplot2::scale_fill_manual(values = farger)   +
    ggplot2::scale_colour_manual(values = farger) +
    ggplot2::labs(x = NULL, y = y_etikett)        +
    ggplot2::theme_bw(base_size = 11)             +
    ggplot2::theme(
      legend.position  = "none",
      axis.text.x      = ggplot2::element_text(angle = 15, hjust = 1)
    )
}

p_metrikker <- (
  .metrikk_boks(fold_data, "MALE", "MALE") |
  .metrikk_boks(fold_data, "MAPE", "MAPE (%)")
) / (
  .metrikk_boks(fold_data, "RMSE", "RMSE") |
  .metrikk_boks(fold_data, "MAE",  "MAE")
) +
  patchwork::plot_annotation(
    title    = "10-fold kryssvalidering \u2014 Laplace-kernel vektet Besag modell (ny) vs. Ikke-intrinsisk Besag modell (gammel)",
    subtitle = "Hvert punkt representerer \u00e9n fold; boks viser median og IKR"
  )

ggsave("results-figures/sammenligning_metrikker_kv.png",
       p_metrikker, width = 12, height = 10, dpi = 150)
message("Lagret: results-figures/sammenligning_metrikker_kv.png")

# ── Trykt oppsummeringstabell ─────────────────────────────────────────────────
.lag_oppsummering <- function(fm, modell_navn) {
  data.frame(
    Modell = modell_navn,
    RMSE   = sprintf("%.0f \u00b1 %.0f", mean(fm$RMSE), sd(fm$RMSE)),
    MAE    = sprintf("%.0f \u00b1 %.0f", mean(fm$MAE),  sd(fm$MAE)),
    MAPE   = sprintf("%.1f \u00b1 %.1f", mean(fm$MAPE), sd(fm$MAPE)),
    MALE   = sprintf("%.3f \u00b1 %.3f", mean(fm$MALE), sd(fm$MALE))
  )
}

sammenligning <- rbind(
  .lag_oppsummering(kv_ny$fold_metrics,    "Laplace-kernel vektet Besag modell (ny)"),
  .lag_oppsummering(fold_metrikker_gammel, "Ikke-intrinsisk Besag modell (gammel)")
)

sep <- paste0(strrep("\u2550", 70), "\n")
cat("\n", sep, sep = "")
cat("  10-fold KV \u2014 Laplace-kernel vektet Besag modell (ny) vs. Ikke-intrinsisk Besag modell (gammel)\n")
cat("  Format: gjennomsnitt \u00b1 standardavvik over 10 folder\n")
cat("  MALE = Median absolutt logaritmisk feil (lavere = bedre)\n")
cat(sep)
print(sammenligning, row.names = FALSE)
cat(sep)
