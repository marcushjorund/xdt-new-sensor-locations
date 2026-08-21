# ──────────────────────────────────────────────────────────────────────────────
# Likelihood Family Comparison: 10-Fold CV — RBF AADT-similarity Model
# ──────────────────────────────────────────────────────────────────────────────
# Compares three likelihood specifications on the same spatial model structure
# (RBF-weighted Besag proper, AADT-only neighbour similarity):
#
#   A. Gaussian on log(AADT)  — numerically stable; log-normal assumption
#   B. Negative Binomial      — count model with estimated overdispersion
#   C. Poisson (stabilised)   — count model; strategy="gaussian" + step.factor
#
# Notes:
#   • Poisson/NB use control.inla = list(strategy="gaussian", step.factor=0.1,
#     diagonal=1e-3) to reduce hyperparameter overshoot during optimisation.
#     Failed folds are excluded from metrics (tryCatch writes NULL).
#   • Fitted values from Gaussian are on the log scale → back-transformed with
#     exp() for all metrics. Poisson/NB fitted values are already E[y]=exp(η).
#   • CI coverage compares posterior credible intervals for E[y] against
#     observed y. For Gaussian this includes residual variance σ² (wider); for
#     Poisson/NB it is the CI for the rate parameter (narrower). Interpret
#     comparatively, not as absolute coverage of observed counts.
# ──────────────────────────────────────────────────────────────────────────────

library(xdtkit)
library(INLA)
library(Matrix)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)   # comma() labels; always available as ggplot2 dependency

# ──────────────────────────────────────────────────────────────────────────────
# 1. DATA PREPARATION
# ──────────────────────────────────────────────────────────────────────────────

year <- 2025
preprocessed_traffic_links <- preprocess_traffic_links(buskerud_directed_traffic_links, year = year)
bus_aadt <- calculate_bus_aadt(stops_on_traffic_links, bus_counts, year = year)

prepared_traffic_links <- fill_missing_values(
  df                    = preprocessed_traffic_links,
  unknown_impute_columns  = c("functionClass", "highestSpeedLimit", "lowestSpeedLimit",
                               "maxLanes", "minLanes"),
  mode_impute_columns     = c("hasOnlyPublicTransportLanes"),
  median_impute_columns   = c("lastYearAadt_aadt", "lastYearAadt_heavyRatio",
                               "lastYearAadt_heavyAadt")
) |>
  remove_negative_aadt() |>
  add_logLastYear() |>
  join_bus_to_traffic(bus_aadt)

adjacency_matrix <- build_adjacency_matrix(
  prepared_traffic_links,
  exclude_public_transport = TRUE
)

n <- nrow(prepared_traffic_links)
prepared_traffic_links$spatial.idx <- seq_len(n)
prepared_traffic_links$log_aadt    <- log(prepared_traffic_links$aadt)
# Poisson/NB require integer counts; AADT is an annual average so we round
prepared_traffic_links$aadt_int    <- as.integer(round(prepared_traffic_links$aadt))

n_measured <- sum(!is.na(prepared_traffic_links$aadt))
cat("n =", n, "links | measured =", n_measured, "\n\n")

# ──────────────────────────────────────────────────────────────────────────────
# 2. AADT-ONLY GRAPH DISTANCES (for RBF kernel)
# ──────────────────────────────────────────────────────────────────────────────
# W_ij = exp(-|AADT_i - AADT_j|^2 / (2*ell^2)), ell estimated from data.

adj_trip   <- Matrix::summary(as(adjacency_matrix, "dgCMatrix"))
mask_upper <- adj_trip$i < adj_trip$j
ui <- adj_trip$i[mask_upper]
uj <- adj_trip$j[mask_upper]

aadt_vals     <- prepared_traffic_links$lastYearAadt_aadt
dist_sq_aadt  <- (aadt_vals[ui] - aadt_vals[uj])^2
ell_init_aadt <- median(sqrt(dist_sq_aadt))

cat("Neighbour pairs (upper triangle):", length(ui), "\n")
cat("Median AADT distance (ell_init):", round(ell_init_aadt, 1), "\n\n")

# ──────────────────────────────────────────────────────────────────────────────
# 3. RGENERIC: RBF-WEIGHTED BESAG PROPER MODEL
# ──────────────────────────────────────────────────────────────────────────────
# Q(tau, d, ell) = tau * (d*I + D_W - W)
# theta[1] = log(tau), theta[2] = log(d_raw), theta[3] = log(ell)
# d = exp(theta[2]) + 1 ensures d >= 1 → Q guaranteed SPD regardless of W

inla.rgeneric.weighted.besag.RBF <- function(cmd, theta) {
  # Injected by inla.rgeneric.define(): n, ui, uj, dist_sq, ell_init

  graph <- function() {
    Matrix::sparseMatrix(
      i    = c(ui, uj, seq_len(n)),
      j    = c(uj, ui, seq_len(n)),
      x    = 1,
      dims = c(n, n)
    )
  }

  Q <- function() {
    tryCatch({
      tau      <- exp(theta[1])
      d        <- exp(theta[2]) + 1.0          # shift ensures d >= 1
      ell      <- pmax(exp(theta[3]), 1e-6)    # pmax handles NaN safely
      exponent <- pmin(dist_sq / (2 * ell^2), 700)
      sim_vals <- pmax(exp(-exponent), 1e-10)
      W  <- Matrix::sparseMatrix(
        i = c(ui, uj), j = c(uj, ui), x = rep(sim_vals, 2), dims = c(n, n)
      )
      Q0 <- Matrix::Diagonal(x = Matrix::rowSums(W)) - W
      tau * (d * Matrix::Diagonal(n) + Q0)
    }, error = function(e) {
      W_fb <- Matrix::sparseMatrix(
        i = c(ui, uj), j = c(uj, ui), x = rep(1, 2 * length(ui)), dims = c(n, n)
      )
      Matrix::Diagonal(x = Matrix::rowSums(W_fb)) - W_fb + Matrix::Diagonal(n)
    })
  }

  mu             <- function() numeric(n)
  log.norm.const <- function() numeric(0)
  quit           <- function() invisible()

  log.prior <- function() {
    lp_gamma <- function(tk, a = 1, b = 5e-5) a * tk - b * exp(tk)
    lp_ell   <- dnorm(theta[3], mean = log(ell_init), sd = 1.0, log = TRUE)
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

# Instantiate with AADT-only distances. The injected name `dist_sq` is what the
# rgeneric function sees internally; `ell_init` likewise.
weighted_besag_RBF_aadt <- inla.rgeneric.define(
  inla.rgeneric.weighted.besag.RBF,
  n        = n,
  ui       = ui,
  uj       = uj,
  dist_sq  = dist_sq_aadt,
  ell_init = ell_init_aadt
)
cat("rgeneric model defined.\n\n")

# ──────────────────────────────────────────────────────────────────────────────
# 4. 10-FOLD CROSS-VALIDATION
# ──────────────────────────────────────────────────────────────────────────────

K        <- 10
set.seed(42)
measured_all <- which(!is.na(prepared_traffic_links$aadt))
# Balanced randomised fold assignment
fold_ids <- sample(rep_len(seq_len(K), length(measured_all)))

# Shared stabilisation settings for Poisson and NB (see script header for rationale)
stable_inla <- list(
  strategy    = "gaussian",  # Gaussian approx to full conditional — avoids exact
  step.factor = 0.1,         # Poisson Hessian evaluation at extreme η
  diagonal    = 1e-3         # ridge on precision matrix; prevents singular Hessian
)

results_list <- vector("list", K * 3)
list_pos     <- 1L

t_total <- proc.time()

for (k in seq_len(K)) {
  t_fold <- proc.time()
  cat(sprintf("── Fold %2d / %d %s\n", k, K, strrep("─", 48)))

  test_idx  <- measured_all[fold_ids == k]
  true_aadt <- prepared_traffic_links$aadt[test_idx]
  cat(sprintf("   train = %d  |  test = %d\n",
              length(measured_all) - length(test_idx), length(test_idx)))

  # Mask test responses; both log_aadt and aadt_int are NA'd
  data_cv <- prepared_traffic_links
  data_cv$log_aadt[test_idx] <- NA_real_
  data_cv$aadt_int[test_idx] <- NA_integer_

  # ── A: Gaussian on log(AADT) ────────────────────────────────────────────────
  cat("   [A] Gaussian (log)........ ")
  t0 <- proc.time()
  fit_a <- tryCatch(
    inla(
      log_aadt ~ lastYearAadt_logAadt +
        f(spatial.idx, model = weighted_besag_RBF_aadt, n = n, constr = FALSE) +
        f(roadSystem, model = "iid"),
      data              = data_cv,
      family            = "gaussian",
      control.predictor = list(compute = TRUE),
      control.compute   = list(dic = TRUE, waic = TRUE),
      num.threads       = "1:1",
      safe              = TRUE
    ),
    error = function(e) {
      cat("FAILED:", conditionMessage(e), "\n")
      NULL
    }
  )
  cat(sprintf("%.0fs\n", (proc.time() - t0)[3]))

  # ── B: Negative Binomial on rounded AADT counts ─────────────────────────────
  # NB allows Var[y] = mu + mu^2/r > mu, absorbing AADT overdispersion.
  # r (size) is estimated; r→∞ recovers Poisson.
  cat("   [B] Neg. Binomial......... ")
  t0 <- proc.time()
  fit_b <- tryCatch(
    inla(
      aadt_int ~ lastYearAadt_logAadt +
        f(spatial.idx, model = weighted_besag_RBF_aadt, n = n, constr = FALSE) +
        f(roadSystem, model = "iid"),
      data              = data_cv,
      family            = "nbinomial",
      control.predictor = list(compute = TRUE),
      control.compute   = list(dic = TRUE, waic = TRUE),
      control.inla      = stable_inla,
      num.threads       = "1:1",
      safe              = TRUE
    ),
    error = function(e) {
      cat("FAILED:", conditionMessage(e), "\n")
      NULL
    }
  )
  cat(sprintf("%.0fs\n", (proc.time() - t0)[3]))

  # ── C: Poisson on rounded AADT counts (stabilised) ──────────────────────────
  # Poisson assumes Var[y] = mu = E[y]. With AADT ~ 10,000+ this is extremely
  # restrictive and the Laplace approximation can overflow during optimisation.
  # Stabilisation: Gaussian inner-loop approx + small step size + diagonal ridge.
  cat("   [C] Poisson (stabilised).. ")
  t0 <- proc.time()
  fit_c <- tryCatch(
    inla(
      aadt_int ~ lastYearAadt_logAadt +
        f(spatial.idx, model = weighted_besag_RBF_aadt, n = n, constr = FALSE) +
        f(roadSystem, model = "iid"),
      data              = data_cv,
      family            = "poisson",
      control.predictor = list(compute = TRUE),
      control.compute   = list(dic = TRUE, waic = TRUE),
      control.inla      = stable_inla,
      num.threads       = "1:1",
      safe              = TRUE
    ),
    error = function(e) {
      cat("FAILED:", conditionMessage(e), "\n")
      NULL
    }
  )
  cat(sprintf("%.0fs\n", (proc.time() - t0)[3]))

  cat(sprintf("   Fold %d done in %.0fs\n", k, (proc.time() - t_fold)[3]))

  # ── Extract predictions ─────────────────────────────────────────────────────
  # Gaussian:     summary.fitted.values gives E[log(AADT)|...] → exp() for AADT
  # Poisson/NB:   summary.fitted.values gives E[y|...] = exp(η) — already counts
  make_df <- function(fit, model_label, log_scale) {
    if (is.null(fit)) return(NULL)
    fv <- fit$summary.fitted.values[test_idx, ]
    if (log_scale) {
      pred  <- exp(fv$mean)
      ci_lo <- exp(fv[["0.025quant"]])
      ci_hi <- exp(fv[["0.975quant"]])
    } else {
      pred  <- fv$mean
      ci_lo <- fv[["0.025quant"]]
      ci_hi <- fv[["0.975quant"]]
    }
    data.frame(
      model     = model_label,
      fold      = k,
      true_aadt = true_aadt,
      pred      = pred,
      ci_lo     = ci_lo,
      ci_hi     = ci_hi
    )
  }

  results_list[[list_pos]]     <- make_df(fit_a, "Gaussian (log)", log_scale = TRUE)
  results_list[[list_pos + 1]] <- make_df(fit_b, "Neg. Binomial",  log_scale = FALSE)
  results_list[[list_pos + 2]] <- make_df(fit_c, "Poisson",        log_scale = FALSE)
  list_pos <- list_pos + 3L
}

cat(sprintf("\nTotal CV time: %.0f minutes\n", (proc.time() - t_total)[3] / 60))

# Combine all fold results; NULL rows (failed fits) are silently dropped
results_df <- do.call(rbind, Filter(Negate(is.null), results_list))
rownames(results_df) <- NULL

model_levels <- c("Gaussian (log)", "Neg. Binomial", "Poisson")
results_df$model <- factor(results_df$model, levels = model_levels)

cat("\nPredictions per model × fold:\n")
print(with(results_df, table(model, fold)))

# Folds that completely failed for a model
succeeded <- with(results_df, tapply(fold, model, function(f) sort(unique(f))))
for (m in model_levels) {
  failed <- setdiff(seq_len(K), succeeded[[m]])
  if (length(failed) > 0) {
    cat(sprintf("\nWARNING: %s had no predictions for fold(s): %s\n",
                m, paste(failed, collapse = ", ")))
  }
}

# ──────────────────────────────────────────────────────────────────────────────
# 5. COMPUTE METRICS
# ──────────────────────────────────────────────────────────────────────────────
# All primary metrics on the AADT (vehicle/day) scale.
# log-RMSE is scale-invariant and comparable across all three families.
# Coverage: see caveat at top of script.

metrics_fold <- results_df |>
  group_by(model, fold) |>
  summarise(
    n_test   = n(),
    RMSE     = sqrt(mean((pred - true_aadt)^2, na.rm = TRUE)),
    MAE      = mean(abs(pred - true_aadt), na.rm = TRUE),
    MAPE     = mean(abs((pred - true_aadt) / true_aadt) * 100, na.rm = TRUE),
    log_RMSE = sqrt(mean((log(pred) - log(true_aadt))^2, na.rm = TRUE)),
    Coverage = mean(true_aadt >= ci_lo & true_aadt <= ci_hi, na.rm = TRUE) * 100,
    .groups  = "drop"
  )

metrics_summary <- metrics_fold |>
  group_by(model) |>
  summarise(
    across(
      c(RMSE, MAE, MAPE, log_RMSE, Coverage),
      list(mean = mean, sd = sd),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) |>
  left_join(
    results_df |>
      group_by(model) |>
      summarise(
        Correlation = cor(pred, true_aadt, use = "complete.obs"),
        .groups     = "drop"
      ),
    by = "model"
  )

cat("\n")
cat(strrep("=", 82), "\n")
cat("10-FOLD CV SUMMARY  (metrics on AADT scale; log-RMSE is scale-invariant)\n")
cat(strrep("=", 82), "\n\n")
cat(sprintf("%-20s %10s %10s %10s %10s %10s %9s\n",
    "Model", "RMSE", "MAE", "MAPE%", "logRMSE", "Coverage%", "Corr"))
cat(strrep("-", 82), "\n")
for (i in seq_len(nrow(metrics_summary))) {
  r <- metrics_summary[i, ]
  cat(sprintf("%-20s %10.1f %10.1f %10.2f %10.4f %10.1f %9.4f\n",
      as.character(r$model),
      r$RMSE_mean, r$MAE_mean, r$MAPE_mean,
      r$log_RMSE_mean, r$Coverage_mean, r$Correlation))
  cat(sprintf("%-20s %10s %10s  (±%.2f) %10s %10s\n",
      "  ±1 SD",
      sprintf("(±%.1f)", r$RMSE_sd),
      sprintf("(±%.1f)", r$MAE_sd),
      r$MAPE_sd, "", ""))
}
cat(strrep("=", 82), "\n\n")

# ──────────────────────────────────────────────────────────────────────────────
# 6. PLOTS
# ──────────────────────────────────────────────────────────────────────────────

model_colours <- c(
  "Gaussian (log)" = "#1F78B4",   # blue
  "Neg. Binomial"  = "#E6550D",   # orange
  "Poisson"        = "#31A354"    # green
)

# ── P1: Predicted vs Observed (log-log axes) ─────────────────────────────────
# Log-log scale is standard for AADT: multiplicative effects dominate, and the
# distribution spans 2+ orders of magnitude.

p1 <- ggplot(results_df, aes(x = true_aadt, y = pred, colour = model)) +
  geom_point(alpha = 0.30, size = 0.9, shape = 16) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "black", linewidth = 0.65) +
  facet_wrap(~model, nrow = 1) +
  scale_x_log10(labels = scales::comma,
                breaks  = 10^(2:5)) +
  scale_y_log10(labels = scales::comma,
                breaks  = 10^(2:5)) +
  scale_colour_manual(values = model_colours, guide = "none") +
  labs(
    title    = "Predicted vs Observed AADT  —  10-Fold CV",
    subtitle = "Log-log axes; dashed = perfect prediction; each point = one held-out link",
    x        = "Observed AADT (vehicles/day, log scale)",
    y        = "Predicted AADT (vehicles/day, log scale)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    panel.grid.minor = element_blank()
  )

# ── P2: Absolute error distributions (violin + box) ──────────────────────────
errors_df <- results_df |>
  mutate(abs_err = abs(pred - true_aadt))

p2_ylim <- quantile(errors_df$abs_err, 0.99, na.rm = TRUE)

p2 <- ggplot(errors_df, aes(x = model, y = abs_err, fill = model)) +
  geom_violin(trim = TRUE, alpha = 0.75, colour = NA) +
  geom_boxplot(width = 0.16, outlier.shape = NA,
               colour = "grey20", fill = "white", linewidth = 0.5) +
  coord_cartesian(ylim = c(0, p2_ylim)) +
  scale_fill_manual(values = model_colours, guide = "none") +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title    = "Absolute Prediction Error Distribution",
    subtitle = "Trimmed at 99th percentile; box shows median and IQR",
    x        = NULL,
    y        = "| Predicted − Observed | AADT (vehicles/day)"
  ) +
  theme_bw(base_size = 11)

# ── P3: RMSE per fold (line plot) ─────────────────────────────────────────────
p3 <- ggplot(metrics_fold, aes(x = fold, y = RMSE, colour = model, group = model)) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = model_colours) +
  scale_x_continuous(breaks = seq_len(K)) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title  = "RMSE per CV Fold",
    x      = "Fold",
    y      = "RMSE (vehicles/day)",
    colour = "Model"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# ── P4: log-RMSE per fold ─────────────────────────────────────────────────────
# log-RMSE = RMSE on log scale; directly comparable across all three families
# because it does not depend on the response scale (AADT vs log-AADT).
# Interpretation: log-RMSE ≈ 0.10 means predictions are off by a factor ~e^0.1 ≈ 1.10.

p4 <- ggplot(metrics_fold, aes(x = fold, y = log_RMSE, colour = model, group = model)) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = model_colours) +
  scale_x_continuous(breaks = seq_len(K)) +
  labs(
    title    = "Log-RMSE per CV Fold  (scale-invariant)",
    subtitle = "log-RMSE ≈ 0.10 → predictions off by factor ~1.10 on average",
    x        = "Fold",
    y        = "log-RMSE (dimensionless)",
    colour   = "Model"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# ── P5: MAPE per fold ─────────────────────────────────────────────────────────
p5 <- ggplot(metrics_fold, aes(x = fold, y = MAPE, colour = model, group = model)) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = model_colours) +
  scale_x_continuous(breaks = seq_len(K)) +
  labs(
    title  = "MAPE per CV Fold",
    x      = "Fold",
    y      = "MAPE (%)",
    colour = "Model"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# ── P6: Summary metrics bar chart (mean ± 1 SD across folds) ─────────────────
metric_long <- metrics_summary |>
  select(model, RMSE_mean, RMSE_sd, MAE_mean, MAE_sd,
         MAPE_mean, MAPE_sd, log_RMSE_mean, log_RMSE_sd) |>
  pivot_longer(
    -model,
    names_to      = c("metric", ".value"),
    names_pattern = "^(.+)_(mean|sd)$"
  ) |>
  mutate(metric = factor(metric,
    levels = c("RMSE", "MAE", "MAPE", "log_RMSE"),
    labels = c("RMSE (veh/day)", "MAE (veh/day)", "MAPE (%)", "log-RMSE")
  ))

p6 <- ggplot(metric_long, aes(x = model, y = mean, fill = model)) +
  geom_col(alpha = 0.85, colour = "white") +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    width     = 0.28,
    colour    = "grey25",
    linewidth = 0.6
  ) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = model_colours, guide = "none") +
  labs(
    title   = "Mean CV Accuracy Metrics  ±  1 SD across 10 Folds",
    x       = NULL,
    y       = "Mean value",
    caption = "Error bars: ±1 SD across folds."
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 22, hjust = 1),
    strip.background = element_rect(fill = "grey92")
  )

# ── P7: 95% CI coverage by fold ─────────────────────────────────────────────
# Gaussian CIs are for E[log(AADT)] → after exp(), they include σ² (wider).
# Poisson/NB CIs are for the rate parameter E[y]=λ — inherently narrower.
# The dashed line shows nominal 95%.

p7 <- ggplot(metrics_fold, aes(x = model, y = Coverage, fill = model)) +
  geom_violin(trim = TRUE, alpha = 0.75, colour = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA,
               colour = "grey20", fill = "white", linewidth = 0.5) +
  geom_hline(yintercept = 95, linetype = "dashed", colour = "black", linewidth = 0.7) +
  scale_fill_manual(values = model_colours, guide = "none") +
  scale_y_continuous(limits = c(0, 100)) +
  labs(
    title    = "95% CI Coverage by Fold",
    subtitle = paste("Dashed = nominal 95%.",
                     "Gaussian CI includes residual variance (wider);",
                     "Poisson/NB CI is for rate parameter only."),
    x        = NULL,
    y        = "Coverage (%)"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.subtitle = element_text(size = 8))

# ── Print all plots ───────────────────────────────────────────────────────────
for (p in list(p1, p2, p3, p4, p5, p6, p7)) print(p)

cat("\nAll plots rendered. See RStudio plot history to cycle through them.\n")
cat("Session info:\n")
print(sessionInfo())
