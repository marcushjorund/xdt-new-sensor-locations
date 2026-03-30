fit_inla_rbf_model <- function(
    data,                                  #prepared_traffic_links
    adjacency_matrix,                      #adjacency matrix of the traffic links
    spatial_term      = "besagproper_rbf", #either "besagproper" or "besagproper_rbf" with default the latter
    similarity_covariates        = NULL,              #IF spatial_term = "besagproper_rbf" include similarity similarity_covariates
    similarity_interaction_pairs = NULL,              #IF spatial_term = "besagproper_rbf" DEPRECATED: no longer used
    ordinal_levels               = list(),            #named list of ordered levels for ordinal covariates, e.g. list(functionalRoadClass = as.character(c(7,6,5,4,3,2,1,0)), roadCategory = c("KOMMUNAL_VEG","FYLKESVEG","RIKSVEG","EUROPAVEG"))
    fixed_effects     = ~ 1,               #formula for fixed effects in inla model
    iid_effects       = "roadSystem",      #iid random effects, default "roadSystem"
    heavy_vehicle     = FALSE,             #when modelling heavy_vehicle set heavy_vehicle = TRUE
    family            = "gaussian",        #either "gaussian", "poisson" or "nbinomial", default "gaussian"
    verbose           = FALSE)             #if returning verbose output from INLA
  {

  # ────────────────────────────────────────────────────────────────────────────
  # Input validation ----
  # ────────────────────────────────────────────────────────────────────────────

  family       <- match.arg(family,       c("gaussian", "poisson", "nbinomial"))
  spatial_term <- match.arg(spatial_term, c("besagproper", "besagproper_rbf"))

  response <- if (heavy_vehicle) "heavyAadt" else "aadt"

  required_cols <- c("id", response)
  missing_cols  <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0)
    stop("Missing required columns in data: ", paste(missing_cols, collapse = ", "))

  if (!is.null(iid_effects)) {
    missing_iid <- setdiff(iid_effects, names(data))
    if (length(missing_iid) > 0)
      stop("IID effect variables not found in data: ", paste(missing_iid, collapse = ", "))
  }

  if (nrow(data) != nrow(adjacency_matrix) || nrow(data) != ncol(adjacency_matrix))
    stop("Dimensions of adjacency_matrix (", nrow(adjacency_matrix), "x",
         ncol(adjacency_matrix), ") do not match number of rows in data (", nrow(data), ")")

  if (!inherits(fixed_effects, "formula"))
    stop("fixed_effects must be a formula (e.g., ~ lastYearAadt_logAadt + isRamp)")

  # ────────────────────────────────────────────────────────────────────────────
  # Prepare data and response ----
  # ────────────────────────────────────────────────────────────────────────────

  message("Preparing data for INLA model...")

  data_model             <- data
  n                      <- nrow(data_model)
  data_model$spatial.idx <- seq_len(n)

  if (family == "gaussian") {
    # Model log(aadt) with Gaussian likelihood; back-transform predictions to AADT
    response_col               <- paste0("log_", response)
    data_model[[response_col]] <- log(data_model[[response]])
  } else {
    response_col <- response
  }

  # ────────────────────────────────────────────────────────────────────────────
  # Build RBF graph structure ----
  # ────────────────────────────────────────────────────────────────────────────

  if (spatial_term == "besagproper_rbf") {
    message("Building RBF feature matrix and graph structure...")

    # ── Nested helper: Gower-style pairwise distance (adjacent pairs only) ────
    # Each covariate contributes equally (weight 1/p) to the mean squared distance.
    # Numeric : ((x_i - x_j) / range(x))^2                   → [0, 1]
    # Binary / nominal categorical : 1 if unequal, 0 if equal → {0, 1}
    # Ordinal : ((|rank_i - rank_j|) / (K-1))^2              → [0, 1]
    # Result dist_sq is always in (0, 1], so ell is on an interpretable scale.
    compute_pairwise_gower_dist_sq <- function(df, ui, uj,
                                               similarity_covariates,
                                               ordinal_levels = list()) {

      if (is.null(similarity_covariates))
        similarity_covariates <- c("maxLanes", "minLanes",
                                   "highestSpeedLimit", "lowestSpeedLimit",
                                   "lastYearAadt_logAadt",
                                   "hasOnlyPublicTransportLanes", "isRamp",
                                   "roadCategory", "functionClass", "roadSystem")

      missing_cols <- setdiff(similarity_covariates, names(df))
      if (length(missing_cols) > 0)
        stop("similarity_covariates not found in data: ", paste(missing_cols, collapse = ", "))

      detect_type <- function(x) {
        if (is.logical(x)) return("binary")
        if (is.factor(x) || is.character(x)) return("categorical")
        if (is.numeric(x) || is.integer(x)) {
          if (all(unique(x[!is.na(x)]) %in% c(0L, 1L))) return("binary")
          return("numeric")
        }
        "numeric"
      }

      p             <- length(similarity_covariates)
      dist_sq_total <- numeric(length(ui))

      for (col in similarity_covariates) {
        xi <- df[[col]][ui]
        xj <- df[[col]][uj]

        if (col %in% names(ordinal_levels)) {
          # Ordinal: normalised squared rank difference
          lvls <- ordinal_levels[[col]]
          K    <- length(lvls)
          ri   <- match(as.character(xi), as.character(lvls)) - 1L
          rj   <- match(as.character(xj), as.character(lvls)) - 1L
          d    <- (abs(ri - rj) / (K - 1L))^2

        } else if (detect_type(df[[col]]) == "numeric") {
          r <- diff(range(df[[col]], na.rm = TRUE))
          if (r == 0) r <- 1
          d <- ((as.numeric(xi) - as.numeric(xj)) / r)^2

        } else {
          # binary or nominal categorical: 0 if equal, 1 if different
          d <- as.numeric(as.character(xi) != as.character(xj))
        }

        d[is.na(d)]   <- 1   # treat NA as maximally dissimilar
        dist_sq_total <- dist_sq_total + d
      }

      dist_sq_total / p   # Gower mean: dist_sq in (0, 1]
    }
    # ────────────────────────────────────────────────────────────────────────

    adj_trip   <- Matrix::summary(as(adjacency_matrix, "dgCMatrix"))
    mask_upper <- adj_trip$i < adj_trip$j
    ui         <- adj_trip$i[mask_upper]
    uj         <- adj_trip$j[mask_upper]

    dist_sq  <- compute_pairwise_gower_dist_sq(
      df                    = data_model,
      ui                    = ui,
      uj                    = uj,
      similarity_covariates = similarity_covariates,
      ordinal_levels        = ordinal_levels
    )
    ell_init <- median(sqrt(dist_sq[dist_sq > 0]))

    message("  Gower dist_sq range: [", round(min(dist_sq), 4), ", ",
            round(max(dist_sq), 4), "] | ell_init = ", round(ell_init, 4))

    # ── Nested rgeneric model: RBF-weighted Besag proper ─────────────────────
    # Q(tau, d, ell) = tau * (d*I + D_W - W)
    # where W_ij = exp(-dist_sq_ij / (2*ell^2)) for adjacent pairs.
    # theta[1] = log(tau), theta[2] = log(d), theta[3] = log(ell).
    # All data injected via inla.rgeneric.define(): n, ui, uj, dist_sq, ell_init.
    inla.rgeneric.weighted.besag.RBF <- function(cmd, theta) {

      graph <- function() {
        Matrix::sparseMatrix(
          i    = c(ui, uj, seq_len(n)),
          j    = c(uj, ui, seq_len(n)),
          x    = 1,
          dims = c(n, n)
        )
      }

      Q <- function() {
        # tryCatch prevents R errors from crashing INLA's C code (segfault).
        # Fallback returns a safe SPD matrix matching graph() sparsity pattern.
        tryCatch({
          tau      <- exp(theta[1])
          d        <- exp(theta[2]) + 1.0       # shift ≥ 1 keeps Q well-conditioned
          ell      <- pmax(exp(theta[3]), 1e-6)  # pmax handles NaN safely (max does not)
          exponent <- pmin(dist_sq / (2 * ell^2), 700)
          sim_vals <- pmax(exp(-exponent), 1e-10)
          W        <- Matrix::sparseMatrix(i = c(ui, uj), j = c(uj, ui),
                                           x = rep(sim_vals, 2), dims = c(n, n))
          Q0_local <- Matrix::Diagonal(x = Matrix::rowSums(W)) - W
          tau * (d * Matrix::Diagonal(n) + Q0_local)
        }, error = function(e) {
          W_fb <- Matrix::sparseMatrix(i = c(ui, uj), j = c(uj, ui),
                                       x = rep(1, 2 * length(ui)), dims = c(n, n))
          Matrix::Diagonal(x = Matrix::rowSums(W_fb)) - W_fb + Matrix::Diagonal(n)
        })
      }

      mu             <- function() numeric(n)
      log.norm.const <- function() numeric(0)
      quit           <- function() invisible()

      log.prior <- function() {
        # log-Gamma(1, 5e-5) on tau and d; Normal on log(ell) centred at log(ell_init)
        lp_gamma <- function(theta_k, a = 1, b = 5e-5) a * theta_k - b * exp(theta_k)
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
    # ─────────────────────────────────────────────────────────────────────────

    rbf_model <- INLA::inla.rgeneric.define(
      inla.rgeneric.weighted.besag.RBF,
      n        = n,
      ui       = ui,
      uj       = uj,
      dist_sq  = dist_sq,
      ell_init = ell_init
    )
  }

  # ────────────────────────────────────────────────────────────────────────────
  # Build formula ----
  # ────────────────────────────────────────────────────────────────────────────

  formula_full <- stats::as.formula(paste(response_col, "~ 1"))

  if (spatial_term == "besagproper") {
    formula_full <- stats::update(
      formula_full,
      ~ . + f(spatial.idx,
              model               = "besagproper",
              graph               = adjacency_matrix,
              adjust.for.con.comp = FALSE,
              constr              = TRUE)
    )
  } else {
    formula_full <- stats::update(
      formula_full,
      ~ . + f(spatial.idx,
              model  = rbf_model,
              n      = n,
              constr = FALSE)
    )
  }

  if (!is.null(iid_effects) && length(iid_effects) > 0) {
    for (var in iid_effects) {
      formula_full <- stats::update(
        formula_full,
        stats::as.formula(paste('~ . + f(', var, ', model = "iid")')))
    }
  }

  if (length(all.vars(fixed_effects)) > 0) {
    fixed_terms  <- as.character(fixed_effects)[2]
    formula_full <- stats::update(
      formula_full,
      stats::as.formula(paste("~ . +", fixed_terms)))
  }

  # ────────────────────────────────────────────────────────────────────────────
  # Fit INLA model ----
  # ────────────────────────────────────────────────────────────────────────────

  message("Fitting INLA model (family = ", family,
          ", spatial = ", spatial_term, ")...")

  model <- INLA::inla(
    formula_full,
    family            = family,
    data              = data_model,
    control.predictor = list(link = 1),
    control.compute   = list(dic = TRUE, waic = TRUE),
    verbose           = verbose
  )

  message("Model fitting complete.")

  # ────────────────────────────────────────────────────────────────────────────
  # Extract predictions ----
  # ────────────────────────────────────────────────────────────────────────────

  fv     <- model$summary.fitted.values
  suffix <- if (heavy_vehicle) "_heavy" else ""

  if (family == "gaussian") {
    # Fitted values are log(aadt); back-transform median to AADT scale.
    # SD approximated on AADT scale via delta method: sd_aadt ≈ exp(mu_log) * sd_log
    pred_median_log <- fv[, "0.5quant"]
    pred_sd_log     <- fv[, "sd"]
    pred_vals       <- round(exp(pred_median_log))
    sd_vals         <- round(sqrt(exp(pred_sd_log^2) - 1) * exp(pred_median_log + pred_sd_log^2 / 2))
  } else {
    pred_vals <- round(fv[, "0.5quant"])
    sd_vals   <- round(fv[, "sd"])
  }

  predictions <- data.frame(
    id        = data_model$id,
    inla_pred = pred_vals,
    inla_sd   = sd_vals
  ) |>
    stats::setNames(c("id", paste0("inla_pred", suffix), paste0("inla_sd", suffix)))

  # ────────────────────────────────────────────────────────────────────────────
  # Return ----
  # ────────────────────────────────────────────────────────────────────────────

  result <- list(
    predictions           = predictions,
    model_summary         = summary(model),
    inla_model            = model,
    formula               = formula_full,
    family                = family,
    spatial_term          = spatial_term,
    similarity_covariates = similarity_covariates,
    ordinal_levels        = ordinal_levels
  )

  class(result) <- "inla_rbf_traffic_model"
  result
}

# ────────────────────────────────────────────────────────────────────────────────────
# S3 print method ----
# ────────────────────────────────────────────────────────────────────────────────────

print.inla_rbf_traffic_model <- function(x, ...) {
  cat("INLA RBF Traffic Model\n")
  cat("======================\n\n")
  cat("Number of predictions:", nrow(x$predictions), "\n")
  cat("Family:              ", x$family, "\n")
  cat("Spatial term:        ", x$spatial_term, "\n")
  if (!is.null(x$similarity_covariates))
    cat("RBF similarity_covariates: ", paste(x$similarity_covariates, collapse = ", "), "\n")
  else
    cat("RBF similarity_covariates:  <defaults>\n")
  if (length(x$ordinal_levels) > 0)
    cat("Ordinal covariates:        ", paste(names(x$ordinal_levels), collapse = ", "), "\n")
  cat("Formula: ")
  print(x$formula, showEnv = FALSE)
  cat("\n")
  cat("Use $model_summary for model details\n")
  cat("Use $predictions to access predictions data frame\n")
  invisible(x)
}

# ════════════════════════════════════════════════════════════════════════════════
# Metric helpers (private) ----
# ════════════════════════════════════════════════════════════════════════════════

.rmse <- function(p, o) sqrt(mean((p - o)^2, na.rm = TRUE))
.mae  <- function(p, o) mean(abs(p - o),      na.rm = TRUE)
.mape <- function(p, o) {
  nz <- o > 0
  mean(abs((p[nz] - o[nz]) / o[nz]), na.rm = TRUE) * 100
}

.metric_label <- function(pred, obs) {
  sprintf("RMSE = %.0f\nMAE  = %.0f\nMAPE = %.1f%%",
          .rmse(pred, obs), .mae(pred, obs), .mape(pred, obs))
}

# ════════════════════════════════════════════════════════════════════════════════
# plot_inla_model() ----
# ════════════════════════════════════════════════════════════════════════════════

#' Diagnostic plot for a fitted inla_rbf_traffic_model
#'
#' @param model   An \code{inla_rbf_traffic_model} object.
#' @param observed Numeric vector of true AADT values (full-length, same row order
#'   as the data passed to \code{fit_inla_rbf_model}).
#' @param test_idx Optional integer index vector.  When supplied, only these
#'   positions are used (e.g. a held-out test set).
#' @param type    One of \code{"pred_vs_obs"}, \code{"residuals_vs_fitted"},
#'   \code{"qq"}, \code{"residuals_hist"}.  Produces a single ggplot.
#' @param title   Optional title string.  Auto-generated if \code{NULL}.
#' @param col     Hex colour for points.  Defaults to a blue.
#' @return A \code{ggplot} object (printed automatically when not assigned).
#'   Also returns \code{invisible(list(fitted, observed, residuals))}.
plot_inla_model <- function(model,
                            observed,
                            test_idx = NULL,
                            type     = c("pred_vs_obs", "residuals_vs_fitted",
                                         "qq", "residuals_hist"),
                            title    = NULL,
                            col      = "#0066CC") {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required for plot_inla_model(). Install it with install.packages('ggplot2').")

  type <- match.arg(type)

  pred_col <- if (!is.null(model$heavy_vehicle) && isTRUE(model$heavy_vehicle))
    "inla_pred_heavy" else "inla_pred"

  fitted <- model$predictions[[pred_col]]

  if (!is.null(test_idx)) {
    fitted   <- fitted[test_idx]
    observed <- observed[test_idx]
  }

  resid <- fitted - observed

  if (is.null(title)) {
    family_lab <- switch(model$family,
                         gaussian  = "Gaussian",
                         nbinomial = "NB",
                         poisson   = "Poisson",
                         model$family)
    spatial_lab <- switch(model$spatial_term,
                          besagproper     = "Besag",
                          besagproper_rbf = "RBF-Besag",
                          model$spatial_term)
    title <- paste(family_lab, "\u2014", spatial_lab)
  }

  pt_col   <- col
  pt_alpha <- 0.45
  pt_size  <- 1.4

  p <- switch(type,

    pred_vs_obs = {
      lims <- range(c(observed, fitted), na.rm = TRUE)
      df   <- data.frame(obs = observed, pred = fitted)
      ggplot2::ggplot(df, ggplot2::aes(x = obs, y = pred)) +
        ggplot2::geom_abline(slope = 1, intercept = 0,
                             colour = "firebrick", linetype = "dashed", linewidth = 0.8) +
        ggplot2::geom_point(colour = pt_col, alpha = pt_alpha, size = pt_size) +
        ggplot2::coord_fixed(xlim = lims, ylim = lims) +
        ggplot2::annotate("label", x = lims[1], y = lims[2],
                          label  = .metric_label(fitted, observed),
                          hjust  = 0, vjust = 1, size = 3,
                          family = "mono", label.size = 0) +
        ggplot2::labs(x = "Observed AADT", y = "Predicted AADT", title = title) +
        ggplot2::theme_bw(base_size = 11)
    },

    residuals_vs_fitted = {
      df <- data.frame(fitted = fitted, resid = resid)
      ggplot2::ggplot(df, ggplot2::aes(x = fitted, y = resid)) +
        ggplot2::geom_hline(yintercept = 0,
                            colour = "firebrick", linetype = "dashed", linewidth = 0.8) +
        ggplot2::geom_point(colour = pt_col, alpha = pt_alpha, size = pt_size) +
        ggplot2::geom_smooth(method = "loess", formula = y ~ x,
                             se = FALSE, colour = "black", linewidth = 0.9) +
        ggplot2::labs(x = "Fitted AADT", y = "Residual (Fitted \u2212 Observed)",
                      title = title) +
        ggplot2::theme_bw(base_size = 11)
    },

    qq = {
      df <- data.frame(resid = resid[!is.na(resid)])
      ggplot2::ggplot(df, ggplot2::aes(sample = resid)) +
        ggplot2::stat_qq(colour = pt_col, alpha = pt_alpha, size = pt_size) +
        ggplot2::stat_qq_line(colour = "firebrick", linetype = "dashed", linewidth = 0.8) +
        ggplot2::labs(x = "Theoretical quantiles", y = "Sample quantiles",
                      title = paste(title, "\u2014 Q-Q plot")) +
        ggplot2::theme_bw(base_size = 11)
    },

    residuals_hist = {
      df <- data.frame(resid = resid[!is.na(resid)])
      ggplot2::ggplot(df, ggplot2::aes(x = resid)) +
        ggplot2::geom_histogram(fill = pt_col, colour = "white", alpha = 0.8,
                                bins = 30) +
        ggplot2::geom_rug(colour = pt_col, alpha = 0.3) +
        ggplot2::geom_vline(xintercept = 0,
                            colour = "firebrick", linetype = "dashed", linewidth = 0.8) +
        ggplot2::labs(x = "Residual (Fitted \u2212 Observed)", y = "Count",
                      title = paste(title, "\u2014 Residual distribution")) +
        ggplot2::theme_bw(base_size = 11)
    }
  )

  print(p)
  invisible(list(fitted = fitted, observed = observed, residuals = resid))
}

# ════════════════════════════════════════════════════════════════════════════════
# kfold_cv_inla() ----
# ════════════════════════════════════════════════════════════════════════════════

#' K-fold cross-validation for \code{fit_inla_rbf_model}
#'
#' @param data            Data frame passed to \code{fit_inla_rbf_model}.
#' @param adjacency_matrix Square adjacency matrix.
#' @param k               Number of folds (default 5).
#' @param seed            Random seed for fold assignment (default 42).
#' @param fixed_effects   Formula of fixed effects.
#' @param spatial_term    Passed to \code{fit_inla_rbf_model}.
#' @param similarity_covariates Character vector; passed through.
#' @param ordinal_levels  Named list; passed through.
#' @param iid_effects     Character; passed through.
#' @param family          Likelihood family; passed through.
#' @param heavy_vehicle   Logical; passed through.
#' @param verbose         Logical; passed to INLA.
#' @return A list with elements \code{fold_metrics}, \code{summary_metrics},
#'   \code{oof_predictions}, \code{k}, \code{seed}, \code{family},
#'   \code{spatial_term}.
kfold_cv_inla <- function(data,
                          adjacency_matrix,
                          k                    = 5,
                          seed                 = 42,
                          fixed_effects        = ~ 1,
                          spatial_term         = "besagproper_rbf",
                          similarity_covariates = NULL,
                          ordinal_levels       = list(),
                          iid_effects          = "roadSystem",
                          family               = "gaussian",
                          heavy_vehicle        = FALSE,
                          verbose              = FALSE) {

  response     <- if (heavy_vehicle) "heavyAadt" else "aadt"
  pred_col     <- if (heavy_vehicle) "inla_pred_heavy" else "inla_pred"
  measured_idx <- which(!is.na(data[[response]]))

  if (k > length(measured_idx))
    stop("k (", k, ") exceeds number of measured links (", length(measured_idx), ")")

  set.seed(seed)
  fold_id <- sample(rep_len(seq_len(k), length(measured_idx)))

  oof_rows   <- vector("list", k)
  fold_stats <- vector("list", k)

  for (fold in seq_len(k)) {
    message("\n── Fold ", fold, " / ", k, " ──────────────────────────────────────────")

    test_pos  <- measured_idx[fold_id == fold]
    data_fold <- data
    data_fold[[response]][test_pos] <- NA
    if (!heavy_vehicle) data_fold[["heavyAadt"]][test_pos] <- NA

    fit <- fit_inla_rbf_model(
      data                  = data_fold,
      adjacency_matrix      = adjacency_matrix,
      fixed_effects         = fixed_effects,
      spatial_term          = spatial_term,
      similarity_covariates = similarity_covariates,
      ordinal_levels        = ordinal_levels,
      iid_effects           = iid_effects,
      family                = family,
      heavy_vehicle         = heavy_vehicle,
      verbose               = verbose
    )

    pred_vals <- fit$predictions[[pred_col]][test_pos]
    true_vals <- data[[response]][test_pos]

    oof_rows[[fold]] <- data.frame(
      link_idx   = test_pos,
      true_aadt  = true_vals,
      pred_aadt  = pred_vals,
      fold       = fold
    )

    fold_stats[[fold]] <- data.frame(
      fold = fold,
      RMSE = .rmse(pred_vals, true_vals),
      MAE  = .mae(pred_vals,  true_vals),
      MAPE = .mape(pred_vals, true_vals)
    )

    message(sprintf("  Fold %d  RMSE = %.0f  MAE = %.0f  MAPE = %.1f%%",
                    fold,
                    fold_stats[[fold]]$RMSE,
                    fold_stats[[fold]]$MAE,
                    fold_stats[[fold]]$MAPE))
  }

  fold_metrics <- do.call(rbind, fold_stats)
  oof_preds    <- do.call(rbind, oof_rows)

  summary_metrics <- data.frame(
    fold = "mean \u00b1 SD",
    RMSE = sprintf("%.0f \u00b1 %.0f", mean(fold_metrics$RMSE), sd(fold_metrics$RMSE)),
    MAE  = sprintf("%.0f \u00b1 %.0f", mean(fold_metrics$MAE),  sd(fold_metrics$MAE)),
    MAPE = sprintf("%.1f \u00b1 %.1f", mean(fold_metrics$MAPE), sd(fold_metrics$MAPE))
  )

  message("\n── CV complete ──────────────────────────────────────────────────────")
  print(fold_metrics,    digits = 4, row.names = FALSE)
  print(summary_metrics, row.names = FALSE)

  structure(
    list(
      fold_metrics     = fold_metrics,
      summary_metrics  = summary_metrics,
      oof_predictions  = oof_preds,
      k                = k,
      seed             = seed,
      family           = family,
      spatial_term     = spatial_term
    ),
    class = "inla_kfold_cv"
  )
}

# ════════════════════════════════════════════════════════════════════════════════
# plot_kfold_cv() ----
# ════════════════════════════════════════════════════════════════════════════════

#' Diagnostic plots for k-fold CV results
#'
#' @param cv_result An \code{inla_kfold_cv} object from \code{kfold_cv_inla()}.
#' @param type      One of \code{"pred_vs_obs"}, \code{"residuals_vs_fitted"},
#'   \code{"qq"}, \code{"residuals_hist"}, \code{"metrics_by_fold"}.
#' @return A \code{ggplot} object (printed automatically).
plot_kfold_cv <- function(cv_result,
                          type = c("pred_vs_obs", "residuals_vs_fitted",
                                   "qq", "residuals_hist", "metrics_by_fold")) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required. Install with install.packages('ggplot2').")

  type <- match.arg(type)

  oof   <- cv_result$oof_predictions
  oof$fold  <- factor(oof$fold)
  oof$resid <- oof$pred_aadt - oof$true_aadt

  family_lab <- switch(cv_result$family,
                       gaussian  = "Gaussian",
                       nbinomial = "NB",
                       poisson   = "Poisson",
                       cv_result$family)
  spatial_lab <- switch(cv_result$spatial_term,
                        besagproper     = "Besag",
                        besagproper_rbf = "RBF-Besag",
                        cv_result$spatial_term)
  base_title <- paste0(family_lab, " \u2014 ", spatial_lab,
                       "  (", cv_result$k, "-fold CV)")

  p <- switch(type,

    pred_vs_obs = {
      lims <- range(c(oof$true_aadt, oof$pred_aadt), na.rm = TRUE)
      ggplot2::ggplot(oof, ggplot2::aes(x = true_aadt, y = pred_aadt,
                                        colour = fold)) +
        ggplot2::geom_abline(slope = 1, intercept = 0,
                             colour = "firebrick", linetype = "dashed", linewidth = 0.8) +
        ggplot2::geom_point(alpha = 0.45, size = 1.4) +
        ggplot2::coord_fixed(xlim = lims, ylim = lims) +
        ggplot2::annotate("label",
                          x     = lims[1],
                          y     = lims[2],
                          label = .metric_label(oof$pred_aadt, oof$true_aadt),
                          hjust = 0, vjust = 1, size = 3,
                          family = "mono", label.size = 0) +
        ggplot2::labs(x = "Observed AADT", y = "Predicted AADT",
                      colour = "Fold", title = base_title) +
        ggplot2::theme_bw(base_size = 11)
    },

    residuals_vs_fitted = {
      ggplot2::ggplot(oof, ggplot2::aes(x = pred_aadt, y = resid,
                                        colour = fold)) +
        ggplot2::geom_hline(yintercept = 0,
                            colour = "firebrick", linetype = "dashed", linewidth = 0.8) +
        ggplot2::geom_point(alpha = 0.45, size = 1.4) +
        ggplot2::geom_smooth(ggplot2::aes(group = 1),
                             method = "loess", formula = y ~ x,
                             se = FALSE, colour = "black", linewidth = 0.9) +
        ggplot2::labs(x = "Fitted AADT", y = "Residual (Fitted \u2212 Observed)",
                      colour = "Fold", title = base_title) +
        ggplot2::theme_bw(base_size = 11)
    },

    qq = {
      oof_clean <- oof[!is.na(oof$resid), ]
      ggplot2::ggplot(oof_clean, ggplot2::aes(sample = resid, colour = fold)) +
        ggplot2::stat_qq(alpha = 0.45, size = 1.4) +
        ggplot2::stat_qq_line(ggplot2::aes(group = 1),
                              colour = "firebrick", linetype = "dashed",
                              linewidth = 0.8) +
        ggplot2::labs(x = "Theoretical quantiles", y = "Sample quantiles",
                      colour = "Fold",
                      title = paste(base_title, "\u2014 Q-Q plot")) +
        ggplot2::theme_bw(base_size = 11)
    },

    residuals_hist = {
      oof_clean <- oof[!is.na(oof$resid), ]
      ggplot2::ggplot(oof_clean, ggplot2::aes(x = resid, fill = fold)) +
        ggplot2::geom_histogram(colour = "white", alpha = 0.7,
                                bins = 30, position = "identity") +
        ggplot2::geom_vline(xintercept = 0,
                            colour = "firebrick", linetype = "dashed",
                            linewidth = 0.8) +
        ggplot2::labs(x = "Residual (Fitted \u2212 Observed)", y = "Count",
                      fill = "Fold",
                      title = paste(base_title, "\u2014 Residual distribution")) +
        ggplot2::theme_bw(base_size = 11)
    },

    metrics_by_fold = {
      fm      <- cv_result$fold_metrics
      fm$fold <- factor(fm$fold)
      # Pivot to long form without tidyr dependency
      long <- rbind(
        data.frame(fold = fm$fold, metric = "RMSE", value = fm$RMSE),
        data.frame(fold = fm$fold, metric = "MAE",  value = fm$MAE),
        data.frame(fold = fm$fold, metric = "MAPE", value = fm$MAPE)
      )
      long$metric <- factor(long$metric, levels = c("RMSE", "MAE", "MAPE"))

      means <- aggregate(value ~ metric, data = long, FUN = mean)

      ggplot2::ggplot(long, ggplot2::aes(x = fold, y = value, fill = fold)) +
        ggplot2::geom_col(width = 0.65, show.legend = FALSE) +
        ggplot2::geom_hline(data = means,
                            ggplot2::aes(yintercept = value),
                            colour = "firebrick", linetype = "dashed",
                            linewidth = 0.8) +
        ggplot2::facet_wrap(~ metric, scales = "free_y") +
        ggplot2::labs(x = "Fold", y = NULL,
                      title = paste(base_title, "\u2014 Metrics by fold")) +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(strip.background = ggplot2::element_blank(),
                       strip.text = ggplot2::element_text(face = "bold"))
    }
  )

  print(p)
  invisible(p)
}


