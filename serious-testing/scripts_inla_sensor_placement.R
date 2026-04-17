fit_inla_rbf_model <- function(
    data,                                  #prepared_traffic_links
    adjacency_matrix,                      #adjacency matrix of the traffic links
    spatial_term      = "besagproper_rbf", #either "besagproper" or "besagproper_rbf" with default the latter
    similarity_covariates        = NULL,              #IF spatial_term = "besagproper_rbf" include similarity similarity_covariates
    ordinal_levels               = list(),            #named list of ordered levels for ordinal covariates, e.g. list(functionalRoadClass = as.character(c(7,6,5,4,3,2,1,0)), roadCategory = c("KOMMUNAL_VEG","FYLKESVEG","RIKSVEG","EUROPAVEG"))
    distance_type                = "gower",             #"gower" (standard absolute Gower, default) or "gower_squared" (squared per-component)
    weight_fn                    = "laplacian",         #"laplacian" exp(-d/sigma) (default), "gaussian" exp(-d/sigma²), or "linear" (1-d, no length-scale)
    fixed_effects     = ~ 1,               #formula for fixed effects in inla model
    iid_effects       = "roadSystem",      #iid random effects, default "roadSystem"
    heavy_vehicle     = FALSE,             #when modelling heavy_vehicle set heavy_vehicle = TRUE
    family            = "gaussian",        #either "gaussian", "poisson" or "nbinomial", default "gaussian"
    verbose           = FALSE)             #if returning verbose output from INLA
{
  
  # ────────────────────────────────────────────────────────────────────────────
  # Input validation ----
  # ────────────────────────────────────────────────────────────────────────────
  
  family        <- match.arg(family,        c("gaussian", "poisson", "nbinomial"))
  spatial_term  <- match.arg(spatial_term,  c("besagproper", "besagproper_rbf"))
  distance_type <- match.arg(distance_type, c("gower_squared", "gower"))
  weight_fn     <- match.arg(weight_fn,     c("gaussian", "laplacian", "linear"))
  
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
    # distance_type = "gower_squared" (default): squared per-component contributions
    #   Numeric : ((x_i - x_j) / range(x))^2    → [0, 1]
    #   Ordinal : ((|rank_i - rank_j|) / (K-1))^2 → [0, 1]
    # distance_type = "gower": standard absolute Gower
    #   Numeric : |x_i - x_j| / range(x)          → [0, 1]
    #   Ordinal : |rank_i - rank_j| / (K-1)        → [0, 1]
    # Binary / nominal: 0 if equal, 1 if different (same for both types)
    # Result is always in [0, 1].
    compute_pairwise_gower_dist <- function(df, ui, uj,
                                            similarity_covariates,
                                            ordinal_levels = list(),
                                            squared = TRUE) {
      
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
      
      p          <- length(similarity_covariates)
      dist_total <- numeric(length(ui))
      
      for (col in similarity_covariates) {
        xi <- df[[col]][ui]
        xj <- df[[col]][uj]
        
        if (col %in% names(ordinal_levels)) {
          # Ordinal: normalised rank difference (absolute, or squared if distance_type = "gower_squared")
          lvls  <- ordinal_levels[[col]]
          K     <- length(lvls)
          ri    <- match(as.character(xi), as.character(lvls)) - 1L
          rj    <- match(as.character(xj), as.character(lvls)) - 1L
          d_abs <- abs(ri - rj) / (K - 1L)
          d     <- if (squared) d_abs^2 else d_abs
          
        } else if (detect_type(df[[col]]) == "numeric") {
          r     <- diff(range(df[[col]], na.rm = TRUE))
          if (r == 0) r <- 1
          d_abs <- abs(as.numeric(xi) - as.numeric(xj)) / r
          d     <- if (squared) d_abs^2 else d_abs
          
        } else {
          # binary or nominal categorical: 0 if equal, 1 if different (same for both)
          d <- as.numeric(as.character(xi) != as.character(xj))
        }
        
        d[is.na(d)] <- 1   # treat NA as maximally dissimilar
        dist_total  <- dist_total + d
      }
      
      dist_total / p   # mean per-component distance in [0, 1]
    }
    # ────────────────────────────────────────────────────────────────────────
    
    adj_trip   <- Matrix::summary(as(adjacency_matrix, "dgCMatrix"))
    mask_upper <- adj_trip$i < adj_trip$j
    ui         <- adj_trip$i[mask_upper]
    uj         <- adj_trip$j[mask_upper]
    
    dist_vals  <- compute_pairwise_gower_dist(
      df                    = data_model,
      ui                    = ui,
      uj                    = uj,
      similarity_covariates = similarity_covariates,
      ordinal_levels        = ordinal_levels,
      squared               = (distance_type == "gower_squared")
    )
    # sigma_init: median non-zero distance — used as length-scale prior centre
    # and initial value for theta[3] (laplacian/gaussian). Not used for "linear".
    # Falls back to 0.5 if all distances are zero (e.g. binary covariate constant
    # across all neighbours), preventing log(NA) from crashing INLA.
    nz_dists   <- dist_vals[dist_vals > 0]
    sigma_init <- if (length(nz_dists) > 0) median(nz_dists) else 0.5
    
    message("  Gower dist range: [", round(min(dist_vals), 4), ", ",
            round(max(dist_vals), 4), "] | sigma_init = ", round(sigma_init, 4),
            " | distance_type = ", distance_type, " | weight_fn = ", weight_fn)
    
    # ── Nested rgeneric model: similarity-weighted Besag proper ─────────────
    # Q(tau, d, [scale]) = tau * (d*I + D_W - W)
    # W_ij depends on weight_fn (all dist_vals in [0,1]):
    #   "gaussian"  : exp(-dist_vals_ij / sigma^2)      theta[3] = log(sigma)
    #   "laplacian" : exp(-dist_vals_ij / sigma)         theta[3] = log(sigma)
    #   "linear"    : max(1 - dist_vals_ij, 0)           (no theta[3])
    # theta[1] = log(tau), theta[2] = log(d-1), [theta[3] = log(scale)].
    # Injected via inla.rgeneric.define(): n, ui, uj, dist_vals, sigma_init, weight_fn.
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
          tau <- exp(theta[1])
          d   <- exp(theta[2]) + 1.0        # shift ≥ 1 keeps Q well-conditioned
          
          # Compute edge weights from dist_vals using the chosen kernel
          sim_vals <- if (weight_fn == "linear") {
            # Parameter-free: w = max(1 - d, 0). Larger distance → less similarity.
            pmax(1 - dist_vals, 0)
          } else if (weight_fn == "laplacian") {
            # Exponential / OU kernel: sharper decay, better for discrete-heavy data.
            sigma <- pmax(exp(theta[3]), 1e-6)
            pmax(exp(-pmin(dist_vals / sigma, 700)), 1e-10)
          } else {
            # Gaussian (RBF) kernel: smooth decay.
            sigma <- pmax(exp(theta[3]), 1e-6)
            pmax(exp(-pmin(dist_vals / sigma^2, 700)), 1e-10)
          }
          
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
        # log-Gamma(1, 5e-5) on tau and d.
        # For laplacian/gaussian: Normal on log(scale) centred at log(sigma_init).
        # For linear: no scale parameter, prior is only on tau and d.
        lp_gamma <- function(theta_k, a = 1, b = 5e-5) a * theta_k - b * exp(theta_k)
        base_lp  <- lp_gamma(theta[1]) + lp_gamma(theta[2])
        if (weight_fn == "linear") base_lp
        else base_lp + dnorm(theta[3], mean = log(sigma_init), sd = 1.0, log = TRUE)
      }
      
      initial <- function() {
        # "linear" has 2 hyperparameters; laplacian/gaussian have 3.
        if (weight_fn == "linear") c(0, 0)
        else c(0, 0, log(sigma_init))
      }
      
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
      n          = n,
      ui         = ui,
      uj         = uj,
      dist_vals  = dist_vals,   # distances in [0,1]; scale depends on distance_type
      sigma_init = sigma_init,  # prior centre for log(scale); unused for "linear"
      weight_fn  = weight_fn    # controls kernel shape inside Q()
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
  # Transform spatial hyperparameters to natural scale ----
  # theta[1] = log(tau)    → tau   = exp(theta[1])
  # theta[2] = log(d - 1)  → d     = exp(theta[2]) + 1
  # theta[3] = log(sigma)  → sigma = exp(theta[3])   [laplacian/gaussian only]
  # Quantiles and mode transform exactly under monotone functions.
  # Mean and sd use the delta method: sd_nat ≈ |f'(mu_log)| * sd_log.
  # ────────────────────────────────────────────────────────────────────────────
  
  spatial_hyperpar <- NULL
  if (spatial_term == "besagproper_rbf") {
    .bt_row <- function(hp, idx, fn, label) {
      mu <- hp[idx, "mean"]
      s  <- hp[idx, "sd"]
      for (col in intersect(c("0.025quant", "0.5quant", "0.975quant", "mode"), names(hp)))
        hp[idx, col] <- fn(hp[idx, col])
      hp[idx, "mean"] <- fn(mu)
      hp[idx, "sd"]   <- abs((fn(mu + 1e-7) - fn(mu - 1e-7)) / 2e-7) * s
      rownames(hp)[idx] <- label
      hp
    }
    
    sp <- model$summary.hyperpar
    rn <- rownames(sp)
    
    i1 <- grep("Theta1 for spatial\\.idx", rn)
    i2 <- grep("Theta2 for spatial\\.idx", rn)
    i3 <- grep("Theta3 for spatial\\.idx", rn)
    
    if (length(i1) == 1)
      sp <- .bt_row(sp, i1, exp,                       "tau (precision) for spatial.idx")
    if (length(i2) == 1)
      sp <- .bt_row(sp, i2, function(x) exp(x) + 1,   "d (diagonal offset) for spatial.idx")
    if (length(i3) == 1)
      sp <- .bt_row(sp, i3, exp,                       "sigma (length-scale) for spatial.idx")
    
    spatial_hyperpar <- sp
  }
  
  # ────────────────────────────────────────────────────────────────────────────
  # Return ----
  # ────────────────────────────────────────────────────────────────────────────
  
  result <- list(
    predictions           = predictions,
    model_summary         = summary(model),
    spatial_hyperpar      = spatial_hyperpar,
    distances             = dist_vals,
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
.male <- function(p, o) {
  ok <- p > 0 & o > 0
  median(abs(log(p[ok]) - log(o[ok])), na.rm = TRUE)
}

.metric_label <- function(pred, obs) {
  sprintf("RMSE = %.0f\nMAE  = %.0f\nMAPE = %.1f%%\nMALE = %.3f",
          .rmse(pred, obs), .mae(pred, obs), .mape(pred, obs), .male(pred, obs))
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
                                    family = "mono") +
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
                          verbose              = FALSE,
                          weight_fn = "laplacian", 
                          distance_type = "gower") {
  
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
      verbose               = verbose,
      weight_fn = weight_fn,
      distance_type = distance_type
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
      MAPE = .mape(pred_vals, true_vals),
      MALE = .male(pred_vals, true_vals)
    )
    
    message(sprintf("  Fold %d  RMSE = %.0f  MAE = %.0f  MAPE = %.1f%%  MALE = %.3f",
                    fold,
                    fold_stats[[fold]]$RMSE,
                    fold_stats[[fold]]$MAE,
                    fold_stats[[fold]]$MAPE,
                    fold_stats[[fold]]$MALE))
  }
  
  fold_metrics <- do.call(rbind, fold_stats)
  oof_preds    <- do.call(rbind, oof_rows)
  
  summary_metrics <- data.frame(
    fold = "mean \u00b1 SD",
    RMSE = sprintf("%.0f \u00b1 %.0f", mean(fold_metrics$RMSE), sd(fold_metrics$RMSE)),
    MAE  = sprintf("%.0f \u00b1 %.0f", mean(fold_metrics$MAE),  sd(fold_metrics$MAE)),
    MAPE = sprintf("%.1f \u00b1 %.1f", mean(fold_metrics$MAPE), sd(fold_metrics$MAPE)),
    MALE = sprintf("%.3f \u00b1 %.3f", mean(fold_metrics$MALE), sd(fold_metrics$MALE))
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
                
                # Per-fold metrics for legend labels
                fold_labs <- vapply(levels(oof$fold), function(f) {
                  sub <- oof[oof$fold == f, ]
                  sprintf("Fold %s  RMSE=%.0f  MAE=%.0f  MAPE=%.1f%%  MALE=%.3f",
                          f,
                          .rmse(sub$pred_aadt, sub$true_aadt),
                          .mae( sub$pred_aadt, sub$true_aadt),
                          .mape(sub$pred_aadt, sub$true_aadt),
                          .male(sub$pred_aadt, sub$true_aadt))
                }, character(1))
                
                # Mean ± SD across folds for annotation
                fm    <- cv_result$fold_metrics
                annot <- sprintf("Mean \u00b1 SD\nRMSE = %.0f \u00b1 %.0f\nMAE  = %.0f \u00b1 %.0f\nMAPE = %.1f \u00b1 %.1f%%\nMALE = %.3f \u00b1 %.3f",
                                 mean(fm$RMSE), sd(fm$RMSE),
                                 mean(fm$MAE),  sd(fm$MAE),
                                 mean(fm$MAPE), sd(fm$MAPE),
                                 mean(fm$MALE), sd(fm$MALE))
                
                ggplot2::ggplot(oof, ggplot2::aes(x = true_aadt, y = pred_aadt,
                                                  colour = fold)) +
                  ggplot2::geom_abline(slope = 1, intercept = 0,
                                       colour = "firebrick", linetype = "dashed", linewidth = 0.8) +
                  ggplot2::geom_point(alpha = 0.45, size = 1.4) +
                  ggplot2::coord_fixed(xlim = lims, ylim = lims) +
                  ggplot2::annotate("label", x = lims[1], y = lims[2],
                                    label = annot, hjust = 0, vjust = 1, size = 3,
                                    family = "mono") +
                  ggplot2::scale_colour_discrete(labels = fold_labs) +
                  ggplot2::labs(x = "Observed AADT", y = "Predicted AADT",
                                colour = NULL, title = base_title) +
                  ggplot2::theme_bw(base_size = 11) +
                  ggplot2::theme(legend.text = ggplot2::element_text(family = "mono", size = 8))
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
                  data.frame(fold = fm$fold, metric = "MAPE", value = fm$MAPE),
                  data.frame(fold = fm$fold, metric = "MALE", value = fm$MALE)
                )
                long$metric <- factor(long$metric, levels = c("RMSE", "MAE", "MAPE", "MALE"))
                
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

# ════════════════════════════════════════════════════════════════════════════════
# select_similarity_covariates() ----
# ════════════════════════════════════════════════════════════════════════════════

#' Forward selection of optimal similarity_covariates for the RBF model.
#'
#' Greedy forward selection: starts from an empty set and at each step adds
#' the candidate that most reduces mean CV-RMSE across folds.  Stops when no
#' addition improves performance.  Uses \code{kfold_cv_inla()} internally.
#' Ordinal levels are automatically subsetted to the active covariate set at
#' each step so irrelevant specs never affect the Gower distance computation.
#'
#' @param data                Data frame passed to \code{fit_inla_rbf_model}.
#' @param adjacency_matrix    Square adjacency matrix.
#' @param candidate_covariates Character vector of all covariates to consider.
#' @param ordinal_levels      Named list of ordered levels for ordinal
#'   covariates; subsetted automatically per step.
#' @param fixed_effects       Formula; held constant throughout the search.
#' @param iid_effects         Character; passed through.
#' @param family              Likelihood family; passed through.
#' @param distance_type       Passed to \code{fit_inla_rbf_model}.
#' @param weight_fn           Passed to \code{fit_inla_rbf_model}.
#' @param k                   Number of CV folds (default 10).
#' @param seed                RNG seed for fold assignment (default 420).
#' @param verbose             Passed to INLA.
#' @return A list: \code{best_covariates}, \code{best_ordinal_levels},
#'   \code{best_male}, \code{selection_history}.
select_similarity_covariates <- function(
    data,
    adjacency_matrix,
    candidate_covariates = c(
      "maxLanes", "minLanes",
      "highestSpeedLimit", "lowestSpeedLimit", "isRamp",
      "roadCategory", "functionClass", "functionalRoadClass",
      "roadSystem", "lastYearAadt_logAadt"
    ),
    ordinal_levels = list(
      functionalRoadClass = as.character(c(7, 6, 5, 4, 3, 2, 1, 0)),
      roadCategory        = c("KOMMUNAL_VEG", "FYLKESVEG", "RIKSVEG", "EUROPAVEG"),
      functionClass       = c("unknown", "E", "D", "C", "B", "A"),
      highestSpeedLimit   = c("unknown", "30", "40", "50", "60", "70",
                              "80", "90", "100", "110"),
      lowestSpeedLimit    = c("unknown", "20", "30", "40", "50", "60",
                              "70", "80", "90", "100")
    ),
    fixed_effects = ~ functionalRoadClass:maxLanes +
      functionalRoadClass:roadCategory +
      minLanes:roadCategory + functionalRoadClass +
      maxLanes + roadCategory +
      hasOnlyPublicTransportLanes +
      functionalRoadClass * isRamp +
      lastYearAadt_logAadt,
    iid_effects   = "roadSystem",
    family        = "nbinomial",
    distance_type = "gower",
    weight_fn     = "laplacian",
    k             = 10,
    seed          = 420,
    verbose       = FALSE) {
  
  selected   <- character(0)
  remaining  <- candidate_covariates
  best_male  <- Inf
  history    <- data.frame(step       = integer(0),
                           added      = character(0),
                           mean_male  = numeric(0),
                           sd_male    = numeric(0),
                           covariates = character(0),
                           stringsAsFactors = FALSE)
  
  message("\n\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550")
  message("  Forward selection of similarity_covariates")
  message("  Candidates : ", paste(candidate_covariates, collapse = ", "))
  message("  Family     : ", family, "  |  k = ", k, "  |  seed = ", seed)
  message("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n")
  
  for (step in seq_along(candidate_covariates)) {
    if (length(remaining) == 0) break
    
    # ── Evaluate adding each remaining candidate ──────────────────────────────
    # Use lapply to capture the full cv object (needed for sd later).
    step_cvs <- lapply(remaining, function(cov) {
      trial_set  <- c(selected, cov)
      active_ord <- ordinal_levels[intersect(names(ordinal_levels), trial_set)]
      message(sprintf("  Step %d | trying +%-30s (set size %d)",
                      step, cov, length(trial_set)))
      kfold_cv_inla(
        data                  = data,
        adjacency_matrix      = adjacency_matrix,
        k                     = k,
        seed                  = seed,
        fixed_effects         = fixed_effects,
        spatial_term          = "besagproper_rbf",
        similarity_covariates = trial_set,
        ordinal_levels        = active_ord,
        iid_effects           = iid_effects,
        family                = family,
        distance_type         = distance_type,
        weight_fn             = weight_fn,
        verbose               = verbose
      )
    })
    names(step_cvs) <- remaining
    
    step_means  <- vapply(step_cvs, function(cv) mean(cv$fold_metrics$MALE), numeric(1))
    best_cov      <- names(step_means)[which.min(step_means)]
    best_male_now <- step_means[[best_cov]]
    best_sd_now   <- sd(step_cvs[[best_cov]]$fold_metrics$MALE)
    
    message(sprintf("\n  Step %d best: +%s  mean MALE = %.4f (sd = %.4f)  previous best = %.4f\n",
                    step, best_cov, best_male_now, best_sd_now, best_male))
    
    if (best_male_now >= best_male) {
      message("  No improvement \u2014 stopping forward selection.")
      break
    }
    
    # ── Accept the addition ───────────────────────────────────────────────────
    selected   <- c(selected, best_cov)
    remaining  <- setdiff(remaining, best_cov)
    best_male  <- best_male_now
    
    history <- rbind(history, data.frame(
      step       = step,
      added      = best_cov,
      mean_male  = best_male_now,
      sd_male    = best_sd_now,
      covariates = paste(selected, collapse = ", "),
      stringsAsFactors = FALSE
    ))
  }
  
  best_ordinal <- ordinal_levels[intersect(names(ordinal_levels), selected)]
  
  message("\n\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550")
  message("  Selection complete")
  message("  Optimal covariates : ", paste(selected, collapse = ", "))
  message(sprintf("  Best mean CV-MALE  : %.4f", best_male))
  message("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n")
  cat("\nSelection history:\n")
  print(history[, c("step", "added", "mean_male", "sd_male")],
        digits = 4, row.names = FALSE)
  
  invisible(list(
    best_covariates     = selected,
    best_ordinal_levels = best_ordinal,
    best_male           = best_male,
    selection_history   = history
  ))
}


#Creating functions/helper-functions for the near-optimal sensor selection algoritm
#The algorithm in Krause et al. (2008) needs a covariance matrix in order to compute the near-optimal sensor locations

create_covariance_and_precision_matrix <- function(
    adjacency_matrix,
    tau,
    d,
    sigma,
    distances,   # dist_vals vector from fit_inla_rbf_model()$distances
    data = NULL,
    with_and_against = FALSE
) {
  n <- nrow(adjacency_matrix)
  
  # Reconstruct upper-triangle edge indices (same logic as fit_inla_rbf_model)
  adj_trip   <- Matrix::summary(as(adjacency_matrix, "dgCMatrix"))
  mask_upper <- adj_trip$i < adj_trip$j
  ui         <- adj_trip$i[mask_upper]
  uj         <- adj_trip$j[mask_upper]
  
  stopifnot(length(distances) == length(ui))
  
  # Symmetric sparse distance matrix: 0 on diagonal, d_ij on edges
  dist_matrix <- Matrix::sparseMatrix(
    i    = c(ui, uj),
    j    = c(uj, ui),
    x    = rep(distances, 2),
    dims = c(n, n)
  )
  
  # Edge weights: exp(-d_ij / sigma), applied only to stored entries
  W      <- dist_matrix
  #laplacian kernel
  W@x    <- exp(-W@x / sigma)
  
  # Precision matrix: Q = tau * (d*I + D_W - W)
  D_W    <- Matrix::Diagonal(x = Matrix::rowSums(W))
  precision_matrix      <- tau * (d * Matrix::Diagonal(n) + D_W - W)
  
  # Covariance matrix: Sigma = Q^{-1}
  # chol2inv() on a sparse positive definite Q returns a dense matrix
  covariance_matrix  <- chol2inv(chol(precision_matrix))
  if(with_and_against){
    #check whether data != NULL
    assertthat::assert_that(!is.null(data), msg = "data must be provided when with_and_against = TRUE")
    with_idx    <- which(grepl("-WITH$",    data$id))
    against_idx <- which(grepl("-AGAINST$", data$id))
    
    
    # Strip the suffix to get base IDs
    base_with    <- sub("-WITH$",    "", data$id[with_idx])
    base_against <- sub("-AGAINST$", "", data$id[against_idx])
    
    # Find paired base IDs
    paired_bases <- intersect(base_with, base_against)
    unpaired_bases <- setdiff(base_with, base_against)
    
    # Get the original row indices in prepared_traffic_links for each group
    paired_with_idx    <- match(paste0(paired_bases,    "-WITH"),    data$id)
    paired_against_idx <- match(paste0(paired_bases,    "-AGAINST"), data$id)
    unpaired_with_idx  <- match(paste0(unpaired_bases,  "-WITH"),    data$id)
    
    n <- nrow(data)
    m <- length(paired_bases) + length(unpaired_bases)
    
    # Build sparse C matrix
    row_idx <- c(
      seq_along(paired_bases),                          
      seq_along(paired_bases),                          
      length(paired_bases) + seq_along(unpaired_bases)  
    )
    col_idx <- c(
      paired_with_idx,
      paired_against_idx,
      unpaired_with_idx
    )
    values <- rep(1, length(row_idx))
    
    C <- sparseMatrix(i = row_idx, j = col_idx, x = values,
                      dims = c(m, n))
    #adding with + against by linear transformation of covariance_matrix by C
    covariance_matrix_sum <- C %*% covariance_matrix %*% t(C)
    return(list(covariance_matrix_sum = covariance_matrix_sum,
                ids = c(paired_bases, unpaired_bases),
                precision_matrix = precision_matrix,
                covariance_matrix = covariance_matrix,
                weights = W))
  }
  else{
    return(list(
      precision_matrix  = precision_matrix,
      covariance_matrix = covariance_matrix,
      weights           = W
    ))}
}

#WITH + AGAINST transformation
"When placing a traffic sensor in the road, the sensors are almost without exception placed such that 
it measures traffic volume going in both directions. This could either be implemented as selecting a pair of
locations, or we can simplify this and say that we are selecting the sum of them. In order to be able
to compare one-way-driven road and those having two"
create_covariance_and_precision_matrix_both_directions <- function(data, covariance_matrix){
  with_idx    <- which(grepl("-WITH$",    data$id))
  against_idx <- which(grepl("-AGAINST$", data$id))
  
  
  # Strip the suffix to get base IDs
  base_with    <- sub("-WITH$",    "", data$id[with_idx])
  base_against <- sub("-AGAINST$", "", data$id[against_idx])
  
  # Find paired base IDs
  paired_bases <- intersect(base_with, base_against)
  unpaired_bases <- setdiff(base_with, base_against)
  
  # Get the original row indices in prepared_traffic_links for each group
  paired_with_idx    <- match(paste0(paired_bases,    "-WITH"),    data$id)
  paired_against_idx <- match(paste0(paired_bases,    "-AGAINST"), data$id)
  unpaired_with_idx  <- match(paste0(unpaired_bases,  "-WITH"),    data$id)
  
  n <- nrow(data)
  m <- length(paired_bases) + length(unpaired_bases)
  
  # Build sparse C matrix
  row_idx <- c(
    seq_along(paired_bases),                          
    seq_along(paired_bases),                          
    length(paired_bases) + seq_along(unpaired_bases)  
  )
  col_idx <- c(
    paired_with_idx,
    paired_against_idx,
    unpaired_with_idx
  )
  values <- rep(1, length(row_idx))
  
  C <- sparseMatrix(i = row_idx, j = col_idx, x = values,
                    dims = c(m, n))
  #adding with + against by linear transformation of covariance_matrix by C
  covariance_matrix_sum <- C %*% covariance_matrix %*% t(C)
  return(list(covariance_matrix_sum = covariance_matrix_sum,
              ids                   = c(paired_bases, unpaired_bases)))
}

# ════════════════════════════════════════════════════════════════════════════════
# greedy_mi_sensor_selection() ----
# ════════════════════════════════════════════════════════════════════════════════

greedy_mi_sensor_selection <- function(data, covariance_matrix, ids, k, adjacency_matrix = NULL,
                                       weighting_bias = NULL) {
  m <- length(ids)
  stopifnot(
    k >= 1,
    nrow(covariance_matrix) == m,
    ncol(covariance_matrix) == m
  )
  
  # Resolve measured/unmeasured in m-space via the -WITH link
  with_row       <- match(paste0(ids, "-WITH"), data$id)
  measured_idx   <- which(!is.na(data$aadt[with_row]))
  unmeasured_idx <- which( is.na(data$aadt[with_row]))
  
  if (length(unmeasured_idx) == 0) stop("No unmeasured locations found.")
  k <- min(k, length(unmeasured_idx))
  
  # ── Warm-start: condition on existing sensors S0 via direct Schur complement ─
  # Sequential left-looking Cholesky over ~300+ sensors accumulates enough
  # float error that pivots go negative → sqrt(NaN) → all cond_var1 = NaN.
  # Instead, solve in one shot: if Sigma = cov(all, all), and S0 = measured,
  #   Var(x_i | S0) = Sigma_ii - Sigma_{i,S0} Sigma_{S0,S0}^{-1} Sigma_{S0,i}
  #                 = diag(Sigma) - colSums(A^2)   where A = L_{S0}^{-T} Sigma_{S0, .}
  # Uses only a single Cholesky of the (n_s0 × n_s0) measured block.
  cond_var1 <- diag(covariance_matrix)
  if (length(measured_idx) > 0) {
    Sigma_mm  <- covariance_matrix[measured_idx, measured_idx]
    Sigma_om  <- covariance_matrix[, measured_idx]               # m × n_s0
    L_mm      <- chol(Sigma_mm)                                  # upper triangular
    # Solve L_mm^T A = Sigma_om^T  →  A is n_s0 × m
    A         <- forwardsolve(t(L_mm), t(Sigma_om))
    cond_var1 <- pmax(cond_var1 - colSums(A^2), 0)
  } else {
    A <- matrix(0, 0, m)   # empty; no warm-start needed
  }
  # Zero out measured locations so they never win the score competition
  cond_var1[measured_idx] <- 0
  
  # ── Initialise precision of unmeasured candidates ───────────────────────────
  Sigma_uu  <- covariance_matrix[unmeasured_idx, unmeasured_idx]
  prec      <- chol2inv(chol(Sigma_uu))   # O(|uu|^3), once
  
  # cond_var2[i] = Var(x_i | all other candidates) = 1 / diag(prec_uu)
  # Use 0 (not -1) as sentinel so score = cond_var1 * cond_var2 is never
  # spuriously positive for excluded locations.
  cond_var2 <- numeric(m)
  cond_var2[unmeasured_idx] <- 1 / diag(prec)
  # Build location weights (1 = no bias)
  w <- rep(1.0, m)
  if (!is.null(weighting_bias)) {
    frc   <- suppressWarnings(as.integer(as.character(data$functionalRoadClass[with_row])))
    w_idx <- frc + 1L
    valid <- !is.na(w_idx) & w_idx >= 1L & w_idx <= length(weighting_bias)
    w[valid] <- weighting_bias[w_idx[valid]]
  }
  # ── Greedy MI loop ───────────────────────────────────────────────────────────
  # Each selected sensor updates cond_var1 via an incremental Cholesky column
  # on the covariance ALREADY conditioned on S0 (using A from warm-start).
  candidates   <- unmeasured_idx
  selected_idx <- integer(k)
  mi_scores    <- numeric(k)
  L2           <- matrix(0, m, k)   # Cholesky columns for greedy selections only
  
  for (i in seq_len(k)) {
    score       <- cond_var1 * cond_var2 * w
    best_global <- which.max(score)
    j           <- which(candidates == best_global)
    
    selected_idx[i] <- best_global
    mi_scores[i]    <- score[best_global]
    
    # Cholesky column for best_global in the conditional covariance given S0
    # col_vec = Sigma_{., s}  −  Sigma_{., S0} Sigma_{S0,S0}^{-1} Sigma_{S0, s}
    #         = raw column  −  A^T A[, s]
    col_vec <- covariance_matrix[, best_global] - crossprod(A, A[, best_global])
    # Subtract contributions from previously greedy-selected sensors
    if (i > 1)
      col_vec <- col_vec - L2[, seq_len(i - 1), drop = FALSE] %*%
      L2[best_global, seq_len(i - 1)]
    pivot       <- pmax(col_vec[best_global], .Machine$double.eps * 100)
    L2[, i]     <- col_vec / sqrt(pivot)
    cond_var1   <- pmax(cond_var1 - L2[, i]^2, 0)
    cond_var1[best_global] <- 0   # selected: exclude from future scores
    
    # Schur-complement deletion of candidate j from precision of remaining set
    if (length(candidates) > 1) {
      prec <- prec[-j, -j, drop = FALSE] -
        tcrossprod(prec[-j, j]) / prec[j, j]
    }
    candidates <- candidates[-j]
    if (length(candidates) > 0)
      cond_var2[candidates] <- 1 / diag(prec)
    cond_var2[best_global] <- 0   # selected: exclude from future scores
  }
  selected_base_ids <- ids[selected_idx]
  selected_data_entries = data[data$id %in% c(paste0(selected_base_ids, "-WITH"),
                                              paste0(selected_base_ids, "-AGAINST")), ]
  selected_data_entries$selected <- rep(TRUE, nrow(selected_data_entries))
  # Map mi_score to both -WITH and -AGAINST rows for each selected base ID
  score_lookup <- setNames(mi_scores, selected_base_ids)
  selected_data_entries$mi_score <- score_lookup[
    sub("-WITH$|-AGAINST$", "", selected_data_entries$id)
  ]
  if (!is.null(adjacency_matrix)) {
    with_rows_selected    <- na.omit(match(paste0(selected_base_ids, "-WITH"),    data$id))
    against_rows_selected <- na.omit(match(paste0(selected_base_ids, "-AGAINST"), data$id))
    selected_dir_rows     <- c(with_rows_selected, against_rows_selected)
    adj_sub               <- adjacency_matrix[selected_dir_rows, , drop = FALSE]
    neigh_data_idx        <- setdiff(which(Matrix::colSums(adj_sub != 0) > 0),
                                     which(data$id %in% selected_data_entries$id))
    if (length(neigh_data_idx) > 0) {
      neigh_rows          <- data[neigh_data_idx, ]
      neigh_rows$selected <- FALSE
      neigh_rows$mi_score <- NA_real_
      selected_data_entries <- rbind(selected_data_entries, neigh_rows)
    }
  }
  list(
    selected_ids           = ids[selected_idx],
    selected_idx           = selected_idx,
    selected_data_entries  = selected_data_entries,
    mi_scores              = mi_scores,
    measured_idx           = measured_idx,
    unmeasured_idx_initial = unmeasured_idx
  )
}

#County partitioning 

expand_border <- function(seed_idx, adjacency_matrix, hops = 2) {
  current <- seed_idx
  for (i in seq_len(hops)) {
    adj_sub <- adjacency_matrix[current, , drop = FALSE]
    new_idx <- which(Matrix::colSums(adj_sub != 0) > 0)
    current <- union(current, new_idx)
  }
  setdiff(current, seed_idx)
}

#' Partition a nationwide traffic-link dataset into per-county sub-datasets,
#' each widened by \code{hops} adjacency-matrix steps across county borders.
#'
#' @param data             Data frame with a \code{county} factor column.
#' @param adjacency_matrix Sparse square adjacency matrix matching \code{nrow(data)}.
#' @param hops             Number of hops to expand beyond the county boundary (default 2).
#' @return Named list (one element per county level), each containing:
#'   \code{$data} (subset of \code{data}) and
#'   \code{$adjacency_matrix} (subsetted and reindexed).
#'   Core rows are identified by \code{$data$county == county_name}.
partition_by_county <- function(data, adjacency_matrix, distances, hops = 2) {
  # Pre-compute upper-triangle edge indices from the full adjacency matrix once.
  # distances[k] corresponds to edge (ui[k], uj[k]); same logic as
  # create_covariance_and_precision_matrix().
  adj_trip   <- Matrix::summary(as(adjacency_matrix, "dgCMatrix"))
  mask_upper <- adj_trip$i < adj_trip$j
  ui         <- adj_trip$i[mask_upper]
  uj         <- adj_trip$j[mask_upper]
  
  counties   <- levels(data$county)
  partitions <- lapply(counties, function(county) {
    county_idx <- which(data$county == county)
    border_idx <- expand_border(county_idx, adjacency_matrix, hops = hops)
    subset_idx <- sort(union(county_idx, border_idx))
    # Keep only edges whose both endpoints are in subset_idx.
    # subset_idx is sorted, so original_i < original_j  ↔  i_sub < j_sub,
    # preserving the column-major order that Matrix::summary uses.
    edge_mask  <- ui %in% subset_idx & uj %in% subset_idx
    list(
      data             = data[subset_idx, ],
      adjacency_matrix = adjacency_matrix[subset_idx, subset_idx],
      distances        = distances[edge_mask]
    )
  })
  setNames(partitions, counties)
}

# ════════════════════════════════════════════════════════════════════════════════
# greedy_mi_sensor_selection_norway() ----
# ════════════════════════════════════════════════════════════════════════════════

#' Greedy MI sensor selection for the full Norwegian network.
#'
#' Partitions by county (with \code{hops}-hop border expansion), builds the
#' with+against covariance per partition, and runs greedy MI selection
#' independently per county.  Border links appearing in multiple partitions are
#' deduplicated by keeping the highest \code{mi_score}.  Optionally prunes to
#' the top-\code{r} locations Norway-wide and appends 1-hop adjacency neighbours
#' of the final selected set.
#'
#' @param data             Data frame with a \code{county} factor column.
#' @param adjacency_matrix Nationwide sparse adjacency matrix (\code{nrow(data)} x \code{nrow(data)}).
#' @param distances        Edge Gower-distance vector (\code{inla_rbf_model$distances}).
#' @param tau,d,sigma      CAR hyperparameters on the natural scale.
#' @param hops             Border expansion hops in \code{partition_by_county()} (default 2).
#' @param k                Sensors selected per county (default 10).
#' @param r                If non-NULL, retain only the top-\code{r} links by \code{mi_score}.
#' @param weighting_bias   Location-weight vector forwarded to \code{greedy_mi_sensor_selection()}.
#' @param include_neighbours Append 1-hop neighbours of the final selected set
#'   (\code{selected = FALSE}, \code{mi_score = NA}).
#' @return List: \code{selected_data_entries} (data frame), \code{n_counties},
#'   \code{k_per_county}, \code{r}.
greedy_mi_sensor_selection_norway <- function(
    data,
    adjacency_matrix,
    distances,
    tau,
    d,
    sigma,
    hops               = 2,
    k                  = 10,
    r                  = NULL,
    weighting_bias     = NULL,
    include_neighbours = FALSE,
    verbose            = FALSE) {
  
  if (!"county" %in% names(data) || !is.factor(data$county))
    stop("'data$county' must be a factor.")
  if (!is.null(r) && (!is.numeric(r) || length(r) != 1 || r < 1))
    stop("'r' must be a positive integer scalar or NULL.")
  
  partition  <- partition_by_county(data, adjacency_matrix, distances, hops = hops)
  n_counties <- length(partition)
  
  # Per-county covariance + greedy selection.
  # adjacency_matrix = NULL: neighbour expansion is deferred so it operates on
  # the post-deduplication / post-r set rather than raw per-county output.
  per_county_results <- vector("list", n_counties)
  for (i in seq_len(n_counties)) {
    message("[", i, "/", n_counties, "] ", names(partition)[i])
    part <- partition[[i]]
    cov  <- create_covariance_and_precision_matrix(
      adjacency_matrix = part$adjacency_matrix,
      tau = tau, d = d, sigma = sigma,
      distances        = part$distances,
      data             = part$data,
      with_and_against = TRUE
    )
    res <- greedy_mi_sensor_selection(
      data              = part$data,
      covariance_matrix = cov$covariance_matrix_sum,
      ids               = cov$ids,
      k                 = k,
      adjacency_matrix  = NULL,
      weighting_bias    = weighting_bias
    )
    per_county_results[[i]] <- res$selected_data_entries
  }
  
  # Combine; border links duplicated across partitions → keep highest mi_score.
  all_selected <- do.call(rbind, per_county_results)
  all_selected <- all_selected[order(all_selected$id, -all_selected$mi_score), ]
  all_selected <- all_selected[!duplicated(all_selected$id), ]
  
  # Optional top-r ranking across Norway.
  # r counts SENSORS (base IDs), not directed links — keep both -WITH and
  # -AGAINST for each of the top-r base IDs ranked by mi_score descending.
  if (!is.null(r)) {
    sel_true     <- all_selected[all_selected$selected %in% TRUE, ]
    base_ids     <- sub("-WITH$|-AGAINST$", "", sel_true$id)
    base_scores  <- tapply(sel_true$mi_score, base_ids, max, na.rm = TRUE)
    top_bases    <- names(sort(base_scores, decreasing = TRUE))[seq_len(min(r, length(base_scores)))]
    keep_ids     <- all_selected$selected %in% FALSE |
                      sub("-WITH$|-AGAINST$", "", all_selected$id) %in% top_bases
    all_selected <- all_selected[keep_ids, ]
  }
  
  # Optional 1-hop neighbour expansion on the final selected set.
  # Uses per-partition adjacency matrices; the full nationwide matrix is never materialised.
  if (include_neighbours) {
    selected_ids_final <- all_selected$id
    neigh_list         <- vector("list", n_counties)
    for (i in seq_len(n_counties)) {
      part            <- partition[[i]]
      sel_in_part_idx <- which(part$data$id %in% selected_ids_final)
      if (length(sel_in_part_idx) == 0L) next
      adj_sub        <- part$adjacency_matrix[sel_in_part_idx, , drop = FALSE]
      neigh_idx      <- setdiff(which(Matrix::colSums(adj_sub != 0) > 0), sel_in_part_idx)
      if (length(neigh_idx) == 0L) next
      nr             <- part$data[neigh_idx, ]
      nr$selected    <- FALSE
      nr$mi_score    <- NA_real_
      neigh_list[[i]] <- nr
    }
    all_neighbours <- do.call(rbind, neigh_list)
    if (!is.null(all_neighbours) && nrow(all_neighbours) > 0) {
      all_neighbours <- all_neighbours[!duplicated(all_neighbours$id), ]
      all_neighbours <- all_neighbours[!all_neighbours$id %in% selected_ids_final, ]
      all_selected   <- rbind(all_selected, all_neighbours)
    }
  }
  
  list(
    selected_data_entries = all_selected,
    n_counties            = n_counties,
    k_per_county          = k,
    r                     = r
  )
}

# ════════════════════════════════════════════════════════════════════════════════
# plot_sensor_selection_map() ----
# ════════════════════════════════════════════════════════════════════════════════

#' Interactive leaflet map of sensor selection results
#'
#' Minimal version: plots only the selected sensors (red) and their neighbour
#' context links (steel blue) from \code{selected_data_entries}.
#' Permanently visible numbered circle markers are placed at the midpoint of
#' each selected link, ranked by \code{mi_score} descending.
#'
#' @param sensor_results  List from \code{greedy_mi_sensor_selection_norway()}
#'   (must contain \code{$selected_data_entries}), or the data frame directly.
#' @param color_selected  Hex colour for selected links. Default \code{"#E63946"}.
#' @param color_neighbour Hex colour for neighbour links. Default \code{"#457B9D"}.
#' @param weight_selected   Line weight for selected links (default 6).
#' @param weight_neighbour  Line weight for neighbour links (default 3).
#' @param opacity_selected  Opacity for selected links (default 0.95).
#' @param opacity_neighbour Opacity for neighbour links (default 0.7).
#'
#' @return A \code{leaflet} map object.
#' @export
plot_sensor_selection_map <- function(
    sensor_results,
    color_selected    = "#E63946",
    color_neighbour   = "#457B9D",
    weight_selected   = 6,
    weight_neighbour  = 3,
    opacity_selected  = 0.95,
    opacity_neighbour = 0.7) {

  # ── 1. Resolve input ──────────────────────────────────────────────────────────
  if (is.list(sensor_results) && !is.data.frame(sensor_results)) {
    if (!"selected_data_entries" %in% names(sensor_results))
      stop("sensor_results must contain '$selected_data_entries'.")
    sel_df <- sensor_results$selected_data_entries
  } else if (is.data.frame(sensor_results)) {
    sel_df <- sensor_results
  } else {
    stop("sensor_results must be a list from greedy_mi_sensor_selection_norway() or a data frame.")
  }
  if (!"selected" %in% names(sel_df))
    stop("sel_df must contain a 'selected' column (TRUE = selected, FALSE = neighbour).")

  sel_df <- sel_df[!duplicated(sel_df$id), , drop = FALSE]

  # ── 2. Geometry: add_geometries() directly on sel_df ─────────────────────────
  # sel_df from greedy_mi_sensor_selection already carries 'selected' and
  # 'mi_score' as plain columns — call add_geometries() once, no joins on sf.
  sel_sf <- add_geometries(sel_df)
  sel_sf <- sf::st_transform(sel_sf, crs = 4326)

  # Drop rows with empty or invalid geometry (prevents lat/lon warnings)
  valid <- !sf::st_is_empty(sel_sf) & sf::st_is_valid(sel_sf)
  sel_sf <- sel_sf[valid, ]

  # ── 3. Rank + popup (before split so subsets inherit both columns) ────────────
  # Rank is per SENSOR (base ID = strip -WITH/-AGAINST), so both directions of
  # the same sensor share the same rank number.  With r=200, ranks run 1..200.
  sel_sf$base_id <- sub("-WITH$|-AGAINST$", "", sel_sf$id)
  sel_sf$rank <- NA_integer_
  true_idx <- which(sel_sf$selected %in% TRUE)
  if (length(true_idx) > 0) {
    true_sf <- sel_sf[true_idx, ]
    if ("mi_score" %in% names(true_sf)) {
      # Unique base IDs ordered by mi_score descending (same score for WITH+AGAINST)
      base_order <- unique(true_sf$base_id[order(-true_sf$mi_score, na.last = TRUE)])
    } else {
      base_order <- unique(true_sf$base_id)
    }
    rank_lookup <- stats::setNames(seq_along(base_order), base_order)
    sel_sf$rank[true_idx] <- rank_lookup[sel_sf$base_id[true_idx]]
  }

  has_mi  <- "mi_score"              %in% names(sel_sf)
  has_ly  <- "lastYearAadt_aadt"     %in% names(sel_sf)
  has_co  <- "county"                %in% names(sel_sf)
  has_src <- "traffic_volume_source" %in% names(sel_sf)

  pop <- paste0("<strong>ID:</strong> ", sel_sf$id,
                "<br><strong>Type:</strong> ",
                ifelse(sel_sf$selected %in% TRUE, "Selected sensor", "Neighbour"))
  if (any(!is.na(sel_sf$rank)))
    pop <- paste0(pop, "<br><strong>Sensor rank:</strong> ",
                  ifelse(is.na(sel_sf$rank), "\u2014", as.character(sel_sf$rank)))
  if (has_mi)
    pop <- paste0(pop, "<br><strong>MI score:</strong> ",
                  ifelse(is.na(sel_sf$mi_score), "\u2014",
                         sprintf("%.4f", sel_sf$mi_score)))
  if (has_co)
    pop <- paste0(pop, "<br><strong>County:</strong> ",
                  ifelse(is.na(sel_sf$county), "\u2014", as.character(sel_sf$county)))
  if (has_src)
    pop <- paste0(pop, "<br><strong>Source:</strong> ",
                  ifelse(is.na(sel_sf$traffic_volume_source), "Unmeasured",
                         sel_sf$traffic_volume_source))
  if (has_ly)
    pop <- paste0(pop, "<br><strong>Last year AADT:</strong> ",
                  format(round(sel_sf$lastYearAadt_aadt),
                         big.mark = "\u00a0", scientific = FALSE))
  sel_sf$popup_text <- pop

  # ── 4. Split (inherits popup_text, rank, base_id) ────────────────────────────
  sel_selected  <- sel_sf[sel_sf$selected %in% TRUE,  ]
  sel_neighbour <- sel_sf[sel_sf$selected %in% FALSE, ]

  # ── 5. One midpoint marker per SENSOR (base ID), not per directed link ────────
  # Both WITH and AGAINST share one numbered circle.
  markers_sf <- NULL
  if (nrow(sel_selected) > 0) {
    # Keep one representative directed link per base_id (first occurrence = -WITH)
    one_per_sensor <- sel_selected[!duplicated(sel_selected$base_id), ]
    suppressWarnings(markers_sf <- sf::st_point_on_surface(one_per_sensor))
    coords <- sf::st_coordinates(markers_sf)
    markers_sf <- markers_sf[!is.na(coords[, 1]) & !is.na(coords[, 2]), ]
    if (nrow(markers_sf) == 0) markers_sf <- NULL
  }

  # ── 6. Build leaflet map ──────────────────────────────────────────────────────
  nvdb <- nvdb_objects()

  m <- leaflet::leaflet(
    sel_sf,
    options = leaflet::leafletOptions(crs = nvdb$nvdb_crs, zoomControl = TRUE)
  ) |>
    leaflet::addTiles(urlTemplate   = nvdb$nvdb_url,
                      attribution   = nvdb$nvdb_attribution)

  # Neighbour links — steel blue
  if (nrow(sel_neighbour) > 0)
    m <- m |> leaflet::addPolylines(
      data             = sel_neighbour,
      color            = color_neighbour,
      weight           = weight_neighbour,
      opacity          = opacity_neighbour,
      popup            = sel_neighbour$popup_text,
      group            = "Neighbours",
      highlightOptions = leaflet::highlightOptions(
        weight = weight_neighbour + 2, color = "white", bringToFront = TRUE))

  # Selected links — red, thick
  if (nrow(sel_selected) > 0)
    m <- m |> leaflet::addPolylines(
      data             = sel_selected,
      color            = color_selected,
      weight           = weight_selected,
      opacity          = opacity_selected,
      popup            = sel_selected$popup_text,
      group            = "Selected sensors",
      highlightOptions = leaflet::highlightOptions(
        weight = weight_selected + 2, color = "white", bringToFront = TRUE))

  # Numbered circle markers — no group → always visible regardless of layer toggles
  if (!is.null(markers_sf))
    m <- m |> leaflet::addCircleMarkers(
      data         = markers_sf,
      radius       = 8,
      color        = "#000000",
      weight       = 1.5,
      fillColor    = color_selected,
      fillOpacity  = 0.95,
      popup        = markers_sf$popup_text,
      label        = as.character(markers_sf$rank),
      labelOptions = leaflet::labelOptions(
        permanent = TRUE,
        textOnly  = TRUE,
        style     = list("font-weight" = "bold",
                         "font-size"   = "10px",
                         "color"       = "white",
                         "text-shadow" = "0 0 3px #000")))

  # Layer controls + legend
  m <- m |>
    leaflet::addLayersControl(
      overlayGroups = c("Neighbours", "Selected sensors"),
      options       = leaflet::layersControlOptions(collapsed = FALSE)) |>
    leaflet::addLegend(
      position = "bottomright",
      colors   = c(color_selected, color_neighbour),
      labels   = c("Selected sensor", "Neighbour link"),
      title    = "Sensor selection",
      opacity  = 0.9)

  m
}
