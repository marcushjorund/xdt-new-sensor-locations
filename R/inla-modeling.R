#' Fit a weighted Besag-proper spatial INLA model for traffic AADT
#'
#' Fits a CAR (conditional autoregressive) spatial model to traffic-link AADT data using
#' INLA. When `spatial_term = "besagproper_rbf"`, edge weights in the CAR precision matrix
#' are modulated by an RBF/Laplacian kernel over Gower distance between adjacent links'
#' road attributes (`similarity_covariates`), allowing dissimilar neighbouring links to be
#' weakly coupled.
#'
#' @param data              Prepared traffic links data frame; must contain `id` and the
#'   response column (`aadt` or `heavyAadt`).
#' @param adjacency_matrix  Square sparse adjacency matrix of the traffic links, matching
#'   `nrow(data)`.
#' @param spatial_term      Either `"besagproper"` or `"besagproper_rbf"` (default).
#' @param similarity_covariates Character vector of column names used to compute Gower
#'   distances when `spatial_term = "besagproper_rbf"`. `NULL` uses built-in defaults.
#' @param ordinal_levels    Named list of ordered levels for ordinal covariates, e.g.
#'   `list(functionalRoadClass = as.character(c(7,6,5,4,3,2,1,0)), roadCategory =
#'   c("KOMMUNAL_VEG","FYLKESVEG","RIKSVEG","EUROPAVEG"))`.
#' @param distance_type     `"gower"` (standard absolute Gower, default) or
#'   `"gower_squared"` (squared per-component).
#' @param weight_fn         `"laplacian"` (default), `"gaussian"`, or `"linear"`.
#' @param fixed_effects     Formula for fixed effects in the INLA model.
#' @param iid_effects       Character vector of IID random effect column names, or `NULL`.
#'   Default `"roadSystem"`.
#' @param heavy_vehicle     Logical; if `TRUE`, models `heavyAadt` instead of `aadt`.
#' @param family            Likelihood family: `"gaussian"`, `"poisson"`, or `"nbinomial"`.
#' @param extraconstr       Optional flow-conservation constraint:
#'   `list(A = <m x n matrix>, e = <m vector>)`.
#' @param verbose           Logical; if `TRUE`, returns verbose INLA output.
#' @return A `weighted_inla_traffic_model` object (list) with elements `predictions`,
#'   `model_summary`, `spatial_hyperpar`, `distances`, `inla_model`, `formula`, `family`,
#'   `spatial_term`, `similarity_covariates`, `ordinal_levels`.
fit_weighted_inla_model <- function(
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
    extraconstr       = NULL,              #optional flow conservation: list(A = <m x n matrix>, e = <m vector>)
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
    
    weighted_model <- INLA::inla.rgeneric.define(
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
              model       = weighted_model,
              n           = n,
              constr      = FALSE,
              extraconstr = extraconstr)
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
  
  class(result) <- "weighted_inla_traffic_model"
  result
}

# ────────────────────────────────────────────────────────────────────────────────────
# S3 print method ----
# ────────────────────────────────────────────────────────────────────────────────────

#' Print method for weighted_inla_traffic_model objects
#'
#' @param x   A `weighted_inla_traffic_model` object.
#' @param ... Unused; present for S3 method consistency.
#' @return `x`, invisibly.
#' @export
print.weighted_inla_traffic_model <- function(x, ...) {
  cat("INLA Weighted Traffic Model\n")
  cat("======================\n\n")
  cat("Number of predictions:", nrow(x$predictions), "\n")
  cat("Family:              ", x$family, "\n")
  cat("Spatial term:        ", x$spatial_term, "\n")
  if (!is.null(x$similarity_covariates))
    cat("Weighted similarity_covariates: ", paste(x$similarity_covariates, collapse = ", "), "\n")
  else
    cat("Weighted similarity_covariates:  <defaults>\n")
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

#' @description Root mean squared error between predicted (p) and observed (o) values.
.rmse <- function(p, o) sqrt(mean((p - o)^2, na.rm = TRUE))
#' @description Mean absolute error between predicted (p) and observed (o) values.
.mae  <- function(p, o) mean(abs(p - o),      na.rm = TRUE)
#' @description Mean absolute percentage error (observed > 0 only).
.mape <- function(p, o) {
  nz <- o > 0
  mean(abs((p[nz] - o[nz]) / o[nz]), na.rm = TRUE) * 100
}
#' @description Median absolute log error (both predicted and observed > 0 only).
.male <- function(p, o) {
  ok <- p > 0 & o > 0
  median(abs(log(p[ok]) - log(o[ok])), na.rm = TRUE)
}

#' @description Formats RMSE/MAE/MAPE/MALE as a multi-line annotation label for plots.
.metric_label <- function(pred, obs) {
  sprintf("RMSE = %.0f\nMAE  = %.0f\nMAPE = %.1f%%\nMALE = %.3f",
          .rmse(pred, obs), .mae(pred, obs), .mape(pred, obs), .male(pred, obs))
}

# ════════════════════════════════════════════════════════════════════════════════
# kfold_cv_inla() ----
# ════════════════════════════════════════════════════════════════════════════════

#' K-fold cross-validation for \code{fit_weighted_inla_model}
#'
#' @param data            Data frame passed to \code{fit_weighted_inla_model}.
#' @param adjacency_matrix Square adjacency matrix.
#' @param k               Number of folds (default 5).
#' @param seed            Random seed for fold assignment (default 42).
#' @param fixed_effects   Formula of fixed effects.
#' @param spatial_term    Passed to \code{fit_weighted_inla_model}.
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
    
    fit <- fit_weighted_inla_model(
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
# select_similarity_covariates() ----
# ════════════════════════════════════════════════════════════════════════════════

#' Forward selection of optimal similarity_covariates for the weighted model.
#'
#' Greedy forward selection: starts from an empty set and at each step adds
#' the candidate that most reduces mean CV-RMSE across folds.  Stops when no
#' addition improves performance.  Uses \code{kfold_cv_inla()} internally.
#' Ordinal levels are automatically subsetted to the active covariate set at
#' each step so irrelevant specs never affect the Gower distance computation.
#'
#' @param data                Data frame passed to \code{fit_weighted_inla_model}.
#' @param adjacency_matrix    Square adjacency matrix.
#' @param candidate_covariates Character vector of all covariates to consider.
#' @param ordinal_levels      Named list of ordered levels for ordinal
#'   covariates; subsetted automatically per step.
#' @param fixed_effects       Formula; held constant throughout the search.
#' @param iid_effects         Character; passed through.
#' @param family              Likelihood family; passed through.
#' @param distance_type       Passed to \code{fit_weighted_inla_model}.
#' @param weight_fn           Passed to \code{fit_weighted_inla_model}.
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


