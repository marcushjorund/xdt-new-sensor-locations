# ────────────────────────────────────────────────────────────────────────────────────
# RBF Kernel INLA Estimation with Comparison to Ordinary Besag
# ────────────────────────────────────────────────────────────────────────────────────
# Purpose: Estimate spatial AADT effects using RBF-weighted Besag CAR model and
#          compare predictive performance across three spatial structures:
#            1. Ordinary besagproper (unweighted adjacency)
#            2. RBF-weighted Besag — AADT-only similarity (no normalisation)
#            3. RBF-weighted Besag — all covariates (normalised feature space)
#
# Sections 1-5:  Data prep, feature matrix, graph structure, pairing, rgeneric
# Sections 6-7:  Full-data fits → covariance matrices for sensor selection
# Section  8:    Hold-out CV comparison (20% test set)
# Section  9:    Predicted vs observed plots
# Section 10:    Greedy sensor selection using RBF covariance
# ────────────────────────────────────────────────────────────────────────────────────

library(xdtkit)
library(INLA)
library(Matrix)

# ────────────────────────────────────────────────────────────────────────────────────
# 1. LOAD & PREPARE DATA
# ────────────────────────────────────────────────────────────────────────────────────

year <- 2025

# Load Buskerud traffic links
preprocessed_traffic_links <- preprocess_traffic_links(buskerud_directed_traffic_links, year = year)
bus_aadt <- calculate_bus_aadt(stops_on_traffic_links, bus_counts, year = year)

# Fill missing values
prepared_traffic_links <- fill_missing_values(
  df = preprocessed_traffic_links,
  unknown_impute_columns = c("functionClass", "highestSpeedLimit", "lowestSpeedLimit", "maxLanes", "minLanes"),
  mode_impute_columns = c("hasOnlyPublicTransportLanes"),
  median_impute_columns = c("lastYearAadt_aadt", "lastYearAadt_heavyRatio", "lastYearAadt_heavyAadt")
) |>
  remove_negative_aadt() |>
  add_logLastYear() |>
  join_bus_to_traffic(bus_aadt)

# Optionally load pre-processed Emma version (ensure consistency)
# prepared_traffic_links_emma <- readRDS("prepared_traffic_links_emma.rds")
# prepared_traffic_links <- subset(prepared_traffic_links_emma, select = -c(spatial.idx, inla_pred))

# Build adjacency matrix
adjacency_matrix <- build_adjacency_matrix(
  prepared_traffic_links,
  exclude_public_transport = TRUE
)

n <- nrow(prepared_traffic_links)
prepared_traffic_links$spatial.idx <- seq_len(n)
prepared_traffic_links$log_aadt   <- log(prepared_traffic_links$aadt)

cat("Data: n =", n, "traffic links\n")
cat("Measured AADT:", sum(!is.na(prepared_traffic_links$aadt)), "links\n")
cat("Adjacency matrix: sparse,", sum(adjacency_matrix != 0), "non-zero entries\n\n")

# ────────────────────────────────────────────────────────────────────────────────────
# 2. BUILD NORMALISED FEATURE MATRIX FOR SIMILARITY
# ────────────────────────────────────────────────────────────────────────────────────
# Column-standardise + L2-normalise all covariates so that squared Euclidean
# distance in feature space is bounded on [0, 4] regardless of dimensionality.
# Used as input to the RBF kernel: W_ij = exp(-||f_i - f_j||^2 / (2*ell^2)).

build_normalised_feature_matrix <- function(
  df,
  numeric_features     = c("maxLanes", "minLanes",
                            "highestSpeedLimit", "lowestSpeedLimit",
                            "lastYearAadt_logAadt"),
  binary_features      = c("hasOnlyPublicTransportLanes", "isRamp"),
  categorical_features = c("roadCategory", "functionClass", "roadSystem"),
  interaction_pairs    = list(
    c("functionClass",        "roadCategory"),
    c("functionClass",        "isRamp"),
    c("maxLanes",             "roadCategory"),
    c("lastYearAadt_logAadt", "roadCategory")
  )) {

  one_hot <- function(x, name, levels) {
    mat <- outer(as.character(x), levels, `==`) * 1L
    colnames(mat) <- paste0(name, "_", levels)
    mat
  }

  cat_levels <- lapply(
    intersect(categorical_features, names(df)),
    function(col) sort(unique(as.character(df[[col]])))
  )
  names(cat_levels) <- intersect(categorical_features, names(df))

  parts <- list()

  for (col in intersect(numeric_features, names(df)))
    parts[[col]] <- as.numeric(df[[col]])

  for (col in intersect(binary_features, names(df)))
    parts[[col]] <- as.numeric(df[[col]])

  for (col in intersect(categorical_features, names(df))) {
    oh <- one_hot(df[[col]], col, cat_levels[[col]])
    for (k in seq_len(ncol(oh))) parts[[colnames(oh)[k]]] <- oh[, k]
  }

  for (pair in interaction_pairs) {
    a <- pair[1]; b <- pair[2]
    if (!a %in% names(df) || !b %in% names(df)) next
    a_is_cat <- a %in% categorical_features
    b_is_cat <- b %in% categorical_features
    a_is_num <- a %in% c(numeric_features, binary_features)
    b_is_num <- b %in% c(numeric_features, binary_features)

    if (a_is_cat && b_is_cat) {
      cross <- paste0(as.character(df[[a]]), ":", as.character(df[[b]]))
      lvls  <- sort(unique(cross))
      oh    <- one_hot(cross, paste0(a, "x", b), lvls)
      for (k in seq_len(ncol(oh))) parts[[colnames(oh)[k]]] <- oh[, k]
    } else if (a_is_cat && b_is_num) {
      oh       <- one_hot(df[[a]], a, cat_levels[[a]])
      num_vals <- as.numeric(df[[b]])
      for (k in seq_len(ncol(oh)))
        parts[[paste0(colnames(oh)[k], "x", b)]] <- oh[, k] * num_vals
    } else if (a_is_num && b_is_cat) {
      oh       <- one_hot(df[[b]], b, cat_levels[[b]])
      num_vals <- as.numeric(df[[a]])
      for (k in seq_len(ncol(oh)))
        parts[[paste0(colnames(oh)[k], "x", a)]] <- oh[, k] * num_vals
    } else if (a_is_num && b_is_num) {
      parts[[paste0(a, "x", b)]] <- as.numeric(df[[a]]) * as.numeric(df[[b]])
    }
  }

  mat <- do.call(cbind, parts)

  # Column-standardise; leave constant columns at 0
  col_means <- colMeans(mat, na.rm = TRUE)
  col_sds   <- apply(mat, 2, sd, na.rm = TRUE)
  non_const <- col_sds > 0
  mat[, non_const]  <- sweep(mat[, non_const, drop = FALSE], 2, col_means[non_const], `-`)
  mat[, non_const]  <- sweep(mat[, non_const, drop = FALSE], 2, col_sds[non_const],   `/`)
  mat[, !non_const] <- 0
  mat[is.na(mat)]   <- 0

  # L2-normalise each row: cosine similarity reduces to dot product
  row_norms <- sqrt(rowSums(mat^2))
  row_norms[row_norms == 0] <- 1
  mat / row_norms
}

cat("Building normalized feature matrix...\n")
feat_norm <- build_normalised_feature_matrix(prepared_traffic_links)
cat("  Features: ", ncol(feat_norm), " columns after normalization\n\n")

# ────────────────────────────────────────────────────────────────────────────────────
# 3. EXTRACT GRAPH STRUCTURE & COVARIATE DISTANCES
# ────────────────────────────────────────────────────────────────────────────────────

# Extract upper triangle of adjacency matrix (to avoid double-counting)
adj_trip <- Matrix::summary(as(adjacency_matrix, "dgCMatrix"))
mask_upper <- adj_trip$i < adj_trip$j
ui <- adj_trip$i[mask_upper]
uj <- adj_trip$j[mask_upper]

# Precompute squared covariate distances for all neighbour pairs
# Full model: all covariates (column-standardised + L2-normalised)
dist_sq  <- rowSums((feat_norm[ui, ] - feat_norm[uj, ])^2)
ell_init <- median(sqrt(dist_sq))   # median heuristic starting value

# AADT-only model: raw lastYearAadt_aadt — no normalisation needed since ell is
# estimated from the data and absorbs whatever scale dist_sq is measured on
aadt_vals     <- prepared_traffic_links$lastYearAadt_aadt
dist_sq_aadt  <- (aadt_vals[ui] - aadt_vals[uj])^2
ell_init_aadt <- median(sqrt(dist_sq_aadt))

cat("Graph structure extracted:\n")
cat("  Upper-triangle neighbour pairs:", length(ui), "\n")
cat("  Median covariate distance, full features (ell_init):", round(ell_init, 4), "\n")
cat("  Median AADT distance (ell_init_aadt):", round(ell_init_aadt, 2), "\n\n")

# ────────────────────────────────────────────────────────────────────────────────────
# 4. BUILD WITH + AGAINST PAIRING MATRIX (for sensor selection)
# ────────────────────────────────────────────────────────────────────────────────────
# Traffic links come in WITH/AGAINST pairs. We combine them via linear transformation.

with_idx_int    <- which(grepl("-WITH$",    prepared_traffic_links$id))
against_idx_int <- which(grepl("-AGAINST$", prepared_traffic_links$id))

base_with    <- sub("-WITH$",    "", prepared_traffic_links$id[with_idx_int])
base_against <- sub("-AGAINST$", "", prepared_traffic_links$id[against_idx_int])

paired_bases <- intersect(base_with, base_against)
unpaired_bases <- setdiff(base_with, base_against)

paired_with_idx    <- match(paste0(paired_bases,    "-WITH"),    prepared_traffic_links$id)
paired_against_idx <- match(paste0(paired_bases,    "-AGAINST"), prepared_traffic_links$id)
unpaired_with_idx  <- match(paste0(unpaired_bases,  "-WITH"),    prepared_traffic_links$id)

m <- length(paired_bases) + length(unpaired_bases)  # combined segments

# Build sparse C matrix: sum WITH + AGAINST
row_idx <- c(
  seq_along(paired_bases),
  seq_along(paired_bases),
  length(paired_bases) + seq_along(unpaired_bases)
)
col_idx <- c(paired_with_idx, paired_against_idx, unpaired_with_idx)
values <- rep(1, length(row_idx))

C <- sparseMatrix(i = row_idx, j = col_idx, x = values, dims = c(m, n))

cat("WITH/AGAINST pairing:\n")
cat("  Paired segments:", length(paired_bases), "\n")
cat("  Unpaired WITH segments:", length(unpaired_bases), "\n")
cat("  Combined dimension m:", m, "\n\n")

# ────────────────────────────────────────────────────────────────────────────────────
# 5. DEFINE RBF-WEIGHTED BESAG PROPER MODEL (rgeneric)
# ────────────────────────────────────────────────────────────────────────────────────
# Q(tau, d, ell) = tau * (d*I + D_W - W)
# where W_ij = exp(-||feat_i - feat_j||^2 / (2*ell^2)) if adjacent, 0 otherwise
# theta[1] = log(tau), theta[2] = log(d), theta[3] = log(ell)

inla.rgeneric.weighted.besag.RBF <- function(cmd, theta) {
  # Injected via inla.rgeneric.define(): n, ui, uj, dist_sq, ell_init

  graph <- function() {
    # Sparsity structure: graph + diagonal
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
      d        <- exp(theta[2]) + 1.0    # shift ≥ 1 keeps Q well-conditioned
      ell      <- pmax(exp(theta[3]), 1e-6)  # pmax handles NaN safely (max does not)
      exponent <- pmin(dist_sq / (2 * ell^2), 700)
      sim_vals <- pmax(exp(-exponent), 1e-10)
      W        <- Matrix::sparseMatrix(i = c(ui, uj), j = c(uj, ui),
                                       x = rep(sim_vals, 2), dims = c(n, n))
      Q0_local <- Matrix::Diagonal(x = Matrix::rowSums(W)) - W
      tau * (d * Matrix::Diagonal(n) + Q0_local)
    }, error = function(e) {
      # Fallback: unweighted Laplacian + I, same sparsity as graph()
      W_fb   <- Matrix::sparseMatrix(i = c(ui, uj), j = c(uj, ui),
                                     x = rep(1, 2 * length(ui)), dims = c(n, n))
      Matrix::Diagonal(x = Matrix::rowSums(W_fb)) - W_fb + Matrix::Diagonal(n)
    })
  }

  mu             <- function() numeric(n)
  log.norm.const <- function() numeric(0)   # INLA computes via Cholesky
  quit           <- function() invisible()

  log.prior <- function() {
    # Log-Gamma(1, 5e-5) prior on tau and d (vague, matching INLA besagproper default)
    lp_gamma <- function(theta_k, a = 1, b = 5e-5) a * theta_k - b * exp(theta_k)
    # Normal prior on log(ell) centered at log(ell_init); sd = 1.0 gives a
    # 95% interval of roughly [ell_init/7, 7*ell_init] — wide enough to let
    # the data determine the bandwidth while keeping optimisation stable
    lp_ell <- dnorm(theta[3], mean = log(ell_init), sd = 1.0, log = TRUE)
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

# Create the rgeneric model objects — same function, different injected distances
weighted_besag_model_RBF <- inla.rgeneric.define(
  inla.rgeneric.weighted.besag.RBF,
  n        = n,
  ui       = ui,
  uj       = uj,
  dist_sq  = dist_sq,
  ell_init = ell_init
)

# AADT-only variant: injects raw AADT distances instead of normalised feature distances
weighted_besag_model_RBF_aadt <- inla.rgeneric.define(
  inla.rgeneric.weighted.besag.RBF,
  n        = n,
  ui       = ui,
  uj       = uj,
  dist_sq  = dist_sq_aadt,
  ell_init = ell_init_aadt
)

cat("RBF-weighted Besag rgeneric models defined (full covariates + AADT-only).\n\n")

# ────────────────────────────────────────────────────────────────────────────────────
# 6. FIT RBF MODELS ON FULL DATA
# ────────────────────────────────────────────────────────────────────────────────────
# These full-data fits produce the posterior hyperparameters used to reconstruct
# precision/covariance matrices (section 7) for downstream sensor selection
# (section 10). The CV comparison in section 8 refits from scratch on masked data.

# ── 6a. Full covariates ──────────────────────────────────────────────────────────

formula_rbf <- log_aadt ~
  lastYearAadt_logAadt +
  f(spatial.idx,
    model  = weighted_besag_model_RBF,
    n      = n,
    constr = FALSE) +
  f(roadSystem, model = "iid")

cat("Fitting RBF-weighted Besag model...\n")
inla_model_rbf <- inla(
  formula_rbf,
  data              = prepared_traffic_links,
  family            = "gaussian",
  control.predictor = list(compute = TRUE),
  control.compute   = list(dic = TRUE, waic = TRUE),
  num.threads       = "1:1",
  safe              = TRUE
)

# Extract hyperparameters (rgeneric names them Theta1, Theta2, Theta3)
hp_rbf <- inla_model_rbf$summary.hyperpar
idx_tau <- grep("Theta1", rownames(hp_rbf))
idx_d   <- grep("Theta2", rownames(hp_rbf))
idx_ell <- grep("Theta3", rownames(hp_rbf))
cat("\nRBF-weighted Besag (all covariates) hyperparameter estimates:\n")
cat("  tau:", round(exp(hp_rbf[idx_tau, "mean"]), 4),
    " 95% CI [", round(exp(hp_rbf[idx_tau, "0.025quant"]), 4), ",",
    round(exp(hp_rbf[idx_tau, "0.975quant"]), 4), "]\n")
cat("  d:  ", round(exp(hp_rbf[idx_d, "mean"]), 4),
    " 95% CI [", round(exp(hp_rbf[idx_d, "0.025quant"]), 4), ",",
    round(exp(hp_rbf[idx_d, "0.975quant"]), 4), "]\n")
cat("  ell:", round(exp(hp_rbf[idx_ell, "mean"]), 4),
    " 95% CI [", round(exp(hp_rbf[idx_ell, "0.025quant"]), 4), ",",
    round(exp(hp_rbf[idx_ell, "0.975quant"]), 4), "]\n\n")

# ── 6b. AADT-only similarity ─────────────────────────────────────────────────────────

formula_rbf_aadt <- log_aadt ~
  lastYearAadt_logAadt +
  f(spatial.idx,
    model  = weighted_besag_model_RBF_aadt,
    n      = n,
    constr = FALSE) +
  f(roadSystem, model = "iid")

cat("Fitting RBF-weighted Besag model (AADT-only similarity)...\n")
inla_model_rbf_aadt <- inla(
  formula_rbf_aadt,
  data              = prepared_traffic_links,
  family            = "gaussian",
  control.predictor = list(compute = TRUE),
  control.compute   = list(dic = TRUE, waic = TRUE),
  num.threads       = "1:1",
  safe              = TRUE
)

hp_rbf_aadt <- inla_model_rbf_aadt$summary.hyperpar
idx_tau_a <- grep("Theta1", rownames(hp_rbf_aadt))
idx_d_a   <- grep("Theta2", rownames(hp_rbf_aadt))
idx_ell_a <- grep("Theta3", rownames(hp_rbf_aadt))
cat("\nRBF Besag (AADT-only) hyperparameter estimates:\n")
cat("  tau:", round(exp(hp_rbf_aadt[idx_tau_a, "mean"]), 4),
    " 95% CI [", round(exp(hp_rbf_aadt[idx_tau_a, "0.025quant"]), 4), ",",
    round(exp(hp_rbf_aadt[idx_tau_a, "0.975quant"]), 4), "]\n")
cat("  d:  ", round(exp(hp_rbf_aadt[idx_d_a, "mean"]), 4),
    " 95% CI [", round(exp(hp_rbf_aadt[idx_d_a, "0.025quant"]), 4), ",",
    round(exp(hp_rbf_aadt[idx_d_a, "0.975quant"]), 4), "]\n")
cat("  ell:", round(exp(hp_rbf_aadt[idx_ell_a, "mean"]), 2), "(vehicles/day)",
    " 95% CI [", round(exp(hp_rbf_aadt[idx_ell_a, "0.025quant"]), 2), ",",
    round(exp(hp_rbf_aadt[idx_ell_a, "0.975quant"]), 2), "]\n\n")

# ────────────────────────────────────────────────────────────────────────────────────
# 7. RECONSTRUCT COVARIANCE MATRICES & APPLY SENSOR PAIRING
# ────────────────────────────────────────────────────────────────────────────────────

tau_hat <- exp(hp_rbf[idx_tau, "mean"])
d_hat   <- exp(hp_rbf[idx_d, "mean"]) + 1.0  # match Q() shift
ell_hat <- exp(hp_rbf[idx_ell, "mean"])

# Rebuild Q at posterior mean hyperparameters
exponent_hat <- pmin(dist_sq / (2 * ell_hat^2), 700)
sim_vals_hat <- pmax(exp(-exponent_hat), 1e-10)
W_hat <- Matrix::sparseMatrix(
  i = c(ui, uj), j = c(uj, ui),
  x = rep(sim_vals_hat, 2), dims = c(n, n)
)
Q0_hat <- Matrix::Diagonal(x = Matrix::rowSums(W_hat)) - W_hat
Q_fitted_rbf <- tau_hat * (d_hat * Matrix::Diagonal(n) + Q0_hat)

# Invert to covariance (dense)
covariance_rbf <- as.matrix(Matrix::solve(Q_fitted_rbf))
# Apply WITH+AGAINST linear combination
covariance_rbf_sum <- as.matrix(C %*% covariance_rbf %*% t(C))

cat("Full-covariate RBF covariance reconstructed.\n")

# AADT-only: reconstruct at posterior mean
tau_hat_aadt <- exp(hp_rbf_aadt[idx_tau_a, "mean"])
d_hat_aadt   <- exp(hp_rbf_aadt[idx_d_a, "mean"]) + 1.0  # match Q() shift
ell_hat_aadt <- exp(hp_rbf_aadt[idx_ell_a, "mean"])
ell_hat
ell_hat_aadt/ell_init_aadt
exponent_aadt <- pmin(dist_sq_aadt / (2 * ell_hat_aadt^2), 700)
sim_vals_aadt <- pmax(exp(-exponent_aadt), 1e-10)
W_hat_aadt <- Matrix::sparseMatrix(
  i = c(ui, uj), j = c(uj, ui),
  x = rep(sim_vals_aadt, 2), dims = c(n, n)
)
Q0_hat_aadt       <- Matrix::Diagonal(x = Matrix::rowSums(W_hat_aadt)) - W_hat_aadt
Q_fitted_rbf_aadt <- tau_hat_aadt * (d_hat_aadt * Matrix::Diagonal(n) + Q0_hat_aadt)

covariance_rbf_aadt     <- as.matrix(Matrix::solve(Q_fitted_rbf_aadt))
covariance_rbf_aadt_sum <- as.matrix(C %*% covariance_rbf_aadt %*% t(C))

cat("AADT-only RBF covariance reconstructed.\n\n")

# ────────────────────────────────────────────────────────────────────────────────────
# 8. TRAIN/TEST COMPARISON: ORDINARY BESAG vs RBF (AADT-only) vs RBF (all covariates)
# ────────────────────────────────────────────────────────────────────────────────────
# Hold out 20% of measured links, refit both models, compare predictive metrics.

set.seed(42)
measured_all <- which(!is.na(prepared_traffic_links$aadt))
n_test <- round(0.2 * length(measured_all))
test_idx <- sample(measured_all, size = n_test)

cat("Cross-validation setup: Hold-out test set\n")
cat("  Measured links:", length(measured_all), "\n")
cat("  Train set:", length(measured_all) - n_test, "\n")
cat("  Test set: ", n_test, "\n\n")

# Create masked copy for cross-validation
data_cv <- prepared_traffic_links
data_cv$log_aadt[test_idx] <- NA
data_cv$spatial.idx         <- seq_len(n)

# ────────────────────────────────────────────────────────────────────────────────────
# 8a. Fit ordinary Besag proper on CV fold
# ────────────────────────────────────────────────────────────────────────────────────

formula_besag_cv <- log_aadt ~
  lastYearAadt_logAadt +
  f(spatial.idx,
    model  = "besagproper",
    graph  = adjacency_matrix,
    constr = FALSE) +
  f(roadSystem, model = "iid")

cat("Fitting ordinary besagproper on CV fold...\n")
inla_cv_besag <- inla(
  formula_besag_cv,
  data              = data_cv,
  family            = "gaussian",
  control.predictor = list(compute = TRUE),
  control.compute   = list(dic = TRUE, waic = TRUE),
  num.threads       = "1:1",
  safe              = TRUE
)

cat("  DIC:", round(inla_cv_besag$dic$dic, 1), "\n")
cat("  WAIC:", round(inla_cv_besag$waic$waic, 1), "\n\n")

# ────────────────────────────────────────────────────────────────────────────────────
# 8b. Fit RBF-weighted Besag (AADT-only) on same CV fold
# ────────────────────────────────────────────────────────────────────────────────────

formula_rbf_aadt_cv <- log_aadt ~
  lastYearAadt_logAadt +
  f(spatial.idx,
    model  = weighted_besag_model_RBF_aadt,
    n      = n,
    constr = FALSE) +
  f(roadSystem, model = "iid")

cat("Fitting RBF-weighted Besag (AADT-only) on CV fold...\n")
inla_cv_rbf_aadt <- inla(
  formula_rbf_aadt_cv,
  data              = data_cv,
  family            = "gaussian",
  control.predictor = list(compute = TRUE),
  control.compute   = list(dic = TRUE, waic = TRUE),
  num.threads       = "1:1",
  safe              = TRUE
)

cat("  DIC:", round(inla_cv_rbf_aadt$dic$dic, 1), "\n")
cat("  WAIC:", round(inla_cv_rbf_aadt$waic$waic, 1), "\n\n")

# ────────────────────────────────────────────────────────────────────────────────────
# 8c. Fit RBF-weighted Besag (all covariates) on same CV fold
# ────────────────────────────────────────────────────────────────────────────────────

formula_rbf_cv <- log_aadt ~
  lastYearAadt_logAadt +
  f(spatial.idx,
    model  = weighted_besag_model_RBF,
    n      = n,
    constr = FALSE) +
  f(roadSystem, model = "iid")

cat("Fitting RBF-weighted Besag (all covariates) on CV fold...\n")
inla_cv_rbf <- inla(
  formula_rbf_cv,
  data              = data_cv,
  family            = "gaussian",
  control.predictor = list(compute = TRUE),
  control.compute   = list(dic = TRUE, waic = TRUE),
  num.threads       = "1:1",
  safe              = TRUE
)

cat("  DIC:", round(inla_cv_rbf$dic$dic, 1), "\n")
cat("  WAIC:", round(inla_cv_rbf$waic$waic, 1), "\n\n")

# ────────────────────────────────────────────────────────────────────────────────────
# 8d. Extract test-set predictions & compute metrics
# ────────────────────────────────────────────────────────────────────────────────────

true_aadt <- prepared_traffic_links$aadt[test_idx]

# Fitted values are on log scale; back-transform to AADT for comparison
fv_besag    <- inla_cv_besag$summary.fitted.values[test_idx, ]
fv_rbf_aadt <- inla_cv_rbf_aadt$summary.fitted.values[test_idx, ]
fv_rbf      <- inla_cv_rbf$summary.fitted.values[test_idx, ]

pred_besag    <- exp(fv_besag$mean)
pred_rbf_aadt <- exp(fv_rbf_aadt$mean)
pred_rbf      <- exp(fv_rbf$mean)

# Utility functions
rmse     <- function(p, o) sqrt(mean((p - o)^2))
mae      <- function(p, o) mean(abs(p - o))
coverage <- function(lo, hi, o) mean(o >= lo & o <= hi)

cat("\n")
cat(strrep("=", 82), "\n")
cat("TEST-SET PERFORMANCE COMPARISON (", n_test, " held-out measured links)\n")
cat(strrep("=", 82), "\n\n")

cat(sprintf("%-32s %15s %16s %20s\n",
    "Metric", "Ordinary Besag", "RBF (AADT only)", "RBF (all covariates)"))
cat(strrep("-", 82), "\n")

cat(sprintf("%-32s %15.2f %16.2f %20.2f\n", "RMSE",
    rmse(pred_besag, true_aadt),
    rmse(pred_rbf_aadt, true_aadt),
    rmse(pred_rbf, true_aadt)))

cat(sprintf("%-32s %15.2f %16.2f %20.2f\n", "MAE",
    mae(pred_besag, true_aadt),
    mae(pred_rbf_aadt, true_aadt),
    mae(pred_rbf, true_aadt)))

cat(sprintf("%-32s %15.4f %16.4f %20.4f\n", "Correlation (pred vs obs)",
    cor(pred_besag, true_aadt),
    cor(pred_rbf_aadt, true_aadt),
    cor(pred_rbf, true_aadt)))

cat(sprintf("%-32s %14.1f%% %15.1f%% %19.1f%%\n", "95% CI Coverage",
    100 * coverage(exp(fv_besag$`0.025quant`),    exp(fv_besag$`0.975quant`),    true_aadt),
    100 * coverage(exp(fv_rbf_aadt$`0.025quant`), exp(fv_rbf_aadt$`0.975quant`), true_aadt),
    100 * coverage(exp(fv_rbf$`0.025quant`),      exp(fv_rbf$`0.975quant`),      true_aadt)))

cat(sprintf("%-32s %15.2f %16.2f %20.2f\n", "Mean 95% PI width",
    mean(exp(fv_besag$`0.975quant`)    - exp(fv_besag$`0.025quant`)),
    mean(exp(fv_rbf_aadt$`0.975quant`) - exp(fv_rbf_aadt$`0.025quant`)),
    mean(exp(fv_rbf$`0.975quant`)      - exp(fv_rbf$`0.025quant`))))

cat(sprintf("%-32s %15.1f %16.1f %20.1f\n", "DIC",
    inla_cv_besag$dic$dic, inla_cv_rbf_aadt$dic$dic, inla_cv_rbf$dic$dic))

cat(sprintf("%-32s %15.1f %16.1f %20.1f\n", "WAIC",
    inla_cv_besag$waic$waic, inla_cv_rbf_aadt$waic$waic, inla_cv_rbf$waic$waic))

cat(strrep("=", 82), "\n\n")

# Differences relative to ordinary Besag (negative = improvement)
models_cv <- list(
  "RBF (AADT only)"      = list(pred = pred_rbf_aadt, cv = inla_cv_rbf_aadt),
  "RBF (all covariates)" = list(pred = pred_rbf,      cv = inla_cv_rbf)
)
cat("Differences vs ordinary Besag (negative = improvement):\n")
for (nm in names(models_cv)) {
  p  <- models_cv[[nm]]$pred
  cv <- models_cv[[nm]]$cv
  cat(sprintf("  %-22s  RMSE %+.2f  MAE %+.2f  DIC %+.1f  WAIC %+.1f\n",
      nm,
      rmse(p, true_aadt) - rmse(pred_besag, true_aadt),
      mae(p, true_aadt)  - mae(pred_besag, true_aadt),
      cv$dic$dic   - inla_cv_besag$dic$dic,
      cv$waic$waic - inla_cv_besag$waic$waic))
}
cat("\n")

# ────────────────────────────────────────────────────────────────────────────────────
# 9. VISUALIZATION: Predicted vs Observed
# ────────────────────────────────────────────────────────────────────────────────────

lims <- range(c(true_aadt, pred_besag, pred_rbf_aadt, pred_rbf), na.rm = TRUE)

par(mfrow = c(1, 3), pty = "s")

plot(true_aadt, pred_besag, pch = 16, col = rgb(0, 0, 1, 0.4),
     xlim = lims, ylim = lims,
     xlab = "Observed AADT", ylab = "Predicted AADT",
     main = "Ordinary Besagproper")
abline(0, 1, col = "red", lty = 2, lwd = 2)

plot(true_aadt, pred_rbf_aadt, pch = 16, col = rgb(0.8, 0.4, 0, 0.4),
     xlim = lims, ylim = lims,
     xlab = "Observed AADT", ylab = "Predicted AADT",
     main = "RBF Besag (AADT only)")
abline(0, 1, col = "red", lty = 2, lwd = 2)

plot(true_aadt, pred_rbf, pch = 16, col = rgb(0, 0.6, 0, 0.4),
     xlim = lims, ylim = lims,
     xlab = "Observed AADT", ylab = "Predicted AADT",
     main = "RBF Besag (all covariates)")
abline(0, 1, col = "red", lty = 2, lwd = 2)

par(mfrow = c(1, 1))

# ────────────────────────────────────────────────────────────────────────────────────
# 10. GREEDY SENSOR SELECTION WITH RBF COVARIANCE
# ────────────────────────────────────────────────────────────────────────────────────
# Greedily select k sensor placements that maximise information gain.
#
# For each candidate y in the unmeasured set U, the selection criterion is:
#   delta(y) = sigma^2(y | A) / sigma^2(y | U \ {y})
# where A is the current measured set and U is the current unmeasured set.
#
# Efficiency:
#   - Numerator: sigma^2(y | A) = Sigma_yy - Sigma_yA * Sigma_AA^{-1} * Sigma_Ay
#     Compute L_A = chol(Sigma_AA) once per iteration; vectorise over all y.
#   - Denominator: by the partitioned-inverse identity,
#       sigma^2(y | U \ {y}) = 1 / (Sigma_UU^{-1})_{yy}
#     Compute Sigma_UU^{-1} once per iteration; read diagonal.
#   - Total cost: O(k * max(|A|^3, |U|^3)) instead of O(k * |U| * max(|A|^3, |U|^3))

Sigma <- covariance_rbf_sum

measured_paired     <- which(!is.na(prepared_traffic_links[paired_with_idx, ]$aadt))
unmeasured_paired   <- which( is.na(prepared_traffic_links[paired_with_idx, ]$aadt))
measured_unpaired   <- which(!is.na(prepared_traffic_links[unpaired_with_idx, ]$aadt))
unmeasured_unpaired <- which( is.na(prepared_traffic_links[unpaired_with_idx, ]$aadt))

measured_idxs   <- c(measured_paired, measured_unpaired + length(paired_bases))
unmeasured_idxs <- c(unmeasured_paired, unmeasured_unpaired + length(paired_bases))

k <- 10
selected_sensors <- integer(k)

cat("Greedy sensor selection using RBF covariance matrix (k =", k, "):\n\n")

for (i in seq_len(k)) {
  U <- unmeasured_idxs
  A <- measured_idxs

  # Numerator: conditional variance of each candidate given measured set
  # sigma^2(y | A) = Sigma_yy - Sigma_yA L_A^{-T} L_A^{-1} Sigma_Ay
  L_A   <- chol(Sigma[A, A])                            # upper Cholesky
  V     <- backsolve(L_A, Sigma[A, U], transpose = TRUE) # L_A^{-T} * Sigma_AU
  numer <- diag(Sigma[U, U]) - colSums(V^2)

  # Denominator: conditional variance of each candidate given all other unmeasured
  # sigma^2(y | U \ {y}) = 1 / (Sigma_UU^{-1})_{yy}
  Prec_UU <- chol2inv(chol(Sigma[U, U]))
  denom   <- 1 / diag(Prec_UU)

  delta <- numer / denom

  best         <- which.max(delta)
  selected_idx <- U[best]
  selected_sensors[i] <- selected_idx

  cat(sprintf("  Iteration %2d: sum_idx %4d  (delta = %.4f)\n",
              i, selected_idx, delta[best]))

  measured_idxs   <- c(measured_idxs, selected_idx)
  unmeasured_idxs <- unmeasured_idxs[-best]
}

# Decode sum-space indices back to road IDs
paired_selected   <- selected_sensors[selected_sensors <= length(paired_bases)]
unpaired_selected <- selected_sensors[selected_sensors >  length(paired_bases)]

cat("\nSelected paired road segment base IDs:\n")
print(paired_bases[paired_selected])
cat("\nSelected unpaired WITH road IDs:\n")
print(unpaired_bases[unpaired_selected - length(paired_bases)])

cat("\nSensor selection complete.\n")
