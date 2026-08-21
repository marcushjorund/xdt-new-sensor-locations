# ════════════════════════════════════════════════════════════════════════════════
# compute_weights() ----
# ════════════════════════════════════════════════════════════════════════════════

#' Compute location importance weights for sensor selection.
#'
#' @param data                Traffic links data frame.
#' @param weighting_vars      Named list; names are column names in \code{data}.
#'   Each element is either:
#'   \describe{
#'     \item{Named numeric vector}{Lookup table by level (as character); e.g.
#'       \code{c("0"=8, "1"=7, ...)} for \code{functionalRoadClass}.
#'       Levels not found in the lookup default to 1.0.}
#'     \item{\code{"proportional"}}{Weights proportional to column values:
#'       \code{pmax(col / mean(col), 1e-6)}.  Column must be numeric.
#'       On high-variance columns (e.g. AADT) this can produce very large weights;
#'       prefer \code{"log_proportional"} in that case.}
#'     \item{\code{"log_proportional"}}{Log-compressed proportional weights:
#'       \code{log(pmax(col / mean(col), 0) + 1) + 1}.  Maps a 5× road to ~1.79
#'       rather than 5.0, keeping the weight range comparable to FRC-based
#'       lookups and making AADT a tiebreaker rather than a dominant driver.}
#'   }
#'   When \code{NULL} (default), all weights are 1.0 (no weighting).
#' @param weighting_vars_alpha Optional numeric vector of length
#'   \code{length(weighting_vars)} with non-negative mixing weights for each
#'   variable.  Values are normalised to sum to 1 internally, so only relative
#'   magnitudes matter.  Each raw weight vector is first min-max normalised to
#'   \eqn{[0, 1]} so that all variables operate on a common scale regardless
#'   of their original units, then combined as a weighted linear sum:
#'   \deqn{w_i = \sum_k \alpha_k \cdot \tilde{w}_{k,i}}{w_i = sum_k alpha_k * w_norm_{k,i}}
#'   where \eqn{\tilde{w}_{k,i} = (w_{k,i} - \min_k) / (\max_k - \min_k)}.
#'   When \code{NULL} (default), equal contribution (\eqn{\alpha_k = 1/K}) is
#'   used.  Because inputs are normalised first, \eqn{\alpha_k} truly reflects
#'   the fraction of score variance attributable to variable \eqn{k} — e.g.
#'   \code{c(0.9, 0.1)} means FRC accounts for 90\% of the spread regardless
#'   of whether the other variable is heavyRatio (0–1) or raw AADT (0–100k).
#' @param row_lookup          Optional integer index into \code{data} rows (e.g.
#'   the \code{-WITH} row indices in the base-ID covariance space).  When
#'   \code{NULL}, all rows of \code{data} are used directly.
#' @return Numeric vector of weights, length = \code{length(row_lookup)} (or
#'   \code{nrow(data)} when \code{row_lookup} is \code{NULL}).
compute_weights <- function(data, weighting_vars, weighting_vars_alpha = NULL,
                            row_lookup = NULL) {
  m <- if (is.null(row_lookup)) nrow(data) else length(row_lookup)
  if (is.null(weighting_vars)) return(rep(1.0, m))
  if (!is.list(weighting_vars) || is.null(names(weighting_vars)))
    stop("'weighting_vars' must be a named list.")

  w_list <- lapply(names(weighting_vars), function(col_name) {
    if (!col_name %in% names(data))
      stop("Column '", col_name, "' not found in data (weighting_vars).")
    spec     <- weighting_vars[[col_name]]
    col_vals <- if (is.null(row_lookup)) data[[col_name]] else data[[col_name]][row_lookup]

    if (identical(spec, "proportional")) {
      vals <- suppressWarnings(as.numeric(col_vals))
      mn   <- mean(vals, na.rm = TRUE)
      if (is.na(mn) || mn <= 0) {
        warning("Column '", col_name, "': mean is NA or <= 0; proportional weighting disabled.")
        w <- rep(1.0, m)
      } else {
        w <- pmax(vals / mn, 1e-6)
      }
      w[is.na(w)] <- 1.0
    } else if (identical(spec, "log_proportional")) {
      vals <- suppressWarnings(as.numeric(col_vals))
      mn   <- mean(vals, na.rm = TRUE)
      if (is.na(mn) || mn <= 0) {
        warning("Column '", col_name, "': mean is NA or <= 0; log_proportional weighting disabled.")
        w <- rep(1.0, m)
      } else {
        w <- log(pmax(vals / mn, 0) + 1) + 1
      }
      w[is.na(w)] <- 1.0
    } else if (identical(spec, "identity")) {
      vals <- suppressWarnings(as.numeric(col_vals))
      w    <- pmax(vals, 1e-6)
      w[is.na(w)] <- 1.0
    }
    else if (is.numeric(spec) && !is.null(names(spec))) {
      w         <- rep(1.0, m)
      lvl_chars <- as.character(col_vals)
      matched   <- spec[lvl_chars]
      found     <- !is.na(matched)
      if (any(matched[found] < 0, na.rm = TRUE))
        warning("weighting_vars[['", col_name, "']]: ",
                sum(matched[found] < 0, na.rm = TRUE),
                " lookup value(s) are negative and will be clamped to 0. ",
                "Negative weights are likely a mistake (e.g. log() applied to values <= 1).")
      w[found]  <- pmax(matched[found], 0)
    } else {
      stop("Each element of 'weighting_vars' must be a named numeric vector, ",
           "\"proportional\", or \"log_proportional\".")
    }
    w
  })

  # Min-max normalise each weight vector to [0, 1] so that variables on
  # different scales (e.g. FRC lookup vs raw AADT vs heavyRatio) all span the
  # same range before mixing.  A constant vector is set to 0.5 (neutral).
  # When only one variable is supplied, normalisation is a no-op (there is
  # nothing to mix) so we skip it and return the raw weights directly.
  K <- length(w_list)
  if (K == 1L) {
    w <- w_list[[1L]]
    return(pmax(w, 1e-8))
  }
  w_list <- lapply(w_list, function(w) {
    lo <- min(w, na.rm = TRUE)
    hi <- max(w, na.rm = TRUE)
    if (!is.finite(lo) || !is.finite(hi) || (hi - lo) < 1e-10)
      return(rep(0.5, length(w)))
    (w - lo) / (hi - lo)
  })

  # Weighted linear combination: w_i = sum_k alpha_k * w_norm_{k,i}.
  # alpha_k is the true fractional contribution of variable k (sums to 1).
  if (is.null(weighting_vars_alpha)) {
    alpha <- rep(1.0 / K, K)
  } else {
    if (length(weighting_vars_alpha) != K)
      stop("'weighting_vars_alpha' must have the same length as 'weighting_vars' (", K, ").")
    if (any(weighting_vars_alpha < 0))
      stop("All elements of 'weighting_vars_alpha' must be non-negative.")
    s <- sum(weighting_vars_alpha)
    if (s <= 0) stop("'weighting_vars_alpha' must not sum to zero.")
    alpha <- weighting_vars_alpha / s
  }
  w <- Reduce("+", mapply(function(wi, ai) ai * wi, w_list, alpha, SIMPLIFY = FALSE))
  pmax(w, 1e-8)
}

# ════════════════════════════════════════════════════════════════════════════════
# compute_neighbourhood_weight() ----
# ════════════════════════════════════════════════════════════════════════════════

#' Compute dynamic neighbourhood-importance weights for greedy MI selection.
#'
#' For each candidate location, its effective weight is
#' \code{f0*own + f1*NEIGH(1-hop) + ... + fH*NEIGH(H-hop)}, where NEIGH
#' aggregates neighbour importance according to \code{neighbourhood_scoring}.
#' Hop shells are computed in the \emph{residual graph}: columns corresponding to
#' \code{exclude_idx} (measured S0 + already-selected sensors) are zeroed in
#' \code{A1} before computing higher hops, so paths that pass through a measured
#' or already-selected node are blocked.  This enforces the Markov conditional
#' independence property — a node reachable only via a measured bridge provides
#' no additional information beyond what that bridge already supplies.
#'
#' @param w_base      Numeric vector of base composite importance weights from
#'   \code{compute_weights()}.
#' @param A1          Binary sparse m×m 1-hop base-space adjacency (diagonal 0).
#' @param decay       Numeric vector of any length \code{c(f0, f1, ..., fH)}.
#'   \code{f0} weights the candidate itself; \code{f1} weights the 1-hop exact
#'   shell; ...; \code{fH} weights the H-hop exact shell.  Must be non-negative
#'   with a positive sum.  Normalised internally to sum 1.
#' @param exclude_idx Integer index of m-space locations to exclude from neighbour
#'   sums and to block as bridges in the residual graph.
#' @param neighbourhood_scoring One of \code{"mean"} (default — per-degree mean),
#'   \code{"sum"} (raw sum), \code{"proportional"} (sum / per-hop baseline),
#'   or \code{"log_proportional"} (log-compressed, min-max normalised within hop).
#' @param neigh_baselines Numeric vector of per-hop baselines for
#'   \code{"proportional"} and \code{"log_proportional"} normalisation (length =
#'   number of hops = \code{length(decay) - 1}).  \code{neigh_baselines[h]} is
#'   the partition-level mean of the h-hop raw sum over unmeasured candidates,
#'   computed once before the greedy loop.  Falls back to \code{mean(s_h)} when
#'   \code{NULL} or shorter than \code{max_hops}.
#' @return Numeric vector of effective weights, same length as \code{w_base},
#'   floored at \code{1e-8}.
compute_neighbourhood_weight <- function(w_base, A1, decay,
                                         exclude_idx           = integer(0),
                                         neighbourhood_scoring = c("mean", "log_proportional", "proportional", "sum"),
                                         neigh_baselines       = NULL) {
  if (!is.numeric(decay) || any(!is.finite(decay)) || any(decay < 0))
    stop("'neighbourhood_decay' must be a non-negative finite numeric vector.")
  if (length(decay) < 1L || sum(decay) <= 0)
    stop("'neighbourhood_decay' must have at least one positive element.")
  f <- decay / sum(decay)

  neighbourhood_scoring <- match.arg(neighbourhood_scoring)

  max_hops <- length(f) - 1L  # f[1] = self; f[2] = 1-hop; ...; f[H+1] = H-hop

  w_known <- w_base
  if (length(exclude_idx) > 0L) w_known[exclude_idx] <- 0

  if (max_hops == 0L) return(pmax(w_base, 1e-8))

  # Residual A1: zero excluded columns so measured/selected nodes cannot act as
  # bridges. This enforces the Markov conditional independence property of the
  # BesagProper CAR model — a node reachable only via a measured location carries
  # no additional information beyond what that measured location already provides.
  A1_res <- A1
  if (length(exclude_idx) > 0L) A1_res[, exclude_idx] <- 0

  if (neighbourhood_scoring == "mean") {
    one_known <- rep(1.0, length(w_base))
    if (length(exclude_idx) > 0L) one_known[exclude_idx] <- 0
  }

  # Iteratively build exact h-hop shell signals in the residual graph.
  # A_power accumulates A1_res^h; A_cumul tracks all nodes reached in <= h hops.
  A_power <- A1_res
  A_cumul <- sign(A1_res); Matrix::diag(A_cumul) <- 0   # exact 1-hop shell

  eff <- f[1L] * w_base

  for (h in seq_len(max_hops)) {
    if (h > 1L) {
      A_power <- A_power %*% A1_res              # A1_res^h
      A_h     <- sign(A_power); Matrix::diag(A_h) <- 0
      A_exact <- A_h; A_exact[A_cumul > 0] <- 0  # exact h-hop shell (no overlap)
      A_cumul <- sign(A_cumul + A_exact)
    } else {
      A_exact <- A_cumul                         # h = 1: exact 1-hop shell
    }

    s_h <- as.numeric(A_exact %*% w_known)

    mu_h <- if (!is.null(neigh_baselines) && h <= length(neigh_baselines))
              pmax(neigh_baselines[[h]], 1e-8)
            else NULL

    neigh_h <- switch(neighbourhood_scoring,
      mean = {
        deg_h <- as.numeric(A_exact %*% one_known)
        s_h / pmax(deg_h, 1)
      },
      sum          = s_h,
      proportional = s_h / pmax(if (!is.null(mu_h)) mu_h else mean(s_h + 1e-8), 1e-8),
      log_proportional = {
        mu   <- pmax(if (!is.null(mu_h)) mu_h else mean(s_h + 1e-8), 1e-8)
        lp_h <- log(s_h / mu + 1) + 1
        lo_h <- min(lp_h); hi_h <- max(lp_h)
        if (hi_h > lo_h) (lp_h - lo_h) / (hi_h - lo_h) else rep(0.5, length(lp_h))
      }
    )

    eff <- eff + f[h + 1L] * neigh_h
  }
  pmax(eff, 1e-8)
}

# ════════════════════════════════════════════════════════════════════════════════
# greedy_mi_sensor_selection() ----
# ════════════════════════════════════════════════════════════════════════════════

#' Greedy mutual-information sensor selection (Krause et al. 2008)
#'
#' Greedily selects `k` unmeasured locations maximising mutual information gain, using
#' incremental Cholesky updates for the conditional variance of already-selected sensors
#' and Schur-complement updates for the conditional variance of remaining candidates, so
#' the covariance/precision matrices are never fully recomputed at each iteration.
#'
#' @param data              Traffic links data frame; `NA` in `aadt` marks unmeasured links.
#' @param covariance_matrix m x m covariance matrix (m = `length(ids)`), e.g. the
#'   `covariance_matrix_sum` from `create_covariance_and_precision_matrix(with_and_against = TRUE)`.
#' @param ids               Character vector of base road IDs corresponding to rows/columns
#'   of `covariance_matrix`.
#' @param k                 Number of sensors to select.
#' @param adjacency_matrix  Optional sparse n x n directed-link adjacency matrix; required
#'   for `neighbourhood_decay` smoothing and `append_neighbours`.
#' @param weighting_vars       Named list forwarded to `compute_weights()`. `NULL` = uniform
#'   weights.
#' @param weighting_vars_alpha Optional numeric mixing weights for `weighting_vars`,
#'   normalised to sum 1. `NULL` = equal contribution.
#' @param neighbourhood_decay  Numeric vector `c(f0, f1, ..., fH)` forwarded to
#'   `compute_neighbourhood_weight()`. `NULL` disables neighbourhood smoothing.
#' @param neighbourhood_scoring One of `"mean"` (default), `"sum"`, `"proportional"`, or
#'   `"log_proportional"`; forwarded to `compute_neighbourhood_weight()`.
#' @param append_neighbours    Logical (default `TRUE`); append 1-hop adjacency neighbours
#'   of selected sensors to the output for context/plotting.
#' @param eligible_ids         Optional character vector restricting candidates to a subset
#'   of `ids` (e.g. from a `filtering` spec). `NULL` = all unmeasured locations eligible.
#' @return A list with `selected_ids`, `selected_idx`, `selected_data_entries` (data frame
#'   with `selected`/`mi_score` columns, WITH + AGAINST rows, plus neighbour rows when
#'   `append_neighbours = TRUE`), `mi_scores`, `measured_idx`, `unmeasured_idx_initial`.
greedy_mi_sensor_selection <- function(data, covariance_matrix, ids, k, adjacency_matrix = NULL,
                                       weighting_vars        = NULL,
                                       weighting_vars_alpha  = NULL,
                                       neighbourhood_decay   = NULL,
                                       neighbourhood_scoring = "mean",
                                       append_neighbours     = TRUE,
                                       eligible_ids          = NULL) {
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

  # ── Restrict candidates to eligible_ids (filtering) ──────────────────────
  if (!is.null(eligible_ids)) {
    eligible_m_idx <- which(ids %in% eligible_ids)
    candidates_eligible <- intersect(unmeasured_idx, eligible_m_idx)
    if (length(candidates_eligible) == 0L)
      stop("No unmeasured eligible locations found after applying filter.")
  } else {
    candidates_eligible <- unmeasured_idx
  }
  k <- min(k, length(candidates_eligible))

  # ── Build location weights ────────────────────────────────────────────────
  w_base <- compute_weights(data, weighting_vars, weighting_vars_alpha, row_lookup = with_row)

  # ── Neighbourhood adjacency matrix (built once; hop shells computed per step) ─
  # neighbourhood_decay = c(f0, f1, ..., fH): f0 weights the candidate itself,
  # f1 the 1-hop shell, ..., fH the H-hop shell. Exact shells and Markov-blocking
  # of already-selected nodes are handled inside compute_neighbourhood_weight().
  A1_neigh <- NULL
  if (!is.null(neighbourhood_decay) && !is.null(adjacency_matrix)) {
    against_row <- match(paste0(ids, "-AGAINST"), data$id)
    row_idx  <- c(seq_len(m)[!is.na(with_row)],  seq_len(m)[!is.na(against_row)])
    col_idx  <- c(with_row[!is.na(with_row)],     against_row[!is.na(against_row)])
    C        <- Matrix::sparseMatrix(i = row_idx, j = col_idx, x = 1.0,
                                     dims = c(m, nrow(adjacency_matrix)))
    A_base_m <- C %*% adjacency_matrix %*% Matrix::t(C)
    A1_neigh <- sign(A_base_m); Matrix::diag(A1_neigh) <- 0
  }

  # ── Neighbourhood baselines for proportional/log_proportional scoring ─────
  # Computed at step 0 using the residual graph (measured_idx excluded as
  # bridges) so baselines reflect the typical within-partition connectivity
  # of unmeasured candidates only, mirroring what compute_neighbourhood_weight()
  # will see at greedy step 1.
  neigh_baselines <- NULL
  if (!is.null(A1_neigh) && neighbourhood_scoring %in% c("log_proportional", "proportional")) {
    max_hops_b  <- length(neighbourhood_decay) - 1L
    A1_init_res <- A1_neigh
    if (length(measured_idx) > 0L) A1_init_res[, measured_idx] <- 0
    A_power_b <- A1_init_res
    A_cumul_b <- sign(A1_init_res); Matrix::diag(A_cumul_b) <- 0
    s_init <- vector("list", max_hops_b)
    s_init[[1L]] <- as.numeric(A_cumul_b %*% w_base)
    for (h in seq_len(max_hops_b - 1L)) {
      A_power_b <- A_power_b %*% A1_init_res
      A_h_b     <- sign(A_power_b); Matrix::diag(A_h_b) <- 0
      A_exact_b <- A_h_b; A_exact_b[A_cumul_b > 0] <- 0
      s_init[[h + 1L]] <- as.numeric(A_exact_b %*% w_base)
      A_cumul_b <- sign(A_cumul_b + A_exact_b)
    }
    neigh_baselines <- sapply(seq_len(max_hops_b),
                              function(h) mean(s_init[[h]][unmeasured_idx]))
  }
  
  # ── Warm-start: Var(x_i | S0) = Sigma_ii - colSums(A^2), A = L_{S0}^{-T} Sigma_{S0,.} ─
  cond_var1 <- diag(covariance_matrix)
  if (length(measured_idx) > 0) {
    Sigma_mm  <- covariance_matrix[measured_idx, measured_idx]
    Sigma_om  <- covariance_matrix[, measured_idx]
    L_mm      <- chol(Sigma_mm)
    A         <- forwardsolve(t(L_mm), t(Sigma_om))   # n_s0 × m
    cond_var1 <- pmax(cond_var1 - colSums(A^2), 0)
  } else {
    A <- matrix(0, 0, m)
  }
  cond_var1[measured_idx] <- 0
  
  # ── Precision of unmeasured candidates conditioned on S0 ─────────────────
  Sigma_uu <- covariance_matrix[unmeasured_idx, unmeasured_idx]
  if (length(measured_idx) > 0) {
    A_uu     <- A[, unmeasured_idx, drop = FALSE]
    Sigma_uu <- Sigma_uu - crossprod(A_uu)
  }
  prec <- chol2inv(chol(Sigma_uu))
  
  # cond_var2[i] = 1/diag(prec_{uu|S0}); Inf sentinel for non-candidates
  cond_var2 <- rep(Inf, m)
  cond_var2[unmeasured_idx] <- 1 / diag(prec)
  # ── Greedy MI loop ──────────────────────────────────────────────────────────
  # candidates is restricted to eligible unmeasured locations; cond_var1/cond_var2
  # are still computed over all unmeasured so spatial correlations remain correct.
  candidates   <- candidates_eligible
  prec_remaining <- unmeasured_idx
  selected_idx <- integer(k)
  mi_scores    <- numeric(k)
  L2           <- matrix(0, m, k)
  
  for (i in seq_len(k)) {
    # MI gain ∝ Var(x_s | A) / Var(x_s | V_u \ {s}, S0)  (Krause et al. 2008)
    # Dynamic neighbourhood weight: exclude measured S0 + sensors selected so far
    # so the bonus decreases as the algorithm commits, distributing selections.
    if (!is.null(A1_neigh)) {
      excl  <- c(measured_idx, selected_idx[seq_len(i - 1L)])
      w_use <- compute_neighbourhood_weight(w_base, A1_neigh,
                                            neighbourhood_decay, excl,
                                            neighbourhood_scoring,
                                            neigh_baselines)
    } else {
      w_use <- w_base
    }
    score       <- (cond_var1 / cond_var2) * w_use
    best_global <- candidates[which.max(score[candidates])]
    if (length(best_global) == 0L || is.na(best_global)) {
      warning("greedy_mi_sensor_selection: no valid candidate at iteration ", i,
              ". All scores are NaN/NA — check weighting_vars for non-positive values.",
              " Stopping early with ", i - 1L, " sensors selected.")
      selected_idx <- selected_idx[seq_len(i - 1L)]
      mi_scores    <- mi_scores[seq_len(i - 1L)]
      break
    }
    j           <- which(candidates == best_global)
    
    selected_idx[i] <- best_global
    mi_scores[i]    <- score[best_global]
    
    # Incremental Cholesky column conditioned on S0 and previously selected
    col_vec <- covariance_matrix[, best_global] - crossprod(A, A[, best_global])
    if (i > 1)
      col_vec <- col_vec - L2[, seq_len(i - 1), drop = FALSE] %*%
      L2[best_global, seq_len(i - 1)]
    pivot       <- pmax(col_vec[best_global], .Machine$double.eps * 100)
    L2[, i]     <- col_vec / sqrt(pivot)
    cond_var1   <- pmax(cond_var1 - L2[, i]^2, 0)
    cond_var1[best_global] <- 0
    if (length(candidates) > 1) {
      j_prec <- which(prec_remaining == best_global)
      prec   <- prec[-j_prec, -j_prec, drop = FALSE] -
        tcrossprod(prec[-j_prec, j_prec]) / prec[j_prec, j_prec]
      prec_remaining <- prec_remaining[-j_prec]
    }
    candidates <- candidates[-j]
    if (length(candidates) > 0) {
      cond_var2[candidates] <- 1 / diag(prec)
    }
    cond_var2[best_global] <- Inf
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
  if (!is.null(adjacency_matrix) && isTRUE(append_neighbours)) {
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

# ════════════════════════════════════════════════════════════════════════════════
# greedy_mi_sensor_selection_norway() ----
# ════════════════════════════════════════════════════════════════════════════════

#' Greedy MI sensor selection for the full Norwegian network.
#'
#' Partitions by county (with \code{hops}-hop border expansion), builds the
#' covariance per partition, and runs greedy MI selection independently per
#' county.  Border links appearing in multiple partitions are deduplicated by
#' keeping the highest \code{mi_score}.  Optionally prunes to the top-\code{r}
#' locations Norway-wide and appends multi-hop adjacency neighbours of the
#' final selected set.
#'
#' Two MI-scoring modes are available via \code{bundle_scoring}:
#' \describe{
#'   \item{\code{bundle_scoring = FALSE} (default)}{Standard mode.  Works in
#'     m×m with+against covariance space.  Each sensor base-ID is scored
#'     individually.  If \code{data} contains the \code{enables_derivable_links}
#'     column (from \code{identify_derivable_nodes()}), flow-conservation derived
#'     links are appended post-selection as informational rows
#'     (\code{is_derived = TRUE}) without affecting scores.}
#'   \item{\code{bundle_scoring = TRUE}}{Bundle mode.  Works in n×n directed
#'     covariance space (no with+against summing).  Each greedy step scores the
#'     full bundle — placed sensor + any WITH/AGAINST + flow-derived links —
#'     jointly.  Requires \code{enables_derivable_links} column.}
#' }
#'
#' @param data             Data frame with a \code{county} factor column.
#'   When \code{bundle_scoring = TRUE}, must also contain
#'   \code{enables_derivable_links} list-column from
#'   \code{identify_derivable_nodes()}.
#' @param adjacency_matrix Nationwide sparse adjacency matrix (\code{nrow(data)} x \code{nrow(data)}).
#' @param distances        Edge Gower-distance vector (\code{fit_weighted_inla_model()$distances}).
#' @param tau,d,sigma      CAR hyperparameters on the natural scale.
#' @param hops             Border expansion hops in \code{partition_by_county()} (default 2).
#' @param k                Sensors (or bundles) selected per county (default 10).
#' @param r                If non-NULL, retain only the top-\code{r} sensors/bundles Norway-wide.
#' @param weighting_vars       Named list forwarded to \code{compute_weights()}.  \code{NULL} = uniform.
#' @param weighting_vars_alpha Optional numeric vector of exponents (one per entry in
#'   \code{weighting_vars}); normalised to sum 1.  \code{NULL} = equal contribution.
#' @param neighbourhood_decay  Numeric vector of any length \code{c(f0, f1, ..., fH)}.
#'   \code{f0} weights the candidate itself; \code{f1} the 1-hop exact shell; ...;
#'   \code{fH} the H-hop exact shell.  Hop shells are computed in the residual graph
#'   (measured and already-selected nodes excluded as bridges) to enforce the Markov
#'   conditional independence property.  \code{NULL} disables neighbourhood
#'   smoothing.
#' @param neighbourhood_scoring One of \code{"mean"} (default, per-degree mean),
#'   \code{"sum"} (raw sum), \code{"proportional"} (sum / county-mean, no
#'   compression), or \code{"log_proportional"} (log-compressed county-normalised
#'   sum; see \code{compute_neighbourhood_weight()}).
#'   Forwarded to both \code{greedy_mi_sensor_selection()} and
#'   \code{greedy_mi_sensor_selection_derivable()}.
#' @param neighbour_hops   Number of context hops to append for plotting (default \code{0L}).
#' @param bundle_scoring   Logical (default \code{FALSE}).  See Details.
#' @return List: \code{selected_data_entries} (data frame), \code{n_counties},
#'   \code{k_per_county}, \code{r}, \code{summary}.
greedy_mi_sensor_selection_norway <- function(
    data,
    adjacency_matrix,
    distances,
    tau,
    d,
    sigma,
    hops                  = 2,
    k                     = 10,
    r                     = NULL,
    weighting_vars        = NULL,
    weighting_vars_alpha  = NULL,
    neighbourhood_decay   = NULL,
    neighbourhood_scoring = "mean",
    neighbour_hops        = 0L,
    bundle_scoring        = FALSE,
    filtering             = NULL,
    verbose               = FALSE) {

  if (!"county" %in% names(data) || !is.factor(data$county))
    stop("'data$county' must be a factor.")
  if (isTRUE(bundle_scoring) && !"enables_derivable_links" %in% names(data))
    stop("'data' must contain 'enables_derivable_links' (from identify_derivable_nodes()) when bundle_scoring = TRUE.")
  if (!is.null(r) && (!is.numeric(r) || length(r) != 1 || r < 1))
    stop("'r' must be a positive integer scalar or NULL.")

  # ── Compute eligible base IDs from filtering spec ─────────────────────────
  # filtering = list(col1 = values1, col2 = values2, ...): a link is eligible
  # when data[[col]] %in% values for every named entry. NULL = all links eligible.
  if (!is.null(filtering)) {
    if (!is.list(filtering) || is.null(names(filtering)))
      stop("'filtering' must be a named list, e.g. list(county = 'Nordland', roadCategory = 'FYLKESVEG').")
    keep <- rep(TRUE, nrow(data))
    for (col in names(filtering)) {
      if (!col %in% names(data))
        stop("filtering column '", col, "' not found in data.")
      keep <- keep & (data[[col]] %in% filtering[[col]])
    }
    eligible_base_ids <- unique(sub("-WITH$|-AGAINST$", "", data$id[keep]))
    if (length(eligible_base_ids) == 0L)
      stop("No links match the supplied filter: ", paste(names(filtering), filtering, sep = " = ", collapse = ", "))
    message("Filtering active: ", sum(keep), " links eligible across ",
            length(eligible_base_ids), " base IDs.")
    base_id_vec       <- sub("-WITH$|-AGAINST$", "", data$id)
    eligible_counties <- unique(as.character(data$county[base_id_vec %in% eligible_base_ids]))
  } else {
    eligible_base_ids <- NULL
    eligible_counties <- NULL
  }

  partition  <- partition_by_county(data, adjacency_matrix, distances, hops = hops)
  n_counties <- length(partition)

  # Only process partitions whose core county has at least one eligible link.
  # When filtering = NULL every partition is active.
  active_indices <- if (!is.null(eligible_counties)) {
    idx <- which(names(partition) %in% eligible_counties)
    message("Active partitions (", length(idx), "/", n_counties, "): ",
            paste(sort(names(partition)[idx]), collapse = ", "))
    idx
  } else {
    seq_len(n_counties)
  }

  per_county_results <- vector("list", n_counties)
  for (i in active_indices) {
    message("[", i, "/", n_counties, "] ", names(partition)[i])
    part <- partition[[i]]

    # ── Determine eligible IDs for this partition ────────────────────────────
    # Eligibility is checked against the partition's CORE links only (those
    # whose county matches the partition name), not the full border-expanded
    # partition.  This prevents a border-adjacent partition (e.g. Trøndelag)
    # from treating Nordland links that were included for spatial context as
    # eligible candidates and re-selecting them in the wrong covariance space.
    if (!is.null(eligible_base_ids)) {
      core_idx          <- which(part$data$county == names(partition)[i])
      core_base_ids     <- unique(sub("-WITH$|-AGAINST$", "", part$data$id[core_idx]))
      part_eligible_ids <- intersect(eligible_base_ids, core_base_ids)
      if (length(part_eligible_ids) == 0L) {
        per_county_results[[i]] <- NULL
        next
      }
    } else {
      part_eligible_ids <- NULL
    }

    if (!bundle_scoring) {
      # ── Normal mode: m×m with+against covariance ──────────────────────────
      cov <- create_covariance_and_precision_matrix(
        adjacency_matrix = part$adjacency_matrix,
        tau = tau, d = d, sigma = sigma,
        distances        = part$distances,
        data             = part$data,
        with_and_against = TRUE
      )
      res <- greedy_mi_sensor_selection(
        data                  = part$data,
        covariance_matrix     = cov$covariance_matrix_sum,
        ids                   = cov$ids,
        k                     = k,
        adjacency_matrix      = part$adjacency_matrix,
        weighting_vars        = weighting_vars,
        weighting_vars_alpha  = weighting_vars_alpha,
        neighbourhood_decay   = neighbourhood_decay,
        neighbourhood_scoring = neighbourhood_scoring,
        append_neighbours     = FALSE,
        eligible_ids          = part_eligible_ids
      )
      per_county_results[[i]] <- res$selected_data_entries

    } else {
      # ── Bundle mode: n×n directed covariance ──────────────────────────────
      cov <- create_covariance_and_precision_matrix(
        adjacency_matrix = part$adjacency_matrix,
        tau = tau, d = d, sigma = sigma,
        distances        = part$distances,
        data             = part$data,
        with_and_against = FALSE
      )
      res <- greedy_mi_sensor_selection_derivable(
        data                  = part$data,
        covariance_matrix     = cov$covariance_matrix,
        k                     = k,
        adjacency_matrix      = part$adjacency_matrix,
        weighting_vars        = weighting_vars,
        weighting_vars_alpha  = weighting_vars_alpha,
        neighbourhood_decay   = neighbourhood_decay,
        neighbourhood_scoring = neighbourhood_scoring,
        append_neighbours     = FALSE,
        eligible_ids          = part_eligible_ids
      )
      # Offset selection_order so bundle IDs are globally unique across counties.
      res_df <- res$selected_data_entries
      ssel   <- !is.na(res_df$selection_order)
      res_df$selection_order[ssel] <- (i - 1L) * k + res_df$selection_order[ssel]
      per_county_results[[i]] <- res_df
    }
  }

  # Combine; NULL entries arise from filtered-out partitions.
  # Border links duplicated across partitions → keep highest mi_score.
  all_selected <- do.call(rbind, Filter(Negate(is.null), per_county_results))
  all_selected <- all_selected[order(all_selected$id, -all_selected$mi_score), ]
  all_selected <- all_selected[!duplicated(all_selected$id), ]

  # Optional top-r ranking across Norway.
  if (!is.null(r)) {
    if (!bundle_scoring) {
      # Normal mode: r counts SENSORS (base IDs).
      sel_true    <- all_selected[all_selected$selected %in% TRUE, ]
      base_ids    <- sub("-WITH$|-AGAINST$", "", sel_true$id)
      base_scores <- tapply(sel_true$mi_score, base_ids, max, na.rm = TRUE)
      top_bases   <- names(sort(base_scores, decreasing = TRUE))[seq_len(min(r, length(base_scores)))]
      keep        <- all_selected$selected %in% FALSE |
        sub("-WITH$|-AGAINST$", "", all_selected$id) %in% top_bases
      all_selected <- all_selected[keep, ]
    } else {
      # Bundle mode: r counts BUNDLES (selection_order).
      sel_true      <- all_selected[all_selected$selected %in% TRUE, ]
      bundle_scores <- tapply(sel_true$mi_score, sel_true$selection_order,
                              max, na.rm = TRUE)
      top_orders    <- names(sort(bundle_scores, decreasing = TRUE))[
        seq_len(min(r, length(bundle_scores)))]
      keep <- all_selected$selected %in% FALSE |
        as.character(all_selected$selection_order) %in% top_orders
      all_selected <- all_selected[keep, ]
    }
  }

  # ── Derivable annotation (normal mode only) ────────────────────────────────
  # If enables_derivable_links exists, append flow-conservation derived links
  # as informational rows after selection.  Operates at data-ID level; no
  # covariance access needed regardless of which space MI was scored in.
  if (!bundle_scoring && "enables_derivable_links" %in% names(data)) {
    # Assign selection_order per SENSOR (base ID), not per directed link.
    # Both -WITH and -AGAINST of the same sensor share the same order number
    # so ranks run 1..r (= number of unique physical sensors selected).
    sel_rows   <- which(all_selected$selected %in% TRUE)
    sel_base   <- sub("-WITH$|-AGAINST$", "", all_selected$id[sel_rows])
    base_scores  <- tapply(all_selected$mi_score[sel_rows], sel_base, max, na.rm = TRUE)
    base_ranked  <- names(sort(base_scores, decreasing = TRUE))
    base_rank_lk <- stats::setNames(seq_along(base_ranked), base_ranked)
    all_selected$selection_order    <- NA_integer_
    all_selected$selection_order[sel_rows] <- base_rank_lk[sel_base]
    all_selected$is_derived         <- FALSE
    all_selected$n_derivable_gained <- NA_integer_

    derived_list <- vector("list", length(sel_rows))
    for (idx in seq_along(sel_rows)) {
      row_idx     <- sel_rows[idx]
      link_id     <- all_selected$id[row_idx]
      data_row    <- match(link_id, data$id)
      if (is.na(data_row)) next
      enabled_ids <- data$enables_derivable_links[[data_row]]
      if (is.null(enabled_ids) || length(enabled_ids) == 0L) next
      enabled_rows <- match(enabled_ids, data$id)
      enabled_rows <- enabled_rows[!is.na(enabled_rows)]
      enabled_rows <- enabled_rows[!data$id[enabled_rows] %in% all_selected$id]
      if (length(enabled_rows) == 0L) next
      drv                    <- data[enabled_rows, ]
      drv$selected           <- TRUE
      drv$is_derived         <- TRUE
      drv$mi_score           <- NA_real_
      drv$selection_order    <- all_selected$selection_order[row_idx]
      drv$n_derivable_gained <- NA_integer_
      all_selected$n_derivable_gained[row_idx] <- nrow(drv)
      derived_list[[idx]] <- drv
    }
    derived_all <- do.call(rbind, derived_list)
    if (!is.null(derived_all) && nrow(derived_all) > 0L)
      all_selected <- rbind(all_selected, derived_all)
  }

  # ── Ensure consistent columns across both modes ────────────────────────────
  if (!"selection_order"    %in% names(all_selected)) all_selected$selection_order    <- NA_integer_
  if (!"is_derived"         %in% names(all_selected)) all_selected$is_derived         <- FALSE
  if (!"n_derivable_gained" %in% names(all_selected)) all_selected$n_derivable_gained <- NA_integer_

  # ── Multi-hop neighbour expansion ─────────────────────────────────────────
  if (neighbour_hops > 0L) {
    all_selected$neighbour_hop <- NA_integer_
    frontier_ids <- all_selected$id
    seen_ids     <- frontier_ids
    for (hop in seq_len(neighbour_hops)) {
      neigh_list <- vector("list", n_counties)
      for (i in active_indices) {
        part          <- partition[[i]]
        front_in_part <- which(part$data$id %in% frontier_ids)
        if (length(front_in_part) == 0L) next
        adj_sub   <- part$adjacency_matrix[front_in_part, , drop = FALSE]
        neigh_idx <- which(Matrix::colSums(adj_sub != 0) > 0)
        neigh_idx <- neigh_idx[!part$data$id[neigh_idx] %in% seen_ids]
        if (length(neigh_idx) == 0L) next
        nr                    <- part$data[neigh_idx, ]
        nr$selected           <- FALSE
        nr$mi_score           <- NA_real_
        nr$selection_order    <- NA_integer_
        nr$is_derived         <- FALSE
        nr$n_derivable_gained <- NA_integer_
        nr$neighbour_hop      <- hop
        if (bundle_scoring) {
          nr$bundle_size   <- NA_integer_
          nr$derived_links <- vector("list", nrow(nr))
        }
        neigh_list[[i]] <- nr
      }
      hop_df <- do.call(rbind, neigh_list)
      if (is.null(hop_df) || nrow(hop_df) == 0L) break
      hop_df       <- hop_df[!duplicated(hop_df$id), ]
      frontier_ids <- hop_df$id
      seen_ids     <- c(seen_ids, hop_df$id)
      all_selected <- rbind(all_selected, hop_df)
    }
  }

  sel_summary <- summarise_sensor_selection(
    selected_data = all_selected,
    params        = list(
      r                    = r,
      k                    = k,
      hops                 = hops,
      weighting_vars       = weighting_vars,
      weighting_vars_alpha = weighting_vars_alpha,
      neighbourhood_decay  = neighbourhood_decay,
      neighbour_hops       = neighbour_hops,
      bundle_scoring       = bundle_scoring,
      filtering            = filtering
    )
  )

  # Build mi_scores: named numeric vector ordered by selection_order.
  # One entry per physical sensor (base ID, i.e. without -WITH/-AGAINST suffix).
  # Derived rows (is_derived == TRUE) are excluded; only placed sensors are included.
  sel_true_final  <- all_selected[
    all_selected$selected %in% TRUE & !all_selected$is_derived %in% TRUE, ]
  sel_base_final  <- sub("-WITH$|-AGAINST$", "", sel_true_final$id)
  base_score_lk   <- tapply(sel_true_final$mi_score,       sel_base_final, max, na.rm = TRUE)
  base_order_lk   <- tapply(sel_true_final$selection_order, sel_base_final, min, na.rm = TRUE)
  sorted_bases    <- names(sort(base_order_lk))
  mi_scores_out   <- unname(base_score_lk[sorted_bases])
  names(mi_scores_out) <- sorted_bases

  list(
    selected_data_entries = all_selected,
    n_counties            = n_counties,
    k_per_county          = k,
    r                     = r,
    summary               = sel_summary,
    mi_scores             = mi_scores_out
  )
}

# ════════════════════════════════════════════════════════════════════════════════
# greedy_mi_sensor_selection_derivable() ----
# ════════════════════════════════════════════════════════════════════════════════

#' Greedy MI sensor selection with flow-conservation derivable-link bonus
#'
#' Operates on the **full directed-link covariance matrix** (n×n, no with+against
#' transformation).  For each candidate sensor placement (base ID), the score is
#' the sum of sequential MI gains from all links that become known:
#'   1. The -WITH direction (if unmeasured).
#'   2. The -AGAINST direction (if it exists and is unmeasured).
#'   3. Any links unlocked via flow conservation (\code{enables_derivable_links}).
#'
#' Score formula:   score(b) = bundle_MI(b) × effective_w(b)
#'   bundle_MI(b)  = Σ_j  (cond_var1[j|prev_j] / cond_var2[j|prev_j])
#'                         computed sequentially inside the bundle
#'   effective_w(b)= f0·mean_w(bundle) + f1·mean_w(hop-1 \ bundle)
#'                         + f2·mean_w(hop-2 \ (bundle ∪ hop-1))
#'   where f0:f1:f2 are the normalised \code{neighbourhood_decay} fractions.
#'   When \code{neighbourhood_decay} is NULL, effective_w(b) = mean_w(bundle).
#'
#' The covariance is updated link-by-link for every member of the selected bundle
#' using the same incremental Cholesky approach as \code{greedy_mi_sensor_selection()}.
#' After each full-bundle selection the live \code{enables_derivable_links} lookup
#' is pruned: any link whose flow-conservation equation is now satisfied (all
#' unknowns in the equation are now measured/derived) is removed from other
#' candidates' bundles.
#'
#' @param data              Data frame produced by \code{identify_derivable_nodes()};
#'   must contain columns \code{id}, \code{aadt}, \code{enables_derivable_links}
#'   (list-column, may be NULL per row if none), and optionally
#'   \code{functionalRoadClass} (for FRC weighting).
#' @param covariance_matrix Dense n×n covariance matrix in directed-link order
#'   matching \code{nrow(data)}.  Produced by
#'   \code{create_covariance_and_precision_matrix(with_and_against = FALSE)$covariance_matrix}.
#' @param k                 Number of sensor LOCATIONS (base IDs) to select.
#' @param adjacency_matrix  Optional sparse n×n adjacency matrix (directed links).
#'   Required when \code{neighbourhood_decay} is non-NULL or
#'   \code{append_neighbours = TRUE}.
#' @param weighting_vars       Named list of weighting specifications; forwarded to
#'   \code{compute_weights()}.  Each entry is either a named numeric vector
#'   (level → weight lookup), \code{"proportional"}, or \code{"log_proportional"}.
#'   \code{NULL} = uniform weights (default).
#' @param weighting_vars_alpha Optional numeric vector of exponents (one per entry in
#'   \code{weighting_vars}); normalised to sum 1.  \code{NULL} = equal contribution.
#' @param neighbourhood_decay  Numeric vector of any length \code{c(f0, f1, ..., fH)},
#'   or NULL (no neighbourhood smoothing).  See
#'   \code{greedy_mi_sensor_selection_norway()} for full details.
#' @param neighbourhood_scoring One of \code{"mean"} (default, per-degree mean),
#'   \code{"sum"} (raw sum), \code{"proportional"} (sum / county-mean, no
#'   compression), or \code{"log_proportional"} (log-compressed county-normalised
#'   sum; see \code{compute_neighbourhood_weight()}).
#' @param append_neighbours    Logical; append 1-hop neighbours to output for
#'   visualisation (default TRUE).
#'
#' @return A list mirroring \code{greedy_mi_sensor_selection()} with elements:
#'   \code{selected_base_ids}, \code{selected_data_entries} (data frame with
#'   extra columns \code{selected}, \code{mi_score}, \code{bundle_size},
#'   \code{n_derivable_gained}, \code{derived_links} list-column,
#'   \code{selection_order}), \code{mi_scores}, \code{bundle_details},
#'   \code{measured_idx}, \code{unmeasured_idx_initial}.
greedy_mi_sensor_selection_derivable <- function(
    data,
    covariance_matrix,
    k,
    adjacency_matrix      = NULL,
    weighting_vars        = NULL,
    weighting_vars_alpha  = NULL,
    neighbourhood_decay   = NULL,
    neighbourhood_scoring = "mean",
    append_neighbours     = TRUE,
    eligible_ids          = NULL) {

  n <- nrow(data)
  stopifnot(k >= 1L, nrow(covariance_matrix) == n, ncol(covariance_matrix) == n)
  if (!"enables_derivable_links" %in% names(data))
    stop("'data' must contain 'enables_derivable_links' (from identify_derivable_nodes()).")

  # ── 1. Measured / unmeasured ──────────────────────────────────────────────
  is_unmeasured <- is.na(data$aadt)

  # ── 2. Candidate base IDs ────────────────────────────────────────────────
  base_ids_all    <- sub("-WITH$|-AGAINST$", "", data$id)
  uniq_bases      <- unique(base_ids_all)
  with_row_all    <- match(paste0(uniq_bases, "-WITH"),    data$id)
  against_row_all <- match(paste0(uniq_bases, "-AGAINST"), data$id)

  has_with        <- !is.na(with_row_all)
  uniq_bases      <- uniq_bases[has_with]
  with_row_all    <- with_row_all[has_with]
  against_row_all <- against_row_all[has_with]

  with_unmeas    <- is_unmeasured[with_row_all]
  against_unmeas <- !is.na(against_row_all) & is_unmeasured[against_row_all]
  is_candidate   <- with_unmeas | against_unmeas

  # ── Restrict candidates to eligible_ids (filtering) ──────────────────────
  if (!is.null(eligible_ids)) {
    is_candidate <- is_candidate & (uniq_bases %in% eligible_ids)
    if (!any(is_candidate))
      stop("No unmeasured eligible locations found after applying filter.")
  }

  # ── 3. Live derivable lookup (pruned each iteration) ──────────────────────
  enables_list <- vector("list", n)
  for (i in seq_len(n)) {
    ids_enabled <- data$enables_derivable_links[[i]]
    if (length(ids_enabled) == 0L || is.null(ids_enabled)) next
    row_idx <- match(ids_enabled, data$id)
    enables_list[[i]] <- row_idx[!is.na(row_idx)]
  }

  # ── 4. Weighting ──────────────────────────────────────────────────────────
  w <- compute_weights(data, weighting_vars, weighting_vars_alpha, row_lookup = NULL)

  # ── 2.5. Hybrid collapse: merge simple WITH+AGAINST pairs into one row ─────
  # A "simple" base: both directions exist, both unmeasured, neither direction
  # enables any derivable link, and neither is a derivable target of another
  # sensor.  Collapsing these pairs reduces the Schur complement matrix from
  # ~|U_n|² to ~|U_mixed|² (~0.6·|U_n|), giving ~3× speedup.
  all_derivable_rows_n <- unique(unlist(enables_list))

    is_simple_base <- logical(length(uniq_bases))
    for (b_idx in seq_along(uniq_bases)) {
      w_row  <- with_row_all[b_idx]
      ag_row <- against_row_all[b_idx]
      is_simple_base[b_idx] <-
        !is.na(ag_row)                                &&
        is_unmeasured[w_row] && is_unmeasured[ag_row] &&
        (is.null(enables_list[[w_row]])  || length(enables_list[[w_row]])  == 0L) &&
        (is.null(enables_list[[ag_row]]) || length(enables_list[[ag_row]]) == 0L) &&
        !w_row  %in% all_derivable_rows_n             &&
        !ag_row %in% all_derivable_rows_n
    }

    # Build n_to_mixed (length n) and mixed_to_n (list, length nw).
    # Simple pairs share one mixed index; all other links get individual indices.
    n_to_mixed <- integer(n)
    mixed_to_n <- list()
    m_counter  <- 0L
    for (b_idx in seq_along(uniq_bases)) {
      if (!is_simple_base[b_idx]) next
      w_row  <- with_row_all[b_idx]
      ag_row <- against_row_all[b_idx]
      m_counter <- m_counter + 1L
      n_to_mixed[w_row]       <- m_counter
      n_to_mixed[ag_row]      <- m_counter
      mixed_to_n[[m_counter]] <- c(w_row, ag_row)
    }
    for (j in seq_len(n)) {
      if (n_to_mixed[j] != 0L) next
      m_counter <- m_counter + 1L
      n_to_mixed[j]           <- m_counter
      mixed_to_n[[m_counter]] <- j
    }
    nw <- m_counter

    C_hybrid <- Matrix::sparseMatrix(
      i = n_to_mixed, j = seq_len(n), x = 1.0, dims = c(nw, n)
    )
    Sigma_mixed <- as.matrix(C_hybrid %*% covariance_matrix %*% t(C_hybrid))

    # Mean weight per mixed row (equal to w[w_row] for simple pairs since both
    # directions of a road share functional road class and AADT attributes).
    w_mixed <- as.numeric((C_hybrid %*% w) / Matrix::rowSums(C_hybrid))

    # A mixed row is unmeasured if any of its constituent n-links is unmeasured.
    is_unmeasured_mixed <- logical(nw)
    for (m in seq_len(nw))
      is_unmeasured_mixed[m] <- any(is_unmeasured[mixed_to_n[[m]]])

    # Map enables_list to mixed space.  Derivable links are individual directed
    # links (never simple pairs), so each maps to a unique mixed index.
    enables_list_mixed <- vector("list", nw)
    for (j in seq_len(n)) {
      en_n <- enables_list[[j]]
      if (is.null(en_n) || length(en_n) == 0L) next
      m <- n_to_mixed[j]
      enables_list_mixed[[m]] <- unique(c(enables_list_mixed[[m]], n_to_mixed[en_n]))
    }

    covariance_matrix <- Sigma_mixed
    w                 <- w_mixed
    is_unmeasured     <- is_unmeasured_mixed
    enables_list      <- enables_list_mixed

    message("  hybrid_collapse: n=", n, " \u2192 nw=", nw, " (",
            sum(is_simple_base), " simple pairs collapsed)")

  # ── 5. Warm-start: condition on existing measurements (S0) ────────────────
  # After hybrid collapse, covariance_matrix is nw×nw and is_unmeasured is
  # length nw; measured_all / unmeasured_all are mixed-space indices.
  measured_all <- which(!is_unmeasured)
  k <- min(k, sum(is_candidate))

  cond_var1 <- diag(covariance_matrix)
  if (length(measured_all) > 0L) {
    Sigma_mm  <- covariance_matrix[measured_all, measured_all, drop = FALSE]
    Sigma_om  <- covariance_matrix[, measured_all, drop = FALSE]
    L_mm      <- chol(Sigma_mm)
    A_warm    <- forwardsolve(t(L_mm), t(Sigma_om))
    cond_var1 <- pmax(cond_var1 - colSums(A_warm^2), 0)
  } else {
    A_warm <- matrix(0, 0, nw)
  }
  cond_var1[measured_all] <- 0

  unmeasured_all <- which(is_unmeasured)
  Sigma_uu <- covariance_matrix[unmeasured_all, unmeasured_all, drop = FALSE]
  if (length(measured_all) > 0L) {
    A_uu     <- A_warm[, unmeasured_all, drop = FALSE]
    Sigma_uu <- Sigma_uu - crossprod(A_uu)
  }
  prec_uu <- chol2inv(chol(Sigma_uu))

  cond_var2 <- rep(Inf, nw)
  cond_var2[unmeasured_all] <- 1 / diag(prec_uu)

  # Incremental Cholesky columns; grows dynamically if needed.
  L2 <- matrix(0, nw, k * 4L)
  n_added <- 0L

  # ── 6. Score one bundle via b×b Cholesky (read-only, no side effects) ────
  # Key identity: diag(chol(S_bb))²[j] equals the conditional variance of
  # x_{bundle[j]} given x_{bundle[1..j-1]}, S0, and all committed links —
  # i.e. exactly local_cv1[bundle[j]] at each step of the old sequential loop.
  # Cost: O(b² × n_added) per call instead of O(b × n × n_added), giving
  # roughly a 1000× reduction at typical county partition sizes.
  score_bundle_fast <- function(bundle_rows) {
    b <- length(bundle_rows)
    if (b == 0L) return(0)
    # Residual covariance submatrix of the bundle given S0 + committed links.
    S_bb <- covariance_matrix[bundle_rows, bundle_rows, drop = FALSE]
    if (nrow(A_warm) > 0L)
      S_bb <- S_bb - crossprod(A_warm[, bundle_rows, drop = FALSE])
    if (n_added > 0L)
      S_bb <- S_bb - tcrossprod(L2[bundle_rows, seq_len(n_added), drop = FALSE])
    R <- tryCatch(chol(S_bb), error = function(e) NULL)
    if (is.null(R)) return(0)
    sum(diag(R)^2 / cond_var2[bundle_rows])
  }

  # ── 7. Commit one link to global state ────────────────────────────────────
  commit_link <- function(j_row) {
    col_vec <- covariance_matrix[, j_row] - crossprod(A_warm, A_warm[, j_row])
    if (n_added > 0L)
      col_vec <- col_vec -
        L2[, seq_len(n_added), drop = FALSE] %*% L2[j_row, seq_len(n_added)]
    pivot   <- pmax(col_vec[j_row], .Machine$double.eps * 100)
    n_added <<- n_added + 1L
    # Grow L2 if we've exhausted pre-allocated columns.
    if (n_added > ncol(L2))
      L2 <<- cbind(L2, matrix(0, nw, ncol(L2)))
    L2[, n_added]    <<- col_vec / sqrt(pivot)
    cond_var1        <<- pmax(cond_var1 - L2[, n_added]^2, 0)
    cond_var1[j_row] <<- 0

    j_pos <- which(unmeasured_candidates == j_row)
    if (length(j_pos) == 1L) {
      if (length(unmeasured_candidates) > 1L) {
        prec_uu <<- prec_uu[-j_pos, -j_pos, drop = FALSE] -
          tcrossprod(prec_uu[-j_pos, j_pos]) / prec_uu[j_pos, j_pos]
      } else {
        prec_uu <<- matrix(numeric(0), 0, 0)
      }
      unmeasured_candidates <<- unmeasured_candidates[-j_pos]
      if (length(unmeasured_candidates) > 0L)
        cond_var2[unmeasured_candidates] <<- 1 / diag(prec_uu)
    }
    cond_var2[j_row] <<- Inf
  }

  # ── 8. Base-ID adjacency (for neighbourhood weighting) ─────────────────────
  A_base <- NULL
  if (!is.null(neighbourhood_decay) && !is.null(adjacency_matrix)) {
    base_of_link <- match(base_ids_all, uniq_bases)
    m_bases   <- length(uniq_bases)
    adj_trip  <- Matrix::summary(as(adjacency_matrix, "dgCMatrix"))
    base_i    <- base_of_link[adj_trip$i]
    base_j    <- base_of_link[adj_trip$j]
    ok        <- !is.na(base_i) & !is.na(base_j) & base_i != base_j
    A_base    <- sign(Matrix::sparseMatrix(
      i = base_i[ok], j = base_j[ok], x = 1,
      dims = c(m_bases, m_bases)
    ))
  }

  # ── 9. Greedy loop with lazy re-scoring ───────────────────────────────────
  # cached_scores[b_idx] is a valid upper bound on the true score by submodularity:
  # committing any bundle can only decrease marginal MI gains for all remaining
  # candidates. At each iteration we sort by cached score descending and evaluate
  # lazily, stopping as soon as the current candidate's fresh score is >= the
  # cached score of every remaining unchecked candidate.
  # Precompute mixed-space row indices for WITH and AGAINST per base.
  # In the identity path (hybrid_collapse = FALSE), these equal with_row_all
  # and against_row_all respectively, so behaviour is unchanged.
  mixed_w_all   <- n_to_mixed[with_row_all]
  mixed_ag_safe <- ifelse(is.na(against_row_all), 1L, n_to_mixed[against_row_all])

  unmeasured_candidates <- unmeasured_all    # shrinks as links are committed
  selected_base_ids  <- character(k)
  mi_scores_out      <- numeric(k)
  bundle_details_out <- vector("list", k)
  cached_scores      <- rep(Inf, length(uniq_bases))  # Inf = not yet evaluated
  committed_bases    <- integer(0)  # base indices committed in previous iterations
  measured_bases     <- which(!is_unmeasured[mixed_w_all])  # bases with a measured WITH link in S0

  # ── Neighbourhood baselines for proportional/log_proportional scoring ─────
  # Uses residual A_base (measured_bases columns zeroed) so baselines reflect
  # the typical connectivity seen by unmeasured candidates at greedy step 1.
  neigh_baselines <- NULL; lp_scales <- NULL
  if (!is.null(A_base) && neighbourhood_scoring %in% c("log_proportional", "proportional")) {
    max_hops_b   <- length(neighbourhood_decay) - 1L
    w_for_bases  <- w[mixed_w_all]
    A_base_res_b <- A_base
    if (length(measured_bases) > 0L) A_base_res_b[, measured_bases] <- 0
    A_power_b <- A_base_res_b
    A_cumul_b <- sign(A_base_res_b); Matrix::diag(A_cumul_b) <- 0
    s_init <- vector("list", max_hops_b)
    s_init[[1L]] <- as.numeric(A_cumul_b %*% w_for_bases)
    for (h in seq_len(max_hops_b - 1L)) {
      A_power_b <- A_power_b %*% A_base_res_b
      A_h_b     <- sign(A_power_b); Matrix::diag(A_h_b) <- 0
      A_exact_b <- A_h_b; A_exact_b[A_cumul_b > 0] <- 0
      s_init[[h + 1L]] <- as.numeric(A_exact_b %*% w_for_bases)
      A_cumul_b <- sign(A_cumul_b + A_exact_b)
    }
    neigh_baselines <- sapply(s_init, mean)
    if (neighbourhood_scoring == "log_proportional") {
      lp_scales <- mapply(function(s, mu)
        max(log(max(s) / pmax(mu, 1e-8) + 1), 1e-8),
        s_init, neigh_baselines)
    }
  }

  for (iter in seq_len(k)) {

    best_score  <- -Inf
    best_base   <- NA_character_
    best_bundle <- integer(0)

    # Rebuild active candidate set for this iteration.
    active_with_unmeas    <- is_unmeasured[mixed_w_all] &
                              !is.nan(cond_var1[mixed_w_all])
    active_against_unmeas <- !is.na(against_row_all) & !is_simple_base &
                              is_unmeasured[mixed_ag_safe] &
                              !is.nan(cond_var1[mixed_ag_safe])
    active_candidate      <- (active_with_unmeas | active_against_unmeas) &
                              is_candidate
    active_idx            <- which(active_candidate)

    if (length(active_idx) == 0L) {
      message("No valid candidate found at iteration ", iter, ". Stopping early.")
      k <- iter - 1L
      break
    }

    # Sort by cached score descending; cached scores are upper bounds so we can
    # stop evaluating the moment the current candidate's fresh score beats all
    # remaining cached values.
    sort_ord          <- order(-cached_scores[active_idx])
    active_idx_sorted <- active_idx[sort_ord]

    for (b_idx in active_idx_sorted) {

      # Lazy upper-bound check: cached score cannot beat current best → done.
      if (cached_scores[b_idx] <= best_score) break

      base     <- uniq_bases[b_idx]
      mixed_w  <- mixed_w_all[b_idx]
      mixed_ag <- mixed_ag_safe[b_idx]

      # Build bundle in mixed space.
      bundle <- integer(0)
      if (is_simple_base[b_idx]) {
        # Collapsed pair: one mixed row represents both directions.
        if (is_unmeasured[mixed_w] && cond_var1[mixed_w] > 0)
          bundle <- mixed_w
      } else {
        if (!is.na(with_row_all[b_idx]) && is_unmeasured[mixed_w] && cond_var1[mixed_w] > 0)
          bundle <- c(bundle, mixed_w)
        if (!is.na(against_row_all[b_idx]) && is_unmeasured[mixed_ag] && cond_var1[mixed_ag] > 0)
          bundle <- c(bundle, mixed_ag)
      }

      # Derivable rows (enables_list is now in mixed space).
      # Simple pairs have empty enables_list by construction; no deriv_rows.
      deriv_rows    <- integer(0)
      src_mixed_rows <- c(
        mixed_w,
        if (!is_simple_base[b_idx] && !is.na(against_row_all[b_idx])) mixed_ag
      )
      for (src_m in src_mixed_rows) {
        en <- enables_list[[src_m]]
        if (is.null(en)) next
        for (dr in en) {
          if (is_unmeasured[dr] && cond_var1[dr] > 0 && !dr %in% bundle)
            deriv_rows <- c(deriv_rows, dr)
        }
      }
      bundle <- c(bundle, deriv_rows)

      if (length(bundle) == 0L) { cached_scores[b_idx] <- 0; next }

      # b×b Cholesky score (fast path).
      mi_b <- score_bundle_fast(bundle)

      w_bundle <- w[bundle]
      mean_w   <- mean(w_bundle)

      if (!is.null(A_base) && !is.null(neighbourhood_decay)) {
        f        <- neighbourhood_decay / sum(neighbourhood_decay)
        max_hops <- length(f) - 1L
        # Exclude current bundle's own bases, committed bases, and initially
        # measured bases. Known bases cannot serve as bridges in the residual
        # graph, enforcing the Markov conditional independence property.
        # Drop NA: derivable links may be -AGAINST-only rows with no entry in
        # uniq_bases, so base_of_link returns NA for them; those links have no
        # corresponding base-adjacency column and must be excluded.
        bundle_bases <- unique(na.omit(base_of_link[unlist(mixed_to_n[bundle])]))
        known_bases  <- unique(c(bundle_bases, committed_bases, measured_bases))
        # Iteratively build exact h-hop shells in the residual graph.
        A_row_res <- A_base[b_idx, , drop = FALSE]
        if (length(known_bases) > 0L) A_row_res[, known_bases] <- 0
        hop_bases_list <- vector("list", max_hops)
        visited        <- known_bases
        A_reach_row    <- A_row_res
        for (h in seq_len(max_hops)) {
          if (h > 1L) {
            A_reach_row <- sign(A_reach_row %*% A_base)
            if (length(known_bases) > 0L) A_reach_row[, known_bases] <- 0
          }
          hop_h               <- setdiff(which(A_reach_row > 0), visited)
          hop_bases_list[[h]] <- hop_h
          visited             <- c(visited, hop_h)
        }
        s_hops <- switch(neighbourhood_scoring,
          mean = sapply(hop_bases_list, function(hb)
            if (length(hb) > 0L) mean(w[mixed_w_all[hb]], na.rm = TRUE) else 0),
          sum  = sapply(hop_bases_list, function(hb)
            if (length(hb) > 0L) sum(w[mixed_w_all[hb]], na.rm = TRUE) else 0),
          proportional = sapply(seq_len(max_hops), function(h) {
            raw <- if (length(hop_bases_list[[h]]) > 0L)
                     sum(w[mixed_w_all[hop_bases_list[[h]]]], na.rm = TRUE) else 0
            nb  <- if (!is.null(neigh_baselines) && h <= length(neigh_baselines))
                     neigh_baselines[[h]] else 1e-8
            raw / pmax(nb, 1e-8)
          }),
          log_proportional = sapply(seq_len(max_hops), function(h) {
            raw <- if (length(hop_bases_list[[h]]) > 0L)
                     sum(w[mixed_w_all[hop_bases_list[[h]]]], na.rm = TRUE) else 0
            nb  <- if (!is.null(neigh_baselines) && h <= length(neigh_baselines))
                     neigh_baselines[[h]] else 1e-8
            sc  <- if (!is.null(lp_scales) && h <= length(lp_scales))
                     lp_scales[[h]] else 1
            log(raw / pmax(nb, 1e-8) + 1) / pmax(sc, 1e-8)
          })
        )
        # Note: neighbourhood scoring is not guaranteed monotone decreasing as
        # committed_bases grows, so cached_scores is a heuristic upper bound.
        eff_w <- max(f[1L] * mean_w + sum(f[-1L] * s_hops), 1e-8)
      } else {
        eff_w <- mean_w
      }

      score_b              <- mi_b * eff_w
      cached_scores[b_idx] <- score_b   # update; by submodularity score only decreases

      if (score_b > best_score) {
        best_score  <- score_b
        best_base   <- base
        best_bundle <- bundle
        best_deriv  <- deriv_rows
        best_b_idx  <- b_idx
      }
    }

    if (is.na(best_base)) {
      message("No valid candidate found at iteration ", iter, ". Stopping early.")
      k <- iter - 1L
      break
    }

    # Commit the entire best bundle link-by-link.
    for (j_row in best_bundle) {
      commit_link(j_row)
      is_unmeasured[j_row] <- FALSE
    }

    # Prune enables_list: remove committed links from all other enables lists.
    for (j_row in best_bundle) {
      for (src in seq_len(nw)) {
        en <- enables_list[[src]]
        if (!is.null(en))
          enables_list[[src]] <- en[en != j_row]
      }
    }

    selected_base_ids[iter] <- best_base
    mi_scores_out[iter]     <- best_score
    committed_bases         <- c(committed_bases, best_b_idx)
    # Convert mixed-space indices back to n-space for output.
    # In the identity path (hybrid_collapse=FALSE) this is a no-op.
    best_bundle_n <- unlist(mixed_to_n[best_bundle])
    best_deriv_n  <- unlist(mixed_to_n[best_deriv])
    bundle_details_out[[iter]] <- list(
      base_id       = best_base,
      bundle_rows   = best_bundle_n,
      deriv_rows    = best_deriv_n,
      bundle_ids    = data$id[best_bundle_n],
      deriv_ids     = data$id[best_deriv_n],
      mi_score      = best_score,
      bundle_size   = length(best_bundle_n),
      n_derivable   = length(best_deriv_n)
    )
  }

  # Trim to actual number selected (may be < k if early stop).
  selected_base_ids  <- selected_base_ids[seq_len(k)]
  mi_scores_out      <- mi_scores_out[seq_len(k)]
  bundle_details_out <- bundle_details_out[seq_len(k)]

  # ── 10. Build selected_data_entries ───────────────────────────────────────
  all_sel_ids <- unlist(lapply(bundle_details_out, `[[`, "bundle_ids"))
  sel_df <- data[data$id %in% all_sel_ids, ]
  sel_df$selected          <- TRUE
  sel_df$mi_score          <- NA_real_
  sel_df$bundle_size       <- NA_integer_
  sel_df$n_derivable_gained <- NA_integer_
  sel_df$selection_order   <- NA_integer_
  sel_df$derived_links     <- vector("list", nrow(sel_df))
  sel_df$is_derived        <- FALSE

  for (iter in seq_len(k)) {
    bd  <- bundle_details_out[[iter]]
    idx <- which(sel_df$id %in% bd$bundle_ids)
    sel_df$mi_score[idx]           <- bd$mi_score
    sel_df$bundle_size[idx]        <- bd$bundle_size
    sel_df$n_derivable_gained[idx] <- bd$n_derivable
    sel_df$selection_order[idx]    <- iter
    for (ii in idx)
      sel_df$derived_links[[ii]]   <- bd$deriv_ids
    sel_df$is_derived[sel_df$id %in% bd$deriv_ids] <- TRUE
  }

  # Optionally append 1-hop neighbours.
  if (!is.null(adjacency_matrix) && isTRUE(append_neighbours)) {
    sel_row_idx    <- which(data$id %in% all_sel_ids)
    adj_sub        <- adjacency_matrix[sel_row_idx, , drop = FALSE]
    neigh_data_idx <- setdiff(which(Matrix::colSums(adj_sub != 0) > 0),
                              which(data$id %in% sel_df$id))
    if (length(neigh_data_idx) > 0L) {
      neigh_rows                   <- data[neigh_data_idx, ]
      neigh_rows$selected          <- FALSE
      neigh_rows$mi_score          <- NA_real_
      neigh_rows$bundle_size       <- NA_integer_
      neigh_rows$n_derivable_gained <- NA_integer_
      neigh_rows$selection_order   <- NA_integer_
      neigh_rows$derived_links     <- vector("list", nrow(neigh_rows))
      neigh_rows$is_derived        <- FALSE
      sel_df <- rbind(sel_df, neigh_rows)
    }
  }

  list(
    selected_base_ids      = selected_base_ids,
    selected_data_entries  = sel_df,
    mi_scores              = mi_scores_out,
    bundle_details         = bundle_details_out,
    measured_idx           = measured_all,
    unmeasured_idx_initial = unmeasured_all
  )
}
