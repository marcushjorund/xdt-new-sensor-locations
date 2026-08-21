# Translation to R of Python implementations for greedy GP sensor placement algorithms.
#
# Reference:
#   Krause, Singh, Guestrin (2008). "Near-Optimal Sensor Placements in Gaussian
#   Processes: Theory, Efficient Algorithms and Empirical Studies."
#   Journal of Machine Learning Research, 9, 235-284.
#   https://www.jmlr.org/papers/volume9/krause08a/krause08a.pdf
#
# Three greedy criteria are implemented, each in multiple variants of increasing
# computational efficiency:
#   1. Maximum entropy (differential entropy of the selected set)       -- Algorithm 1
#   2. Maximum mutual information (MI between selected and unselected)  -- Improvement to Alg. 1
#   3. Minimum Nystrom KL divergence                                    -- greedy column sparsity
#
# IMPORTANT - index convention:
#   All functions return 1-based integer indices (standard R convention).
#   The equivalent Python code returns 0-based indices.
#
# KERNEL INTERFACE:
#   Functions expect two kernel arguments instead of a single sklearn Kernel object:
#     kernel(X1, X2) -- returns an nrow(X1) x nrow(X2) covariance matrix
#     kernel_diag(X) -- returns a length-nrow(X) vector of marginal variances (diagonal of kernel(X, X))
#   Any kernel that satisfies this contract can be used.
#
# CHOLESKY CONVENTION:
#   R's chol() returns an UPPER triangular U such that  t(U) %*% U = Sigma.
#   NumPy's np.linalg.cholesky() returns a LOWER triangular L such that L %*% t(L) = Sigma.
#   Where the code manually builds or updates a Cholesky factor it uses LOWER triangular
#   matrices (matching the Python convention) via forwardsolve / backsolve with transpose.

# =============================================================================
# Section 1 -- Helper / utility functions
# =============================================================================

# Invert a symmetric positive definite (SPD) matrix m.
# chol2inv(chol(m)) is faster than solve(m) for SPD matrices: chol() costs O(n^3/3)
# while a general LU inversion costs O(n^3). chol2inv then completes the inverse
# from the Cholesky factor in O(n^3/3) additional work, roughly halving the total cost.
mat_inv <- function(m) {
  chol2inv(chol(m))
}

# Solve the linear system A x = b for SPD matrix A.
# Uses the Cholesky factorisation A = t(U) U to solve via two triangular systems:
#   t(U) y = b  (forward substitution)  then  U x = y  (back substitution).
# This avoids forming the full inverse A^{-1}, which would require an additional
# O(n^3) multiply step and accumulates more floating-point error.
solve_spd <- function(A, b) {
  U <- chol(A)
  backsolve(U, forwardsolve(t(U), b))
}

# Solve the triangular system L x = b.
#   lower = TRUE  (default) : L is lower triangular -> use forwardsolve
#   lower = FALSE           : L is upper triangular -> use backsolve
# The 'transpose' argument allows solving L^T x = b when set to TRUE.
solve_tri <- function(L, b, lower = TRUE, transpose = FALSE) {
  if (lower) {
    forwardsolve(L, b, transpose = transpose)
  } else {
    backsolve(L, b, transpose = transpose)
  }
}

# Compute the quadratic form  x^T Theta^{-1} x  column-wise.
#
# Arguments:
#   theta   -- p x p SPD matrix (or, when factor = FALSE, its upper Cholesky factor U)
#   x       -- p x q matrix; each column is a vector to form x_i^T Theta^{-1} x_i
#   factor  -- if TRUE (default), Theta is given and chol() is called internally;
#              if FALSE, Theta is already an upper Cholesky factor U (t(U) %*% U = Theta)
#
# Returns a length-q vector of scalar quadratic forms.
#
# Algorithm: write Theta = t(U) U  (U upper triangular from chol()).
#   x^T Theta^{-1} x = x^T U^{-1} (t(U))^{-1} x = || (t(U))^{-1} x ||^2
#   Solve  t(U) temp = x  by forward substitution, then return colSums(temp^2).
quad_form <- function(theta, x, factor = TRUE) {
  if (nrow(x) == 0) {
    # empty selection set: conditional == marginal variance
    return(numeric(ncol(x)))
  }
  U <- if (factor) chol(theta) else theta
  # backsolve with transpose=TRUE solves t(U) temp = x  (forward substitution)
  temp <- backsolve(U, x, transpose = TRUE)
  colSums(temp * temp)
}

# Rank-one Cholesky DOWNDATE in-place: L -> chol_lower(L L^T - u u^T).
#
# Given the lower Cholesky factor L of some matrix A (A = L L^T), and a vector u,
# this function returns the lower Cholesky factor of  A - u u^T.
#
# Arguments:
#   L -- n x n lower triangular Cholesky factor
#   u -- length-n vector to remove
#   j -- number of columns to process (defaults to ncol(L)); set to a smaller
#        integer to process only the first j columns (used in mi_chol for submatrices)
#
# Algorithm: Givens rotations applied column by column.
#   At each column i the rotation parameters are:
#     dp = sqrt(L[i,i]^2 - u[i]^2)     (new diagonal)
#     c1 = L[i,i] / dp ,  c2 = u[i] / dp
#   The column is updated:   L[:,i] <- c1 * L[:,i] - c2 * u
#   The carry vector is updated: u  <- u / c1 - (c2 / c1) * L[:,i]_new
#
# Note: R copies on modify, so L and u are NOT modified in the caller's scope.
chol_downdate <- function(L, u, j = NULL) {
  n <- if (is.null(j)) ncol(L) else j
  for (i in seq_len(n)) {
    c1 <- L[i, i]
    c2 <- u[i]
    # new diagonal entry after the downdate.
    # max(0, ...) prevents NaN from floating-point rounding when c1 ≈ c2.
    dp <- sqrt(max(0, c1 * c1 - c2 * c2))
    # normalised Givens rotation coefficients
    c1 <- c1 / dp
    c2 <- c2 / dp
    # update i-th column of L (must come before updating u, as u update uses new L[:,i])
    L[, i] <- c1 * L[, i] - c2 * u
    # update carry vector using the already-updated L[:,i]
    u <- u / c1 - (c2 / c1) * L[, i]
  }
  L
}

# =============================================================================
# Section 2 -- Entropy-based greedy sensor placement (four variants)
#
# All four functions implement the same greedy criterion:
#   At each step, add the point x_k that maximises the conditional variance
#   Var(x_k | x_{S}) where S is the already-selected set.
#
# This is equivalent to greedily maximising the differential entropy of the
# selected set (Krause et al., Algorithm 1).
#
# Computational complexities:
#   entropy_naive   : O(s*(s^3 + n*s^2)) = O(n s^3 + s^4)
#   entropy_prec    : O(s*(n*s + s^2))   = O(n s^2)
#   entropy_prechol : O(s*(n*s + s^2))   = O(n s^2)
#   entropy_chol    : O(s*(n*s + s^2))   = O(n s^2)
# =============================================================================

# Greedy max-entropy selection -- naive O(n s^3 + s^4) implementation.
#
# At each iteration the full submatrix theta[S, S] is passed to quad_form,
# which re-factors it from scratch. This is correct but wasteful compared to
# the incremental Cholesky variants below.
#
# Arguments:
#   X           -- n x d matrix of candidate locations (one location per row)
#   kernel      -- function(X1, X2) returning the nrow(X1) x nrow(X2) covariance matrix
#   kernel_diag -- function(X) returning marginal variances (diagonal of kernel(X, X))
#   s           -- number of sensors to select
#
# Returns:
#   Integer vector of length min(s, n) with the selected 1-based row indices of X.
entropy_naive <- function(X, kernel, kernel_diag, s) {
  n <- nrow(X)
  s <- min(s, n)

  indexes <- integer(s)        # selected indices (1-based), filled iteratively
  theta   <- matrix(0, s, s)   # s x s submatrix of K for the selected set
  cov     <- matrix(0, s, n)   # cov[i, ] = K(X, X[indexes[i], ]) -- cross-covariances
  var     <- kernel_diag(X)    # marginal variances (updated by setting selected to -1)

  for (i in seq_len(s)) {
    # Conditional variance of every point given the current selected set S.
    # When i = 1 the set is empty and cond_var = marginal variance.
    # quad_form(theta[1:(i-1), 1:(i-1)], cov[1:(i-1), ]) computes:
    #   col_k of  cov[S, ] * K[S,S]^{-1} * cov[S, ]^T  -- the reduction in variance.
    if (i == 1) {
      cond_var <- var
    } else {
      theta_sel <- theta[1:(i - 1), 1:(i - 1), drop = FALSE]
      cov_sel   <- cov[1:(i - 1), , drop = FALSE]
      cond_var  <- var - quad_form(theta_sel, cov_sel)
    }

    # Pick the point with the largest conditional variance.
    k <- which.max(cond_var)
    indexes[i] <- k
    var[k] <- -1  # mark as selected so it won't be picked again

    # Store cross-covariance between X[k,] and every other point.
    cov[i, ] <- as.vector(kernel(X, X[k, , drop = FALSE]))

    # Update theta: fill row and column i with K values at the selected indices.
    theta[1:i, i] <- cov[i, indexes[1:i]]
    theta[i, 1:i] <- cov[i, indexes[1:i]]
  }

  indexes
}

# Greedy max-entropy selection -- precision-matrix O(n s^2) implementation.
#
# Instead of re-factoring K[S,S] every iteration, we maintain the PRECISION
# matrix P = K[S,S]^{-1} and update it with a rank-1 Schur complement formula.
# The conditional variance update uses a rank-1 column update.
#
# See entropy_naive for argument and return documentation.
entropy_prec <- function(X, kernel, kernel_diag, s) {
  n <- nrow(X)
  s <- min(s, n)

  indexes  <- integer(s)       # selected 1-based indices
  prec     <- matrix(0, s, s)  # growing precision matrix P = K[S,S]^{-1}
  cov      <- matrix(0, n, s)  # cov[, i] = K omit (X, X[indexes[i], ])
  cond_var <- kernel_diag(X)   # conditional variances (updated each step)

  for (i in seq_len(s)) {
    # Select the point with the greatest conditional variance.
    k <- which.max(cond_var)
    indexes[i] <- k

    # Store the cross-covariance vector K(:, k).
    cov[, i] <- as.vector(kernel(X, X[k, , drop = FALSE]))

    # --- Rank-1 update of the precision matrix ---
    # New precision after adding point k uses the Schur complement:
    #   [K_S  k]^{-1}  block inversion formula.
    #   [k^T  s_k]
    # Where s_k = cond_var[k] = K[k,k] - K[k,S] K[S,S]^{-1} K[S,k].
    w <- 1 / cond_var[k]  # = 1 / s_k

    if (i > 1) {
      # v = P K[S, k] -- precision applied to cross-covariance.
      # as.vector() prevents %*% returning a 1x1 matrix, which would make outer()
      # produce a 4D array instead of a 2D matrix when i-1 = 1.
      v <- as.vector(prec[1:(i - 1), 1:(i - 1), drop = FALSE] %*% cov[indexes[1:(i - 1)], i])
      prec[1:(i - 1), 1:(i - 1)] <- prec[1:(i - 1), 1:(i - 1)] + w * outer(v, v)
      prec[1:(i - 1), i] <- -w * v
      prec[i, 1:(i - 1)] <- -w * v
    } else {
      v <- numeric(0)  # empty; the update blocks are vacuous for i = 1
    }
    prec[i, i] <- w

    # --- Rank-1 update of conditional variances ---
    # The conditional covariance Cov(X, x_k | S) = K(:, k) - K(:,S) K[S,S]^{-1} K[S,k]
    # is computed as  cov(:, i) - cov(:, 1:(i-1)) * v  (where v = K[S,S]^{-1} K[S,k]).
    # Normalised by sqrt(cond_var[k]) it gives the vector by which cond_var decreases.
    cond_cov_k <- cov[, i]
    if (i > 1) {
      cond_cov_k <- cond_cov_k - as.vector(cov[, 1:(i - 1), drop = FALSE] %*% v)
    }
    cond_cov_k <- cond_cov_k / sqrt(cond_var[k])
    cond_var   <- cond_var - cond_cov_k^2
    cond_var[k] <- -1  # mark as selected
  }

  indexes
}

# Greedy max-entropy selection -- pre-Cholesky O(n s^2) implementation.
#
# Maintains a growing LOWER triangular Cholesky factor L of K[S,S] using the
# "up-looking" incremental Cholesky algorithm (fill the new row, then forward-solve).
# This avoids forming the precision matrix explicitly.
#
# See entropy_naive for argument and return documentation.
entropy_prechol <- function(X, kernel, kernel_diag, s) {
  n <- nrow(X)
  s <- min(s, n)

  indexes  <- integer(s)       # selected 1-based indices
  L        <- matrix(0, s, s)  # lower triangular Cholesky factor of K[S,S]
  cov      <- matrix(0, n, s)  # cov[, i] = K(X, X[indexes[i], ])
  cond_var <- kernel_diag(X)

  for (i in seq_len(s)) {
    k <- which.max(cond_var)
    indexes[i] <- k
    cov[, i] <- as.vector(kernel(X, X[k, , drop = FALSE]))

    # --- "Up-looking" incremental Cholesky update ---
    # Fill the new row of L: L[i, 1:(i-1)] = K[selected_prev, k].
    # Then forward-solve with the existing (i-1)×(i-1) lower block to get
    # the off-diagonal part of the new Cholesky row.
    if (i > 1) {
      Lp  <- L[1:(i - 1), 1:(i - 1), drop = FALSE]  # existing lower Cholesky block
      row <- cov[indexes[1:(i - 1)], i]               # K[S_prev, k]

      # Solve Lp * t = row  (forward substitution); result overwrites row.
      # After this: row = Lp^{-1} K[S_prev, k].
      row <- forwardsolve(Lp, row)
      L[i, 1:(i - 1)] <- row  # store solved row back into L

      # New diagonal: sqrt of conditional variance of k given S_prev.
      # cond_var[k] = K[k,k] - K[k,S_prev] K[S_prev,S_prev]^{-1} K[S_prev,k]
      #             = K[k,k] - ||row||^2
      L[i, i] <- sqrt(cov[k, i] - sum(row^2))

      # --- Conditional covariance vector update ---
      # Cov(X, x_k | S_prev) = K(:,k) - K(:,S_prev) * K[S_prev,S_prev]^{-1} * K[S_prev,k]
      # K[S_prev,S_prev]^{-1} K[S_prev,k] = (Lp Lp^T)^{-1} K[S_prev,k]
      #                                    = Lp^{-T} Lp^{-1} K[S_prev,k]
      #                                    = Lp^{-T} row   (row already = Lp^{-1} K[S_prev,k])
      # Solve Lp^T z = row  (back-substitution on lower-triangular Lp^T = upper t(Lp)).
      z          <- forwardsolve(Lp, row, transpose = TRUE)
      cond_cov_k <- cov[, i] - cov[, 1:(i - 1), drop = FALSE] %*% z
    } else {
      # First point: no conditioning, diagonal entry is just sqrt(K[k,k]).
      L[i, i]    <- sqrt(cov[k, i])
      cond_cov_k <- cov[, i]
    }

    # Normalise to obtain the "direction" that reduces conditional variance.
    cond_cov_k <- cond_cov_k / sqrt(cond_var[k])
    cond_var   <- cond_var - cond_cov_k^2
    cond_var[k] <- -1
  }

  indexes
}

# Greedy max-entropy selection -- left-looking Cholesky O(n s^2) implementation.
#
# Maintains an n x s lower-triangular factor L where column i encodes the
# Cholesky update for the i-th selected point across ALL n candidate locations.
# The "left-looking" update fills column i using all previously filled columns.
# This is the most cache-friendly variant for large n.
#
# See entropy_naive for argument and return documentation.
entropy_chol <- function(X, kernel, kernel_diag, s) {
  n <- nrow(X)
  s <- min(s, n)

  indexes  <- integer(s)       # selected 1-based indices
  L        <- matrix(0, n, s)  # n x s partial lower Cholesky factor
  cond_var <- kernel_diag(X)

  for (i in seq_len(s)) {
    k <- which.max(cond_var)
    indexes[i] <- k

    # --- Left-looking Cholesky column update ---
    # Column i of L stores the i-th column of  chol_lower(K)  for ALL n points.
    # Initialise with K(:, k) then subtract the contribution of previous columns.
    L[, i] <- as.vector(kernel(X, X[k, , drop = FALSE]))
    if (i > 1) {
      # Remove the projection onto the span of the first (i-1) columns.
      # L[k, 1:(i-1)] is the k-th row of the Cholesky factor up to column i-1.
      L[, i] <- L[, i] - L[, 1:(i - 1), drop = FALSE] %*% L[k, 1:(i - 1)]
    }
    # Normalise so that L[k, i]^2 = cond_var[k] (checked: equals sqrt of conditional var).
    L[, i] <- L[, i] / sqrt(L[k, i])

    # The squared i-th column gives the reduction in conditional variance for all points.
    cond_var   <- cond_var - L[, i]^2
    cond_var[k] <- -1  # mark selected
  }

  indexes
}

# =============================================================================
# Section 3 -- Mutual-information greedy sensor placement (three variants)
#
# Criterion: at each step, add the point x_k that maximises the mutual
# information  I(x_k ; X_V\S)  between x_k and the UNSELECTED set V\S,
# conditioned on the already-selected set S.
#
# This factorises as:
#   I(x_k ; X_{V\S} | X_S) ∝ cond_var_S(k) * cond_var_{V\(S∪{k})}(k)
#
# where:
#   cond_var_S(k)        = Var(x_k | X_S)         -- marginal cond. var. given selected
#   cond_var_{V\S}(k)    = 1 / [K^{-1}_{V\S}]_{kk} -- full conditional var. given unselected
#
# The full conditional var. of k given the rest equals the reciprocal of the k-th
# diagonal entry of the precision matrix of the unselected set (Schur complement property).
#
# Complexities:
#   mi_naive : O(s * (s^3 + n^3)) = O(n^4 s) -- recomputes inverses each step
#   mi_prec  : O(n^3 + s*n^2)     = O(n^3)   -- maintains precision, updates rank-1
#   mi_chol  : O(n^3 + s*n^2)     = O(n^3)   -- same but avoids explicit precision via Cholesky
# =============================================================================

# Greedy max-MI selection -- naive O(n^4 s) implementation.
#
# Computes cond_var_S from scratch via quad_form each iteration (O(n s^3) per step).
# Computes the precision of X_{V\S} by inverting the candidate submatrix (O(n^3) per step).
#
# See entropy_naive for argument documentation.
mi_naive <- function(X, kernel, kernel_diag, s) {
  n <- nrow(X)
  s <- min(s, n)

  indexes    <- integer(0)      # grows from empty
  candidates <- seq_len(n)      # indices not yet selected
  cov_full   <- kernel(X, X)    # full n x n kernel matrix (computed once)
  var        <- diag(cov_full)  # marginal variances

  for (iter in seq_len(s)) {
    # cond_var1[k] = Var(x_k | X_S) for k = 1..n.
    # quad_form expects x as (p x n): rows = selected points, cols = all candidates.
    theta1    <- cov_full[indexes, indexes, drop = FALSE]
    cond_var1 <- var - quad_form(theta1, cov_full[indexes, , drop = FALSE])

    # cond_var2[k] = full conditional variance of k given ALL other unselected points.
    # This equals 1 / [K_{V\S}^{-1}]_{kk} = 1 / diag(inv(K_{V\S})).
    # Equivalently: quad_form applied to the identity returns the diagonal of the inverse,
    # so cond_var2[candidates] = diag(K_{V\S}^{-1}) for the candidate subblock.
    # (This is the expensive O(n^3) step per iteration.)
    theta2    <- cov_full[candidates, candidates, drop = FALSE]
    cond_var2 <- rep(-1, n)
    # quad_form(theta2, I) = diag(theta2^{-1}), returned as length-nc vector
    cond_var2[candidates] <- quad_form(theta2, diag(length(candidates)))

    # Score = product of the two conditional variances (proportional to MI).
    k <- which.max(cond_var1 * cond_var2)
    indexes    <- c(indexes, k)
    candidates <- candidates[candidates != k]
  }

  indexes
}

# Greedy max-MI selection -- precision-matrix O(n^3) implementation.
#
# Initialises the full precision P = K^{-1} (O(n^3) once), then maintains
# the precision of the UNSELECTED set via Schur-complement row/column deletion.
# Marginalising out a variable in covariance is equivalent to conditioning in precision.
#
# See entropy_naive for argument documentation.
mi_prec <- function(X, kernel, kernel_diag, s) {
  n <- nrow(X)
  s <- min(s, n)

  indexes    <- integer(s)    # selected 1-based indices
  candidates <- seq_len(n)    # 1-based indices of unselected points
  L1         <- matrix(0, n, s)   # n x s left-looking Cholesky for cond_var1
  # Full precision matrix of all candidates (initialised to K^{-1}).
  prec    <- mat_inv(kernel(X, X))
  cond_var1 <- kernel_diag(X)
  # The full conditional variance of point i given all other points equals
  # 1 / P[i, i]  (diagonal of the precision = inverse of full-conditional variance).
  cond_var2 <- diag(prec)  # length-n; updated as candidates shrink
  for (i in seq_len(s)) {
    # Score = cond_var1 (entropy criterion) * cond_var2 (full conditional variance)
    k <- which.max(cond_var1 * cond_var2)
    indexes[i] <- k

    # Position of k within the current candidate vector (1-based within candidates).
    j <- which(candidates == k)
    candidates <- candidates[-j]

    # --- Update left-looking Cholesky for cond_var1 (same as entropy_chol) ---
    L1[, i] <- as.vector(kernel(X, X[k, , drop = FALSE]))
    if (i > 1) {
      L1[, i] <- L1[, i] - L1[, 1:(i - 1), drop = FALSE] %*% L1[k, 1:(i - 1)]
    }
    L1[, i]   <- L1[, i] / sqrt(L1[k, i])
    cond_var1 <- cond_var1 - L1[, i]^2
    cond_var1[k] <- -1

    # --- Schur-complement deletion from precision of candidates ---
    # Removing point k from the candidate set is equivalent to marginalising it
    # out of the precision (conditioning in covariance = marginalising in precision).
    # P_new = P[-j, -j] - outer(P[-j, j], P[j, -j]) / P[j, j]
    # This one-liner combines the rank-1 update and the row/column deletion.
    prec <- prec[-j, -j, drop = FALSE] -
      outer(prec[-j, j], prec[j, -j]) / prec[j, j]

    # Update cond_var2 for remaining candidates using the new diagonal.
    cond_var2[candidates] <- diag(prec)
  }

  indexes
}

# Greedy max-MI selection -- Cholesky-based O(n^3) implementation.
#
# Like mi_prec, but maintains a Cholesky factorisation of the precision matrix
# for the unselected set, avoiding explicit matrix inversions. The precision
# Cholesky L2 is initialised once via a "flipped Cholesky" trick, then updated
# by rank-1 DOWNDATES as points are removed from the candidate set.
#
# Initialisation of L2 (the lower Cholesky factor of K^{-1}):
#   Let X_rev be X with rows reversed. Compute U = chol(kernel(X_rev, X_rev)).
#   Because reversing permutes the kernel matrix (K_rev = J K J where J is the
#   flip permutation), the inverse is  K^{-1} = J (t(U) U)^{-1} J = L2^T L2
#   where L2 is obtained by flipping U^{-T} in both dimensions and transposing.
#   This gives a valid lower Cholesky factor of K^{-1} without explicitly inverting K.
#
# See entropy_naive for argument documentation.
mi_chol <- function(X, kernel, kernel_diag, s) {
  n <- nrow(X)
  s <- min(s, n)

  indexes    <- integer(s)
  candidates <- seq_len(n)
  L1         <- matrix(0, n, s)   # left-looking Cholesky for cond_var1 (Var given S)
  cond_var1  <- kernel_diag(X)

  # --- Initialise L2: lower Cholesky factor of K^{-1} ---
  # Reverse the rows of X so that the kernel matrix evaluates as K_rev = J K J.
  X_rev     <- X[n:1, , drop = FALSE]
  K_rev     <- kernel(X_rev, X_rev)
  # chol() returns upper U such that K_rev = t(U) U.
  # Lower Cholesky:  L_rev = t(U)  satisfies  L_rev L_rev^T = K_rev.
  L_rev_lower <- t(chol(K_rev))
  # Solve L_rev Z = I by forward substitution: Z = L_rev^{-1}  (n x n matrix).
  Z <- forwardsolve(L_rev_lower, diag(n))
  # Flip Z in both dimensions and transpose to get L2 = (J Z J)^T.
  # Proof that L2^T L2 = K^{-1}:
  #   L2^T L2 = J Z J (J Z^T J) = J (Z Z^T) J = J K_rev^{-1} J
  #           = J (J K J)^{-1} J = K^{-1}.
  Z_flipped <- Z[n:1, n:1]      # flip rows and columns
  L2        <- t(Z_flipped)     # L2 is lower triangular

  # The full conditional variance of each unselected point = row sum of squares of L2,
  # since diag(K^{-1}) = diag(L2^T L2) = rowSums(L2 * L2).
  cond_var2 <- rowSums(L2 * L2)

  for (i in seq_len(s)) {
    k <- which.max(cond_var1 * cond_var2)
    indexes[i] <- k
    j <- which(candidates == k)  # position of k within current candidates
    candidates <- candidates[-j]

    # --- Update L1 and cond_var1 (same as entropy_chol) ---
    L1[, i] <- as.vector(kernel(X, X[k, , drop = FALSE]))
    if (i > 1) {
      L1[, i] <- L1[, i] - L1[, 1:(i - 1), drop = FALSE] %*% L1[k, 1:(i - 1)]
    }
    L1[, i]   <- L1[, i] / sqrt(L1[k, i])
    cond_var1 <- cond_var1 - L1[, i]^2
    cond_var1[k] <- -1

    # --- Rank-1 downdate of L2 (Cholesky of precision of remaining candidates) ---
    # Removing candidate j is a rank-1 downdate:  L2 L2^T - u u^T
    # where u = L2 L2^T[:,j] / sqrt([L2 L2^T]_{jj}) = (K^{-1}[:,j]) / sqrt(K^{-1}_{jj}).
    # We compute u directly from L2 without forming K^{-1} explicitly.
    u <- drop(L2 %*% L2[j, ])       # = K^{-1} e_j  (j-th column of K^{-1})
    u <- u / sqrt(u[j])              # normalise by sqrt of j-th diagonal
    L2 <- chol_downdate(L2, u, j)   # rank-1 downdate
    L2 <- L2[-j, -j, drop = FALSE]  # remove row and column j

    # Update cond_var2 for remaining candidates.
    cond_var2[candidates] <- rowSums(L2 * L2)
  }

  indexes
}

# =============================================================================
# Section 4 -- Nystrom ordering (greedy KL divergence minimisation)
#
# Chooses a subset S to minimise the KL divergence between the full GP and its
# Nystrom approximation induced by S:
#   KL(GP_full || GP_Nystrom(S)) = (1/2) * sum_{j not in S} log(sigma_j^2 / sigma_hat_j^2)
# where sigma_j^2 is the full marginal variance and sigma_hat_j^2 is the Nystrom
# approximation's variance at j given S.
#
# Complexities:
#   nystrom_order_naive : O(s * n * n s)   = O(n^2 s^2)  -- inner loop over candidates
#   nystrom_order       : O(s * n^2)       = O(n^2 s)    -- vectorised over candidates
# =============================================================================

# Greedy Nystrom ordering -- naive O(n^2 s^2) implementation.
#
# For each candidate j at each step, tentatively adds j and computes the
# resulting KL score; picks the j that reduces it most. This inner loop is
# O(n * n) per outer step, giving O(n^2 s^2) overall.
#
# See entropy_naive for argument documentation.
nystrom_order_naive <- function(X, kernel, kernel_diag, s) {
  n <- nrow(X)
  s <- min(s, n)

  indexes    <- integer(s)
  candidates <- seq_len(n)
  L          <- matrix(0, n, s)   # left-looking Cholesky factor (same as entropy_chol)
  cond_var   <- kernel_diag(X)

  for (i in seq_len(s)) {
    # For each candidate, try adding it and compute the KL score = sum of log-ratios
    # over the remaining candidates.  Pick the one that minimises the KL increase.
    best_k     <- -1L
    best_score <- -Inf

    for (jj in seq_along(candidates)) {
      j_try         <- candidates[jj]
      other_cands   <- candidates[-jj]

      # Tentative conditional variance update for point j_try.
      L_try         <- as.vector(kernel(X, X[j_try, , drop = FALSE]))
      if (i > 1) {
        L_try <- L_try - L[, 1:(i - 1), drop = FALSE] %*% L[j_try, 1:(i - 1)]
      }
      L_try         <- L_try / sqrt(L_try[j_try])
      ncond_var     <- cond_var - L_try^2

      # KL contribution: sum over remaining candidates of  log(old_var / new_var).
      score         <- sum(log(cond_var[other_cands] / ncond_var[other_cands]))
      if (score > best_score) {
        best_score <- score
        best_k     <- j_try
      }
    }

    indexes[i] <- best_k
    jj         <- which(candidates == best_k)
    candidates <- candidates[-jj]

    # Commit the Cholesky update for the selected point.
    L[, i] <- as.vector(kernel(X, X[best_k, , drop = FALSE]))
    if (i > 1) {
      L[, i] <- L[, i] - L[, 1:(i - 1), drop = FALSE] %*% L[best_k, 1:(i - 1)]
    }
    L[, i]   <- L[, i] / sqrt(L[best_k, i])
    cond_var <- cond_var - L[, i]^2
    cond_var[best_k] <- -1
  }

  indexes
}

# Greedy Nystrom ordering -- efficient O(n^2 s) implementation.
#
# Instead of re-computing each tentative update independently, maintains the
# full CONDITIONAL COVARIANCE matrix  Cov_S(X, X) = K - K[:,S] K[S,S]^{-1} K[S,:]
# and derives all candidate scores simultaneously using a vectorised formula.
#
# Score of candidate j (after conditioning on the current S):
#   improve(j) = sum_{k not in S, k ≠ j}  log( cond_var[k] / (cond_var[k] - cond_cov[j,k]^2 / cond_var[j]) )
#
# The candidate minimising the KL divergence is picked (which.min, not which.max).
#
# The conditional covariance matrix is updated by a rank-1 downdate after each selection.
#
# See entropy_naive for argument documentation.
nystrom_order <- function(X, kernel, kernel_diag, s) {
  n <- nrow(X)
  s <- min(s, n)

  indexes    <- integer(s)
  candidates <- seq_len(n)       # shrinks by one each iteration
  cond_cov   <- kernel(X, X)    # starts as full K; updated by rank-1 downdates

  for (i in seq_len(s)) {
    # Current conditional variances along the diagonal.
    cond_var_vec <- diag(cond_cov)  # length = current n_remaining

    # Build the score matrix: scores[a, b] = cond_var_vec[a] - cond_cov[a,b]^2 / cond_var_vec[b]
    # This is the conditional variance of point a if b were added next (contribution to KL score).
    # Broadcasting: cond_var_vec as a column (row-varying) minus scaled cond_cov^2
    scores <- outer(cond_var_vec, rep(1, length(cond_var_vec))) -
      t(t(cond_cov^2) / cond_var_vec)

    # Set the diagonal to cond_var_vec (a point doesn't affect its own variance).
    diag(scores) <- cond_var_vec

    # KL score for candidate j = sum_a log(sigma_a^2 / scores[a, j]).
    # Equivalently: -sum_a log(scores[a, j] / sigma_a^2) = -colSums(log(scores / outer(cond_var_vec, rep(1,...)))).
    # Using log(ifelse(scores > 0, scores, 1)) avoids log of zero on the diagonal.
    log_scores <- log(ifelse(scores > 0, scores, 1))
    kl_scores  <- colSums(log_scores)

    # Pick the candidate (position within current candidates) that minimises KL.
    j           <- which.min(kl_scores)
    k           <- candidates[j]
    indexes[i]  <- k
    candidates  <- candidates[-j]

    # --- Rank-1 downdate of conditional covariance ---
    # After adding point k, cond_cov shrinks:
    #   cond_cov_new = cond_cov[-j, -j] - outer(u, u)
    # where u = cond_cov[:, j] / sqrt(cond_cov[j, j]).
    u        <- cond_cov[, j] / sqrt(cond_cov[j, j])
    cond_cov <- cond_cov - outer(u, u)
    cond_cov <- cond_cov[-j, -j, drop = FALSE]
  }

  indexes
}
