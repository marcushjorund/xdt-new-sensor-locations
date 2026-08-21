#Creating functions/helper-functions for the near-optimal sensor selection algoritm
#The algorithm in Krause et al. (2008) needs a covariance matrix in order to compute the near-optimal sensor locations

#' Build the CAR precision and covariance matrices for greedy MI sensor selection
#'
#' Reconstructs the weighted Besag-proper precision matrix `Q = tau * (d*I + D_W - W)` from
#' the natural-scale hyperparameters and edge distances, then inverts it (via Cholesky) to
#' get the covariance matrix used by `greedy_mi_sensor_selection()`. Optionally collapses the
#' directed WITH/AGAINST link pairs into a single summed-covariance space via `with_and_against`.
#'
#' @param adjacency_matrix Square sparse adjacency matrix of the traffic links.
#' @param tau       CAR precision hyperparameter (natural scale).
#' @param d         CAR diagonal-offset hyperparameter (natural scale).
#' @param sigma     RBF/Laplacian kernel length-scale (natural scale).
#' @param distances Numeric vector of edge Gower distances (upper-triangle order), as
#'   returned by `fit_weighted_inla_model()$distances`.
#' @param data      Traffic links data frame with an `id` column ending in `-WITH`/`-AGAINST`.
#'   Required when `with_and_against = TRUE`.
#' @param with_and_against Logical (default `FALSE`). When `TRUE`, sums each base link's
#'   WITH and AGAINST covariance via a linear transform `C`, returning `covariance_matrix_sum`
#'   and `ids` (base IDs) instead of the raw n x n matrices.
#' @return A list. When `with_and_against = FALSE`: `precision_matrix`, `covariance_matrix`,
#'   `weights`. When `with_and_against = TRUE`: `covariance_matrix_sum`, `ids`,
#'   `precision_matrix`, `covariance_matrix`, `weights`.
create_covariance_and_precision_matrix <- function(
    adjacency_matrix,
    tau,
    d,
    sigma,
    distances,   # dist_vals vector from fit_weighted_inla_model()$distances
    data = NULL,
    with_and_against = FALSE
) {
  n <- nrow(adjacency_matrix)
  
  # Reconstruct upper-triangle edge indices (same logic as fit_weighted_inla_model)
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

