# Internal Helper Functions for acor()
# These functions are not exported and are used internally by acor() and acor.test()

# ============================================================================
# INPUT VALIDATION
# ============================================================================

#' Validate and standardize X, Y inputs for acor functions
#' @return List with validated X (matrix), Y, n, m
#' @keywords internal
#' @noRd
validate_acor_inputs <- function(X, Y) {
  if (is.data.frame(Y)) Y <- as.numeric(Y[[1]])
  if (is.data.frame(X)) X <- as.matrix(X)
  
  if (!is.numeric(Y) || !is.vector(Y)) {
    stop("Y must be a numeric vector")
  }

  if (is.vector(X)) {
    if (!is.numeric(X)) stop("X must be a numeric vector or matrix")
    X <- matrix(X, ncol = 1)
  } else if (!is.matrix(X) || !is.numeric(X)) {
    stop("X must be a numeric vector or matrix")
  }
  
  n <- length(Y)
  if (nrow(X) != n) {
    stop("X and Y must have the same number of observations")
  }
  if (any(is.na(X)) || any(is.na(Y))) stop("NA values not supported; remove NAs first")
  if (length(unique(Y)) < 2) stop("Y must have at least 2 distinct values")
  m <- ncol(X)
  list(X = X, Y = Y, n = n, m = m)
}


# ============================================================================
# KERNEL VERSION SELECTION
# ============================================================================

#' Determine optimal kernel version based on data characteristics
#' 
#' @param Y Numeric outcome vector
#' @param X Numeric predictor vector or matrix
#' @return Character: "v1" or "v2" 
#' @keywords internal
#' @noRd
select_kernel_version <- function(Y, X) {
  n <- length(Y)
  M <- length(unique(Y))
  
  if (M == 2) {
    return("v1")
  } else {
    return("v2")
  }
}


# ============================================================================
# INTERNAL DISPATCHER FUNCTIONS
# ============================================================================

#' Compute univariate AKC with variance using specified kernel version
#' 
#' @param X Numeric predictor vector
#' @param Y Numeric outcome vector
#' @param IID Logical; if TRUE uses IID variance, if FALSE uses HAC
#' @param version Character: "v1" or "v2"
#' @return List with akc, var, and var_ind
#' @keywords internal
#' @noRd
compute_akc_variance_auto <- function(X, Y, IID = TRUE, version = "v2") {
  
  # For v1 and v2, compute tau values first (shared across var and var_ind)
  if (is_binary(Y)) {
    tau_Y_result <- tau_Y_func_binary(Y)
    akc_result   <- kendall_tau_sign_binary(X, Y)
  } else {
    tau_Y_result <- tau_Y_func(Y)
    akc_result   <- kendall_tau_sign_cpp(X, Y)
  }
  
  tau_Y  <- tau_Y_result$expectation
  p_Y    <- tau_Y_result$p_tie_y
  tau_XY <- akc_result$expectation
  
  if (1 - p_Y < 1e-10) {
    stop("Y has near-total ties (nearly constant); ",
         "AKC variance is undefined")
  }
  
  if (IID) {
    if (version == "v1") {
      return(Sigma_akc_v1(X, Y, tau_XY, tau_Y, p_Y))
    } else {  # v2
      return(Sigma_akc_v2(X, Y, tau_XY, tau_Y, p_Y))
    }
  } else {
    if (version == "v1") {
      return(Sigma_akc_ts_v1(X, Y, tau_XY, tau_Y, p_Y))
    } else {  # v2
      return(Sigma_akc_ts_v2(X, Y, tau_XY, tau_Y, p_Y))
    }
  }
}

#' Compute multivariate AKC with covariance matrix using specified kernel version
#' 
#' @param X Numeric matrix of predictors (n x m)
#' @param Y Numeric outcome vector
#' @param IID Logical; if TRUE uses IID variance, if FALSE uses HAC
#' @param version Character: "v1"or "v2"
#' @return List with akc_vector, Sigma, and Sigma_ind
#' @keywords internal
#' @noRd
compute_akc_multivariate_variance_auto <- function(X, Y, IID = TRUE, version = "v2") {
  
  if (IID) {
    if (version == "v1") {
      return(Sigma_akc_multivariate_v1(X, Y))
    } else {  # v2
      return(Sigma_akc_multivariate_v2(X, Y))
    }
  } else {
    if (version == "v1") {
      return(Sigma_akc_multivariate_ts_v1(X, Y))
    } else {  # v2
      return(Sigma_akc_multivariate_ts_v2(X, Y))
    }
  }
}


# ============================================================================
# AGC INTERNAL DISPATCHER FUNCTIONS
# ============================================================================

#' Compute univariate AGC with variance using specified kernel version
#' 
#' @param y_rank Numeric vector of average ranks for Y
#' @param x_rank Numeric vector of average ranks for X
#' @param IID Logical; if TRUE uses IID variance, if FALSE uses HAC
#' @param version Character: "v2" or "binary"
#' @return List with agc, var, and var_ind
#' @keywords internal
#' @noRd
compute_agc_variance_auto <- function(y_rank, x_rank, IID = TRUE) {
  binary <- is_binary(y_rank)
  if (IID) {
    if (binary) {
      return(Sigma_agc_binary(y_rank, x_rank))
    } else {  # v2
      return(Sigma_agc_v2(y_rank, x_rank))
    }
  } else {
    if (binary) {
      return(Sigma_agc_ts_binary(y_rank, x_rank))
    } else {  # v2
      return(Sigma_agc_ts_v2(y_rank, x_rank))
    }
  }
}

#' Compute multivariate AGC with covariance matrix using specified kernel version
#' 
#' @param y_rank Numeric vector of average ranks for Y
#' @param xarray_ranks Matrix of ranks (predictors x n)
#' @param IID Logical; if TRUE uses IID variance, if FALSE uses HAC
#' @param version Character: "v2" or "binary"
#' @return List with agc_vector, Sigma, and Sigma_ind
#' @keywords internal
#' @noRd
compute_agc_multivariate_variance_auto <- function(y_rank, xarray_ranks, IID = TRUE) {
  binary <- is_binary(y_rank)
  if (IID) {
    if (binary) {
      return(Sigma_agc_multivariate_binary(y_rank, xarray_ranks))
    } else {  # v2
      return(Sigma_agc_multivariate_v2(y_rank, xarray_ranks))
    }
  } else {
    if (binary) {
      return(Sigma_agc_multivariate_ts_binary(y_rank, xarray_ranks))
    } else {  # v2
      return(Sigma_agc_multivariate_ts_v2(y_rank, xarray_ranks))
    }
  }
}

#' @keywords internal
#' @noRd
is_binary <- function(Y) {
  unique_vals <- unique(Y)
  length(unique_vals) == 2
}

#' Univariate HAC long-run variance of grade functions
#'
#' Computes the Bartlett-kernel weighted long-run variance
#' of centered grade functions for X and Y. All univariate
#' HAC independence variance functions are multiples of this.
#'
#' @param x_grade_centered Numeric vector: (rank - 0.5)/N - 0.5
#' @param y_grade_centered Numeric vector: (rank - 0.5)/N - 0.5
#' @param N Integer sample size.
#' @param b Integer bandwidth (Bartlett kernel).
#' @return Scalar long-run variance (unscaled).
#' @keywords internal
#' @noRd
ind_lrv_univariate <- function(x_grade_centered, y_grade_centered, N, b) {
  h_vec <- 1:(N - 1)
  w <- pmax(1 - abs(h_vec) / (b + 1), 0)
  
  x_autoc <- stats::acf(x_grade_centered, plot = FALSE, type = "covariance",
                        demean = FALSE, lag.max = N - 1)$acf
  y_autoc <- stats::acf(y_grade_centered, plot = FALSE, type = "covariance",
                        demean = FALSE, lag.max = N - 1)$acf
  
  sum(x_autoc[1] * y_autoc[1], 2 * (w * x_autoc[-1] * y_autoc[-1]))
}


#' Multivariate HAC long-run covariance of grade functions
#'
#' Computes the Bartlett-kernel weighted long-run covariance matrix
#' for multiple X predictors against a single Y. All multivariate
#' HAC independence covariance functions are multiples of this.
#'
#' @param x_grades_centered Matrix of centered grades.
#'   Rows = predictors (k), columns = observations (N) for AGC/rho_a convention.
#'   Rows = observations (N), columns = predictors (m) for AKC/tau_a convention.
#' @param y_grade_centered Numeric vector of length N.
#' @param N Integer sample size.
#' @param b Integer bandwidth.
#' @param x_by_row Logical; TRUE if x_grades_centered is k x N (AGC/rho_a),
#'   FALSE if N x m (AKC/tau_a).
#' @return k x k (or m x m) unscaled long-run covariance matrix.
#' @keywords internal
#' @noRd
ind_lrv_multivariate <- function(x_grades_centered, y_grade_centered, N, b,
                                 x_by_row = FALSE) {
  # Normalize to k x N layout
  if (!x_by_row) {
    # Input is N x m, transpose to m x N
    x_grades_centered <- t(x_grades_centered)
  }
  k <- nrow(x_grades_centered)
  
  h_vec <- 1:(N - 1)
  w <- pmax(1 - abs(h_vec) / (b + 1), 0)
  
  y_autoc <- stats::acf(y_grade_centered, plot = FALSE, type = "covariance",
                        demean = FALSE, lag.max = N - 1)$acf
  
  Sigma <- matrix(0, nrow = k, ncol = k)
  
  for (j in 1:k) {
    for (l in j:k) {
      x_grade_j <- x_grades_centered[j, ]
      x_grade_l <- x_grades_centered[l, ]
      
      xcov_0 <- mean(x_grade_j * x_grade_l)
      hac_sum <- xcov_0 * y_autoc[1]
      
      for (h in seq_len(min(b, N - 1))) {
        xcov_h <- mean(x_grade_j[1:(N - h)] * x_grade_l[(h + 1):N] +
                         x_grade_j[(h + 1):N] * x_grade_l[1:(N - h)]) / 2
        hac_sum <- hac_sum + 2 * w[h] * xcov_h * y_autoc[h + 1]
      }
      
      Sigma[j, l] <- hac_sum
      if (j != l) Sigma[l, j] <- hac_sum
    }
  }
  
  Sigma
}
