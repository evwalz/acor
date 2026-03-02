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
#' @return Character: "v1", "v2", or "original"
#' @keywords internal
#' @noRd
select_kernel_version <- function(Y, X) {
  n <- length(Y)
  M <- length(unique(Y))
  
  # Handle X: compute R for each column, take minimum
  if (is.vector(X)) {
    R <- length(unique(X))
    m <- 1
  } else {
    R <- min(apply(X, 2, function(col) length(unique(col))))
    m <- ncol(X)
  }
  
  if (M == 2) {
    return("v1")
  } else if (M / n < 0.25) {
    # More than 25% ties in Y → use V1
    return("v1")
  } else if (R / n < 0.25) {
    # More than 25% ties in X → use V1
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
#' @param version Character: "v1", "v2", or "original"
#' @return List with akc, var, and var_ind
#' @keywords internal
#' @noRd
compute_akc_variance_auto <- function(X, Y, IID = TRUE, version = "original") {
  
  if (version == "original") {
    if (IID) {
      return(Sigma_akc(X, Y))
    } else {
      return(Sigma_akc_ts(X, Y))
    }
  }
  
  # For v1 and v2, compute tau values first (shared across var and var_ind)
  if (is_binary(Y)) {
    tau_Y_result <- tau_Y_func_binary(Y)
    akc_result   <- kendall_tau_sign_binary(X, Y)
  } else {
    tau_Y_result <- tau_Y_func(Y)
    akc_result   <- kendall_tau_sign(X, Y)
  }
  
  tau_Y  <- tau_Y_result$expectation
  p_Y    <- tau_Y_result$p_tie_y
  tau_XY <- akc_result$expectation
  
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
#' @param version Character: "v1", "v2", or "original"
#' @return List with akc_vector, Sigma, and Sigma_ind
#' @keywords internal
#' @noRd
compute_akc_multivariate_variance_auto <- function(X, Y, IID = TRUE, version = "original") {
  
  if (IID) {
    if (version == "original") {
      return(Sigma_akc_multivariate(X, Y))
    } else if (version == "v1") {
      return(Sigma_akc_multivariate_v1(X, Y))
    } else {  # v2
      return(Sigma_akc_multivariate_v2(X, Y))
    }
  } else {
    if (version == "original") {
      return(Sigma_akc_multivariate_ts(X, Y))
    } else if (version == "v1") {
      return(Sigma_akc_multivariate_ts_v1(X, Y))
    } else {  # v2
      return(Sigma_akc_multivariate_ts_v2(X, Y))
    }
  }
}


# ============================================================================
# AGC KERNEL VERSION SELECTION
# ============================================================================

#' Determine optimal AGC kernel version based on data characteristics
#' 
#' AGC has two implementations:
#' - "original": kernel_ties_optim2 using sign matrices — O(R*M + R*N + M*N)
#'   where R = unique X values, M = unique Y values.  Fast when R or M are
#'   small (ties/discrete data), but O(n^2) when data is continuous.
#' - "v2": Fenwick tree based kernel — O(n log n).  Faster for large n
#'   continuous data.
#' 
#' @param Y Numeric outcome vector
#' @param X Numeric predictor vector or matrix
#' @return Character: "v2" or "original"
#' @keywords internal
#' @noRd
select_agc_kernel_version <- function(Y, X) {
  n <- length(Y)
  M <- length(unique(Y))
  
  # Binary Y gets its own specialization
  if (M == 2) {
    return("binary")
  } else {
    return("v2")
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
#' @param version Character: "v2" or "original"
#' @return List with agc, var, and var_ind
#' @keywords internal
#' @noRd
compute_agc_variance_auto <- function(y_rank, x_rank, IID = TRUE, version = "original") {
  
  if (IID) {
    if (version == "binary") {
      return(Sigma_agc_binary(y_rank, x_rank))
    } else if (version == "original") {
      return(Sigma_agc(y_rank, x_rank))
    } else {  # v2
      return(Sigma_agc_v2(y_rank, x_rank))
    }
  } else {
    if (version == "binary") {
      return(Sigma_agc_ts_binary(y_rank, x_rank))
    } else if (version == "original") {
      return(Sigma_agc_ts(y_rank, x_rank))
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
#' @param version Character: "v2" or "original"
#' @return List with agc_vector, Sigma, and Sigma_ind
#' @keywords internal
#' @noRd
compute_agc_multivariate_variance_auto <- function(y_rank, xarray_ranks, IID = TRUE, version = "original") {
  
  if (IID) {
    if (version == "binary") {
      return(Sigma_agc_multivariate_binary(y_rank, xarray_ranks))
    } else if (version == "original") {
      return(Sigma_agc_multivariate(y_rank, xarray_ranks))
    } else {  # v2
      return(Sigma_agc_multivariate_v2(y_rank, xarray_ranks))
    }
  } else {
    if (version == "binary") {
      return(Sigma_agc_multivariate_ts_binary(y_rank, xarray_ranks))
    } else if (version == "original") {
      return(Sigma_agc_multivariate_ts(y_rank, xarray_ranks))
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