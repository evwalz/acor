# Asymmetric Kendall Correlation (AKC) Functions
# Refactored to eliminate duplicated logic via shared helpers.
# All Sigma_* functions return both main variance and independence variance
# to maintain compatibility with acor_internal.R dispatchers.


# ============================================================================
# CORE HELPERS
# ============================================================================

#' Dispatch tau_Y computation (binary-optimized or general)
#'
#' @param Y Numeric outcome vector.
#' @return List with elements \code{expectation} and \code{p_tie_y}.
#' @keywords internal
#' @noRd
compute_tau_Y <- function(Y) {
  if (is_binary(Y)) tau_Y_func_binary(Y) else tau_Y_func(Y)
}

#' Dispatch Kendall tau sign computation (binary-optimized or general)
#'
#' @param X Numeric predictor vector.
#' @param Y Numeric outcome vector.
#' @return List with elements \code{tau} and \code{expectation}.
#' @keywords internal
#' @noRd
compute_kendall <- function(X, Y) {
  if (is_binary(Y)) kendall_tau_sign_binary(X, Y) else kendall_tau_sign_cpp(X, Y)
}

#' Ensure X is a matrix
#' @param X Numeric vector or matrix.
#' @return Numeric matrix.
#' @keywords internal
#' @noRd
ensure_matrix <- function(X) {
  if (is.vector(X)) matrix(X, ncol = 1) else X
}


# ============================================================================
# SHARED VARIANCE COMPUTATION HELPERS
# ============================================================================

#' Bartlett-kernel HAC correction for a univariate adjusted-kernel series
#'
#' @param adjusted_K Numeric vector of adjusted kernel values (length n).
#' @return Scalar HAC correction term.
#' @keywords internal
#' @noRd
hac_correction_univariate <- function(adjusted_K) {
  n <- length(adjusted_K)
  b <- floor(2 * n^(1 / 3))
  correction <- 0
  
  for (h in seq_len(b)) {
    omega <- 1 - h / (b + 1)
    autocov_h <- (2 / n) * sum(adjusted_K[1:(n - h)] * adjusted_K[(h + 1):n])
    correction <- correction + omega * autocov_h
  }
  
  correction
}

#' Bartlett-kernel HAC correction for a multivariate adjusted-kernel matrix
#'
#' @param adjusted_K_tau Numeric matrix (n x m) of adjusted kernel values.
#' @return m x m HAC correction matrix.
#' @keywords internal
#' @noRd
hac_correction_multivariate <- function(adjusted_K_tau) {
  n <- nrow(adjusted_K_tau)
  m <- ncol(adjusted_K_tau)
  b <- floor(2 * n^(1 / 3))
  Sigma_hac <- matrix(0, nrow = m, ncol = m)
  
  for (h in seq_len(b)) {
    omega <- 1 - h / (b + 1)
    K_lag  <- adjusted_K_tau[1:(n - h), , drop = FALSE]
    K_lead <- adjusted_K_tau[(h + 1):n, , drop = FALSE]
    autocov_h <- (t(K_lag) %*% K_lead + t(K_lead) %*% K_lag) / n
    Sigma_hac <- Sigma_hac + omega * autocov_h
  }
  
  Sigma_hac
}

#' IID variance from a univariate adjusted-kernel vector
#'
#' @param adjusted_K Numeric vector.
#' @param scale_factor Scalar.
#' @return Scalar variance estimate.
#' @keywords internal
#' @noRd
iid_variance_univariate <- function(adjusted_K, scale_factor) {
  scale_factor * mean(adjusted_K^2)
}

#' HAC variance from a univariate adjusted-kernel vector
#'
#' @param adjusted_K Numeric vector.
#' @param scale_factor Scalar.
#' @return Scalar variance estimate (IID + HAC).
#' @keywords internal
#' @noRd
hac_variance_univariate <- function(adjusted_K, scale_factor) {
  iid <- mean(adjusted_K^2)
  hac <- hac_correction_univariate(adjusted_K)
  scale_factor * (iid + hac)
}

#' IID covariance matrix from a multivariate adjusted-kernel matrix
#'
#' @param adjusted_K_tau Numeric matrix (n x m).
#' @param scale_factor Scalar.
#' @return m x m covariance matrix.
#' @keywords internal
#' @noRd
iid_covariance_multivariate <- function(adjusted_K_tau, scale_factor) {
  n <- nrow(adjusted_K_tau)
  scale_factor * (t(adjusted_K_tau) %*% adjusted_K_tau) / n
}

#' HAC covariance matrix from a multivariate adjusted-kernel matrix
#'
#' @param adjusted_K_tau Numeric matrix (n x m).
#' @param scale_factor Scalar.
#' @return m x m covariance matrix (IID + HAC).
#' @keywords internal
#' @noRd
hac_covariance_multivariate <- function(adjusted_K_tau, scale_factor) {
  n <- nrow(adjusted_K_tau)
  Sigma_iid <- (t(adjusted_K_tau) %*% adjusted_K_tau) / n
  Sigma_hac <- hac_correction_multivariate(adjusted_K_tau)
  scale_factor * (Sigma_iid + Sigma_hac)
}

#' Compute adjusted kernel values for a single predictor
#'
#' @param K_tau_values Numeric vector of K_tau values.
#' @param K_p_values Numeric vector of K_p values.
#' @param tau_XY Scalar expectation from Kendall tau.
#' @param p_Y Scalar probability of ties in Y.
#' @return Numeric vector of adjusted kernel values.
#' @keywords internal
#' @noRd
compute_adjusted_K <- function(K_tau_values, K_p_values, tau_XY, p_Y) {
  K_tau_values - (tau_XY / (1 - p_Y)) * K_p_values
}

#' Build adjusted kernel matrix for multiple predictors
#'
#' @param K_tau_matrix n x m matrix of K_tau values.
#' @param K_p_values Length-n vector of K_p values.
#' @param tau_vector Length-m vector of tau expectations.
#' @param p_Y Scalar.
#' @return n x m matrix of adjusted kernel values.
#' @keywords internal
#' @noRd
compute_adjusted_K_matrix <- function(K_tau_matrix, K_p_values, tau_vector, p_Y) {
  m <- ncol(K_tau_matrix)
  adjusted <- K_tau_matrix
  for (k in seq_len(m)) {
    adjusted[, k] <- K_tau_matrix[, k] - (tau_vector[k] / (1 - p_Y)) * K_p_values
  }
  adjusted
}


# ============================================================================
# INDEPENDENCE VARIANCE HELPERS
# ============================================================================
# These compute the variance under the null hypothesis of independence.
# Used by all Sigma_* functions to return var_ind / Sigma_ind.

#' Univariate IID independence variance for AKC
#'
#' @param X Numeric predictor vector.
#' @param Y Numeric outcome vector.
#' @param p_Y Scalar probability of ties in Y (pre-computed).
#' @return Scalar independence variance.
#' @keywords internal
#' @noRd
ind_variance_akc_iid <- function(X, Y, p_Y) {
  N <- length(Y)
  y_rank <- rank(Y, ties.method = "average")
  x_rank <- rank(X, ties.method = "average")
  
  var_y_rank <- sum((y_rank - mean(y_rank))^2) / N
  zeta_3Y <- 1 - (12 / N^2) * var_y_rank
  
  var_x_rank <- sum((x_rank - mean(x_rank))^2) / N
  zeta_3X <- 1 - (12 / N^2) * var_x_rank
  
  (4 / 9) * ((1 - zeta_3X) * (1 - zeta_3Y)) / (1 - p_Y)^2
}

#' Univariate HAC independence variance for AKC
#'
#' @param X Numeric predictor vector.
#' @param Y Numeric outcome vector.
#' @param p_Y Scalar probability of ties in Y (pre-computed).
#' @return Scalar independence variance with HAC correction.
#' @keywords internal
#' @noRd
ind_variance_akc_hac <- function(X, Y, p_Y) {
  N <- length(Y)

  b <- floor(2 * N^(1 / 3))
  h_vec <- 1:(N - 1)
  w <- pmax(1 - abs(h_vec) / (b + 1), 0)

  x_grade_centered <- (rank(X, ties.method = "average") - 0.5) / N - 0.5
  y_grade_centered <- (rank(Y, ties.method = "average") - 0.5) / N - 0.5

  x_autoc <- stats::acf(x_grade_centered, plot = FALSE, type = "covariance",
                        demean = FALSE, lag.max = N - 1)$acf
  y_autoc <- stats::acf(y_grade_centered, plot = FALSE, type = "covariance",
                        demean = FALSE, lag.max = N - 1)$acf

  tau_ind_lrv <- 64 * sum(x_autoc[1] * y_autoc[1],
                          2 * (w * x_autoc[-1] * y_autoc[-1]))

  tau_ind_lrv / (1 - p_Y)^2
}

#' Multivariate IID independence covariance for AKC
#'
#' @param X Numeric matrix (n x m).
#' @param Y Numeric outcome vector.
#' @param p_Y Scalar.
#' @return m x m independence covariance matrix.
#' @keywords internal
#' @noRd
ind_covariance_akc_iid <- function(X, Y, p_Y) {
  N <- length(Y)
  m <- ncol(X)
  
  y_rank <- rank(Y, ties.method = "average")
  var_y_rank <- sum((y_rank - mean(y_rank))^2) / N
  zeta_3Y <- 1 - (12 / N^2) * var_y_rank
  
  scale_factor <- 4 / (1 - p_Y)^2
  zeta_3X <- numeric(m)
  x_ranks <- matrix(0, nrow = N, ncol = m)
  
  for (k in seq_len(m)) {
    x_ranks[, k] <- rank(X[, k], ties.method = "average")
    var_x_rank <- sum((x_ranks[, k] - mean(x_ranks[, k]))^2) / N
    zeta_3X[k] <- 1 - (12 / N^2) * var_x_rank
  }
  
  Sigma_ind <- matrix(0, nrow = m, ncol = m)
  for (k in seq_len(m)) {
    for (l in k:m) {
      if (k == l) {
        Sigma_ind[k, k] <- (4 / 9) * (1 - zeta_3X[k]) * (1 - zeta_3Y) / (1 - p_Y)^2
      } else {
        x_grade_k <- (x_ranks[, k] - 0.5) / N - 0.5
        x_grade_l <- (x_ranks[, l] - 0.5) / N - 0.5
        rho_kl <- 12 * mean(x_grade_k * x_grade_l)
        Sigma_ind[k, l] <- scale_factor * rho_kl * (1 - zeta_3Y) / 9
        Sigma_ind[l, k] <- Sigma_ind[k, l]
      }
    }
  }
  
  Sigma_ind
}

#' Multivariate HAC independence covariance for AKC
#'
#' @param X Numeric matrix (n x m).
#' @param Y Numeric outcome vector.
#' @param p_Y Scalar.
#' @return m x m independence covariance matrix with HAC correction.
#' @keywords internal
#' @noRd
ind_covariance_akc_hac <- function(X, Y, p_Y) {
  N <- length(Y)
  m <- ncol(X)

  b <- floor(2 * N^(1 / 3))
  h_vec <- 1:(N - 1)
  w <- pmax(1 - abs(h_vec) / (b + 1), 0)

  x_grades_centered <- matrix(0, nrow = N, ncol = m)
  for (k in seq_len(m)) {
    x_grades_centered[, k] <- (rank(X[, k], ties.method = "average") - 0.5) / N - 0.5
  }
  y_grade_centered <- (rank(Y, ties.method = "average") - 0.5) / N - 0.5

  y_autoc <- stats::acf(y_grade_centered, plot = FALSE, type = "covariance",
                        demean = FALSE, lag.max = N - 1)$acf

  Sigma_ind <- matrix(0, nrow = m, ncol = m)

  for (k in seq_len(m)) {
    for (l in k:m) {
      x_grade_k <- x_grades_centered[, k]
      x_grade_l <- x_grades_centered[, l]

      xcov_0 <- mean(x_grade_k * x_grade_l)
      hac_sum_xy <- xcov_0 * y_autoc[1]

      for (h in seq_len(min(b, N - 1))) {
        xcov_h <- mean(x_grade_k[1:(N - h)] * x_grade_l[(h + 1):N] +
                         x_grade_k[(h + 1):N] * x_grade_l[1:(N - h)]) / 2
        hac_sum_xy <- hac_sum_xy + 2 * w[h] * xcov_h * y_autoc[h + 1]
      }

      Sigma_ind[k, l] <- 64 * hac_sum_xy / (1 - p_Y)^2
      if (k != l) Sigma_ind[l, k] <- Sigma_ind[k, l]
    }
  }

  Sigma_ind
}


# ============================================================================
# BINARY-OPTIMIZED CORE FUNCTIONS
# ============================================================================

#' Binary-optimized tau_Y computation -- O(1) instead of O(n^2)
#' @param Y Numeric vector with exactly 2 unique values.
#' @return List with \code{expectation} and \code{p_tie_y}.
#' @keywords internal
#' @noRd
tau_Y_func_binary <- function(Y) {
  n <- length(Y)
  unique_vals <- unique(Y)
  n0 <- sum(Y == unique_vals[1])
  n1 <- sum(Y == unique_vals[2])
  
  num_pairs <- n * (n - 1) / 2
  n_ties_y <- n0 * (n0 - 1) / 2 + n1 * (n1 - 1) / 2
  p_tie_y <- n_ties_y / num_pairs
  expectation <- 1 - p_tie_y
  
  list(expectation = expectation, p_tie_y = p_tie_y)
}

#' Binary-optimized Kendall tau sign
#' @param X Numeric predictor vector.
#' @param Y Numeric vector with exactly 2 unique values.
#' @return List with \code{tau} and \code{expectation}.
#' @keywords internal
#' @noRd
kendall_tau_sign_binary <- function(X, Y) {
  n <- length(X)
  unique_vals <- sort(unique(Y))
  idx0 <- which(Y == unique_vals[1])
  idx1 <- which(Y == unique_vals[2])
  
  X0 <- X[idx0]
  X1_sorted <- sort(X[idx1])
  n0 <- length(X0)
  n1 <- length(X1_sorted)
  
  # For each x0: count how many in X1 are strictly less, and strictly greater
  count_leq <- findInterval(X0, X1_sorted)           # X1 <= x0
  count_less <- findInterval(X0, X1_sorted, left.open = TRUE)  # X1 < x0
  count_greater <- n1 - count_leq
  
  concordant <- sum(count_greater)
  discordant <- sum(count_less)
  
  num_pairs <- n * (n - 1) / 2
  n_ties_y <- n0 * (n0 - 1) / 2 + n1 * (n1 - 1) / 2
  expectation <- (concordant - discordant) / num_pairs
  p_tie_y <- n_ties_y / num_pairs
  scale_factor <- 1 - p_tie_y
  
  tau <- if (scale_factor > 1e-10) expectation / scale_factor else 0
  
  list(tau = tau, expectation = expectation)
}


# ============================================================================
# GENERAL (NON-BINARY) CORE FUNCTIONS
# ============================================================================

#' General tau_Y computation for non-binary Y
#' @param Y Numeric outcome vector.
#' @return List with \code{expectation} and \code{p_tie_y}.
#' @keywords internal
#' @noRd
tau_Y_func <- function(Y) {
  n <- length(Y)
  num_pairs <- n * (n - 1) / 2
  freq <- as.numeric(table(Y))
  n_ties_y <- sum(freq * (freq - 1) / 2)
  p_tie_y <- n_ties_y / num_pairs
  expectation <- (num_pairs - n_ties_y) / num_pairs
  list(expectation = expectation, p_tie_y = p_tie_y)
}

#' General Kendall tau sign for non-binary Y
#' @param X Numeric predictor vector.
#' @param Y Numeric outcome vector.
#' @return List with \code{tau} and \code{expectation}.
#' @keywords internal
#' @noRd
kendall_tau_sign <- function(X, Y) {
  n <- length(X)
  num_pairs <- n * (n - 1) / 2
  
  sign_products <- numeric(num_pairs)
  n_ties_y <- 0
  idx <- 1
  
  for (i in 1:(n - 1)) {
    j_vals <- (i + 1):n
    sgn_x <- sign(X[j_vals] - X[i])
    sgn_y <- sign(Y[j_vals] - Y[i])
    
    n_j <- length(j_vals)
    sign_products[idx:(idx + n_j - 1)] <- sgn_x * sgn_y
    n_ties_y <- n_ties_y + sum(Y[j_vals] == Y[i])
    idx <- idx + n_j
  }
  
  expectation <- mean(sign_products)
  p_tie_y <- n_ties_y / num_pairs
  scale_factor <- 1 - p_tie_y
  
  tau <- if (scale_factor > 1e-10) expectation / scale_factor else 0
  
  list(tau = tau, expectation = expectation)
}


# ============================================================================
# KERNEL FUNCTIONS (point-wise, O(n) per observation)
# ============================================================================

#' @keywords internal
#' @noRd
F_bar <- function(x, X) {
  mean(X < x) + 0.5 * mean(X == x)
}

#' @keywords internal
#' @noRd
G_bar <- function(y, Y) {
  mean(Y < y) + 0.5 * mean(Y == y)
}

#' @keywords internal
#' @noRd
H_bar <- function(x, y, X, Y) {
  p_both_less <- mean((X < x) & (Y < y))
  p_x_equal_y_less <- mean((X == x) & (Y < y))
  p_x_less_y_equal <- mean((X < x) & (Y == y))
  p_both_equal <- mean((X == x) & (Y == y))
  
  p_both_less + 0.5 * p_x_equal_y_less + 0.5 * p_x_less_y_equal + 0.25 * p_both_equal
}

#' @keywords internal
#' @noRd
K_tau <- function(x, y, X, Y, tau_XY) {
  h_bar <- H_bar(x, y, X, Y)
  f_bar <- F_bar(x, X)
  g_bar <- G_bar(y, Y)
  4 * h_bar - 2 * (f_bar + g_bar) + 1 - tau_XY
}

#' @keywords internal
#' @noRd
K_p <- function(y, Y, tau_y) {
  p_y_neq_y <- mean(Y != y)
  tau_y - p_y_neq_y
}

#' @keywords internal
#' @noRd
kernel_expectation <- function(X, Y, tau_XY, tau_y, p_Y) {
  n <- length(X)
  squared_diffs <- numeric(n)
  
  for (i in 1:n) {
    k_tau_i <- K_tau(X[i], Y[i], X, Y, tau_XY)
    k_p_i <- K_p(Y[i], Y, tau_y)
    diff <- k_tau_i - (tau_XY / (1 - p_Y)) * k_p_i
    squared_diffs[i] <- diff^2
  }
  
  mean(squared_diffs)
}


# ============================================================================
# UNIVARIATE AKC (no variance)
# ============================================================================

#' Compute AKC for a single predictor (no variance)
#'
#' @param X Numeric predictor vector.
#' @param Y Numeric outcome vector.
#' @return AKC value (scalar).
#' @keywords internal
#' @noRd
compute_akc <- function(X, Y) {
  compute_kendall(X, Y)$tau
}

#' Compute AKC for multiple predictors (no variance)
#'
#' @param X Numeric matrix of predictors (n x m), or a vector.
#' @param Y Numeric outcome vector.
#' @return Numeric vector of AKC values (length m).
#' @keywords internal
#' @noRd
compute_akc_multivariate <- function(X, Y) {
  X <- ensure_matrix(X)
  m <- ncol(X)
  vapply(seq_len(m), function(k) compute_kendall(X[, k], Y)$tau, numeric(1))
}


# ============================================================================
# UNIVARIATE VARIANCE ESTIMATORS
# ============================================================================
# All return list(akc, var, var_ind) for compatibility with acor_internal.R

#' AKC variance (IID, point-wise kernels)
#'
#' @param X Numeric predictor vector.
#' @param Y Numeric outcome vector.
#' @return List with \code{akc}, \code{var}, and \code{var_ind}.
#' @keywords internal
#' @noRd
Sigma_akc <- function(X, Y) {
  tau_Y_result <- compute_tau_Y(Y)
  tau_Y <- tau_Y_result$expectation
  p_Y   <- tau_Y_result$p_tie_y
  
  akc_result <- compute_kendall(X, Y)
  akc    <- akc_result$tau
  tau_XY <- akc_result$expectation
  
  scale_factor <- 4 / (1 - p_Y)^2
  expected_val <- kernel_expectation(X, Y, tau_XY, tau_Y, p_Y)
  
  list(akc = akc,
       var = scale_factor * expected_val,
       var_ind = ind_variance_akc_iid(X, Y, p_Y))
}

#' AKC variance for time series (HAC, point-wise kernels)
#'
#' @param X Numeric predictor vector.
#' @param Y Numeric outcome vector.
#' @return List with \code{akc}, \code{var}, and \code{var_ind}.
#' @keywords internal
#' @noRd
Sigma_akc_ts <- function(X, Y) {
  n <- length(Y)
  
  tau_Y_result <- compute_tau_Y(Y)
  tau_Y <- tau_Y_result$expectation
  p_Y   <- tau_Y_result$p_tie_y
  
  akc_result <- compute_kendall(X, Y)
  akc    <- akc_result$tau
  tau_XY <- akc_result$expectation
  
  scale_factor <- 4 / (1 - p_Y)^2
  
  K_tau_values <- numeric(n)
  K_p_values   <- numeric(n)
  for (i in seq_len(n)) {
    K_tau_values[i] <- K_tau(X[i], Y[i], X, Y, tau_XY)
    K_p_values[i]   <- K_p(Y[i], Y, tau_Y)
  }
  
  adjusted_K <- compute_adjusted_K(K_tau_values, K_p_values, tau_XY, p_Y)
  
  list(akc = akc,
       var = hac_variance_univariate(adjusted_K, scale_factor),
       var_ind = ind_variance_akc_hac(X, Y, p_Y))
}


# ============================================================================
# MULTIVARIATE VARIANCE ESTIMATORS (original point-wise kernels)
# ============================================================================
# All return list(akc_vector, Sigma, Sigma_ind)

#' Multivariate AKC covariance (IID, point-wise kernels)
#'
#' @param X Numeric matrix of predictors (n x m).
#' @param Y Numeric outcome vector.
#' @return List with \code{akc_vector}, \code{Sigma}, and \code{Sigma_ind}.
#' @keywords internal
#' @noRd
Sigma_akc_multivariate <- function(X, Y) {
  X <- ensure_matrix(X)
  n <- length(Y)
  m <- ncol(X)
  
  tau_Y_result <- compute_tau_Y(Y)
  tau_Y <- tau_Y_result$expectation
  p_Y   <- tau_Y_result$p_tie_y
  
  akc_vector <- numeric(m)
  tau_vector <- numeric(m)
  K_tau_values <- matrix(0, nrow = n, ncol = m)
  K_p_values <- numeric(n)
  
  for (i in seq_len(n)) {
    K_p_values[i] <- K_p(Y[i], Y, tau_Y)
  }
  
  for (k in seq_len(m)) {
    X_k <- X[, k]
    akc_result <- compute_kendall(X_k, Y)
    akc_vector[k] <- akc_result$tau
    tau_vector[k] <- akc_result$expectation
    
    for (i in seq_len(n)) {
      K_tau_values[i, k] <- K_tau(X_k[i], Y[i], X_k, Y, tau_vector[k])
    }
  }
  
  scale_factor <- 4 / ((1 - p_Y)^2)
  adjusted_K_tau <- compute_adjusted_K_matrix(K_tau_values, K_p_values,
                                              tau_vector, p_Y)
  
  Sigma <- matrix(0, nrow = m, ncol = m)
  for (k in seq_len(m)) {
    for (l in k:m) {
      Sigma[k, l] <- scale_factor * mean(adjusted_K_tau[, k] * adjusted_K_tau[, l])
      if (k != l) Sigma[l, k] <- Sigma[k, l]
    }
  }
  
  list(akc_vector = akc_vector,
       Sigma = Sigma,
       Sigma_ind = ind_covariance_akc_iid(X, Y, p_Y))
}

#' Multivariate AKC covariance for time series (HAC, point-wise kernels)
#'
#' @param X Numeric matrix of predictors (n x m).
#' @param Y Numeric outcome vector.
#' @return List with \code{akc_vector}, \code{Sigma}, and \code{Sigma_ind}.
#' @keywords internal
#' @noRd
Sigma_akc_multivariate_ts <- function(X, Y) {
  X <- ensure_matrix(X)
  n <- length(Y)
  m <- ncol(X)
  
  tau_Y_result <- compute_tau_Y(Y)
  tau_Y <- tau_Y_result$expectation
  p_Y   <- tau_Y_result$p_tie_y
  
  akc_vector <- numeric(m)
  tau_vector <- numeric(m)
  K_tau_values <- matrix(0, nrow = n, ncol = m)
  K_p_values <- numeric(n)
  
  for (i in seq_len(n)) {
    K_p_values[i] <- K_p(Y[i], Y, tau_Y)
  }
  
  for (k in seq_len(m)) {
    X_k <- X[, k]
    akc_result <- compute_kendall(X_k, Y)
    akc_vector[k] <- akc_result$tau
    tau_vector[k] <- akc_result$expectation
    
    for (i in seq_len(n)) {
      K_tau_values[i, k] <- K_tau(X_k[i], Y[i], X_k, Y, tau_vector[k])
    }
  }
  
  scale_factor <- 4 / ((1 - p_Y)^2)
  adjusted_K_tau <- compute_adjusted_K_matrix(K_tau_values, K_p_values,
                                              tau_vector, p_Y)
  
  # Element-wise HAC covariance (for original kernel compatibility)
  b <- floor(2 * n^(1 / 3))
  Sigma <- matrix(0, nrow = m, ncol = m)
  
  for (k in seq_len(m)) {
    for (l in k:m) {
      iid_component <- mean(adjusted_K_tau[, k] * adjusted_K_tau[, l])
      
      K_k <- adjusted_K_tau[, k]
      K_l <- adjusted_K_tau[, l]
      hac_corr <- 0
      
      for (h in seq_len(b)) {
        omega <- 1 - h / (b + 1)
        autocov_h <- sum(K_k[1:(n - h)] * K_l[(h + 1):n] +
                           K_k[(h + 1):n] * K_l[1:(n - h)]) / n
        hac_corr <- hac_corr + omega * autocov_h
      }
      
      Sigma[k, l] <- scale_factor * (iid_component + hac_corr)
      if (k != l) Sigma[l, k] <- Sigma[k, l]
    }
  }
  
  list(akc_vector = akc_vector,
       Sigma = Sigma,
       Sigma_ind = ind_covariance_akc_hac(X, Y, p_Y))
}


# ============================================================================
# OPTIMIZED KERNEL COMPUTATIONS
# ============================================================================
#
# Two optimized implementations:
# 1. V1: O(R x M + n) using unique values -- faster when many ties
# 2. V2: O(n log n) using Fenwick trees -- fastest for continuous data
#
# R = number of unique X values, M = number of unique Y values


# ---------- V1: unique-values approach ----------

#' F_bar for all observations using ranks -- O(n log n)
#' @keywords internal
#' @noRd
F_bar_vec_v1 <- function(X) {
  n <- length(X)
  (rank(X, ties.method = "average") - 0.5) / n
}

#' @keywords internal
#' @noRd
G_bar_vec_v1 <- function(Y) F_bar_vec_v1(Y)

#' K_p for all observations -- O(n)
#' @keywords internal
#' @noRd
K_p_vec_v1 <- function(Y, tau_y) {
  n <- length(Y)
  Y_factor <- factor(Y)
  Y_counts <- tabulate(Y_factor)
  Y_count_map <- Y_counts[as.integer(Y_factor)]
  tau_y - (1 - Y_count_map / n)
}

#' H_bar for all observations using unique-value grid -- O(R x M + n)
#' @keywords internal
#' @noRd
H_bar_vec_v1 <- function(X, Y) {
  n <- length(X)
  
  X_unique <- sort(unique(X))
  Y_unique <- sort(unique(Y))
  R <- length(X_unique)
  M <- length(Y_unique)
  
  X_idx <- match(X, X_unique)
  Y_idx <- match(Y, Y_unique)
  
  cell_idx <- (X_idx - 1) * M + Y_idx
  counts_vec <- tabulate(cell_idx, nbins = R * M)
  counts <- matrix(counts_vec, nrow = R, ncol = M, byrow = TRUE)
  
  cum_leq <- t(apply(apply(counts, 2, cumsum), 1, cumsum))
  
  p_both_less <- matrix(0, nrow = R, ncol = M)
  if (R > 1 && M > 1) {
    p_both_less[2:R, 2:M] <- cum_leq[1:(R - 1), 1:(M - 1)]
  }
  
  p_x_eq_y_less <- matrix(0, nrow = R, ncol = M)
  if (M > 1) {
    for (i in seq_len(R)) {
      p_x_eq_y_less[i, 2:M] <- cumsum(counts[i, 1:(M - 1)])
    }
  }
  
  p_x_less_y_eq <- matrix(0, nrow = R, ncol = M)
  if (R > 1) {
    for (j in seq_len(M)) {
      p_x_less_y_eq[2:R, j] <- cumsum(counts[1:(R - 1), j])
    }
  }
  
  p_both_eq <- counts
  
  H_bar_all <- numeric(n)
  for (i in seq_len(n)) {
    xi <- X_idx[i]
    yi <- Y_idx[i]
    H_bar_all[i] <- (p_both_less[xi, yi] +
                       0.5 * p_x_eq_y_less[xi, yi] +
                       0.5 * p_x_less_y_eq[xi, yi] +
                       0.25 * p_both_eq[xi, yi]) / n
  }
  
  H_bar_all
}

#' K_tau for all observations -- O(R x M + n)
#' @keywords internal
#' @noRd
K_tau_vec_v1 <- function(X, Y, tau_XY) {
  4 * H_bar_vec_v1(X, Y) - 2 * (F_bar_vec_v1(X) + G_bar_vec_v1(Y)) + 1 - tau_XY
}


#' @keywords internal
#' @noRd
F_bar_vec_v2 <- function(X) F_bar_vec_v1(X)

#' @keywords internal
#' @noRd
G_bar_vec_v2 <- function(Y) F_bar_vec_v1(Y)

#' @keywords internal
#' @noRd
K_p_vec_v2 <- function(Y, tau_y) K_p_vec_v1(Y, tau_y)

#' K_tau for all observations -- O(n log n)
#' @keywords internal
#' @noRd
K_tau_vec_v2 <- function(X, Y, tau_XY) {
  4 * H_bar_vec_v2_cpp(X, Y) - 2 * (F_bar_vec_v2(X) + G_bar_vec_v2(Y)) + 1 - tau_XY
}


# ============================================================================
# OPTIMIZED UNIVARIATE VARIANCE ESTIMATORS (v1 / v2 kernels)
# ============================================================================
# All return list(akc, var, var_ind)

#' IID variance using v1 (unique-values) kernels
#' @keywords internal
#' @noRd
Sigma_akc_v1 <- function(X, Y, tau_XY, tau_Y, p_Y) {
  scale_factor <- 4 / (1 - p_Y)^2
  adjusted_K <- compute_adjusted_K(K_tau_vec_v1(X, Y, tau_XY),
                                   K_p_vec_v1(Y, tau_Y), tau_XY, p_Y)
  list(akc = tau_XY / (1 - p_Y),
       var = scale_factor * mean(adjusted_K^2),
       var_ind = ind_variance_akc_iid(X, Y, p_Y))
}

#' IID variance using v2 (Fenwick tree) kernels
#' @keywords internal
#' @noRd
Sigma_akc_v2 <- function(X, Y, tau_XY, tau_Y, p_Y) {
  scale_factor <- 4 / (1 - p_Y)^2
  adjusted_K <- compute_adjusted_K(K_tau_vec_v2(X, Y, tau_XY),
                                   K_p_vec_v2(Y, tau_Y), tau_XY, p_Y)
  list(akc = tau_XY / (1 - p_Y),
       var = scale_factor * mean(adjusted_K^2),
       var_ind = ind_variance_akc_iid(X, Y, p_Y))
}

#' HAC variance using v1 kernels
#' @keywords internal
#' @noRd
Sigma_akc_ts_v1 <- function(X, Y, tau_XY, tau_Y, p_Y) {
  scale_factor <- 4 / (1 - p_Y)^2
  adjusted_K <- compute_adjusted_K(K_tau_vec_v1(X, Y, tau_XY),
                                   K_p_vec_v1(Y, tau_Y), tau_XY, p_Y)
  list(akc = tau_XY / (1 - p_Y),
       var = hac_variance_univariate(adjusted_K, scale_factor),
       var_ind = ind_variance_akc_hac(X, Y, p_Y))
}

#' HAC variance using v2 kernels
#' @keywords internal
#' @noRd
Sigma_akc_ts_v2 <- function(X, Y, tau_XY, tau_Y, p_Y) {
  scale_factor <- 4 / (1 - p_Y)^2
  adjusted_K <- compute_adjusted_K(K_tau_vec_v2(X, Y, tau_XY),
                                   K_p_vec_v2(Y, tau_Y), tau_XY, p_Y)
  list(akc = tau_XY / (1 - p_Y),
       var = hac_variance_univariate(adjusted_K, scale_factor),
       var_ind = ind_variance_akc_hac(X, Y, p_Y))
}


# ============================================================================
# OPTIMIZED MULTIVARIATE ESTIMATORS (v1 / v2 kernels)
# ============================================================================
# All return list(akc_vector, Sigma, Sigma_ind)

#' Shared preamble for optimized multivariate functions
#'
#' @param X Matrix (n x m).
#' @param Y Numeric vector.
#' @param K_tau_vec_fn Function to compute vectorised K_tau (v1 or v2).
#' @param K_p_vec_fn Function to compute vectorised K_p (v1 or v2).
#' @return List with shared intermediate results.
#' @keywords internal
#' @noRd
multivariate_kernel_preamble <- function(X, Y, K_tau_vec_fn, K_p_vec_fn) {
  X <- ensure_matrix(X)
  n <- length(Y)
  m <- ncol(X)
  
  tau_Y_result <- compute_tau_Y(Y)
  tau_Y <- tau_Y_result$expectation
  p_Y   <- tau_Y_result$p_tie_y
  
  akc_vector <- numeric(m)
  tau_vector <- numeric(m)
  K_tau_values <- matrix(0, nrow = n, ncol = m)
  
  K_p_values <- K_p_vec_fn(Y, tau_Y)
  
  for (k in seq_len(m)) {
    X_k <- X[, k]
    akc_result <- compute_kendall(X_k, Y)
    akc_vector[k] <- akc_result$tau
    tau_vector[k] <- akc_result$expectation
    K_tau_values[, k] <- K_tau_vec_fn(X_k, Y, tau_vector[k])
  }
  
  scale_factor <- 4 / ((1 - p_Y)^2)
  adjusted_K_tau <- compute_adjusted_K_matrix(K_tau_values, K_p_values,
                                              tau_vector, p_Y)
  
  list(n = n, m = m, p_Y = p_Y, akc_vector = akc_vector,
       tau_vector = tau_vector, adjusted_K_tau = adjusted_K_tau,
       scale_factor = scale_factor)
}

#' Multivariate IID covariance using v1 kernels
#' @keywords internal
#' @noRd
Sigma_akc_multivariate_v1 <- function(X, Y) {
  X <- ensure_matrix(X)
  pre <- multivariate_kernel_preamble(X, Y, K_tau_vec_v1, K_p_vec_v1)
  Sigma <- iid_covariance_multivariate(pre$adjusted_K_tau, pre$scale_factor)
  list(akc_vector = pre$akc_vector,
       Sigma = Sigma,
       Sigma_ind = ind_covariance_akc_iid(X, Y, pre$p_Y))
}

#' Multivariate IID covariance using v2 kernels
#' @keywords internal
#' @noRd
Sigma_akc_multivariate_v2 <- function(X, Y) {
  X <- ensure_matrix(X)
  pre <- multivariate_kernel_preamble(X, Y, K_tau_vec_v2, K_p_vec_v2)
  Sigma <- iid_covariance_multivariate(pre$adjusted_K_tau, pre$scale_factor)
  list(akc_vector = pre$akc_vector,
       Sigma = Sigma,
       Sigma_ind = ind_covariance_akc_iid(X, Y, pre$p_Y))
}

#' Multivariate HAC covariance using v1 kernels
#' @keywords internal
#' @noRd
Sigma_akc_multivariate_ts_v1 <- function(X, Y) {
  X <- ensure_matrix(X)
  pre <- multivariate_kernel_preamble(X, Y, K_tau_vec_v1, K_p_vec_v1)
  Sigma <- hac_covariance_multivariate(pre$adjusted_K_tau, pre$scale_factor)
  list(akc_vector = pre$akc_vector,
       Sigma = Sigma,
       Sigma_ind = ind_covariance_akc_hac(X, Y, pre$p_Y))
}

#' Multivariate HAC covariance using v2 kernels
#' @keywords internal
#' @noRd
Sigma_akc_multivariate_ts_v2 <- function(X, Y) {
  X <- ensure_matrix(X)
  pre <- multivariate_kernel_preamble(X, Y, K_tau_vec_v2, K_p_vec_v2)
  Sigma <- hac_covariance_multivariate(pre$adjusted_K_tau, pre$scale_factor)
  list(akc_vector = pre$akc_vector,
       Sigma = Sigma,
       Sigma_ind = ind_covariance_akc_hac(X, Y, pre$p_Y))
}
