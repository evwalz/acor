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


##############################################################
#### AGC original version
##############################################################

# ============================================================================
# KERNEL: original (kernel_ties_optim2)
# ============================================================================

#' Optimized kernel computation with ties (vectorized version)
#' @keywords internal
#' @noRd
kernel_ties_optim2 <- function(x_rank_sort, y_rank_sort, rho,
                               y_rank_unique, y_num, pos_y, ix_y,
                               x_rank_unique) {
  N <- length(x_rank_sort)
  R <- length(x_rank_unique)
  M <- length(y_rank_unique)
  
  G_x_unique <- (x_rank_unique - 0.5) / N
  G_y <- (y_rank_unique - 0.5) / N
  
  x_rank_indices <- match(x_rank_sort, x_rank_unique)
  G_x <- G_x_unique[x_rank_indices]
  
  # Vectorized sign matrices
  sign_x <- outer(x_rank_unique, x_rank_sort, function(a, b) sign(b - a))
  sign_y <- outer(y_rank_unique, y_rank_sort, function(a, b) sign(b - a))
  mean_sign_x <- (sign_x %*% t(sign_y)) / N
  
  mean_sign_x_indexed <- mean_sign_x[x_rank_indices, , drop = FALSE]
  
  all_exp <- mean_sign_x_indexed +
    2 * outer(G_x, rep(1, M)) +
    outer(rep(1, N), 2 * G_y) - 1
  
  # Vectorized g_1
  col_sums <- colSums(all_exp)
  g_1 <- numeric(N)
  for (i in 1:M) {
    start_idx <- pos_y[i] + 1
    end_idx <- pos_y[i + 1]
    g_1[start_idx:end_idx] <- col_sums[i] / N
  }
  
  # Vectorized g_2
  g_2 <- (all_exp %*% y_num) / N
  
  g_1 <- g_1 / 4
  g_2 <- as.vector(g_2) / 4
  
  G_y_full <- (y_rank_sort - 0.5) / N
  
  k_p <- 4 * (g_1 + g_2 + G_x * G_y_full - G_y_full - G_x) + 1 - rho
  
  # Reverse sorting
  rev_order <- order(ix_y)
  k_p[rev_order]
}

#' Prepare kernel_ties_optim2 arguments from ranks
#' @keywords internal
#' @noRd
prepare_original_kernel_args <- function(y_rank, x_rank) {
  y_rank_unique <- sort(unique(y_rank))
  y_num <- as.numeric(table(factor(y_rank, levels = y_rank_unique)))
  pos_y <- c(0, cumsum(y_num))
  ix_y <- order(y_rank)
  y_rank_sort <- y_rank[ix_y]
  x_rank_sort <- x_rank[ix_y]
  x_rank_unique <- sort(unique(x_rank_sort))
  
  list(x_rank_sort = x_rank_sort, y_rank_sort = y_rank_sort,
       y_rank_unique = y_rank_unique, y_num = y_num, pos_y = pos_y,
       ix_y = ix_y, x_rank_unique = x_rank_unique)
}

#' Compute kernel using original method
#' @keywords internal
#' @noRd
compute_kp_original <- function(y_rank, x_rank, rho) {
  args <- prepare_original_kernel_args(y_rank, x_rank)
  kernel_ties_optim2(args$x_rank_sort, args$y_rank_sort, rho,
                     args$y_rank_unique, args$y_num, args$pos_y,
                     args$ix_y, args$x_rank_unique)
}

#' @keywords internal
#' @noRd
kfn_original <- function(y_rank, x_rank, rho) {
  compute_kp_original(y_rank, x_rank, rho)
}



#' AGC variance, univariate IID (original kernel)
#' @keywords internal
#' @noRd
Sigma_agc <- function(y_rank, x_rank) {
  agc_sigma_univariate_iid(y_rank, x_rank, kfn_original)
}

#' AGC variance, univariate HAC (original kernel)
#' @keywords internal
#' @noRd
Sigma_agc_ts <- function(y_rank, x_rank) {
  agc_sigma_univariate_hac(y_rank, x_rank, kfn_original)
}

#' AGC covariance, multivariate IID (original kernel)
#' @keywords internal
#' @noRd
Sigma_agc_multivariate <- function(y_rank, xarray_ranks) {
  agc_sigma_multivariate_iid(y_rank, xarray_ranks, kfn_original)
}

#' AGC covariance, multivariate HAC (original kernel)
#' @keywords internal
#' @noRd
Sigma_agc_multivariate_ts <- function(y_rank, xarray_ranks) {
  agc_sigma_multivariate_hac(y_rank, xarray_ranks, kfn_original)
}

