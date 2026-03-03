# Asymmetric Grade Correlation (AGC) Functions
# Refactored to eliminate duplicated variance/HAC logic via shared helpers.
# Three kernel implementations: original (kernel_ties_optim2), v2 (Fenwick tree),
# binary (specialized for binary Y).
# All Sigma_* functions return both main variance and independence variance.


# ============================================================================
# CORE HELPERS
# ============================================================================

#' Compute probability for each value in y based on frequencies
#' @param y numeric vector
#' @return Named numeric vector of probabilities
#' @keywords internal
#' @noRd
prob_y <- function(y) {
  freq_table <- table(y)
  probs <- as.vector(freq_table) / length(y)
  names(probs) <- names(freq_table)
  probs[as.character(y)]
}

#' Compute rho (Spearman-style) and AGC
#' @param y_rank numeric vector of ranks
#' @param x_rank numeric vector of ranks
#' @return List with \code{rho} and \code{agc}.
#' @keywords internal
#' @noRd
comp_rho_agc <- function(y_rank, x_rank) {
  N <- length(y_rank)
  mean_rank <- (N + 1) / 2
  var_y <- sum((y_rank - mean(y_rank))^2) * (1 / (N - 1))
  rho_val <- (12 / (N^3)) * sum((x_rank - mean_rank) * (y_rank - mean_rank))
  agc_val <- cov(y_rank, x_rank) / var_y
  
  list(rho = rho_val, agc = agc_val)
}

#' Compute AGC for a single predictor (no variance)
#' @param y_rank numeric ranks for Y
#' @param x_rank numeric ranks for X
#' @return AGC value
#' @keywords internal
#' @noRd
compute_agc <- function(y_rank, x_rank) {
  N <- length(y_rank)
  var_y <- sum((y_rank - mean(y_rank))^2) * (1 / (N - 1))
  cov(y_rank, x_rank) / var_y
}

#' Compute AGC for multiple predictors (no variance)
#' @param y_rank numeric ranks for Y
#' @param xarray_ranks matrix of ranks (predictors x n)
#' @return Vector of AGC values
#' @keywords internal
#' @noRd
compute_agc_multivariate <- function(y_rank, xarray_ranks) {
  k <- nrow(xarray_ranks)
  agcs <- numeric(k)
  for (j in 1:k) {
    agcs[j] <- compute_agc(y_rank, xarray_ranks[j, ])
  }
  agcs
}


# ============================================================================
# SHARED VARIANCE COMPUTATION HELPERS
# ============================================================================

#' Compute Y-side preamble shared across all AGC variance functions
#'
#' @param y_rank Numeric vector of average ranks for Y.
#' @return List with \code{N}, \code{zeta_3Y}, \code{k_zeta}, \code{sigma_zeta}.
#' @keywords internal
#' @noRd
agc_y_preamble <- function(y_rank) {
  N <- length(y_rank)
  var_y_rank <- sum((y_rank - mean(y_rank))^2) / N  # biased variance (n denominator)
  zeta_3Y <- 1 - (12 / N^2) * var_y_rank
  k_zeta <- prob_y(y_rank)^2 - zeta_3Y
  sigma_zeta <- 9 * mean(k_zeta^2)
  
  list(N = N, zeta_3Y = zeta_3Y, k_zeta = k_zeta, sigma_zeta = sigma_zeta)
}

#' Compute IID variance for a single predictor from kernel values
#'
#' @param k_p Numeric vector of kernel influence function values.
#' @param k_zeta Numeric vector of zeta kernel values.
#' @param rho Scalar Spearman-style correlation.
#' @param zeta_3Y Scalar.
#' @param sigma_zeta Scalar.
#' @return Scalar variance estimate.
#' @keywords internal
#' @noRd
agc_iid_variance <- function(k_p, k_zeta, rho, zeta_3Y, sigma_zeta) {
  factor_val <- 1 / ((1 - zeta_3Y)^2)
  sigma_rho <- 9 * mean(k_p^2)
  sigma_pz <- 9 * mean(k_p * k_zeta)
  
  factor_val * (sigma_rho + (2 * rho * sigma_pz) / (1 - zeta_3Y) +
                  (rho^2 * sigma_zeta) / ((1 - zeta_3Y)^2))
}

#' Compute HAC correction for C_zeta (autocovariance of k_zeta)
#'
#' @param k_zeta Numeric vector.
#' @param b Integer bandwidth.
#' @param N Integer sample size.
#' @return Scalar HAC correction.
#' @keywords internal
#' @noRd
agc_hac_C_zeta <- function(k_zeta, b, N) {
  C_zeta_hac <- 0
  for (h in 1:b) {
    omega <- 1 - h / (b + 1)
    c_k <- sum(k_zeta[1:(N - h)] * k_zeta[(h + 1):N])
    C_zeta_hac <- C_zeta_hac + omega * (2 / N) * c_k
  }
  C_zeta_hac
}

#' Compute univariate HAC correction terms (A_spear, B_both)
#'
#' @param k_p Numeric vector of kernel values.
#' @param k_zeta Numeric vector of zeta kernel values.
#' @param b Integer bandwidth.
#' @param N Integer sample size.
#' @return List with \code{A_spear} and \code{B_both}.
#' @keywords internal
#' @noRd
agc_hac_univariate_terms <- function(k_p, k_zeta, b, N) {
  A_spear <- 0
  B_both <- 0
  
  for (h in 1:b) {
    omega <- 1 - h / (b + 1)
    K1_lag  <- k_p[1:(N - h)]
    K1_lead <- k_p[(h + 1):N]
    K0_lag  <- k_zeta[1:(N - h)]
    K0_lead <- k_zeta[(h + 1):N]
    
    a_k <- sum(K1_lag * K1_lead)
    b_k <- sum(K1_lag * K0_lead + K0_lag * K1_lead)
    
    A_spear <- A_spear + omega * (2 / N) * a_k
    B_both  <- B_both  + omega * (1 / N) * b_k
  }
  
  list(A_spear = A_spear, B_both = B_both)
}

#' Compute off-diagonal HAC correction terms for two predictors
#'
#' @param k_p1 Numeric vector of kernel values for predictor 1.
#' @param k_p2 Numeric vector of kernel values for predictor 2.
#' @param k_zeta Numeric vector of zeta kernel values.
#' @param b Integer bandwidth.
#' @param N Integer sample size.
#' @return List with \code{A2_spear}, \code{B_one}, \code{B_two}.
#' @keywords internal
#' @noRd
agc_hac_offdiag_terms <- function(k_p1, k_p2, k_zeta, b, N) {
  A2_spear <- 0
  B_one <- 0
  B_two <- 0
  
  for (h in 1:b) {
    omega <- 1 - h / (b + 1)
    kp1_lag  <- k_p1[1:(N - h)]
    kp1_lead <- k_p1[(h + 1):N]
    kp2_lag  <- k_p2[1:(N - h)]
    kp2_lead <- k_p2[(h + 1):N]
    K0_lag   <- k_zeta[1:(N - h)]
    K0_lead  <- k_zeta[(h + 1):N]
    
    a_k  <- sum(kp1_lag * kp2_lead + kp1_lead * kp2_lag)
    b_k1 <- sum(kp1_lag * K0_lead + kp1_lead * K0_lag)
    b_k2 <- sum(kp2_lag * K0_lead + kp2_lead * K0_lag)
    
    A2_spear <- A2_spear + omega * (1 / N) * a_k
    B_one    <- B_one    + omega * (1 / N) * b_k1
    B_two    <- B_two    + omega * (1 / N) * b_k2
  }
  
  list(A2_spear = A2_spear, B_one = B_one, B_two = B_two)
}

#' Assemble univariate HAC variance from pre-computed terms
#'
#' @param rho Scalar.
#' @param zeta_3Y Scalar.
#' @param sigma_zeta Scalar.
#' @param A_spear Scalar from HAC correction.
#' @param B_both Scalar from HAC correction.
#' @param C_zeta_hac Scalar from HAC correction.
#' @return Scalar HAC variance correction.
#' @keywords internal
#' @noRd
agc_assemble_hac_variance <- function(rho, zeta_3Y, sigma_zeta,
                                      A_spear, B_both, C_zeta_hac) {
  factor_val <- 1 / ((1 - zeta_3Y)^2)
  factor_val * 9 * (A_spear + (2 * rho * B_both) / (1 - zeta_3Y) +
                  (rho^2 * C_zeta_hac) / ((1 - zeta_3Y)^2))
}

#' Assemble off-diagonal HAC covariance from pre-computed terms
#'
#' @param rho_j Scalar rho for predictor j.
#' @param rho_i Scalar rho for predictor i.
#' @param zeta_3Y Scalar.
#' @param sigma_zeta Scalar.
#' @param A2_spear Scalar.
#' @param B_one Scalar.
#' @param B_two Scalar.
#' @param C_zeta_hac Scalar.
#' @return Scalar HAC covariance correction.
#' @keywords internal
#' @noRd
agc_assemble_hac_covariance <- function(rho_j, rho_i, zeta_3Y, sigma_zeta,
                                        A2_spear, B_one, B_two, C_zeta_hac) {
  factor_val <- 1 / ((1 - zeta_3Y)^2)
  factor_val * 9 * (A2_spear + (rho_j * B_one) / (1 - zeta_3Y) +
                  (rho_i * B_two) / (1 - zeta_3Y) +
                  (rho_j * rho_i * C_zeta_hac) / ((1 - zeta_3Y)^2))
}


# ============================================================================
# INDEPENDENCE VARIANCE HELPERS
# ============================================================================

#' Univariate IID independence variance for AGC
#'
#' @param x_rank Numeric vector of ranks for X.
#' @param N Integer sample size.
#' @param zeta_3Y Scalar.
#' @return Scalar independence variance.
#' @keywords internal
#' @noRd
ind_variance_agc_iid <- function(x_rank, N, zeta_3Y) {
  var_x_rank <- sum((x_rank - mean(x_rank))^2) / N
  zeta_3X <- 1 - (12 / N^2) * var_x_rank
  (1 - zeta_3X) / (1 - zeta_3Y)
}

#' Univariate HAC independence variance for AGC
#'
#' @param x_rank Numeric vector of ranks for X.
#' @param y_rank Numeric vector of ranks for Y.
#' @param N Integer sample size.
#' @param zeta_3Y Scalar.
#' @param b Integer bandwidth.
#' @return Scalar independence variance with HAC correction.
#' @keywords internal
#' @noRd
ind_variance_agc_hac <- function(x_rank, y_rank, N, zeta_3Y, b) {
 h_vec <- 1:(N - 1)
 w <- pmax(1 - abs(h_vec) / (b + 1), 0)

 x_grade_centered <- (x_rank - 0.5) / N - 0.5
 y_grade_centered <- (y_rank - 0.5) / N - 0.5

 x_autoc <- stats::acf(x_grade_centered, plot = FALSE, type = "covariance",
                       demean = FALSE, lag.max = N - 1)$acf
 y_autoc <- stats::acf(y_grade_centered, plot = FALSE, type = "covariance",
                       demean = FALSE, lag.max = N - 1)$acf

 rho_ind_lrv <- 144 * sum(x_autoc[1] * y_autoc[1],
                          2 * (w * x_autoc[-1] * y_autoc[-1]))
 rho_ind_lrv / ((1 - zeta_3Y)^2)
}


#' Multivariate IID independence covariance for AGC
#'
#' @param xarray_ranks Matrix of ranks (k x N).
#' @param N Integer sample size.
#' @param zeta_3Y Scalar.
#' @return k x k independence covariance matrix.
#' @keywords internal
#' @noRd
ind_covariance_agc_iid <- function(xarray_ranks, N, zeta_3Y) {
  k <- nrow(xarray_ranks)
  
  zeta_3X <- numeric(k)
  for (j in 1:k) {
    var_x_rank <- sum((xarray_ranks[j, ] - mean(xarray_ranks[j, ]))^2) / N
    zeta_3X[j] <- 1 - (12 / N^2) * var_x_rank
  }
  
  Sigma_ind <- matrix(0, nrow = k, ncol = k)
  for (j in 1:k) {
    Sigma_ind[j, j] <- (1 - zeta_3X[j]) / (1 - zeta_3Y)
    if (j < k) {
      for (i in (j + 1):k) {
        x_grade_j <- (xarray_ranks[j, ] - 0.5) / N - 0.5
        x_grade_i <- (xarray_ranks[i, ] - 0.5) / N - 0.5
        rho_ji <- 12 * mean(x_grade_j * x_grade_i)
        Sigma_ind[j, i] <- Sigma_ind[i, j] <- rho_ji / (1 - zeta_3Y)
      }
    }
  }
  Sigma_ind
}

#' Multivariate HAC independence covariance for AGC
#'
#' @param xarray_ranks Matrix of ranks (k x N).
#' @param y_rank Numeric vector of ranks for Y.
#' @param N Integer sample size.
#' @param zeta_3Y Scalar.
#' @param b Integer bandwidth.
#' @return k x k independence covariance matrix with HAC correction.
#' @keywords internal
#' @noRd
ind_covariance_agc_hac <- function(xarray_ranks, y_rank, N, zeta_3Y, b) {
  k <- nrow(xarray_ranks)

  h_vec <- 1:(N - 1)
  w <- pmax(1 - abs(h_vec) / (b + 1), 0)

  x_grades_centered <- matrix(0, nrow = k, ncol = N)
  for (j in 1:k) {
    x_grades_centered[j, ] <- (xarray_ranks[j, ] - 0.5) / N - 0.5
  }
  y_grade_centered <- (y_rank - 0.5) / N - 0.5

  y_autoc <- stats::acf(y_grade_centered, plot = FALSE, type = "covariance",
                        demean = FALSE, lag.max = N - 1)$acf

  Sigma_ind <- matrix(0, nrow = k, ncol = k)

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

      Sigma_ind[j, l] <- 144 * hac_sum / ((1 - zeta_3Y)^2)
      if (j != l) Sigma_ind[l, j] <- Sigma_ind[j, l]
    }
  }
  Sigma_ind
}




# ============================================================================
# IID COVARIANCE MATRIX BUILDER (shared by multivariate IID functions)
# ============================================================================

#' Build IID covariance matrix from kernel values
#'
#' @param kps Matrix (k x N) of kernel values per predictor.
#' @param k_zeta Numeric vector.
#' @param rhos Numeric vector of rho values per predictor.
#' @param zeta_3Y Scalar.
#' @param sigma_zeta Scalar.
#' @return k x k covariance matrix.
#' @keywords internal
#' @noRd
agc_build_iid_covariance <- function(kps, k_zeta, rhos, zeta_3Y, sigma_zeta) {
  k <- nrow(kps)
  factor_val <- 1 / ((1 - zeta_3Y)^2)
  Sigma <- matrix(0, nrow = k, ncol = k)
  
  # Pre-compute sigma_pz for each predictor (reused in off-diagonal)
  sigma_pz_vec <- numeric(k)
  for (j in 1:k) {
    sigma_pz_vec[j] <- 9 * mean(kps[j, ] * k_zeta)
  }
  
  for (j in 1:k) {
    k_p <- kps[j, ]
    sigma_rho <- 9 * mean(k_p^2)
    
    Sigma[j, j] <- factor_val * (
      sigma_rho + (2 * rhos[j] * sigma_pz_vec[j]) / (1 - zeta_3Y) +
        (rhos[j]^2 * sigma_zeta) / ((1 - zeta_3Y)^2))
    
    if (j < k) {
      for (i in (j + 1):k) {
        sigma_rho2 <- 9 * mean(k_p * kps[i, ])
        
        cov_agc <- factor_val * (
          sigma_rho2 + (rhos[j] * sigma_pz_vec[j]) / (1 - zeta_3Y) +
            (rhos[i] * sigma_pz_vec[i]) / (1 - zeta_3Y) +
            (rhos[j] * rhos[i] * sigma_zeta) / ((1 - zeta_3Y)^2))
        Sigma[j, i] <- Sigma[i, j] <- cov_agc
      }
    }
  }
  Sigma
}


# ============================================================================
# HAC COVARIANCE MATRIX BUILDER (shared by multivariate HAC functions)
# ============================================================================

#' Build HAC covariance matrix from kernel values
#'
#' @param kps Matrix (k x N) of kernel values per predictor.
#' @param k_zeta Numeric vector.
#' @param rhos Numeric vector of rho values per predictor.
#' @param zeta_3Y Scalar.
#' @param sigma_zeta Scalar.
#' @param N Integer sample size.
#' @param b Integer bandwidth.
#' @return k x k covariance matrix (IID + HAC).
#' @keywords internal
#' @noRd
agc_build_hac_covariance <- function(kps, k_zeta, rhos, zeta_3Y, sigma_zeta, N, b) {
  k <- nrow(kps)
  factor_val <- 1 / ((1 - zeta_3Y)^2)
  
  # Pre-compute C_zeta_hac (same for all predictors)
  C_zeta_hac <- agc_hac_C_zeta(k_zeta, b, N)
  
  # Pre-compute sigma_pz for each predictor
  sigma_pz_vec <- numeric(k)
  for (j in 1:k) {
    sigma_pz_vec[j] <- 9 * mean(kps[j, ] * k_zeta)
  }
  
  Sigma <- matrix(0, nrow = k, ncol = k)
  
  for (j in 1:k) {
    k_p <- kps[j, ]
    sigma_rho <- 9 * mean(k_p^2)
    
    # IID component
    var_iid <- factor_val * (
      sigma_rho + (2 * rhos[j] * sigma_pz_vec[j]) / (1 - zeta_3Y) +
        (rhos[j]^2 * sigma_zeta) / ((1 - zeta_3Y)^2))
    
    # HAC correction for diagonal
    hac_terms <- agc_hac_univariate_terms(k_p, k_zeta, b, N)
    var_hac <- agc_assemble_hac_variance(rhos[j], zeta_3Y, sigma_zeta,
                                         hac_terms$A_spear, hac_terms$B_both,
                                         C_zeta_hac)
    Sigma[j, j] <- var_iid + var_hac
    
    # Off-diagonal elements
    if (j < k) {
      for (i in (j + 1):k) {
        k_p2 <- kps[i, ]
        sigma_rho2 <- 9 * mean(k_p * k_p2)
        
        # IID component for off-diagonal
        cov_iid <- factor_val * (
          sigma_rho2 + (rhos[j] * sigma_pz_vec[j]) / (1 - zeta_3Y) +
            (rhos[i] * sigma_pz_vec[i]) / (1 - zeta_3Y) +
            (rhos[j] * rhos[i] * sigma_zeta) / ((1 - zeta_3Y)^2))
        
        # HAC correction for off-diagonal
        offdiag <- agc_hac_offdiag_terms(k_p, k_p2, k_zeta, b, N)
        cov_hac <- agc_assemble_hac_covariance(rhos[j], rhos[i], zeta_3Y,
                                               sigma_zeta, offdiag$A2_spear,
                                               offdiag$B_one, offdiag$B_two,
                                               C_zeta_hac)
        Sigma[j, i] <- Sigma[i, j] <- cov_iid + cov_hac
      }
    }
  }
  Sigma
}


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


# ============================================================================
# KERNEL: Binary Y specialization -- O(n log n)
# ============================================================================

#' Compute AGC kernel influence function values for binary Y
#'
#' @param x_rank numeric vector of average ranks for X
#' @param y_rank numeric vector of average ranks for Y (exactly 2 unique values)
#' @param rho Spearman-style correlation
#' @return numeric vector k_p of length N
#' @keywords internal
#' @noRd
kernel_agc_binary <- function(x_rank, y_rank, rho) {
  N <- length(x_rank)
  
  y_unique <- sort(unique(y_rank))
  y0 <- y_unique[1]
  y1 <- y_unique[2]
  
  mask0 <- y_rank == y0
  mask1 <- y_rank == y1
  n0 <- sum(mask0)
  n1 <- sum(mask1)
  
  G_x <- (x_rank - 0.5) / N
  G_y0 <- (y0 - 0.5) / N
  G_y1 <- (y1 - 0.5) / N
  G_y <- numeric(N)
  G_y[mask0] <- G_y0
  G_y[mask1] <- G_y1
  
  # Sort x_ranks within each Y group
  x_rank_g0 <- sort(x_rank[mask0])
  x_rank_g1 <- sort(x_rank[mask1])
  
  # Count sign sums using findInterval
  pos_g1_strict <- findInterval(x_rank, x_rank_g1, left.open = TRUE)
  count_less_g1 <- pos_g1_strict
  count_leq_g1 <- findInterval(x_rank, x_rank_g1)
  count_greater_g1 <- n1 - count_leq_g1
  sum_sign_g1 <- count_greater_g1 - count_less_g1
  
  pos_g0_strict <- findInterval(x_rank, x_rank_g0, left.open = TRUE)
  count_less_g0 <- pos_g0_strict
  count_leq_g0 <- findInterval(x_rank, x_rank_g0)
  count_greater_g0 <- n0 - count_leq_g0
  sum_sign_g0 <- count_greater_g0 - count_less_g0
  
  mean_sign_x_0 <- sum_sign_g1 / N
  mean_sign_x_1 <- -sum_sign_g0 / N
  
  all_exp_0 <- mean_sign_x_0 + 2 * G_x + 2 * G_y0 - 1
  all_exp_1 <- mean_sign_x_1 + 2 * G_x + 2 * G_y1 - 1
  
  col_sum_0 <- sum(all_exp_0)
  col_sum_1 <- sum(all_exp_1)
  
  g_1 <- numeric(N)
  g_1[mask0] <- col_sum_0 / (4 * N)
  g_1[mask1] <- col_sum_1 / (4 * N)
  
  g_2 <- (all_exp_0 * n0 + all_exp_1 * n1) / (4 * N)
  
  4 * (g_1 + g_2 + G_x * G_y - G_y - G_x) + 1 - rho
}


# ============================================================================
# GENERIC SIGMA FUNCTIONS (dispatch to kernel, then shared variance logic)
# ============================================================================
# These are the main entry points. Each computes kernels using the appropriate
# method, then delegates to the shared variance/HAC helpers.
# All return list(agc, var, var_ind) for univariate
#   or list(agc_vector, Sigma, Sigma_ind) for multivariate.

# ---------- Univariate IID ----------

#' @keywords internal
#' @noRd
agc_sigma_univariate_iid <- function(y_rank, x_rank, kernel_fn) {
  pre <- agc_y_preamble(y_rank)
  result <- comp_rho_agc(y_rank, x_rank)
  
  k_p <- kernel_fn(y_rank, x_rank, result$rho)
  
  var_agc <- agc_iid_variance(k_p, pre$k_zeta, result$rho,
                              pre$zeta_3Y, pre$sigma_zeta)
  var_ind <- ind_variance_agc_iid(x_rank, pre$N, pre$zeta_3Y)
  
  list(agc = result$agc, var = var_agc, var_ind = var_ind)
}

# ---------- Univariate HAC ----------

#' @keywords internal
#' @noRd
agc_sigma_univariate_hac <- function(y_rank, x_rank, kernel_fn) {
  pre <- agc_y_preamble(y_rank)
  result <- comp_rho_agc(y_rank, x_rank)
  
  k_p <- kernel_fn(y_rank, x_rank, result$rho)
  
  var_iid <- agc_iid_variance(k_p, pre$k_zeta, result$rho,
                              pre$zeta_3Y, pre$sigma_zeta)
  
  b <- floor(2 * pre$N^(1 / 3))
  C_zeta_hac <- agc_hac_C_zeta(pre$k_zeta, b, pre$N)
  hac_terms <- agc_hac_univariate_terms(k_p, pre$k_zeta, b, pre$N)
  var_hac <- agc_assemble_hac_variance(result$rho, pre$zeta_3Y,
                                       pre$sigma_zeta,
                                       hac_terms$A_spear, hac_terms$B_both,
                                       C_zeta_hac)
  
  var_ind <- ind_variance_agc_hac(x_rank, y_rank, pre$N, pre$zeta_3Y, b)
  
  list(agc = result$agc, var = var_iid + var_hac, var_ind = var_ind)
}

# ---------- Multivariate IID ----------

#' @keywords internal
#' @noRd
agc_sigma_multivariate_iid <- function(y_rank, xarray_ranks, kernel_fn) {
  pre <- agc_y_preamble(y_rank)
  k <- nrow(xarray_ranks)
  
  rhos <- numeric(k)
  agcs <- numeric(k)
  kps <- matrix(0, nrow = k, ncol = pre$N)
  
  for (j in 1:k) {
    result <- comp_rho_agc(y_rank, xarray_ranks[j, ])
    rhos[j] <- result$rho
    agcs[j] <- result$agc
    kps[j, ] <- kernel_fn(y_rank, xarray_ranks[j, ], rhos[j])
  }
  
  Sigma <- agc_build_iid_covariance(kps, pre$k_zeta, rhos,
                                    pre$zeta_3Y, pre$sigma_zeta)
  Sigma_ind <- ind_covariance_agc_iid(xarray_ranks, pre$N, pre$zeta_3Y)
  
  list(agc_vector = agcs, Sigma = Sigma, Sigma_ind = Sigma_ind)
}

# ---------- Multivariate HAC ----------

#' @keywords internal
#' @noRd
agc_sigma_multivariate_hac <- function(y_rank, xarray_ranks, kernel_fn) {
  pre <- agc_y_preamble(y_rank)
  k <- nrow(xarray_ranks)
  
  rhos <- numeric(k)
  agcs <- numeric(k)
  kps <- matrix(0, nrow = k, ncol = pre$N)
  
  for (j in 1:k) {
    result <- comp_rho_agc(y_rank, xarray_ranks[j, ])
    rhos[j] <- result$rho
    agcs[j] <- result$agc
    kps[j, ] <- kernel_fn(y_rank, xarray_ranks[j, ], rhos[j])
  }
  
  b <- floor(2 * pre$N^(1 / 3))
  Sigma <- agc_build_hac_covariance(kps, pre$k_zeta, rhos,
                                    pre$zeta_3Y, pre$sigma_zeta,
                                    pre$N, b)
  Sigma_ind <- ind_covariance_agc_hac(xarray_ranks, y_rank,
                                      pre$N, pre$zeta_3Y, b)
  
  if (k == 1) {
    return(list(agc = agcs[1], var = Sigma[1, 1], var_ind = Sigma_ind[1, 1]))
  }
  
  list(agc_vector = agcs, Sigma = Sigma, Sigma_ind = Sigma_ind)
}


# ============================================================================
# KERNEL WRAPPER FUNCTIONS
# ============================================================================
# Each kernel type needs a uniform interface: fn(y_rank, x_rank, rho) -> k_p
# The original kernel requires extra preparation, so we wrap it.

#' @keywords internal
#' @noRd
kfn_original <- function(y_rank, x_rank, rho) {
  compute_kp_original(y_rank, x_rank, rho)
}

#' @keywords internal
#' @noRd
kfn_v2 <- function(y_rank, x_rank, rho) {
  kernel_agc_v2_cpp(x_rank, y_rank, rho)
}

#' @keywords internal
#' @noRd
kfn_binary <- function(y_rank, x_rank, rho) {
  kernel_agc_binary(x_rank, y_rank, rho)
}


# ============================================================================
# PUBLIC SIGMA FUNCTIONS (original kernel)
# ============================================================================

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


# ============================================================================
# PUBLIC SIGMA FUNCTIONS (v2 / Fenwick tree kernel)
# ============================================================================

#' AGC variance, univariate IID (Fenwick tree kernel)
#' @keywords internal
#' @noRd
Sigma_agc_v2 <- function(y_rank, x_rank) {
  agc_sigma_univariate_iid(y_rank, x_rank, kfn_v2)
}

#' AGC variance, univariate HAC (Fenwick tree kernel)
#' @keywords internal
#' @noRd
Sigma_agc_ts_v2 <- function(y_rank, x_rank) {
  agc_sigma_univariate_hac(y_rank, x_rank, kfn_v2)
}

#' AGC covariance, multivariate IID (Fenwick tree kernel)
#' @keywords internal
#' @noRd
Sigma_agc_multivariate_v2 <- function(y_rank, xarray_ranks) {
  agc_sigma_multivariate_iid(y_rank, xarray_ranks, kfn_v2)
}

#' AGC covariance, multivariate HAC (Fenwick tree kernel)
#' @keywords internal
#' @noRd
Sigma_agc_multivariate_ts_v2 <- function(y_rank, xarray_ranks) {
  agc_sigma_multivariate_hac(y_rank, xarray_ranks, kfn_v2)
}


# ============================================================================
# PUBLIC SIGMA FUNCTIONS (binary Y kernel)
# ============================================================================

#' AGC variance, univariate IID (binary Y kernel)
#' @keywords internal
#' @noRd
Sigma_agc_binary <- function(y_rank, x_rank) {
  agc_sigma_univariate_iid(y_rank, x_rank, kfn_binary)
}

#' AGC variance, univariate HAC (binary Y kernel)
#' @keywords internal
#' @noRd
Sigma_agc_ts_binary <- function(y_rank, x_rank) {
  agc_sigma_univariate_hac(y_rank, x_rank, kfn_binary)
}

#' AGC covariance, multivariate IID (binary Y kernel)
#' @keywords internal
#' @noRd
Sigma_agc_multivariate_binary <- function(y_rank, xarray_ranks) {
  agc_sigma_multivariate_iid(y_rank, xarray_ranks, kfn_binary)
}

#' AGC covariance, multivariate HAC (binary Y kernel)
#' @keywords internal
#' @noRd
Sigma_agc_multivariate_ts_binary <- function(y_rank, xarray_ranks) {
  agc_sigma_multivariate_hac(y_rank, xarray_ranks, kfn_binary)
}