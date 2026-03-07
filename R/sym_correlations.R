#' Compute Spearman rho with no tie correction
#' @param y_rank numeric vector of ranks
#' @param x_rank numeric vector of ranks
#' @return List with \code{rho} and \code{agc}.
#' @keywords internal
#' @noRd
comp_spearman_rho_a <- function(y_rank, x_rank) {
  N <- length(y_rank)
  mean_rank <- (N + 1) / 2
  rho_val <- (12 / (N^3)) * sum((x_rank - mean_rank) * (y_rank - mean_rank))
  return(rho_val)
}


#' Compute Spearman rho with tie correction
#' @param y_rank numeric vector of ranks
#' @param x_rank numeric vector of ranks
#' @return List with \code{rho} and \code{agc}.
#' @keywords internal
#' @noRd
comp_spearman_rho_b <- function(y_rank, x_rank) {
  rho_val_b <- comp_spearman_rho_a(y_rank, x_rank) / sqrt(comp_spearman_rho_a(x_rank, x_rank)*comp_spearman_rho_a(y_rank, y_rank) )
  return(rho_val_b)
}


#' Compute Pearson correlation
#' @param x numeric vector
#' @param y numeric vector
#' @return Pearson correlation coefficient.
#' @keywords internal
#' @noRd
comp_pearson_rho <- function(x, y) {
  n <- length(x)
  mx <- mean(x)
  my <- mean(y)
  sum((x - mx) * (y - my)) / sqrt(sum((x - mx)^2) * sum((y - my)^2))
}


#' IID independence variance for tau-a
#' @keywords internal
#' @noRd
ind_variance_tau_a_iid <- function(X, Y) {
  N <- length(Y)
  var_y_rank <- sum((rank(Y) - mean(rank(Y)))^2) / N
  var_x_rank <- sum((rank(X) - mean(rank(X)))^2) / N
  zeta_3Y <- 1 - (12 / N^2) * var_y_rank
  zeta_3X <- 1 - (12 / N^2) * var_x_rank
  (4 / 9) * (1 - zeta_3X) * (1 - zeta_3Y)
}


#' Multivariate IID independence covariance for tau-a
#' @keywords internal
#' @noRd
ind_covariance_tau_a_iid <- function(X, Y) {
  N <- length(Y)
  m <- ncol(X)
  
  y_rank <- rank(Y, ties.method = "average")
  var_y_rank <- sum((y_rank - mean(y_rank))^2) / N
  zeta_3Y <- 1 - (12 / N^2) * var_y_rank
  
  x_ranks <- matrix(0, nrow = N, ncol = m)
  zeta_3X <- numeric(m)
  for (k in seq_len(m)) {
    x_ranks[, k] <- rank(X[, k], ties.method = "average")
    var_x_rank <- sum((x_ranks[, k] - mean(x_ranks[, k]))^2) / N
    zeta_3X[k] <- 1 - (12 / N^2) * var_x_rank
  }
  
  Sigma_ind <- matrix(0, nrow = m, ncol = m)
  for (k in seq_len(m)) {
    for (l in k:m) {
      if (k == l) {
        Sigma_ind[k, k] <- (4 / 9) * (1 - zeta_3X[k]) * (1 - zeta_3Y)
      } else {
        x_grade_k <- (x_ranks[, k] - 0.5) / N - 0.5
        x_grade_l <- (x_ranks[, l] - 0.5) / N - 0.5
        rho_kl <- 12 * mean(x_grade_k * x_grade_l)
        Sigma_ind[k, l] <- (4 / 9) * rho_kl * (1 - zeta_3Y)
        Sigma_ind[l, k] <- Sigma_ind[k, l]
      }
    }
  }
  
  Sigma_ind
}

#' Compute multivariate tau-a with covariance matrix
#' @param X Numeric matrix (n x m).
#' @param Y Numeric outcome vector.
#' @param IID Logical; if TRUE use IID variance, otherwise HAC.
#' @param version Character; "v1" or "v2".
#' @return List with \code{tau_a_vector}, \code{Sigma}, \code{Sigma_ind}.
#' @keywords internal
#' @noRd
compute_tau_a_multivariate_variance <- function(X, Y, IID = TRUE, version = "v2") {
  X <- ensure_matrix(X)
  n <- length(Y)
  m <- ncol(X)
  
  K_tau_fn <- if (version == "v1") K_tau_vec_v1 else K_tau_vec_v2
  
  tau_vector <- numeric(m)
  K_tau_matrix <- matrix(0, nrow = n, ncol = m)
  
  for (k in seq_len(m)) {
    tau_result <- compute_kendall(X[, k], Y)
    tau_vector[k] <- tau_result$expectation
    K_tau_matrix[, k] <- K_tau_fn(X[, k], Y, tau_vector[k])
  }
  
  if (IID) {
    Sigma <- 4 * (t(K_tau_matrix) %*% K_tau_matrix) / n
  } else {
    Sigma <- 4 * (t(K_tau_matrix) %*% K_tau_matrix) / n +
      4 * hac_correction_multivariate(K_tau_matrix)
  }
  
  Sigma_ind <- ind_covariance_tau_a_iid(X, Y)
  
  list(tau_a_vector = tau_vector, Sigma = Sigma, Sigma_ind = Sigma_ind)
}

#' HAC independence variance for tau-a
#' @keywords internal
#' @noRd
ind_variance_tau_a_hac <- function(X, Y) {
  N <- length(Y)
  b <- floor(2 * N^(1/3))
  h_vec <- 1:(N - 1)
  w <- pmax(1 - abs(h_vec) / (b + 1), 0)
  
  x_grade_centered <- (rank(X, ties.method = "average") - 0.5) / N - 0.5
  y_grade_centered <- (rank(Y, ties.method = "average") - 0.5) / N - 0.5
  
  x_autoc <- stats::acf(x_grade_centered, plot = FALSE, type = "covariance",
                        demean = FALSE, lag.max = N - 1)$acf
  y_autoc <- stats::acf(y_grade_centered, plot = FALSE, type = "covariance",
                        demean = FALSE, lag.max = N - 1)$acf
  
  64 * sum(x_autoc[1] * y_autoc[1], 2 * (w * x_autoc[-1] * y_autoc[-1]))
}

#' Compute tau-a with variance
#' @param X Numeric predictor vector.
#' @param Y Numeric outcome vector.
#' @param IID Logical; if TRUE use IID variance, otherwise HAC.
#' @param version Character; "v2" for Fenwick kernel.
#' @return List with \code{tau_a}, \code{var}, \code{var_ind}.
#' @keywords internal
#' @noRd
compute_tau_a_variance <- function(X, Y, IID = TRUE, version = "v2") {
  tau_result <- compute_kendall(X, Y)
  tau_XY <- tau_result$expectation
  
  K_tau_fn <- if (version == "v1") K_tau_vec_v1 else K_tau_vec_v2
  K_tau_values <- K_tau_fn(X, Y, tau_XY)
  
  var_iid <- 4 * mean(K_tau_values^2)
  
  if (IID) {
    var_est <- var_iid
    var_ind <- ind_variance_tau_a_iid(X, Y)
  } else {
    var_est <- var_iid + 4 * hac_correction_univariate(K_tau_values)
    var_ind <- ind_variance_tau_a_hac(X, Y)
  }
  
  list(tau_a = tau_XY, var = var_est, var_ind = var_ind)
}



#' IID independence variance for rho_a
#' @keywords internal
#' @noRd
ind_variance_rho_a_iid <- function(x_rank, N, zeta_3Y) {
  var_x_rank <- sum((x_rank - mean(x_rank))^2) / N
  zeta_3X <- 1 - (12 / N^2) * var_x_rank
  (1 - zeta_3X) * (1 - zeta_3Y)
}

#' HAC independence variance for rho_a
#' @keywords internal
#' @noRd
ind_variance_rho_a_hac <- function(x_rank, y_rank, N, b) {
  h_vec <- 1:(N - 1)
  w <- pmax(1 - abs(h_vec) / (b + 1), 0)
  
  x_grade_centered <- (x_rank - 0.5) / N - 0.5
  y_grade_centered <- (y_rank - 0.5) / N - 0.5
  
  x_autoc <- stats::acf(x_grade_centered, plot = FALSE, type = "covariance",
                        demean = FALSE, lag.max = N - 1)$acf
  y_autoc <- stats::acf(y_grade_centered, plot = FALSE, type = "covariance",
                        demean = FALSE, lag.max = N - 1)$acf
  
  144 * sum(x_autoc[1] * y_autoc[1], 2 * (w * x_autoc[-1] * y_autoc[-1]))
}

#' Compute rho_a with variance
#' @keywords internal
#' @noRd
compute_rho_a_variance <- function(x_rank, y_rank, IID = TRUE) {
  N <- length(y_rank)
  
  rho_val <- comp_spearman_rho_a(y_rank, x_rank)
  
  kernel_fn <- if (is_binary(y_rank)) kfn_binary else kfn_v2
  k_p <- kernel_fn(y_rank, x_rank, rho_val)
  
  var_iid <- 9 * mean(k_p^2)
  
  if (IID) {
    pre <- agc_y_preamble(y_rank)
    var_est <- var_iid
    var_ind <- ind_variance_rho_a_iid(x_rank, N, pre$zeta_3Y)
  } else {
    b <- floor(2 * N^(1/3))
    k_p_autoc <- stats::acf(k_p, plot = FALSE, type = "covariance",
                            demean = FALSE, lag.max = N - 1)$acf
    h_vec <- 1:(N - 1)
    w <- pmax(1 - abs(h_vec) / (b + 1), 0)
    var_hac <- 9 * (2 * sum(w * k_p_autoc[-1]))
    var_est <- var_iid + var_hac
    var_ind <- ind_variance_rho_a_hac(x_rank, y_rank, N, b)
  }
  
  list(rho_a = rho_val, var = var_est, var_ind = var_ind)
}

#' Multivariate independence covariance for rho_a (IID)
#' @keywords internal
#' @noRd
ind_covariance_rho_a_iid <- function(xarray_ranks, N, zeta_3Y) {
  k <- nrow(xarray_ranks)
  
  zeta_3X <- numeric(k)
  for (j in 1:k) {
    var_x_rank <- sum((xarray_ranks[j, ] - mean(xarray_ranks[j, ]))^2) / N
    zeta_3X[j] <- 1 - (12 / N^2) * var_x_rank
  }
  
  Sigma_ind <- matrix(0, nrow = k, ncol = k)
  for (j in 1:k) {
    for (l in j:k) {
      if (j == l) {
        Sigma_ind[j, j] <- (1 - zeta_3X[j]) * (1 - zeta_3Y)
      } else {
        x_grade_j <- (xarray_ranks[j, ] - 0.5) / N - 0.5
        x_grade_l <- (xarray_ranks[l, ] - 0.5) / N - 0.5
        rho_jl <- 12 * mean(x_grade_j * x_grade_l)
        Sigma_ind[j, l] <- rho_jl * (1 - zeta_3Y)
        Sigma_ind[l, j] <- Sigma_ind[j, l]
      }
    }
  }
  
  Sigma_ind
}

#' Compute multivariate rho_a with covariance matrix
#' @keywords internal
#' @noRd
compute_rho_a_multivariate_variance <- function(y_rank, xarray_ranks, IID = TRUE) {
  k <- nrow(xarray_ranks)
  N <- length(y_rank)
  
  kernel_fn <- if (is_binary(y_rank)) kfn_binary else kfn_v2
  
  rho_vector <- numeric(k)
  kps <- matrix(0, nrow = k, ncol = N)
  
  for (j in 1:k) {
    rho_vector[j] <- comp_spearman_rho_a(y_rank, xarray_ranks[j, ])
    kps[j, ] <- kernel_fn(y_rank, xarray_ranks[j, ], rho_vector[j])
  }
  
  Sigma_iid <- 9 * (kps %*% t(kps)) / N
  
  if (IID) {
    Sigma <- Sigma_iid
  } else {
    b <- floor(2 * N^(1/3))
    Sigma_hac <- matrix(0, nrow = k, ncol = k)
    for (h in 1:b) {
      omega <- 1 - h / (b + 1)
      K_lag <- t(kps[, 1:(N - h), drop = FALSE])
      K_lead <- t(kps[, (h + 1):N, drop = FALSE])
      autocov_h <- (t(K_lag) %*% K_lead + t(K_lead) %*% K_lag) / N
      Sigma_hac <- Sigma_hac + omega * autocov_h
    }
    Sigma <- Sigma_iid + 9 * Sigma_hac
  }
  
  pre <- agc_y_preamble(y_rank)
  if (IID) {
    Sigma_ind <- ind_covariance_rho_a_iid(xarray_ranks, N, pre$zeta_3Y)
  } else {
    b <- floor(2 * N^(1/3))
    Sigma_ind <- ind_covariance_rho_a_hac(xarray_ranks, y_rank, N, b)
  }
  
  list(rho_a_vector = rho_vector, Sigma = Sigma, Sigma_ind = Sigma_ind)
}

#' Multivariate independence covariance for rho_a (HAC)
#' @keywords internal
#' @noRd
ind_covariance_rho_a_hac <- function(xarray_ranks, y_rank, N, b) {
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
      
      Sigma_ind[j, l] <- 144 * hac_sum
      if (j != l) Sigma_ind[l, j] <- Sigma_ind[j, l]
    }
  }
  
  Sigma_ind
}