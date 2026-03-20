# Spearman-family Symmetric Correlations: rho_a, rho_b
# Point estimates, influence functions, and variance (IID / HAC)
# Used by acor() and acor.test()


#' Compute Spearman rho with no tie correction
#' @param y_rank numeric vector of ranks
#' @param x_rank numeric vector of ranks
#' @return Spearman correlation without tie correction
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
#' @return Spearman correlation with tie correction
#' @keywords internal
#' @noRd
comp_spearman_rho_b <- function(y_rank, x_rank) {
  rho_val_b <- comp_spearman_rho_a(y_rank, x_rank) / sqrt(comp_spearman_rho_a(x_rank, x_rank)*comp_spearman_rho_a(y_rank, y_rank) )
  return(rho_val_b)
}


# ============================================================================
# rho_a
# ============================================================================

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
  x_grade <- (x_rank - 0.5) / N - 0.5
  y_grade <- (y_rank - 0.5) / N - 0.5
  144 * ind_lrv_univariate(x_grade, y_grade, N, b)
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

#' Multivariate independence covariance for rho_a (HAC)
#' @keywords internal
#' @noRd
ind_covariance_rho_a_hac <- function(xarray_ranks, y_rank, N, b) {
  k <- nrow(xarray_ranks)

  x_grades <- matrix(0, nrow = k, ncol = N)
  for (j in 1:k) {
    x_grades[j, ] <- (xarray_ranks[j, ] - 0.5) / N - 0.5
  }
  y_grade <- (y_rank - 0.5) / N - 0.5

  144 * ind_lrv_multivariate(x_grades, y_grade, N, b, x_by_row = TRUE)
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
  
  pre <- agc_y_preamble(y_rank)
  
  if (IID) {
    Sigma <- Sigma_iid
    Sigma_ind <- ind_covariance_rho_a_iid(xarray_ranks, N, pre$zeta_3Y)
  } else {
    b <- floor(2 * N^(1/3))
    Sigma_hac <- matrix(0, nrow = k, ncol = k)
    max_lag <- min(b, N - 1)
    for (h in seq_len(max_lag)) {
      omega <- 1 - h / (b + 1)
      K_lag <- t(kps[, 1:(N - h), drop = FALSE])
      K_lead <- t(kps[, (h + 1):N, drop = FALSE])
      autocov_h <- (t(K_lag) %*% K_lead + t(K_lead) %*% K_lag) / N
      Sigma_hac <- Sigma_hac + omega * autocov_h
    }
    Sigma <- Sigma_iid + 9 * Sigma_hac
    Sigma_ind <- ind_covariance_rho_a_hac(xarray_ranks, y_rank, N, b)
  }
  
  list(rho_a_vector = rho_vector, Sigma = Sigma, Sigma_ind = Sigma_ind)
}


# ============================================================================
# rho_b
# ============================================================================

#' IID independence variance for rho_b
#' @keywords internal
#' @noRd
ind_variance_rho_b_iid <- function() {
  1
}

#' HAC independence variance for rho_b
#' @keywords internal
#' @noRd
ind_variance_rho_b_hac <- function(X, Y) {
  n <- length(X)
  b <- floor(2 * n^(1 / 3))
  h <- 1:(n - 1)
  w <- pmax(1 - abs(h) / (b + 1), 0)

  x_grade <- (rank(X, ties.method = "average") - 0.5) / n - 0.5
  y_grade <- (rank(Y, ties.method = "average") - 0.5) / n - 0.5

  x_autoc <- n / (n - 1) *
    stats::acf(x_grade, plot = FALSE, type = "covariance",
               demean = FALSE, lag.max = n - 1)$acf / stats::var(x_grade)
  y_autoc <- n / (n - 1) *
    stats::acf(y_grade, plot = FALSE, type = "covariance",
               demean = FALSE, lag.max = n - 1)$acf / stats::var(y_grade)

  sum(x_autoc[1] * y_autoc[1], 2 * (w * x_autoc[-1] * y_autoc[-1]))
}

#' rho_b tie-normalization preamble for one rank vector
#' @keywords internal
#' @noRd
rho_b_rank_preamble <- function(rank_vec) {
  N <- length(rank_vec)
  var_rank <- sum((rank_vec - mean(rank_vec))^2) / N
  zeta_3 <- 1 - (12 / N^2) * var_rank
  k_zeta <- prob_y(rank_vec)^2 - zeta_3

  list(
    zeta_3 = zeta_3,
    rho_self = 1 - zeta_3,
    k_self = -k_zeta
  )
}

#' Build rho_b kernel components
#' @keywords internal
#' @noRd
rho_b_kernel_components <- function(X, Y) {
  x_rank <- rank(X, ties.method = "average")
  y_rank <- rank(Y, ties.method = "average")

  rho <- comp_spearman_rho_a(y_rank, x_rank)
  rho_b <- comp_spearman_rho_b(y_rank, x_rank)
  pre_x <- rho_b_rank_preamble(x_rank)
  pre_y <- rho_b_rank_preamble(y_rank)

  kernel_fn <- if (is_binary(y_rank)) kfn_binary else kfn_v2

  list(
    rho_b = rho_b,
    rho = rho,
    rho_x = pre_x$rho_self,
    rho_y = pre_y$rho_self,
    k_rho = kernel_fn(y_rank, x_rank, rho),
    k_rho_x = pre_x$k_self,
    k_rho_y = pre_y$k_self
  )
}

#' IID covariance matrix for rho_b delta method
#' @keywords internal
#' @noRd
rho_b_iid_covariance <- function(k_rho, k_rho_x, k_rho_y) {
  Sigma <- matrix(0, nrow = 3, ncol = 3)
  Sigma[1, 1] <- 9 * mean(k_rho^2)
  Sigma[2, 2] <- 4 * mean(k_rho_x^2)
  Sigma[3, 3] <- 4 * mean(k_rho_y^2)
  Sigma[1, 2] <- Sigma[2, 1] <- 9 * mean(k_rho * k_rho_x)
  Sigma[1, 3] <- Sigma[3, 1] <- 9 * mean(k_rho * k_rho_y)
  Sigma[2, 3] <- Sigma[3, 2] <- 9 * mean(k_rho_x * k_rho_y)
  Sigma
}

#' HAC covariance matrix for rho_b delta method
#' @keywords internal
#' @noRd
rho_b_hac_covariance <- function(k_rho, k_rho_x, k_rho_y, N, b) {
  h <- 1:(N - 1)
  w <- pmax(1 - abs(h) / (b + 1), 0)
  cross_weights <- c(sort(w), 1, w)

  k_rho_autoc <- stats::acf(k_rho, plot = FALSE, type = "covariance",
                            demean = FALSE, lag.max = N - 1)$acf
  k_rho_x_autoc <- stats::acf(k_rho_x, plot = FALSE, type = "covariance",
                              demean = FALSE, lag.max = N - 1)$acf
  k_rho_y_autoc <- stats::acf(k_rho_y, plot = FALSE, type = "covariance",
                              demean = FALSE, lag.max = N - 1)$acf
  k_rho_rho_x <- stats::ccf(k_rho, k_rho_x, plot = FALSE, type = "covariance",
                            demean = TRUE, lag.max = N - 1)$acf
  k_rho_rho_y <- stats::ccf(k_rho, k_rho_y, plot = FALSE, type = "covariance",
                            demean = TRUE, lag.max = N - 1)$acf
  k_rho_xrho_y <- stats::ccf(k_rho_x, k_rho_y, plot = FALSE,
                             type = "covariance", demean = TRUE,
                             lag.max = N - 1)$acf

  Sigma <- matrix(0, nrow = 3, ncol = 3)
  Sigma[1, 1] <- 9 * sum(k_rho_autoc[1], 2 * (w * k_rho_autoc[-1]))
  Sigma[2, 2] <- 9 * sum(k_rho_x_autoc[1], 2 * (w * k_rho_x_autoc[-1]))
  Sigma[3, 3] <- 9 * sum(k_rho_y_autoc[1], 2 * (w * k_rho_y_autoc[-1]))
  Sigma[1, 2] <- Sigma[2, 1] <- 9 * sum(cross_weights * k_rho_rho_x)
  Sigma[1, 3] <- Sigma[3, 1] <- 9 * sum(cross_weights * k_rho_rho_y)
  Sigma[2, 3] <- Sigma[3, 2] <- 9 * sum(cross_weights * k_rho_xrho_y)
  Sigma
}

#' Delta-method gradient for rho_b = rho / sqrt(rho_x * rho_y)
#' @keywords internal
#' @noRd
rho_b_gradient <- function(rho, rho_x, rho_y) {
  c(
    1 / sqrt(rho_x * rho_y),
    -rho / (2 * rho_x^(3 / 2) * sqrt(rho_y)),
    -rho / (2 * rho_y^(3 / 2) * sqrt(rho_x))
  )
}

#' Compute rho_b with variance
#' @keywords internal
#' @noRd
compute_rho_b_variance <- function(X, Y, IID = TRUE) {
  n <- length(X)
  comps <- rho_b_kernel_components(X, Y)
  grad <- rho_b_gradient(comps$rho, comps$rho_x, comps$rho_y)

  if (IID) {
    Sigma <- rho_b_iid_covariance(comps$k_rho, comps$k_rho_x, comps$k_rho_y)
    var_est <- as.numeric(t(grad) %*% Sigma %*% grad)
    var_ind <- ind_variance_rho_b_iid()
  } else {
    b <- floor(2 * n^(1 / 3))
    Sigma <- rho_b_hac_covariance(comps$k_rho, comps$k_rho_x, comps$k_rho_y, n, b)
    var_est <- as.numeric(t(grad) %*% Sigma %*% grad)
    var_ind <- ind_variance_rho_b_hac(X, Y)
  }

  list(rho_b = comps$rho_b, var = var_est, var_ind = var_ind)
}
