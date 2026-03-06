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