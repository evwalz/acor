# Kendall-family Symmetric Correlations: tau_a, tau_b, gamma
# Point estimates, influence functions, and variance (IID / HAC)
# Used by acor() and acor.test()
#
# Shared helpers: prob_xy() is also used by rho_b (spearman.R)


#' Observation-wise joint empirical probabilities
#' @keywords internal
#' @noRd
prob_xy <- function(X, Y) {
  freq_table <- table(X, Y)
  probs <- as.vector(freq_table) / length(X)
  pair_keys <- paste(X, Y, sep = "\r")
  level_grid <- expand.grid(dimnames(freq_table), KEEP.OUT.ATTRS = FALSE,
                            stringsAsFactors = FALSE)
  names(probs) <- apply(level_grid, 1, paste, collapse = "\r")
  probs[pair_keys]
}


# ============================================================================
# tau_a
# ============================================================================

#' IID independence variance for tau-a
#' @keywords internal
#' @noRd
ind_variance_tau_a_iid <- function(X, Y) {
  N <- length(Y)
  var_y_rank <- sum((rank(Y,ties.method = "average") - mean(rank(Y, ties.method = "average")))^2) / N
  var_x_rank <- sum((rank(X, ties.method = "average") - mean(rank(X, ties.method = "average")))^2) / N
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


#' Multivariate HAC independence covariance for tau-a
#'
#' @param X Numeric matrix (n x m).
#' @param Y Numeric outcome vector.
#' @return m x m independence covariance matrix with HAC correction.
#' @keywords internal
#' @noRd
ind_covariance_tau_a_hac <- function(X, Y) {
  N <- length(Y)
  m <- ncol(X)
  b <- floor(2 * N^(1/3))
  
  x_grades <- matrix(0, nrow = N, ncol = m)
  for (k in seq_len(m)) {
    x_grades[, k] <- (rank(X[, k], ties.method = "average") - 0.5) / N - 0.5
  }
  y_grade <- (rank(Y, ties.method = "average") - 0.5) / N - 0.5
  
  64 * ind_lrv_multivariate(x_grades, y_grade, N, b, x_by_row = FALSE)
}

#' HAC independence variance for tau-a
#' @keywords internal
#' @noRd
ind_variance_tau_a_hac <- function(X, Y) {
  N <- length(Y)
  b <- floor(2 * N^(1/3))
  x_grade <- (rank(X, ties.method = "average") - 0.5) / N - 0.5
  y_grade <- (rank(Y, ties.method = "average") - 0.5) / N - 0.5
  64 * ind_lrv_univariate(x_grade, y_grade, N, b)
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
  
  if (IID) {
    Sigma_ind <- ind_covariance_tau_a_iid(X, Y)
  } else {
    Sigma_ind <- ind_covariance_tau_a_hac(X, Y)
  }
  
  list(tau_a_vector = tau_vector, Sigma = Sigma, Sigma_ind = Sigma_ind)
}


# ============================================================================
# tau_b
# ============================================================================

#' tau_b tie-normalization preamble
#' @keywords internal
#' @noRd
tau_b_tie_preamble <- function(V) {
  N <- length(V)
  probs <- as.vector(table(V)) / length(V)
  tie_prob2 <- sum(probs^2)
  tie_prob3 <- sum(probs^3)
  tau_self <- (N / (N - 1)) * (1 - tie_prob2)
  denom_ind <- 1 - tie_prob2
  k_self <- 1 - prob_y(V) - tau_self

  list(
    tau_self = tau_self,
    denom_ind = denom_ind,
    tie_prob3 = tie_prob3,
    k_self = k_self
  )
}

#' Univariate HAC long-run variance helper for tau_b components
#' @keywords internal
#' @noRd
tau_b_hac_lrv <- function(k_vec, N, b) {
  h <- 1:(N - 1)
  w <- pmax(1 - abs(h) / (b + 1), 0)
  k_autoc <- stats::acf(k_vec, plot = FALSE, type = "covariance",
                        demean = FALSE, lag.max = N - 1)$acf
  4 * sum(k_autoc[1], 2 * (w * k_autoc[-1]))
}

#' Univariate HAC long-run cross-covariance helper for tau_b components
#' @keywords internal
#' @noRd
tau_b_hac_cross_lrv <- function(k1, k2, N, b) {
  h <- 1:(N - 1)
  w <- pmax(1 - abs(h) / (b + 1), 0)
  k_cross <- stats::ccf(k1, k2, plot = FALSE, type = "covariance",
                        lag.max = N - 1)$acf
  4 * sum(c(sort(w), 1, w) * k_cross)
}

#' Compute tau_b with variance
#' @keywords internal
#' @noRd
compute_tau_b_variance <- function(X, Y, IID = TRUE, version = "v2") {
  N <- length(Y)
  tau_result <- compute_kendall(X, Y)
  tau <- tau_result$expectation
  tau_b <- kendall_tau_b(X, Y)

  pre_x <- tau_b_tie_preamble(X)
  pre_y <- tau_b_tie_preamble(Y)

  K_tau_fn <- if (version == "v1") K_tau_vec_v1 else K_tau_vec_v2
  k_tau <- K_tau_fn(X, Y, tau)
  k_tau_x <- pre_x$k_self
  k_tau_y <- pre_y$k_self

  if (IID) {
    var_tau <- 4 * mean(k_tau^2)
    var_tau_x <- 4 * mean(k_tau_x^2)
    var_tau_y <- 4 * mean(k_tau_y^2)
    var_tautau_x <- 4 * mean(k_tau * k_tau_x)
    var_tautau_y <- 4 * mean(k_tau * k_tau_y)
    var_tau_xtau_y <- 4 * mean(k_tau_x * k_tau_y)

    var_est <- (
      var_tau -
        tau * (var_tautau_x / pre_x$tau_self - var_tautau_y / pre_y$tau_self) +
        tau^2 / 4 * (
          var_tau_x / pre_x$tau_self^2 +
            var_tau_y / pre_y$tau_self^2 +
            (2 * var_tau_xtau_y) / (pre_y$tau_self * pre_x$tau_self)
        )
    ) / (pre_x$tau_self * pre_y$tau_self)

    var_ind <- 4 / 9 * (1 - pre_x$tie_prob3) * (1 - pre_y$tie_prob3) /
      (pre_x$denom_ind * pre_y$denom_ind)
  } else {
    b <- floor(2 * N^(1 / 3))
    sigma_tau_sq <- tau_b_hac_lrv(k_tau, N, b)
    sigma_tau_x_sq <- tau_b_hac_lrv(k_tau_x, N, b)
    sigma_tau_y_sq <- tau_b_hac_lrv(k_tau_y, N, b)
    sigma_tautau_x <- tau_b_hac_cross_lrv(k_tau, k_tau_x, N, b)
    sigma_tautau_y <- tau_b_hac_cross_lrv(k_tau, k_tau_y, N, b)
    sigma_tau_xtau_y <- tau_b_hac_cross_lrv(k_tau_x, k_tau_y, N, b)

    var_est <- (
      sigma_tau_sq -
        tau * (sigma_tautau_x / pre_x$tau_self - sigma_tautau_y / pre_y$tau_self) +
        tau^2 / 4 * (
          sigma_tau_x_sq / pre_x$tau_self^2 +
            sigma_tau_y_sq / pre_y$tau_self^2 +
            (2 * sigma_tau_xtau_y) / (pre_x$tau_self * pre_y$tau_self)
        )
    ) / (pre_x$tau_self * pre_y$tau_self)

    var_ind <- ind_variance_tau_a_hac(X, Y) / (pre_x$denom_ind * pre_y$denom_ind)
  }

  list(tau_b = tau_b, var = var_est, var_ind = var_ind)
}


# ============================================================================
# Goodman-Kruskal gamma
# ============================================================================

#' Gamma tie-probability preamble
#' @keywords internal
#' @noRd
gamma_tie_preamble <- function(X, Y) {
  x_probs <- as.vector(table(X)) / length(X)
  y_probs <- as.vector(table(Y)) / length(Y)
  xy_probs <- as.vector(table(X, Y)) / length(X)

  x_tie_prob <- sum(x_probs^2)
  y_tie_prob <- sum(y_probs^2)
  xy_tie_prob <- sum(xy_probs^2)
  tie_prob <- x_tie_prob + y_tie_prob - xy_tie_prob

  list(
    x_tie_prob = x_tie_prob,
    y_tie_prob = y_tie_prob,
    tie_prob = tie_prob,
    k_tie = prob_y(X) + prob_y(Y) - prob_xy(X, Y) - tie_prob
  )
}

#' Compute gamma with variance
#' @keywords internal
#' @noRd
compute_gamma_variance <- function(X, Y, IID = TRUE, version = "v2") {
  N <- length(Y)
  tau_result <- compute_kendall(X, Y)
  tau <- tau_result$expectation
  gamma <- goodman_kruskal_gamma(X, Y)
  tie_pre <- gamma_tie_preamble(X, Y)
  gamma_scale <- tau / (1 - tie_pre$tie_prob)

  K_tau_fn <- if (version == "v1") K_tau_vec_v1 else K_tau_vec_v2
  k_tau <- K_tau_fn(X, Y, tau)
  k_nu <- tie_pre$k_tie

  if (IID) {
    var_tau <- 4 * mean(k_tau^2)
    var_nu <- 4 * mean(k_nu^2)
    var_taunu <- 4 * mean(k_tau * k_nu)

    var_est <- (var_tau + gamma_scale^2 * var_nu + 2 * gamma_scale * var_taunu) /
      (1 - tie_pre$tie_prob)^2

    x_probs <- as.vector(table(X)) / N
    y_probs <- as.vector(table(Y)) / N
    var_ind <- 4 / 9 * (1 - sum(x_probs^3)) * (1 - sum(y_probs^3)) /
      (1 - tie_pre$x_tie_prob)^2 / (1 - tie_pre$y_tie_prob)^2
  } else {
    b <- floor(2 * N^(1 / 3))
    sigma_tau_sq <- tau_b_hac_lrv(k_tau, N, b)
    sigma_nu_sq <- tau_b_hac_lrv(k_nu, N, b)
    sigma_taunu <- tau_b_hac_cross_lrv(k_tau, k_nu, N, b)

    var_est <- (sigma_tau_sq + gamma_scale^2 * sigma_nu_sq + 2 * gamma_scale * sigma_taunu) /
      (1 - tie_pre$tie_prob)^2

    var_ind <- ind_variance_tau_a_hac(X, Y) /
      (1 - tie_pre$x_tie_prob)^2 / (1 - tie_pre$y_tie_prob)^2
  }

  list(gamma = gamma, var = var_est, var_ind = var_ind)
}
