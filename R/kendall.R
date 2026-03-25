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


#' Delta-method gradient for tau_b = tau / sqrt(sx * sy)
#'
#' Matches \code{kendall_tau_b}: \code{sx}, \code{sy} are \code{tau_self} margins
#' (\eqn{\widehat{\tau}(X,X)}, \eqn{\widehat{\tau}(Y,Y)} in the pair-tie sense).
#'
#' @param tau Kendall tau-a expectation \eqn{\tau(X,Y)} (plug-in).
#' @param sx,sy Marginal \code{tau_self} values.
#' @return Length-3 vector \eqn{(\partial h/\partial\tau, \partial h/\partial s_x,
#'   \partial h/\partial s_y)} for \eqn{h = \tau/\sqrt{s_x s_y}}.
#' @keywords internal
#' @noRd
tau_b_gradient <- function(tau, sx, sy) {
  c(
    1 / sqrt(sx * sy),
    -tau / (2 * sx^(3 / 2) * sqrt(sy)),
    -tau / (2 * sy^(3 / 2) * sqrt(sx))
  )
}


#' IID asymptotic covariance for (tau, sx, sy) influence components
#' @keywords internal
#' @noRd
tau_b_iid_covariance <- function(k_tau, k_sx, k_sy) {
  S <- matrix(0, 3, 3)
  S[1, 1] <- 4 * mean(k_tau^2)
  S[2, 2] <- 4 * mean(k_sx^2)
  S[3, 3] <- 4 * mean(k_sy^2)
  S[1, 2] <- S[2, 1] <- 4 * mean(k_tau * k_sx)
  S[1, 3] <- S[3, 1] <- 4 * mean(k_tau * k_sy)
  S[2, 3] <- S[3, 2] <- 4 * mean(k_sx * k_sy)
  S
}


#' HAC long-run covariance matrix for (tau, sx, sy) components
#' @keywords internal
#' @noRd
tau_b_hac_covariance <- function(k_tau, k_sx, k_sy, N, b) {
  S <- matrix(0, 3, 3)
  S[1, 1] <- tau_b_hac_lrv(k_tau, N, b)
  S[2, 2] <- tau_b_hac_lrv(k_sx, N, b)
  S[3, 3] <- tau_b_hac_lrv(k_sy, N, b)
  S[1, 2] <- S[2, 1] <- tau_b_hac_cross_lrv(k_tau, k_sx, N, b)
  S[1, 3] <- S[3, 1] <- tau_b_hac_cross_lrv(k_tau, k_sy, N, b)
  S[2, 3] <- S[3, 2] <- tau_b_hac_cross_lrv(k_sx, k_sy, N, b)
  S
}


#' Univariate independence variance for tau_b (tie preambles only; no kernels)
#' @keywords internal
#' @noRd
ind_variance_tau_b_univariate <- function(X, Y, IID = TRUE) {
  pre_x <- tau_b_tie_preamble(X)
  pre_y <- tau_b_tie_preamble(Y)
  if (IID) {
    4 / 9 * (1 - pre_x$tie_prob3) * (1 - pre_y$tie_prob3) /
      (pre_x$denom_ind * pre_y$denom_ind)
  } else {
    ind_variance_tau_a_hac(X, Y) / (pre_x$denom_ind * pre_y$denom_ind)
  }
}


#' Compute tau_b with variance
#'
#' Point estimate is \code{kendall_tau_b} = \eqn{\tau(X,Y)/\sqrt{s_x s_y}} with
#' \code{sx, sy = tau_self} from \code{tau_b_tie_preamble}.  Main variance uses
#' the delta method for that ratio (same structure as \code{rho_b}).
#'
#' @keywords internal
#' @noRd
compute_tau_b_variance <- function(X, Y, IID = TRUE, version = "v2") {
  N <- length(Y)
  tau_result <- compute_kendall(X, Y)
  tau <- tau_result$expectation
  tau_b <- kendall_tau_b(X, Y)

  pre_x <- tau_b_tie_preamble(X)
  pre_y <- tau_b_tie_preamble(Y)
  sx <- pre_x$tau_self
  sy <- pre_y$tau_self

  K_tau_fn <- if (version == "v1") K_tau_vec_v1 else K_tau_vec_v2
  k_tau <- K_tau_fn(X, Y, tau)
  k_tau_x <- pre_x$k_self
  k_tau_y <- pre_y$k_self

  grad <- tau_b_gradient(tau, sx, sy)

  if (IID) {
    S <- tau_b_iid_covariance(k_tau, k_tau_x, k_tau_y)
    var_est <- as.numeric(t(grad) %*% S %*% grad)
  } else {
    b <- floor(2 * N^(1 / 3))
    S <- tau_b_hac_covariance(k_tau, k_tau_x, k_tau_y, N, b)
    var_est <- as.numeric(t(grad) %*% S %*% grad)
  }

  list(tau_b = tau_b, var = var_est,
       var_ind = ind_variance_tau_b_univariate(X, Y, IID))
}


#' Independence covariance for multivariate tau-b (IID)
#'
#' Under independence, off-diagonals scale with \eqn{\sqrt{\mathrm{denom}_{x,k}
#' \mathrm{denom}_{x,l}}\,\mathrm{denom}_y}; diagonals match univariate
#' \code{compute_tau_b_variance()$var_ind}.
#'
#' @param X Numeric matrix \code{n x m}.
#' @param Y Numeric vector length \code{n}.
#' @return \code{m x m} matrix.
#' @keywords internal
#' @noRd
ind_covariance_tau_b_iid <- function(X, Y) {
  X <- ensure_matrix(X)
  m <- ncol(X)
  denom_x <- numeric(m)
  pre_x_list <- vector("list", m)
  for (k in seq_len(m)) {
    pre_x_list[[k]] <- tau_b_tie_preamble(X[, k])
    denom_x[k] <- pre_x_list[[k]]$denom_ind
  }
  pre_y <- tau_b_tie_preamble(Y)
  denom_y <- pre_y$denom_ind
  sta <- ind_covariance_tau_a_iid(X, Y)
  Sigma <- matrix(0, nrow = m, ncol = m)
  # Diagonal: same closed form as univariate compute_tau_b_variance()$var_ind (IID)
  for (k in seq_len(m)) {
    px <- pre_x_list[[k]]
    Sigma[k, k] <- (4 / 9) * (1 - px$tie_prob3) * (1 - pre_y$tie_prob3) /
      (px$denom_ind * pre_y$denom_ind)
  }
  for (k in seq_len(m)) {
    for (l in seq_len(m)) {
      if (k != l) {
        Sigma[k, l] <- sta[k, l] / (sqrt(denom_x[k] * denom_x[l]) * denom_y)
      }
    }
  }
  Sigma
}


#' Independence covariance for multivariate tau-b (HAC)
#' @keywords internal
#' @noRd
ind_covariance_tau_b_hac <- function(X, Y) {
  X <- ensure_matrix(X)
  m <- ncol(X)
  denom_x <- numeric(m)
  for (k in seq_len(m)) {
    denom_x[k] <- tau_b_tie_preamble(X[, k])$denom_ind
  }
  denom_y <- tau_b_tie_preamble(Y)$denom_ind
  sta <- ind_covariance_tau_a_hac(X, Y)
  Sigma <- matrix(0, nrow = m, ncol = m)
  for (k in seq_len(m)) {
    Sigma[k, k] <- ind_variance_tau_a_hac(X[, k], Y) / (denom_x[k] * denom_y)
  }
  for (k in seq_len(m)) {
    for (l in seq_len(m)) {
      if (k != l) {
        Sigma[k, l] <- sta[k, l] / (sqrt(denom_x[k] * denom_x[l]) * denom_y)
      }
    }
  }
  Sigma
}


#' Multivariate Kendall tau-b with covariance matrix
#'
#' Point estimates are \code{kendall_tau_b} per column of \code{X}.  Main
#' covariance is \eqn{J \Sigma_\theta J^\top} with \eqn{\theta = (\tau^{(1)},
#' s_x^{(1)}, \ldots, \tau^{(m)}, s_x^{(m)}, s_y)} and one shared \eqn{s_y}
#' column; \eqn{\Sigma_\theta} from stacked influences (IID or HAC).  Same
#' sqrt-denominator delta method as \code{compute_tau_b_variance()}.
#'
#' @param X Numeric matrix \code{n x m}.
#' @param Y Numeric outcome vector.
#' @param IID Logical; IID or HAC for \code{Sigma}.
#' @param version Kernel \code{"v1"} or \code{"v2"}.
#' @return List with \code{tau_b_vector}, \code{Sigma}, \code{Sigma_ind}.
#' @keywords internal
#' @noRd
compute_tau_b_multivariate_variance <- function(X, Y, IID = TRUE, version = "v2") {
  X <- ensure_matrix(X)
  n <- nrow(X)
  m <- ncol(X)

  K_tau_fn <- if (version == "v1") K_tau_vec_v1 else K_tau_vec_v2

  pre_y <- tau_b_tie_preamble(Y)
  sy <- pre_y$tau_self
  ky <- pre_y$k_self

  tau_vec <- numeric(m)
  sx_vec <- numeric(m)
  Psi <- matrix(0, nrow = n, ncol = 2 * m + 1)

  for (k in seq_len(m)) {
    tau_result <- compute_kendall(X[, k], Y)
    tau_vec[k] <- tau_result$expectation
    pre_x <- tau_b_tie_preamble(X[, k])
    sx_vec[k] <- pre_x$tau_self
    Psi[, 2 * k - 1] <- 2 * K_tau_fn(X[, k], Y, tau_vec[k])
    Psi[, 2 * k] <- 2 * pre_x$k_self
  }
  Psi[, 2 * m + 1] <- 2 * ky

  if (IID) {
    Sigma_theta <- crossprod(Psi) / n
  } else {
    Sigma_theta <- crossprod(Psi) / n + hac_correction_multivariate(Psi)
  }

  J <- matrix(0, nrow = m, ncol = 2 * m + 1)
  for (k in seq_len(m)) {
    g <- tau_b_gradient(tau_vec[k], sx_vec[k], sy)
    J[k, 2 * (k - 1) + 1] <- g[1]
    J[k, 2 * k] <- g[2]
    J[k, 2 * m + 1] <- g[3]
  }

  Sigma <- J %*% Sigma_theta %*% t(J)

  tau_b_vector <- vapply(seq_len(m), function(k) kendall_tau_b(X[, k], Y), numeric(1))

  if (IID) {
    Sigma_ind <- ind_covariance_tau_b_iid(X, Y)
  } else {
    Sigma_ind <- ind_covariance_tau_b_hac(X, Y)
  }

  list(tau_b_vector = tau_b_vector, Sigma = Sigma, Sigma_ind = Sigma_ind)
}


# ============================================================================
# Goodman-Kruskal gamma
# ============================================================================
#
# gamma = tau_a(X,Y) / (1 - nu(X,Y)), where nu is the pair-based (without-
# replacement) proportion of pairs with (X_i - X_j)(Y_i - Y_j) = 0, matching
# goodman_kruskal_gamma()'s denominator (C + D) / binom(n,2).

#' Gamma tie-probability preamble (pair-based nu)
#'
#' Computes pair-based tie proportion \code{nu} and its Hoeffding projection
#' kernel \code{k_nu}.  With-replacement quantities \code{x_tie_prob} and
#' \code{y_tie_prob} (\eqn{\sum p_k^2}) are retained for independence variance.
#'
#' @param X Numeric predictor vector.
#' @param Y Numeric outcome vector.
#' @return List with \code{x_tie_prob}, \code{y_tie_prob}, \code{nu},
#'   \code{k_nu}.
#' @keywords internal
#' @noRd
gamma_tie_preamble <- function(X, Y) {
  N <- length(X)

  x_freq <- as.vector(table(X))
  y_freq <- as.vector(table(Y))
  xy_freq <- as.vector(table(X, Y))

  x_probs <- x_freq / N
  y_probs <- y_freq / N
  x_tie_prob <- sum(x_probs^2)
  y_tie_prob <- sum(y_probs^2)

  num_pairs <- N * (N - 1) / 2
  x_tied_pairs <- sum(x_freq * (x_freq - 1)) / 2
  y_tied_pairs <- sum(y_freq * (y_freq - 1)) / 2
  xy_tied_pairs <- sum(xy_freq * (xy_freq - 1)) / 2
  nu <- (x_tied_pairs + y_tied_pairs - xy_tied_pairs) / num_pairs

  freq_x_i <- prob_y(X) * N
  freq_y_i <- prob_y(Y) * N
  freq_xy_i <- prob_xy(X, Y) * N

  H_nu_i <- (freq_x_i - 1) + (freq_y_i - 1) - (freq_xy_i - 1)
  k_nu <- H_nu_i / (N - 1) - nu

  list(
    x_tie_prob = x_tie_prob,
    y_tie_prob = y_tie_prob,
    nu = nu,
    k_nu = k_nu
  )
}


#' Univariate independence variance for gamma (\code{gamma_tie_preamble} only)
#' @keywords internal
#' @noRd
ind_variance_gamma_univariate <- function(X, Y, IID = TRUE) {
  N <- length(Y)
  tie_pre <- gamma_tie_preamble(X, Y)
  if (IID) {
    x_probs <- as.vector(table(X)) / N
    y_probs <- as.vector(table(Y)) / N
    4 / 9 * (1 - sum(x_probs^3)) * (1 - sum(y_probs^3)) /
      (1 - tie_pre$x_tie_prob)^2 / (1 - tie_pre$y_tie_prob)^2
  } else {
    ind_variance_tau_a_hac(X, Y) /
      (1 - tie_pre$x_tie_prob)^2 / (1 - tie_pre$y_tie_prob)^2
  }
}


#' Multivariate independence covariance for gamma (IID)
#' @keywords internal
#' @noRd
ind_covariance_gamma_iid <- function(X, Y) {
  X <- ensure_matrix(X)
  N <- length(Y)
  m <- ncol(X)
  x_tie_probs <- numeric(m)
  y_probs <- as.vector(table(Y)) / N
  y_tie_prob <- sum(y_probs^2)
  for (k in seq_len(m)) {
    x_tie_probs[k] <- gamma_tie_preamble(X[, k], Y)$x_tie_prob
  }
  Sigma_ind_tau <- ind_covariance_tau_a_iid(X, Y)
  denom <- (1 - x_tie_probs) * (1 - y_tie_prob)
  Sigma_ind_tau / (denom %o% denom)
}


#' Multivariate independence covariance for gamma (HAC)
#' @keywords internal
#' @noRd
ind_covariance_gamma_hac <- function(X, Y) {
  X <- ensure_matrix(X)
  N <- length(Y)
  m <- ncol(X)
  x_tie_probs <- numeric(m)
  y_probs <- as.vector(table(Y)) / N
  y_tie_prob <- sum(y_probs^2)
  for (k in seq_len(m)) {
    x_tie_probs[k] <- gamma_tie_preamble(X[, k], Y)$x_tie_prob
  }
  Sigma_ind_tau <- ind_covariance_tau_a_hac(X, Y)
  denom <- (1 - x_tie_probs) * (1 - y_tie_prob)
  Sigma_ind_tau / (denom %o% denom)
}


#' Compute gamma with variance
#' @keywords internal
#' @noRd
compute_gamma_variance <- function(X, Y, IID = TRUE, version = "v2") {
  N <- length(Y)
  tau_result <- compute_kendall(X, Y)
  tau <- tau_result$expectation
  gamma_val <- goodman_kruskal_gamma(X, Y)
  tie_pre <- gamma_tie_preamble(X, Y)

  D <- 1 - tie_pre$nu
  gamma_scale <- tau / D

  K_tau_fn <- if (version == "v1") K_tau_vec_v1 else K_tau_vec_v2
  k_tau <- K_tau_fn(X, Y, tau)
  k_nu <- tie_pre$k_nu

  if (IID) {
    var_tau <- 4 * mean(k_tau^2)
    var_nu <- 4 * mean(k_nu^2)
    var_taunu <- 4 * mean(k_tau * k_nu)

    var_est <- (var_tau + gamma_scale^2 * var_nu + 2 * gamma_scale * var_taunu) /
      D^2
  } else {
    b <- floor(2 * N^(1 / 3))
    sigma_tau_sq <- tau_b_hac_lrv(k_tau, N, b)
    sigma_nu_sq <- tau_b_hac_lrv(k_nu, N, b)
    sigma_taunu <- tau_b_hac_cross_lrv(k_tau, k_nu, N, b)

    var_est <- (sigma_tau_sq + gamma_scale^2 * sigma_nu_sq +
                  2 * gamma_scale * sigma_taunu) / D^2
  }

  list(gamma = gamma_val, var = var_est,
       var_ind = ind_variance_gamma_univariate(X, Y, IID))
}


#' Multivariate Goodman--Kruskal gamma with covariance matrix
#'
#' Same estimand and variance construction as \code{compute_gamma_variance()},
#' extended to several predictors (stacked influence functions, then IID or HAC
#' covariance).
#'
#' @param X Numeric matrix (n x m).
#' @param Y Numeric outcome vector.
#' @param IID Logical; if TRUE use IID variance, otherwise HAC.
#' @param version Character; "v1" or "v2".
#' @return List with \code{gamma_vector}, \code{Sigma}, \code{Sigma_ind}.
#' @keywords internal
#' @noRd
compute_gamma_multivariate_variance <- function(X, Y, IID = TRUE, version = "v2") {
  X <- ensure_matrix(X)
  N <- length(Y)
  m <- ncol(X)

  K_tau_fn <- if (version == "v1") K_tau_vec_v1 else K_tau_vec_v2

  gamma_vector <- numeric(m)
  IC_matrix <- matrix(0, nrow = N, ncol = m)

  for (k in seq_len(m)) {
    tau_result <- compute_kendall(X[, k], Y)
    tau_k <- tau_result$expectation
    gamma_vector[k] <- goodman_kruskal_gamma(X[, k], Y)

    tie_pre_k <- gamma_tie_preamble(X[, k], Y)
    D_k <- 1 - tie_pre_k$nu
    gamma_scale_k <- tau_k / D_k

    k_tau_k <- K_tau_fn(X[, k], Y, tau_k)
    k_nu_k <- tie_pre_k$k_nu

    IC_matrix[, k] <- 2 * (k_tau_k + gamma_scale_k * k_nu_k) / D_k
  }

  if (IID) {
    Sigma <- (t(IC_matrix) %*% IC_matrix) / N
    Sigma_ind <- ind_covariance_gamma_iid(X, Y)
  } else {
    Sigma <- (t(IC_matrix) %*% IC_matrix) / N +
      hac_correction_multivariate(IC_matrix)
    Sigma_ind <- ind_covariance_gamma_hac(X, Y)
  }

  list(gamma_vector = gamma_vector, Sigma = Sigma, Sigma_ind = Sigma_ind)
}
