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


#' Compute Pearson correlation
#' @param x numeric vector
#' @param y numeric vector
#' @return Pearson correlation coefficient.
#' @keywords internal
#' @noRd
comp_pearson_rho <- function(x, y) {
  mx <- mean(x)
  my <- mean(y)
  sum((x - mx) * (y - my)) / sqrt(sum((x - mx)^2) * sum((y - my)^2))
}


#' Compute Pearson influence function values
#'
#' @param X Numeric vector.
#' @param Y Numeric vector.
#' @param r Pearson correlation.
#' @return Numeric vector of influence function values.
#' @keywords internal
#' @noRd
pearson_kernel <- function(X, Y, r) {
  # Use n-denominator scaling so the influence-function HAC variance matches
  # the delta-method formulation used by the legacy Pearson implementation.
  X_centered <- X - mean(X)
  Y_centered <- Y - mean(Y)
  X_std <- X_centered / sqrt(mean(X_centered^2))
  Y_std <- Y_centered / sqrt(mean(Y_centered^2))
  X_std * Y_std - 0.5 * r * (X_std^2 + Y_std^2)
}


#' IID independence variance for Pearson
#'
#' @return Scalar independence variance.
#' @keywords internal
#' @noRd
ind_variance_pearson_iid <- function() {
  1
}


#' HAC independence variance for Pearson
#'
#' @param X Numeric vector.
#' @param Y Numeric vector.
#' @return Scalar independence long-run variance.
#' @keywords internal
#' @noRd
ind_variance_pearson_hac <- function(X, Y) {
  n <- length(X)
  b <- floor(2 * n^(1 / 3))
  h_vec <- 1:(n - 1)
  w <- pmax(1 - abs(h_vec) / (b + 1), 0)
  
  x_autoc <- stats::acf(X, plot = FALSE, type = "correlation",
                        demean = TRUE, lag.max = n - 1)$acf
  y_autoc <- stats::acf(Y, plot = FALSE, type = "correlation",
                        demean = TRUE, lag.max = n - 1)$acf
  
  sum(x_autoc[1] * y_autoc[1], 2 * (w * x_autoc[-1] * y_autoc[-1]))
}


#' Compute Pearson correlation with variance
#'
#' @param X Numeric predictor vector.
#' @param Y Numeric outcome vector.
#' @param IID Logical; if TRUE use IID variance, otherwise HAC.
#' @return List with \code{estimate}, \code{var}, and \code{var_ind}.
#' @keywords internal
#' @noRd
compute_pearson_variance <- function(X, Y, IID = TRUE) {
  r <- comp_pearson_rho(X, Y)
  K_r <- pearson_kernel(X, Y, r)
  var_iid <- mean(K_r^2)
  
  if (IID) {
    var_est <- var_iid
    var_ind <- ind_variance_pearson_iid()
  } else {
    n <- length(X)
    b <- floor(2 * n^(1 / 3))
    k_autoc <- stats::acf(K_r, plot = FALSE, type = "covariance",
                          demean = FALSE, lag.max = n - 1)$acf
    h_vec <- 1:(n - 1)
    w <- pmax(1 - abs(h_vec) / (b + 1), 0)
    var_hac <- 2 * sum(w * k_autoc[-1])
    var_est <- var_iid + var_hac
    var_ind <- ind_variance_pearson_hac(X, Y)
  }
  
  list(estimate = r, var = var_est, var_ind = var_ind)
}


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
  
  # FIX: dispatch based on IID flag
  if (IID) {
    Sigma_ind <- ind_covariance_tau_a_iid(X, Y)
  } else {
    Sigma_ind <- ind_covariance_tau_a_hac(X, Y)
  }
  
  list(tau_a_vector = tau_vector, Sigma = Sigma, Sigma_ind = Sigma_ind)
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



