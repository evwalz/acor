# Pearson Correlation: point estimate, influence function, and variance
# Used by acor() and acor.test() for method = "pearson"


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


#' IID independence covariance for multivariate Pearson
#'
#' Under independence of X and Y, the influence function reduces to
#' X_k_std * Y_std (since r_k = 0), so Cov(r_k, r_l) = E[X_k_std * X_l_std].
#' This is simply the correlation matrix of the predictors.
#'
#' @param X Numeric matrix (n x m).
#' @return m x m independence covariance matrix.
#' @keywords internal
#' @noRd
ind_covariance_pearson_iid <- function(X) {
  n <- nrow(X)
  m <- ncol(X)
  Sigma_ind <- matrix(0, nrow = m, ncol = m)
  X_std <- scale(X, center = TRUE, scale = FALSE)
  X_std <- X_std / rep(sqrt(colMeans(X_std^2)), each = n)
  for (k in seq_len(m)) {
    for (l in k:m) {
      Sigma_ind[k, l] <- mean(X_std[, k] * X_std[, l])
      Sigma_ind[l, k] <- Sigma_ind[k, l]
    }
  }
  Sigma_ind
}


#' HAC independence covariance for multivariate Pearson
#'
#' Under independence, IC_k(t) = X_k_std(t) * Y_std(t), and since X and Y are
#' independent processes the lag-h cross-covariance factors:
#' Cov(IC_k(t), IC_l(t+h)) = Cov(X_k_std(t), X_l_std(t+h)) * Cov(Y_std(t), Y_std(t+h)).
#'
#' @param X Numeric matrix (n x m).
#' @param Y Numeric vector.
#' @return m x m independence covariance matrix with HAC correction.
#' @keywords internal
#' @noRd
ind_covariance_pearson_hac <- function(X, Y) {
  n <- nrow(X)
  m <- ncol(X)
  b <- floor(2 * n^(1 / 3))
  max_lag <- min(b, n - 1)

  y_autoc <- stats::acf(Y, plot = FALSE, type = "correlation",
                        demean = TRUE, lag.max = n - 1)$acf

  Sigma_ind <- matrix(0, nrow = m, ncol = m)
  for (k in seq_len(m)) {
    for (l in k:m) {
      x_cross <- stats::ccf(X[, k], X[, l], plot = FALSE,
                             type = "correlation", lag.max = n - 1)$acf
      lag0_idx <- n
      val <- x_cross[lag0_idx] * y_autoc[1]
      for (h in seq_len(max_lag)) {
        omega <- 1 - h / (b + 1)
        val <- val + omega * (x_cross[lag0_idx + h] +
                                x_cross[lag0_idx - h]) * y_autoc[h + 1]
      }
      Sigma_ind[k, l] <- val
      Sigma_ind[l, k] <- val
    }
  }
  Sigma_ind
}


#' Compute multivariate Pearson correlation with covariance matrix
#'
#' @param X Numeric matrix (n x m).
#' @param Y Numeric outcome vector.
#' @param IID Logical; if TRUE use IID variance, otherwise HAC.
#' @return List with \code{estimate_vector}, \code{Sigma}, \code{Sigma_ind}.
#' @keywords internal
#' @noRd
compute_pearson_multivariate_variance <- function(X, Y, IID = TRUE) {
  n <- nrow(X)
  m <- ncol(X)

  estimate_vector <- numeric(m)
  IC_matrix <- matrix(0, nrow = n, ncol = m)

  for (k in seq_len(m)) {
    r_k <- comp_pearson_rho(X[, k], Y)
    estimate_vector[k] <- r_k
    IC_matrix[, k] <- pearson_kernel(X[, k], Y, r_k)
  }

  if (IID) {
    Sigma <- (t(IC_matrix) %*% IC_matrix) / n
    Sigma_ind <- ind_covariance_pearson_iid(X)
  } else {
    Sigma <- (t(IC_matrix) %*% IC_matrix) / n +
      hac_correction_multivariate(IC_matrix)
    Sigma_ind <- ind_covariance_pearson_hac(X, Y)
  }

  list(estimate_vector = estimate_vector, Sigma = Sigma, Sigma_ind = Sigma_ind)
}
