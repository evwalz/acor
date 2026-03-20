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
