# Asymmetric Correlation Functions
# These functions provide a unified interface for computing and testing
# AKC (Asymmetric Kendall Correlation) and AGC (Asymmetric grade correlation)
#
# Internal helper functions are in acor_internal.R

#' Compute Correlation Coefficients
#' 
#' @param X Predictor variable, supplied either as a numeric vector or as a
#'   numeric matrix with one predictor per column.
#' @param Y Numeric outcome vector.
#' @param method Character string specifying the method:
#'   "akc", "agc", "cid", "cma",
#'   "tau_a", "tau_b",
#'   "gamma",
#'   "rho_a", "rho_b",
#'   "pearson"
#'   
#' @return An object of class `"acor"` with components:
#'   \item{estimate}{A numeric vector of correlation estimates, one per column
#'   of `X`.}
#'   \item{method}{The method used.}
#'   
#' @details
#' Asymmetric measures (directional, Y is the outcome):
#'   - AKC: Asymmetric Kendall Correlation = 2*CID - 1
#'   - CID: Concordance-Discordance Index (base measure from Kendall framework).
#'     Matches `survival::concordance(Y ~ x, timefix = FALSE)$concordance` for a
#'     single numeric predictor (`timefix = FALSE` for exact Y ties, as in this package).
#'   - AGC: Asymmetric Grade Correlation = 2*CMA - 1
#'   - CMA: Coefficient of Monotone Association, equal to `(AGC + 1) / 2`
#' 
#' Kendall rank correlations (symmetric):
#'   - tau_a: Kendall's tau-a (no tie correction)
#'   - tau_b: Kendall's tau-b (pair-based tie correction)
#' 
#' Goodman-Kruskal measure (symmetric):
#'   - gamma: Goodman-Kruskal gamma (concordant-discordant pairs only)
#' 
#' Spearman rank correlations (symmetric):
#'   - rho_a: Spearman's rho without tie correction
#'   - rho_b: Spearman's rho with tie correction
#' 
#' Pearson correlation (symmetric):
#'   - pearson: Pearson product-moment correlation
#' 
#' CID and CMA range from 0 to 1 (with 0.5 = independence).
#' All other measures range from -1 to 1 (with 0 = independence).
#' 
#' For multiple predictors, `X` should be a matrix with predictors in columns.
#' The returned `estimate` vector follows the column order of `X`.
#'
#' @examples
#' # Single predictor
#' x <- rnorm(100)
#' y <- rnorm(100)
#' result_akc <- acor(x, y, method = "akc")
#' result_cma <- acor(x, y, method = "cma")
#' 
#' # Multiple predictors
#' X <- matrix(rnorm(300), ncol = 3)
#' y <- rnorm(100)
#' result <- acor(X, y, method = "agc")
#' 
#' @export
acor <- function(X, Y, method = c("pearson", "akc", "agc", "cid", "cma",
                                  "tau_a", "tau_b",
                                  "gamma",
                                  "rho_a", "rho_b"
                                  )) {
  method <- match.arg(method)
  
  validated <- validate_acor_inputs(X, Y)
  X <- validated$X
  Y <- validated$Y
  m <- validated$m

  # Pre-compute Y ranks once for methods that need them
  y_ranks <- NULL
  rank_methods <- c("agc", "cma", "rho_a", "rho_b")
  if (method %in% rank_methods) {
    y_ranks <- rank(Y, ties.method = "average")
  }
  
  compute_one <- function(xk, Y, method) {
    switch(method,
           akc = compute_akc(xk, Y),
           cid = (compute_akc(xk, Y) + 1) / 2,
           agc = compute_agc(y_ranks, rank(xk, ties.method = "average")),
           cma = (compute_agc(y_ranks, rank(xk, ties.method = "average")) + 1) / 2,
           tau_a = kendall_tau_a(xk, Y),
           tau_b = kendall_tau_b(xk, Y),
           gamma = goodman_kruskal_gamma(xk, Y),
           rho_a = comp_spearman_rho_a(y_ranks, rank(xk, ties.method = "average")),
           rho_b = comp_spearman_rho_b(y_ranks, rank(xk, ties.method = "average")),
           pearson = comp_pearson_rho(xk, Y)
    )
  }
  
  estimates <- vapply(seq_len(m), function(k) compute_one(X[, k], Y, method), numeric(1))
  
  structure(
    list(estimate = estimates, method = method),
    class = "acor"
  )
}


#' @method print acor
#' @export
print.acor <- function(x, ...) {
  cat("\n\t", toupper(x$method), "\n\n")
  cat("estimate:", format(x$estimate, digits = 4), "\n\n")
  invisible(x)
}


#' Statistical Test for Asymmetric Correlation
#' 
#' @param X Predictor variable, supplied either as a numeric vector or as a
#'   numeric matrix with one or more predictor columns.
#' @param Y Numeric outcome vector.
#' @param method Character string specifying the method: `"akc"`, `"agc"`,
#'   `"cid"`, `"cma"`, `"tau_a"`, `"tau_b"`, `"rho_a"`,
#'   `"rho_b"`, `"gamma"`, or `"pearson"`.
#' @param alternative Character string specifying the alternative hypothesis:
#'   * `"two.sided"` (default): tests if correlation differs from null value
#'   * `"greater"`: tests if correlation is greater than null value
#'   * `"less"`: tests if correlation is less than null value
#' @param conf.level Confidence level for confidence intervals. Must be a single
#'   number strictly between 0 and 1.
#' @param fisher Logical; if `TRUE`, uses a Fisher transformation when
#'   constructing confidence intervals.
#' @param IID Logical; if `FALSE`, inference is performed under time-series
#'   assumptions and a HAC variance estimator is used.
#' @param variance Character string specifying the variance estimation method:
#'   * `"ij"` (default): infinitesimal jackknife variance (supported for
#'     `"akc"`, `"agc"`, `"cid"`, `"cma"`, `"tau_b"`, `"gamma"`, and `"rho_b"`).
#'   * `"plugin"`: plug-in asymptotic variance
#'
#'   When \code{variance = "ij"} and \code{IID = FALSE}, the Bartlett-kernel HAC
#'   estimator is applied to the IJ influence function values.
#'   The independence variance uses the same closed-form formula regardless
#'   of the variance method.
#' 
#' @return An object of class `"acor_htest"`. For a single predictor, the
#'   result also inherits from `"htest"` and contains the estimate, its
#'   asymptotic variance, z-statistics and p-values under both the main and
#'   independence variance formulas, a confidence interval, the null value, and
#'   metadata such as `alternative`, `method`, `data.name`, `Fisher`, and
#'   `IID`.
#'
#'   For multiple predictors, the result contains a global chi-squared test of
#'   equality across predictors together with predictor-specific results in
#'   `results`, pairwise differences in `pairwise_results`, the covariance
#'   matrices `variance` and `variance_ind`, the contrast matrix used for the
#'   comparisons, and metadata such as `alternative`, `method`, `conf.level`,
#'   and `IID`.
#'   
#' @details
#' For a single predictor `X`, the null hypothesis is `H0: correlation =
#' null.value`.
#'
#' For multiple predictors (matrix `X`), the global null hypothesis is that all
#' predictor-specific correlations with `Y` are equal. The returned object also
#' includes individual predictor tests against the method-specific null value and
#' all pairwise differences between predictors.
#'
#' Multivariate inference (matrix `X`) is available for all `method` values.
#'
#' In the multivariate case, `alternative` affects the individual predictor
#' tests and pairwise differences. The top-level equality test remains a
#' chi-squared test of equality across predictors.
#' 
#' The test uses asymptotic normality of the correlation estimators.
#' At least 3 observations are required.
#' 
#' Independence null values:
#' - AKC, AGC: H0: correlation = 0
#' - CID, CMA: H0: correlation = 0.5
#'
#' For `method = "cid"`, compare to `survival::concordance(Y ~ x, timefix = FALSE)`
#' (exact Y ties; survival's default merges nearly-equal `Y`).
#'
#' @examples
#' # Test if AKC differs from 0 (independence test)
#' x <- rnorm(100)
#' y <- rnorm(100)
#' test_result <- acor.test(x, y, method = "akc", alternative = "two.sided")
#' 
#' # Test if CMA differs from 0.5 (independence test)
#' test_result <- acor.test(x, y, method = "cma", alternative = "two.sided")
#'
#' # Pearson inference (plugin variance; no IJ for pearson)
#' test_result <- acor.test(x, y, method = "pearson", variance = "plugin")
#'
#' # rho_b inference
#' test_result <- acor.test(x, y, method = "rho_b")
#'
#' # tau_b inference
#' test_result <- acor.test(x, y, method = "tau_b")
#'
#' # gamma inference
#' test_result <- acor.test(x, y, method = "gamma")
#' 
#' # Default IJ variance for AKC (same as omitting variance=)
#' test_ij <- acor.test(x, y, method = "akc")
#'
#' # Plug-in variance instead of default IJ
#' acor.test(x, y, method = "tau_b", variance = "plugin")
#'
#' # Compare multiple predictors
#' x1 <- rnorm(100)
#' x2 <- rnorm(100)
#' x3 <- rnorm(100)
#' y <- rnorm(100)
#' X <- cbind(x1, x2, x3)
#' test_result <- acor.test(X, y, method = "akc")
#' 
#' @export
acor.test <- function(X, Y, 
                      method = c("akc", "agc", "cid", "cma", "tau_a", "tau_b", "rho_a", "rho_b", "gamma", "pearson"),
                      alternative = c("two.sided", "less", "greater"),
                      conf.level = 0.95,
                      fisher = FALSE, 
                      IID = TRUE,
                      variance = c("ij", "plugin")) {
  
  dname <- paste(deparse(substitute(X)), "and", deparse(substitute(Y)))
  method <- match.arg(method)
  alternative <- match.arg(alternative)
  variance_method <- match.arg(variance, c("ij", "plugin"))
  if (!is.numeric(conf.level) || length(conf.level) != 1 || is.na(conf.level) ||
      conf.level <= 0 || conf.level >= 1) {
    stop("'conf.level' must be a single number between 0 and 1")
  }

  validated <- validate_acor_inputs(X, Y)
  X <- validated$X
  Y <- validated$Y
  n <- validated$n
  m <- validated$m
  if (n < 3) {
    stop("At least 3 observations are required")
  }
  
  if (variance_method == "ij") {
    ij_methods <- c("akc", "agc", "cid", "cma", "tau_b", "gamma", "rho_b")
    if (!(method %in% ij_methods)) {
      stop("variance = \"ij\" is only supported for methods: ",
           paste(ij_methods, collapse = ", "))
    }
  }
  
  # Pre-compute ranks for methods that need them
  xarray_ranks <- NULL
  x_ranks <- NULL
  y_ranks <- NULL
  
  rank_methods <- c("agc", "cma", "rho_a")
  if (method %in% rank_methods) {
    y_ranks <- rank(Y, ties.method = "average")
    if (m == 1) {
      x_ranks <- rank(X[, 1], ties.method = "average")
    } else {
      xarray_ranks <- matrix(0, nrow = m, ncol = n)
      for (i in 1:m) {
        xarray_ranks[i, ] <- rank(X[, i], ties.method = "average")
      }
    }
  }
  
  if (method %in% c("akc", "agc", "tau_a", "tau_b", "rho_a", "rho_b", "gamma", "pearson")) {
    null.value <- 0  # Independence for [-1, 1] scale
  } else {  # cid or cma
    null.value <- 0.5  # Independence for [0, 1] scale
  }
  
  # Compute correlation(s) - this already handles variance correctly
  # Compute estimates and variance based on method
  if (method %in% c("akc", "cid")) {
    if (variance_method == "ij") {
      # Default IJ: akc_ij_cpp (sample-ratio).
      tau_Y_result <- compute_tau_Y(Y)
      p_Y <- tau_Y_result$p_tie_y
      p_Y_plugin <- ((n - 1) / n) * p_Y + 1 / n

      if (m == 1) {
        result <- akc_ij_cpp(X[, 1], Y)
        akc_vals <- result$akc
        ic <- result$ic
        if (IID) {
          Sigma_akc <- result$var_ij
          Sigma_akc_ind <- ind_variance_akc_iid(X[, 1], Y, p_Y_plugin)
        } else {
          Sigma_akc <- hac_variance_univariate(ic, scale_factor = 1)
          Sigma_akc_ind <- ind_variance_akc_hac(X[, 1], Y, p_Y_plugin)
        }
      } else {
        IC_mat <- matrix(0, nrow = n, ncol = m)
        akc_vals <- numeric(m)
        for (k in seq_len(m)) {
          res_k <- akc_ij_cpp(X[, k], Y)
          akc_vals[k] <- res_k$akc
          IC_mat[, k] <- res_k$ic
        }
        if (IID) {
          Sigma_akc <- crossprod(IC_mat) / n
          Sigma_akc_ind <- ind_covariance_akc_iid(X, Y, p_Y_plugin)
        } else {
          Sigma_akc <- hac_covariance_multivariate(IC_mat, scale_factor = 1)
          Sigma_akc_ind <- ind_covariance_akc_hac(X, Y, p_Y_plugin)
        }
      }
    } else {
      if (m == 1) {
        version <- select_kernel_version(Y, X)
        result <- compute_akc_variance_auto(X[, 1], Y, IID = IID, version = version)
        akc_vals <- result$akc
        Sigma_akc <- result$var
        Sigma_akc_ind <- result$var_ind
      } else {
        version <- select_kernel_version(Y, X)
        result <- compute_akc_multivariate_variance_auto(X, Y, IID = IID, version = version)
        akc_vals <- result$akc_vector
        Sigma_akc <- result$Sigma
        Sigma_akc_ind <- result$Sigma_ind
      }
    }

    if (method == "cid") {
      estimates <- (akc_vals + 1) / 2
      variance <- Sigma_akc / 4
      variance_ind <- Sigma_akc_ind / 4
    } else {
      estimates <- akc_vals
      variance <- Sigma_akc
      variance_ind <- Sigma_akc_ind
    }
  } else if (method == "tau_a") {
    version <- select_kernel_version(Y, X)
    
    if (m == 1) {
      result <- compute_tau_a_variance(X[, 1], Y, IID = IID, version = version)
      estimates <- result$tau_a
      variance <- result$var
      variance_ind <- result$var_ind
    } else {
      result <- compute_tau_a_multivariate_variance(X, Y, IID = IID, version = version)
      estimates <- result$tau_a_vector
      variance <- result$Sigma
      variance_ind <- result$Sigma_ind
    } 
  } else if (method == "tau_b") {
    if (variance_method == "ij") {
      if (m == 1) {
        result <- tau_b_ij_cpp(X[, 1], Y)
        estimates <- result$tau_b
        ic <- result$ic
        variance_ind <- ind_variance_tau_b_univariate(X[, 1], Y, IID = IID)
        if (IID) {
          variance <- result$var_ij
        } else {
          variance <- hac_variance_univariate(ic, scale_factor = 1)
        }
      } else {
        IC_mat <- matrix(0, nrow = n, ncol = m)
        estimates <- numeric(m)
        for (k in seq_len(m)) {
          res_k <- tau_b_ij_cpp(X[, k], Y)
          estimates[k] <- res_k$tau_b
          IC_mat[, k] <- res_k$ic
        }
        if (IID) {
          variance <- crossprod(IC_mat) / n
          variance_ind <- ind_covariance_tau_b_iid(X, Y)
        } else {
          variance <- hac_covariance_multivariate(IC_mat, scale_factor = 1)
          variance_ind <- ind_covariance_tau_b_hac(X, Y)
        }
      }
    } else {
      version <- select_kernel_version(Y, X)
      if (m == 1) {
        result <- compute_tau_b_variance(X[, 1], Y, IID = IID, version = version)
        estimates <- result$tau_b
        variance <- result$var
        variance_ind <- result$var_ind
      } else {
        result <- compute_tau_b_multivariate_variance(X, Y, IID = IID, version = version)
        estimates <- result$tau_b_vector
        variance <- result$Sigma
        variance_ind <- result$Sigma_ind
      }
    }
  } else if (method == "gamma") {
    if (variance_method == "ij") {
      if (m == 1) {
        result <- gamma_ij_cpp(X[, 1], Y)
        estimates <- result$gamma
        ic <- result$ic
        variance_ind <- ind_variance_gamma_univariate(X[, 1], Y, IID = IID)
        if (IID) {
          variance <- result$var_ij
        } else {
          variance <- hac_variance_univariate(ic, scale_factor = 1)
        }
      } else {
        IC_mat <- matrix(0, nrow = n, ncol = m)
        estimates <- numeric(m)
        for (k in seq_len(m)) {
          res_k <- gamma_ij_cpp(X[, k], Y)
          estimates[k] <- res_k$gamma
          IC_mat[, k] <- res_k$ic
        }
        if (IID) {
          variance <- crossprod(IC_mat) / n
        } else {
          variance <- hac_covariance_multivariate(IC_mat, scale_factor = 1)
        }
        variance_ind <- if (IID) {
          ind_covariance_gamma_iid(X, Y)
        } else {
          ind_covariance_gamma_hac(X, Y)
        }
      }
    } else {
      version <- select_kernel_version(Y, X)
      if (m == 1) {
        result <- compute_gamma_variance(X[, 1], Y, IID = IID, version = version)
        estimates <- result$gamma
        variance <- result$var
        variance_ind <- result$var_ind
      } else {
        result <- compute_gamma_multivariate_variance(X, Y, IID = IID, version = version)
        estimates <- result$gamma_vector
        variance <- result$Sigma
        variance_ind <- result$Sigma_ind
      }
    }
  } else if (method == "rho_a") {
    
    if (m == 1) {
      result <- compute_rho_a_variance(x_ranks, y_ranks, IID = IID)
      estimates <- result$rho_a
      variance <- result$var
      variance_ind <- result$var_ind
    } else {
      result <- compute_rho_a_multivariate_variance(y_ranks, xarray_ranks, IID = IID)
      estimates <- result$rho_a_vector
      variance <- result$Sigma
      variance_ind <- result$Sigma_ind
    }
    
  } else if (method == "rho_b") {
    if (variance_method == "ij") {
      if (m == 1) {
        result <- rho_b_ij_cpp(X[, 1], Y)
        estimates <- result$rho_b
        ic <- result$ic
        variance_ind <- if (IID) {
          ind_variance_rho_b_iid()
        } else {
          ind_variance_rho_b_hac(X[, 1], Y)
        }
        if (IID) {
          variance <- result$var_ij
        } else {
          variance <- hac_variance_univariate(ic, scale_factor = 1)
        }
      } else {
        IC_mat <- matrix(0, nrow = n, ncol = m)
        estimates <- numeric(m)
        for (k in seq_len(m)) {
          res_k <- rho_b_ij_cpp(X[, k], Y)
          estimates[k] <- res_k$rho_b
          IC_mat[, k] <- res_k$ic
        }
        if (IID) {
          variance <- crossprod(IC_mat) / n
          variance_ind <- ind_covariance_rho_b_iid(X, Y)
        } else {
          variance <- hac_covariance_multivariate(IC_mat, scale_factor = 1)
          variance_ind <- ind_covariance_rho_b_hac(X, Y)
        }
      }
    } else if (m == 1) {
      result <- compute_rho_b_variance(X[, 1], Y, IID = IID)
      estimates <- result$rho_b
      variance <- result$var
      variance_ind <- result$var_ind
    } else {
      result <- compute_rho_b_multivariate_variance(X, Y, IID = IID)
      estimates <- result$rho_b_vector
      variance <- result$Sigma
      variance_ind <- result$Sigma_ind
    }
  } else if (method == "pearson") {
    
    if (m == 1) {
      result <- compute_pearson_variance(X[, 1], Y, IID = IID)
      estimates <- result$estimate
      variance <- result$var
      variance_ind <- result$var_ind
    } else {
      result <- compute_pearson_multivariate_variance(X, Y, IID = IID)
      estimates <- result$estimate_vector
      variance <- result$Sigma
      variance_ind <- result$Sigma_ind
    }
    
  } else {
    
    if (variance_method == "ij") {
      # Reuse the y_ranks / x_ranks / xarray_ranks pre-computed at the top of
      # acor.test (for method %in% c("agc", "cma", "rho_a")) -- avoid extra
      # rank() passes that were previously costing ~2 ms at n=5000.
      pre_ij <- agc_y_preamble(y_ranks)
      b_ij <- floor(2 * n^(1 / 3))

      if (m == 1) {
        result <- agc_ij_cpp(X[, 1], Y)
        agc_vals <- result$agc
        ic <- result$ic
        if (IID) {
          Sigma_agc <- result$var_ij
          Sigma_agc_ind <- ind_variance_agc_iid(x_ranks, pre_ij$N, pre_ij$zeta_3Y)
        } else {
          Sigma_agc <- hac_variance_univariate(ic, scale_factor = 1)
          Sigma_agc_ind <- ind_variance_agc_hac(x_ranks, y_ranks,
                                                pre_ij$N, pre_ij$zeta_3Y, b_ij)
        }
      } else {
        IC_mat <- matrix(0, nrow = n, ncol = m)
        agc_vals <- numeric(m)
        for (k in seq_len(m)) {
          res_k <- agc_ij_cpp(X[, k], Y)
          agc_vals[k] <- res_k$agc
          IC_mat[, k] <- res_k$ic
        }
        if (IID) {
          Sigma_agc <- crossprod(IC_mat) / n
          Sigma_agc_ind <- ind_covariance_agc_iid(xarray_ranks, pre_ij$N, pre_ij$zeta_3Y)
        } else {
          Sigma_agc <- hac_covariance_multivariate(IC_mat, scale_factor = 1)
          Sigma_agc_ind <- ind_covariance_agc_hac(xarray_ranks, y_ranks,
                                                  pre_ij$N, pre_ij$zeta_3Y, b_ij)
        }
      }
    } else {
      if (m == 1) {
        result <- compute_agc_variance_auto(y_ranks, x_ranks, IID = IID)
        agc_vals <- result$agc
        Sigma_agc <- result$var
        Sigma_agc_ind <- result$var_ind
      } else {
        result <- compute_agc_multivariate_variance_auto(y_ranks, xarray_ranks, IID = IID)
        agc_vals <- result$agc_vector
        Sigma_agc <- result$Sigma
        Sigma_agc_ind <- result$Sigma_ind
      }
    }
    
    if (method == "cma") {
      estimates <- (agc_vals + 1) / 2
      variance <- Sigma_agc / 4
      variance_ind <- Sigma_agc_ind / 4
    } else {  # agc
      estimates <- agc_vals
      variance <- Sigma_agc
      variance_ind <- Sigma_agc_ind
    }
  }
  
  
  # Compute z critical value for confidence interval
  alpha <- 1 - conf.level
  z_alpha <- stats::qnorm(1 - alpha / 2)
  
  # Perform test
  if (m == 1) {
    # Single predictor: test H0: correlation = null.value
    # The variance from acor() is already correct for the given method
    
    se <- sqrt(variance / n)
    test_stat <- (estimates - null.value) / se
    
    # Independence variance test (NULL when variance = "ij")
    if (!is.null(variance_ind)) {
      se_ind <- sqrt(variance_ind / n)
      test_stat_ind <- (estimates - null.value) / se_ind
    } else {
      test_stat_ind <- NULL
    }
    
    # Compute p-values based on alternative
    if (alternative == "two.sided") {
      p_value <- 2 * (1 - stats::pnorm(abs(test_stat)))
      p_value_ind <- if (!is.null(test_stat_ind))
        2 * (1 - stats::pnorm(abs(test_stat_ind))) else NULL
    } else if (alternative == "greater") {
      p_value <- 1 - stats::pnorm(test_stat)
      p_value_ind <- if (!is.null(test_stat_ind))
        1 - stats::pnorm(test_stat_ind) else NULL
    } else {  # less
      p_value <- stats::pnorm(test_stat)
      p_value_ind <- if (!is.null(test_stat_ind))
        stats::pnorm(test_stat_ind) else NULL
    }
    
    # Confidence interval
    if (fisher) {
      # Fisher transformation - always compute on AKC/AGC scale first
      if (method %in% c("cid", "cma")) {
        # Transform estimate back to [-1, 1] scale for Fisher transformation
        est_transformed <- 2 * estimates - 1  # CID/CMA -> AKC/AGC
        se_transformed <- 2 * se  # SE also scales by 2
      } else {
        est_transformed <- estimates
        se_transformed <- se
      }
      # Clamp to avoid Inf from atanh at boundaries
      est_transformed <- pmax(pmin(est_transformed, 1 - 1e-10), -1 + 1e-10)
      
      # Fisher transformation on [-1, 1] scale
      z_est <- atanh(est_transformed)
      deriv <- 1 / (1 - est_transformed^2)
      se_z <- se_transformed * abs(deriv)
      
      z_lower <- z_est - z_alpha * se_z
      z_upper <- z_est + z_alpha * se_z
      CI_transformed <- c(tanh(z_lower), tanh(z_upper))
      
      if (method %in% c("cid", "cma")) {
        # Transform CI back to [0, 1] scale
        CI <- (CI_transformed + 1) / 2
      } else {
        CI <- CI_transformed
      }
    } else {
      CI <- c(estimates - z_alpha * se, estimates + z_alpha * se)
    }
    
    stat_ind_named <- if (!is.null(test_stat_ind)) c(z_ind = test_stat_ind) else NULL
    
    out <- structure(list(
     statistic = c(z = test_stat),
     statistic_ind = stat_ind_named,
     p.value = p_value,
     p.value_ind = p_value_ind,
     estimate = stats::setNames(estimates, method),
     variance = variance,
     variance_ind = variance_ind,
     Fisher = fisher,
     null.value = c(correlation = null.value),
     alternative = alternative,
     method = paste(toupper(method), "test"),
     data.name = dname,
     conf.int = structure(CI, conf.level = conf.level),
     IID = IID,
     variance_method = variance_method
    ), class = c("acor_htest", "htest"))
    
  } else {
    # Multiple predictors (m >= 2)
    
    # --- Build contrast matrix L for pairwise differences ---
    n_pairs <- m * (m - 1) / 2
    L <- matrix(0, nrow = n_pairs, ncol = m)
    
    row_idx <- 1
    for (i in 1:(m - 1)) {
      for (j in (i + 1):m) {
        L[row_idx, i] <- 1
        L[row_idx, j] <- -1
        row_idx <- row_idx + 1
      }
    }
    
    # --- Chi-square test for H0: all correlations are equal ---
    # Compute pairwise differences
    est_diff <- as.vector(L %*% estimates)
  
    # Variance of differences (scaled by n)
    L_S_Lt <- L %*% (variance / n) %*% t(L)
    L_S_Lt_ind <- L %*% (variance_ind / n) %*% t(L)
    
    # Use pseudoinverse for potential singularity
    L_S_Lt_inv <- MASS::ginv(L_S_Lt)
    L_S_Lt_inv_ind <- MASS::ginv(L_S_Lt_ind)
    
    # Chi-square statistic
    chi_sq_stat <- as.numeric(t(est_diff) %*% L_S_Lt_inv %*% est_diff)
    chi_sq_stat_ind <- as.numeric(t(est_diff) %*% L_S_Lt_inv_ind %*% est_diff)
    
    # Degrees of freedom (rank of L_S_Lt)
    qr_decomp <- qr(L_S_Lt)
    df <- qr_decomp$rank
    
    qr_decomp_ind <- qr(L_S_Lt_ind)
    df_ind <- qr_decomp_ind$rank
    
    # P-values from chi-square distribution
    p_value <- stats::pchisq(chi_sq_stat, df = df, lower.tail = FALSE)
    p_value_ind <- stats::pchisq(chi_sq_stat_ind, df = df_ind, lower.tail = FALSE)
    
    # --- Individual predictor tests ---
    
    # Standard errors for individual predictors
    se_individual <- sqrt(diag(variance) / n)
    se_individual_ind <- sqrt(diag(variance_ind) / n)
    
    # Test statistics for individual predictors
    test_stat_individual <- (estimates - null.value) / se_individual
    test_stat_individual_ind <- (estimates - null.value) / se_individual_ind
    
    # P-values for individual predictors
    if (alternative == "two.sided") {
      p_value_individual <- 2 * (1 - stats::pnorm(abs(test_stat_individual)))
      p_value_individual_ind <- 2 * (1 - stats::pnorm(abs(test_stat_individual_ind)))
    } else if (alternative == "greater") {
      p_value_individual <- 1 - stats::pnorm(test_stat_individual)
      p_value_individual_ind <- 1 - stats::pnorm(test_stat_individual_ind)
    } else {  # less
      p_value_individual <- stats::pnorm(test_stat_individual)
      p_value_individual_ind <- stats::pnorm(test_stat_individual_ind)
    }
    
    # Confidence intervals for individual predictors
    if (fisher) {
      # Fisher transformation for individual CIs
      if (method %in% c("cid", "cma")) {
        est_transformed <- 2 * estimates - 1  # CID/CMA -> AKC/AGC scale
        se_transformed <- 2 * se_individual
      } else {
        est_transformed <- estimates
        se_transformed <- se_individual
      }
      
      # Clamp to avoid Inf from atanh at boundaries
      est_transformed <- pmax(pmin(est_transformed, 1 - 1e-10), -1 + 1e-10)
      
      z_est <- atanh(est_transformed)
      deriv <- 1 / (1 - est_transformed^2)
      se_z <- se_transformed * abs(deriv)
      
      z_lower <- z_est - z_alpha * se_z
      z_upper <- z_est + z_alpha * se_z
      
      if (method %in% c("cid", "cma")) {
        CI_lower_individual <- (tanh(z_lower) + 1) / 2
        CI_upper_individual <- (tanh(z_upper) + 1) / 2
      } else {
        CI_lower_individual <- tanh(z_lower)
        CI_upper_individual <- tanh(z_upper)
      }
    } else {
      CI_lower_individual <- estimates - z_alpha * se_individual
      CI_upper_individual <- estimates + z_alpha * se_individual
    }
    
    # Build results table for individual predictors
    results_table <- data.frame(
      predictor = paste0("X", 1:m),
      estimate = estimates,
      statistic = test_stat_individual,
      statistic_ind = test_stat_individual_ind,
      p.value = p_value_individual,
      p.value_ind = p_value_individual_ind,
      CI_lower = CI_lower_individual,
      CI_upper = CI_upper_individual
    )
    
    # --- Pairwise difference CIs ---
    se_diff <- sqrt(diag(L_S_Lt))
    
    pair_labels <- character(n_pairs)
    row_idx <- 1
    for (i in 1:(m - 1)) {
      for (j in (i + 1):m) {
        pair_labels[row_idx] <- paste0("X", i, " - X", j)
        row_idx <- row_idx + 1
      }
    }
    
    zero_tol <- sqrt(.Machine$double.eps)
    z_diff <- est_diff / se_diff
    zero_var <- se_diff <= zero_tol
    if (any(zero_var)) {
      z_diff[zero_var & abs(est_diff) <= zero_tol] <- 0
      z_diff[zero_var & est_diff > zero_tol] <- Inf
      z_diff[zero_var & est_diff < -zero_tol] <- -Inf
    }
    
    if (alternative == "two.sided") {
      p_diff <- 2 * (1 - stats::pnorm(abs(z_diff)))
    } else if (alternative == "greater") {
      p_diff <- 1 - stats::pnorm(z_diff)
    } else {
      p_diff <- stats::pnorm(z_diff)
    }
    
    pairwise_table <- data.frame(
      pair = pair_labels,
      difference = est_diff,
      statistic = z_diff,
      p.value = p_diff,
      CI_lower = est_diff - z_alpha * se_diff,
      CI_upper = est_diff + z_alpha * se_diff
    )
    
    # Create result object for m >= 2
    out <- list(
      statistic = chi_sq_stat,
      statistic_ind = chi_sq_stat_ind,
      df = df,
      df_ind = df_ind,
      p.value = p_value,
      p.value_ind = p_value_ind,
      estimate = estimates,
      variance = variance,
      variance_ind = variance_ind,
      Fisher = fisher,
      conf.level = conf.level,
      alternative = alternative,
      method = paste(toupper(method), "test"),
      IID = IID,
      results = results_table,
      pairwise_results = pairwise_table, 
      pairwise_differences = est_diff,
      contrast_matrix = L, 
      null.value = c(difference = 0),
      data.name = dname
    )
    class(out) <- "acor_htest"
  }
  
  
  return(out)
}


#' Print method for acor_htest objects
#'
#' @param x An object of class acor_htest
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns x
#' @method print acor_htest
#' @export
print.acor_htest <- function(x, ...) {
  cat("\n\t", x$method, "\n\n")
  cat("data: ", x$data.name, "\n")
  if (!is.null(x$variance_method) && x$variance_method == "ij") {
    cat("variance: infinitesimal jackknife\n")
  }
  
  if (length(x$estimate) == 1) {
    # --- Single predictor ---
    cat("statistic =", format(x$statistic, digits = 4),
        ", p-value =", format.pval(x$p.value), "\n")
    if (!is.null(x$statistic_ind)) {
      cat("statistic (ind) =", format(x$statistic_ind, digits = 4),
          ", p-value (ind) =", format.pval(x$p.value_ind), "\n")
    }
    display_name <- switch(names(x$estimate),
                           akc = "AKC", agc = "AGC", cid = "CID", cma = "CMA",
                           tau_a = "tau-a", tau_b = "tau-b",
                           gamma = "gamma",
                           rho_a = "rho-a", rho_b = "rho-b",
                           pearson = "Pearson r",
                           names(x$estimate)
    )
    
    alt_text <- switch(x$alternative,
                       two.sided = "is not equal to",
                       greater = "is greater than",
                       less = "is less than")
    cat("alternative hypothesis: true",
        display_name, alt_text, x$null.value, "\n")
    if (!is.null(x$conf.int)) {
      cl <- attr(x$conf.int, "conf.level")
      cat(format(100 * cl), "percent confidence interval:\n")
      cat(" ", format(x$conf.int[1], digits = 4),
          format(x$conf.int[2], digits = 4), "\n")
    }
    cat("sample estimates:\n")
    cat(" ", names(x$estimate), "\n")
    cat(" ", format(x$estimate, digits = 4), "\n")
    
  } else {
    # --- Multiple predictors ---
    cat("Equality test (H0: all correlations equal):\n")
    cat("Chi-squared =", format(x$statistic, digits = 4),
        ", df =", x$df,
        ", p-value =", format.pval(x$p.value), "\n")
    if (!is.null(x$statistic_ind)) {
      cat("Chi-squared (ind) =", format(x$statistic_ind, digits = 4),
          ", df =", x$df_ind,
          ", p-value (ind) =", format.pval(x$p.value_ind), "\n")
    }
    
    cat("\nIndividual predictors:\n")
    print(x$results, row.names = FALSE)
    
    if (!is.null(x$pairwise_results)) {
      cl <- x$conf.level
      cat("\nPairwise differences")
      if (!is.null(cl)) cat(" (", format(100 * cl), "% CI)", sep = "")
      cat(":\n")
      print(x$pairwise_results, row.names = FALSE)
    }
    
    cat("\nalternative hypothesis:", x$alternative, "\n")
  }
  
  cat("\n")
  invisible(x)
}