# Asymmetric Correlation Functions
# These functions provide a unified interface for computing and testing
# AKC (Asymmetric Kendall Correlation) and AGC (Asymmetric grade correlation)
#
# Internal helper functions are in acor_internal.R

#' Compute Correlation Coefficients
#' 
#' @param X Predictor variable (vector or matrix for multiple predictors)
#' @param Y Outcome variable (vector)
#' @param method Character string specifying the method:
#'   "akc", "agc", "cid", "cma",
#'   "tau_a", "tau_b", "tau_b_mod",
#'   "gamma",
#'   "rho_a", "rho_b",
#'   "pearson"
#'   
#' @return A list containing:
#'   \item{estimate}{Vector of correlation estimates}
#'   \item{method}{The method used}
#'   
#' @details
#' Asymmetric measures (directional, Y is the outcome):
#'   - AKC: Asymmetric Kendall Correlation = 2*CID - 1
#'   - CID: Concordance-Discordance Index (base measure from Kendall framework)
#'   - AGC: Asymmetric Grade Correlation = 2*CMA - 1
#'   - CMA: Coefficient of Monotone Association
#' 
#' Kendall rank correlations (symmetric):
#'   - tau_a: Kendall's tau-a (no tie correction)
#'   - tau_b: Kendall's tau-b (pair-based tie correction)
#'   - tau_b_mod: Modified tau-b (triple-based tie correction)
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
#' For multiple predictors, X should be a matrix with predictors as columns
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
                                  "tau_a", "tau_b", "tau_b_mod",
                                  "gamma",
                                  "rho_a", "rho_b"
                                  )) {
  method <- match.arg(method)
  
  validated <- validate_acor_inputs(X, Y)
  X <- validated$X
  Y <- validated$Y
  n <- validated$n
  m <- validated$m
  
  # Pre-compute Y ranks once for methods that need them
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
           tau_b_mod = kendall_tau_b_mod(xk, Y),
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
#' @param X Predictor variable (vector) or two predictors (matrix with 2 columns) to compare
#' @param Y Outcome variable (vector)
#' @param method Character string specifying the method: "akc", "agc", "cid", "cma" or "tau_a"
#' @param alternative Character string specifying the alternative hypothesis:
#'   * `"two.sided"` (default): tests if correlation differs from null value
#'   * `"greater"`: tests if correlation is greater than null value
#'   * `"less"`: tests if correlation is less than null value
#' @param conf.level Confidence level (default 0.95)
#' @param fisher Logical; if TRUE, uses Fisher transformation for confidence interval.
#' @param IID Logical; if FALSE inference performed under time series assumptions and thus HAC variance estimator is computed 
#' 
#' @return A list containing:
#'   \item{statistic}{The test statistic (z-score)}
#'   \item{p.value}{The p-value for the test}
#'   \item{estimate}{The correlation estimate(s)}
#'   \item{variance}{The variance of the estimate(s)}
#'   \item{conf.level}{Confidence level}
#'   \item{alternative}{The alternative hypothesis}
#'   \item{method}{The method used}
#'   \item{conf.int}{Confidence interval}
#'   \item{fisher}{Logical indicating if Fisher transformation was used}
#'   \item{IID}{Logical indicating if IID or time series assumptions were used}
#'   
#' @details
#' For a single predictor X, tests H0: correlation = null.value
#' For two predictors (X with 2 columns), tests H0: correlation(X1,Y) = correlation(X2,Y)
#' 
#' The test uses asymptotic normality of the correlation estimators.
#' 
#' Independence null values:
#' - AKC, AGC: H0: correlation = 0
#' - CID, CMA: H0: correlation = 0.5
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
#' # Compare two predictors
#' x1 <- rnorm(100)
#' x2 <- rnorm(100)
#' y <- rnorm(100)
#' X <- cbind(x1, x2)
#' test_result <- acor.test(X, y, method = "akc")
#' 
#' @importFrom stats setNames qnorm pnorm pchisq acf cov
#' @export
acor.test <- function(X, Y, 
                      method = c("akc", "agc", "cid", "cma", "tau_a"),
                      alternative = c("two.sided", "less", "greater"),
                      conf.level = 0.95,
                      fisher = FALSE, 
                      IID = TRUE
) {
  
  dname <- paste(deparse(substitute(X)), "and", deparse(substitute(Y)))
  method <- match.arg(method)
  alternative <- match.arg(alternative)

  validated <- validate_acor_inputs(X, Y)
  X <- validated$X
  Y <- validated$Y
  n <- validated$n
  m <- validated$m
  
  if (method %in% c("akc", "agc", "tau_a")) {
    null.value <- 0  # Independence for [-1, 1] scale
  } else {  # cid or cma
    null.value <- 0.5  # Independence for [0, 1] scale
  }
  
  # Fisher transformation only valid for AKC/AGC with single predictor
  
  # Compute correlation(s) - this already handles variance correctly
  # Compute estimates and variance based on method
  if (method %in% c("akc", "cid")) {
    # Use Kendall framework
    version <- select_kernel_version(Y, X)
    
    if (m == 1) {
      result <- compute_akc_variance_auto(X[, 1], Y, IID = IID, version = version)
      akc_val <- result$akc
      var_akc <- result$var
      var_akc_ind <- result$var_ind  # Independence variance from main function
      
      if (method == "cid") {
        estimates <- (akc_val + 1) / 2
        variance <- var_akc / 4  # Var((X+1)/2) = Var(X)/4
        variance_ind <- var_akc_ind / 4
      } else {
        estimates <- akc_val
        variance <- var_akc
        variance_ind <- var_akc_ind
      }
    } else {
      result <- compute_akc_multivariate_variance_auto(X, Y, IID = IID, version = version)
      akcs <- result$akc_vector
      Sigma <- result$Sigma
      Sigma_ind <- result$Sigma_ind  # Independence covariance from main function
      
      if (method == "cid") {
        estimates <- (akcs + 1) / 2
        variance <- Sigma / 4
        variance_ind <- Sigma_ind / 4
      } else {
        estimates <- akcs
        variance <- Sigma
        variance_ind <- Sigma_ind
      }
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
  } else {

    y_ranks <- rank(Y, ties.method = "average")
    version_agc <- select_agc_kernel_version(Y, X)
    
    if (m == 1) {
      x_ranks <- rank(X[, 1], ties.method = "average")
      result <- compute_agc_variance_auto(y_ranks, x_ranks, IID = IID, version = version_agc)
      agc_val <- result$agc
      var_agc <- result$var
      var_agc_ind <- result$var_ind  # Asymptotic variance for AGC
      
      if (method == "cma") {
        estimates <- (agc_val + 1) / 2
        variance <- var_agc / 4  # Var((X+1)/2) = Var(X)/4
        variance_ind <- var_agc_ind / 4
      } else {  # agc
        estimates <- agc_val
        variance <- var_agc
        variance_ind <- var_agc_ind 
      }
    } else {  # m >= 2
      xarray_ranks <- matrix(0, nrow = m, ncol = n)
      for (i in 1:m) {
        xarray_ranks[i, ] <- rank(X[, i], ties.method = "average")
      }
      
      result <- compute_agc_multivariate_variance_auto(y_ranks, xarray_ranks, IID = IID, version = version_agc)
      agcs <- result$agc_vector
      Sigma_agc_mat <- result$Sigma
      Sigma_agc_ind <- result$Sigma_ind
      
      if (method == "cma") {
        estimates <- (agcs + 1) / 2
        variance <- Sigma_agc_mat / 4
        variance_ind <- Sigma_agc_ind / 4
      } else {  # agc
        estimates <- agcs
        variance <- Sigma_agc_mat
        variance_ind <- Sigma_agc_ind
      }
    }
  }
  
  
  # Compute z critical value for confidence interval
  alpha <- 1 - conf.level
  z_alpha <- qnorm(1 - alpha / 2)
  
  # Perform test
  if (m == 1) {
    # Single predictor: test H0: correlation = null.value
    # The variance from acor() is already correct for the given method
    
    se <- sqrt(variance / n)
    test_stat <- (estimates - null.value) / se
    
    # Independence variance test
    se_ind <- sqrt(variance_ind / n)
    test_stat_ind <- (estimates - null.value) / se_ind
    
    # Compute p-values based on alternative
    if (alternative == "two.sided") {
      p_value <- 2 * (1 - pnorm(abs(test_stat)))
      p_value_ind <- 2 * (1 - pnorm(abs(test_stat_ind)))
    } else if (alternative == "greater") {
      p_value <- 1 - pnorm(test_stat)
      p_value_ind <- 1 - pnorm(test_stat_ind)
    } else {  # less
      p_value <- pnorm(test_stat)
      p_value_ind <- pnorm(test_stat_ind)
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
    
    ## versus htest object:
    out <- structure(list(
     statistic = c(z = test_stat),
     statistic_ind = c(z_ind = test_stat_ind),
     p.value = p_value,
     p.value_ind = p_value_ind,
     estimate = setNames(estimates, method),  # or whichever method
     variance = variance,
     variance_ind = variance_ind,
     Fisher = fisher,
     null.value = c(correlation = null.value),
     alternative = alternative,
     method = paste(toupper(method), "test"),
     data.name = dname,
     conf.int = structure(CI, conf.level = conf.level),
     IID = IID
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
    p_value <- pchisq(chi_sq_stat, df = df, lower.tail = FALSE)
    p_value_ind <- pchisq(chi_sq_stat_ind, df = df_ind, lower.tail = FALSE)
    
    # --- Individual predictor tests ---
    
    # Standard errors for individual predictors
    se_individual <- sqrt(diag(variance) / n)
    se_individual_ind <- sqrt(diag(variance_ind) / n)
    
    # Test statistics for individual predictors
    test_stat_individual <- (estimates - null.value) / se_individual
    test_stat_individual_ind <- (estimates - null.value) / se_individual_ind
    
    # P-values for individual predictors
    if (alternative == "two.sided") {
      p_value_individual <- 2 * (1 - pnorm(abs(test_stat_individual)))
      p_value_individual_ind <- 2 * (1 - pnorm(abs(test_stat_individual_ind)))
    } else if (alternative == "greater") {
      p_value_individual <- 1 - pnorm(test_stat_individual)
      p_value_individual_ind <- 1 - pnorm(test_stat_individual_ind)
    } else {  # less
      p_value_individual <- pnorm(test_stat_individual)
      p_value_individual_ind <- pnorm(test_stat_individual_ind)
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
  
  if (length(x$estimate) == 1) {
    cat("estimate =", format(x$estimate, digits = 4), "\n")
    cat("statistic =", format(x$statistic, digits = 4), ", p-value =", format.pval(x$p.value), "\n")
    if (!is.null(x$statistic_ind)) {
      cat("statistic (ind) =", format(x$statistic_ind, digits = 4), ", p-value (ind) =", format.pval(x$p.value_ind), "\n")
    }
  } else {
    cat("Chi-squared =", format(x$statistic, digits = 4), ", df =", x$df, ", p-value =", format.pval(x$p.value), "\n")
    if (!is.null(x$statistic_ind)) {
      cat("Chi-squared (ind) =", format(x$statistic_ind, digits = 4), ", df =", x$df_ind, ", p-value (ind) =", format.pval(x$p.value_ind), "\n")
    }
    cat("\nIndividual predictors:\n")
    print(x$results)
  }
  
  cat("alternative hypothesis:", x$alternative, "\n")
  if (!is.null(x$conf.int)) {
    cat(format(100 * x$conf.level), " percent confidence interval:\n")
    cat(" ", format(x$conf.int[1], digits = 4), format(x$conf.int[2], digits = 4), "\n")
  }
  cat("\n")
  invisible(x)
}
