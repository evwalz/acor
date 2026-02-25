# ============================================================================
# Tests comparing acor.test() variance to external packages
# ============================================================================
# 
# These tests validate the variance estimators by comparing to:
# 1. survival::concordance - for CID variance
# 2. pROC::var (DeLong) - for binary Y (AUC) variance
# 3. cor.test - for Kendall/Spearman under independence
#
# Note: Some tests require additional packages (pROC)
# ============================================================================

library(testthat)
library(survival)
library(acor)  # Replace with your actual package name

# ============================================================================
# Test A: Compare CID variance to survival::concordance
# ============================================================================

test_that("CID variance matches survival::concordance variance (no ties)", {
  set.seed(2001)
  n <- 200
  
  # Ensure no ties by adding small noise (matching test_akc_fast.R pattern)
  X <- rnorm(n) + rnorm(n, sd = 0.0001)
  Y <- rnorm(n) + rnorm(n, sd = 0.0001)
  
  # Our CID result
  result_cid <- acor.test(X, Y, method = "cid")
  
  # survival::concordance
  conc <- concordance(Y ~ X)
  
  # Compare estimates
  expect_equal(unname(result_cid$estimate), unname(conc$concordance), 
               tolerance = 1e-3,
               info = "CID estimate should match concordance")
  
  # Compare variance (survival gives variance, we give asymptotic variance / n)
  # Our variance is the asymptotic variance (not divided by n in storage)
  # SE = sqrt(variance / n), so variance / n = conc$var
  our_var_scaled <- result_cid$variance / n
  
  expect_equal(our_var_scaled, unname(conc$var), 
               tolerance = 0.05,  # Allow 15% relative difference (different estimators)
               info = "CID variance/n should approximately match concordance variance")
})

test_that("CID variance matches survival::concordance variance (with ties in X)", {
  set.seed(2002)
  n <- 200
  X <- rnorm(n)
  X[seq(1, n, by = 5)] <- X[1]  # Create ties
  Y <- rnorm(n) + rnorm(n, sd = 0.0001)  # No ties in Y
  
  result_cid <- acor.test(X, Y, method = "cid")
  conc <- concordance(Y ~ X)
  
  # Compare estimates
  expect_equal(unname(result_cid$estimate), unname(conc$concordance), 
               tolerance = 1e-4,
               info = "CID estimate should match concordance with ties")
  
  # Compare variance
  our_var_scaled <- result_cid$variance / n
  
  expect_equal(our_var_scaled, unname(conc$var), 
               tolerance = 0.05,  # Allow 20% for tied case
               info = "CID variance should approximately match concordance variance with ties")
})

test_that("CID variance matches survival::concordance for discrete X (10 levels)", {
  set.seed(2005)
  n <- 200
  X <- sample(1:10, n, replace = TRUE)
  Y <- rnorm(n) + rnorm(n, sd = 0.0001)
  
  result_cid <- acor.test(X, Y, method = "cid")
  conc <- concordance(Y ~ X)
  
  # Compare estimates
  expect_equal(unname(result_cid$estimate), unname(conc$concordance), 
               tolerance = 1e-4,
               info = "CID should match concordance for discrete X (10 levels)")
  
  # Compare variance
  our_var_scaled <- result_cid$variance / n
  
  expect_equal(our_var_scaled, unname(conc$var), 
               tolerance = 0.05,  # Allow 25% for discrete case
               info = "CID variance should match concordance for discrete X")
})

test_that("CID variance matches survival::concordance for discrete Y (5 levels)", {
  set.seed(2006)
  n <- 200
  X <- rnorm(n)
  Y <- sample(1:5, n, replace = TRUE)
  
  result_cid <- acor.test(X, Y, method = "cid")
  conc <- concordance(Y ~ X)
  
  # Compare estimates
  expect_equal(unname(result_cid$estimate), unname(conc$concordance), 
               tolerance = 1e-4,
               info = "CID should match concordance for discrete Y (5 levels)")
  
  # Compare variance
  our_var_scaled <- result_cid$variance / n
  
  expect_equal(our_var_scaled, unname(conc$var), 
               tolerance = 0.05,  # Allow 25% for discrete case
               info = "CID variance should match concordance for discrete Y")
})

test_that("CID variance matches survival::concordance for discrete X and Y", {
  set.seed(2007)
  n <- 300
  X <- sample(1:10, n, replace = TRUE)
  Y <- sample(1:5, n, replace = TRUE)
  
  result_cid <- acor.test(X, Y, method = "cid")
  conc <- concordance(Y ~ X)
  
  # Compare estimates
  expect_equal(unname(result_cid$estimate), unname(conc$concordance), 
               tolerance = 1e-4,
               info = "CID should match concordance for discrete X and Y")
  
  # Compare variance
  our_var_scaled <- result_cid$variance / n
  
  expect_equal(our_var_scaled, unname(conc$var), 
               tolerance = 0.05,  # Allow 30% for highly discrete case
               info = "CID variance should match concordance for discrete X and Y")
})

test_that("CID variance matches survival::concordance across discrete levels", {
  # Test multiple discrete configurations
  test_cases <- list(
    list(name = "3 levels in Y", n = 200, X_levels = NA, Y_levels = 3, seed = 2008),
    list(name = "20 levels in X", n = 200, X_levels = 20, Y_levels = NA, seed = 2009),
    list(name = "5x5 discrete", n = 250, X_levels = 5, Y_levels = 5, seed = 2010),
    list(name = "10x3 discrete", n = 250, X_levels = 10, Y_levels = 3, seed = 2011),
    list(name = "Many ties (2x2)", n = 200, X_levels = 2, Y_levels = 2, seed = 2012)
  )
  
  for (tc in test_cases) {
    set.seed(tc$seed)
    
    # Generate X
    if (is.na(tc$X_levels)) {
      X <- rnorm(tc$n)
    } else {
      X <- sample(1:tc$X_levels, tc$n, replace = TRUE)
    }
    
    # Generate Y
    if (is.na(tc$Y_levels)) {
      Y <- rnorm(tc$n) + rnorm(tc$n, sd = 0.0001)
    } else {
      Y <- sample(1:tc$Y_levels, tc$n, replace = TRUE)
    }
    
    result_cid <- acor.test(X, Y, method = "cid")
    conc <- concordance(Y ~ X)
    
    # Compare estimates
    expect_equal(unname(result_cid$estimate), unname(conc$concordance), 
                 tolerance = 1e-4,
                 info = sprintf("%s: CID should match concordance", tc$name))
    
    # Compare variance (with generous tolerance for discrete cases)
    our_var_scaled <- result_cid$variance / tc$n
    
    expect_equal(our_var_scaled, unname(conc$var), 
                 tolerance = 0.1,  # Allow 35% for various discrete cases
                 info = sprintf("%s: CID variance should match concordance", tc$name))
  }
})

test_that("CID variance matches survival::concordance for binary Y", {
  set.seed(2003)
  n <- 200
  X <- rnorm(n)
  Y <- rbinom(n, 1, 0.6)
  
  result_cid <- acor.test(X, Y, method = "cid")
  conc <- concordance(Y ~ X)
  
  # Compare estimates (should be exact for binary Y)
  expect_equal(unname(result_cid$estimate), unname(conc$concordance), 
               tolerance = 1e-4,
               info = "CID should equal concordance for binary Y")
  
  # Variance comparison
  our_var_scaled <- result_cid$variance / n
  
  expect_equal(our_var_scaled, unname(conc$var), 
               tolerance = 0.05,  # Allow 25% difference between estimators
               info = "CID variance should match concordance variance for binary Y")
})

test_that("AKC variance matches survival::concordance variance (scaled)", {
  # AKC = 2*CID - 1, so Var(AKC) = 4 * Var(CID)
  set.seed(2004)
  n <- 200
  X <- rnorm(n) + rnorm(n, sd = 0.0001)
  Y <- rnorm(n) + rnorm(n, sd = 0.0001)
  
  result_akc <- acor.test(X, Y, method = "akc")
  result_cid <- acor.test(X, Y, method = "cid")
  conc <- concordance(Y ~ X)
  
  # AKC variance should be 4x CID variance
  expect_equal(result_akc$variance, 4 * result_cid$variance, 
               tolerance = 1e-5,
               info = "AKC variance = 4 * CID variance")
  
  # And should match 4x concordance variance
  conc_var_scaled <- conc$var * n  # Scale up to asymptotic variance
  expect_equal(result_akc$variance, 4 * conc_var_scaled, 
               tolerance = 0.05,  # Allow 20% difference
               info = "AKC variance should be ~4x concordance variance")
})


# ============================================================================
# Test B: Compare AGC/AKC variance to cor.test under independence
# ============================================================================

# Under independence with no ties:
# - AGC = Spearman's rho, so variances should match
# - AKC = Kendall's tau, so variances should match
#
# cor.test gives z-statistic for Kendall and S-statistic for Spearman
# We can back-calculate the implied variance from the test statistic

test_that("AGC variance matches Spearman variance under independence (no ties)", {
  set.seed(3001)
  n <- 200
  
  # Independent data with no ties (matching test_akc_fast.R pattern)
  X <- rnorm(n) + rnorm(n, sd = 0.0001)
  Y <- rnorm(n) + rnorm(n, sd = 0.0001)
  
  # Our AGC result with independence variance
  result_agc <- acor.test(X, Y, method = "agc")
  
  # cor.test for Spearman
  spearman_test <- cor.test(X, Y, method = "spearman")
  
  # The estimates should match (AGC = Spearman when no ties)
  expect_equal(unname(result_agc$estimate), unname(spearman_test$estimate), 
               tolerance = 1e-5,
               info = "AGC should equal Spearman's rho with no ties")
  
  # Under independence, the variance can be computed from the test statistic
  # For Spearman: statistic S, and p-value computed from t-distribution or exact
  # The asymptotic variance under H0 is 1/(n-1) for Spearman
  # Our variance_ind should be close to 1 (the asymptotic variance before /n scaling)
  
  # Simple check: independence variance should give similar p-value
  our_z <- result_agc$statistic_ind
  our_p <- result_agc$p.value_ind
  
  # P-values should be similar (both testing independence)
  expect_equal(our_p, spearman_test$p.value, 
               tolerance = 0.1,  # Allow some difference due to exact vs asymptotic
               info = "AGC independence p-value should be similar to Spearman p-value")
})

test_that("AKC variance matches Kendall variance under independence (no ties)", {
  set.seed(3002)
  n <- 200
  
  # Independent data with no ties (matching test_akc_fast.R pattern)
  X <- rnorm(n) + rnorm(n, sd = 0.0001)
  Y <- rnorm(n) + rnorm(n, sd = 0.0001)
  
  # Our AKC result
  result_akc <- acor.test(X, Y, method = "akc")
  
  # cor.test for Kendall (uses z-statistic asymptotically)
  kendall_test <- cor.test(X, Y, method = "kendall")
  
  # Estimates should match
  expect_equal(unname(result_akc$estimate), unname(kendall_test$estimate), 
               tolerance = 1e-6,
               info = "AKC should equal Kendall's tau with no ties")
  
  # Kendall's test uses z-statistic: z = tau / sqrt(var)
  # where var = 2(2n+5) / (9n(n-1)) under independence
  # So we can back-calculate the variance used by cor.test
  kendall_var_theory <- 2 * (2 * n + 5) / (9 * n * (n - 1))
  
  # Our variance_ind / n should be close to this
  our_var_scaled <- result_akc$variance_ind / n
  
  expect_equal(our_var_scaled, kendall_var_theory, 
               tolerance = 0.05,
               info = "AKC independence variance should match Kendall theoretical variance")
  
  # P-values should be similar
  expect_equal(result_akc$p.value_ind, kendall_test$p.value, 
               tolerance = 0.1,
               info = "AKC independence p-value should be similar to Kendall p-value")
})

test_that("Multiple samples: AGC/AKC p-values track Spearman/Kendall p-values", {
  # Test across multiple random samples to ensure consistency
  set.seed(3003)
  n_tests <- 20
  n <- 150
  
  agc_p <- numeric(n_tests)
  spearman_p <- numeric(n_tests)
  akc_p <- numeric(n_tests)
  kendall_p <- numeric(n_tests)
  
  for (i in 1:n_tests) {
    X <- rnorm(n) + rnorm(n, sd = 0.0001)
    Y <- rnorm(n) + rnorm(n, sd = 0.0001)
    
    agc_p[i] <- acor.test(X, Y, method = "agc")$p.value_ind
    spearman_p[i] <- cor.test(X, Y, method = "spearman")$p.value
    
    akc_p[i] <- acor.test(X, Y, method = "akc")$p.value_ind
    kendall_p[i] <- cor.test(X, Y, method = "kendall")$p.value
  }
  
  # AGC and Spearman p-values should be highly correlated
  cor_agc_spearman <- cor(agc_p, spearman_p)
  expect_gt(cor_agc_spearman, 0.95)
  
  # AKC and Kendall p-values should be highly correlated
  cor_akc_kendall <- cor(akc_p, kendall_p)
  expect_gt(cor_akc_kendall, 0.95)
})

test_that("AGC/AKC handle discrete data correctly", {
  # Test with discrete X and Y (different number of levels)
  test_cases <- list(
    list(name = "Discrete X (10 levels), continuous Y", seed = 3004,
         gen_X = function(n) sample(1:10, n, replace = TRUE),
         gen_Y = function(n) rnorm(n)),
    list(name = "Continuous X, discrete Y (5 levels)", seed = 3005,
         gen_X = function(n) rnorm(n),
         gen_Y = function(n) sample(1:5, n, replace = TRUE)),
    list(name = "Discrete X (5 levels), discrete Y (3 levels)", seed = 3006,
         gen_X = function(n) sample(1:5, n, replace = TRUE),
         gen_Y = function(n) sample(1:3, n, replace = TRUE)),
    list(name = "Discrete X (20 levels), discrete Y (10 levels)", seed = 3007,
         gen_X = function(n) sample(1:20, n, replace = TRUE),
         gen_Y = function(n) sample(1:10, n, replace = TRUE)),
    list(name = "Many ties in X", seed = 3008,
         gen_X = function(n) { x <- rnorm(n); x[seq(1, n, by = 5)] <- x[1]; x },
         gen_Y = function(n) rnorm(n) + rnorm(n, sd = 0.0001))
  )
  
  n <- 200
  
  for (tc in test_cases) {
    set.seed(tc$seed)
    X <- tc$gen_X(n)
    Y <- tc$gen_Y(n)
    
    # Both methods should run without error
    result_agc <- acor.test(X, Y, method = "agc")
    result_akc <- acor.test(X, Y, method = "akc")
    
    # Results should be valid
    expect_true(!is.na(result_agc$estimate), 
                info = sprintf("%s: AGC estimate should be valid", tc$name))
    expect_true(!is.na(result_akc$estimate), 
                info = sprintf("%s: AKC estimate should be valid", tc$name))
    expect_true(!is.na(result_agc$variance) && result_agc$variance > 0, 
                info = sprintf("%s: AGC variance should be positive", tc$name))
    expect_true(!is.na(result_akc$variance) && result_akc$variance > 0, 
                info = sprintf("%s: AKC variance should be positive", tc$name))
    expect_true(!is.na(result_agc$variance_ind) && result_agc$variance_ind > 0, 
                info = sprintf("%s: AGC independence variance should be positive", tc$name))
    expect_true(!is.na(result_akc$variance_ind) && result_akc$variance_ind > 0, 
                info = sprintf("%s: AKC independence variance should be positive", tc$name))
    
    # P-values should be in [0, 1]
    expect_true(result_agc$p.value >= 0 && result_agc$p.value <= 1,
                info = sprintf("%s: AGC p-value should be in [0,1]", tc$name))
    expect_true(result_akc$p.value >= 0 && result_akc$p.value <= 1,
                info = sprintf("%s: AKC p-value should be in [0,1]", tc$name))
  }
})

test_that("AGC and AKC estimates match Spearman/Kendall for discrete data with few ties", {
  # When discrete data has enough levels relative to n, estimates should still be close
  set.seed(3009)
  n <- 100
  
  # Use enough levels that ties don't dominate
  X <- sample(1:50, n, replace = TRUE) + runif(n, -0.001, 0.001)  # Break ties
  Y <- sample(1:30, n, replace = TRUE) + runif(n, -0.001, 0.001)  # Break ties
  
  result_agc <- acor.test(X, Y, method = "agc")
  result_akc <- acor.test(X, Y, method = "akc")
  
  spearman <- cor(X, Y, method = "spearman")
  kendall <- cor(X, Y, method = "kendall")
  
  # Should be close when ties are minimal
  expect_equal(unname(result_agc$estimate), spearman, tolerance = 0.01,
               info = "AGC should be close to Spearman when ties are minimal")
  expect_equal(unname(result_akc$estimate), kendall, tolerance = 0.01,
               info = "AKC should be close to Kendall when ties are minimal")
})


# ============================================================================
# Test C: Compare binary Y variance to pROC DeLong (if available)
# ============================================================================

# Skip these tests if pROC is not installed
skip_if_no_pROC <- function() {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    skip("pROC package not available")
  }
}

test_that("Binary Y: CMA/CID variance matches DeLong variance (single predictor)", {
  skip_if_no_pROC()
  
  set.seed(4001)
  n <- 200
  X <- rnorm(n)
  Y <- rbinom(n, 1, 0.6)
  
  # Our result
  result_cma <- acor.test(X, Y, method = "cma")
  
  # pROC DeLong
  roc_obj <- pROC::roc(Y, X, direction = '<', quiet = TRUE)
  delong_var <- pROC::var(roc_obj)
  
  # Compare estimates (CMA = AUC for binary Y)
  expect_equal(unname(result_cma$estimate), as.numeric(pROC::auc(roc_obj)), 
               tolerance = 1e-8,
               info = "CMA should equal AUC for binary Y")
  
  # Compare variance
  our_var_scaled <- result_cma$variance / n
  
  expect_equal(our_var_scaled, delong_var, 
               tolerance = 1e-2,
               info = "CMA variance should approximately match DeLong variance")
})

test_that("Binary Y: Two-predictor comparison matches DeLong test", {
  skip_if_no_pROC()
  
  set.seed(4002)
  n <- 200
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  Y <- rbinom(n, 1, 0.6)
  
  # Add some signal to X1
  Y_signal <- Y
  X1 <- X1 + 0.5 * Y_signal
  
  # Our result
  X <- cbind(X1, X2)
  result <- acor.test(X, Y, method = "cma")
  
  # pROC comparison
  roc1 <- pROC::roc(Y, X1,direction = "<", quiet = TRUE)
  roc2 <- pROC::roc(Y, X2,direction = "<", quiet = TRUE)
  delong_test <- pROC::roc.test(roc1, roc2, method = "delong")
  
  # P-values should be similar
  expect_equal(result$p.value, delong_test$p.value, 
               tolerance = 0.2,  # Allow some difference
               info = "Two-predictor test p-value should be similar to DeLong")
})

test_that("Binary Y: Variance estimates are consistent across methods", {
  skip_if_no_pROC()
  
  set.seed(4003)
  n <- 300
  X <- rnorm(n)
  Y <- rbinom(n, 1, 0.5)
  
  # Get all variance estimates
  result_cma <- acor.test(X, Y, method = "cma")
  result_cid <- acor.test(X, Y, method = "cid")
  
  conc <- concordance(Y ~ X)
  roc_obj <- pROC::roc(Y, X, direction = '<', quiet = TRUE)
  delong_var <- pROC::var(roc_obj)
  
  # All estimates should be equal
  expect_equal(result_cma$estimate, result_cid$estimate, tolerance = 1e-10)
  expect_equal(unname(result_cma$estimate), unname(conc$concordance), tolerance = 1e-8)
  expect_equal(unname(result_cma$estimate), as.numeric(pROC::auc(roc_obj)), tolerance = 1e-8)
  
  # Variances should be similar (concordance, DeLong, ours)
  our_var <- result_cma$variance / n
  
  # All three should be within reasonable tolerance of each other
  vars <- c(our_var, conc$var, delong_var)
  max_var <- max(vars)
  min_var <- min(vars)
  relative_range <- (max_var - min_var) / mean(vars)
  
  expect_lt(relative_range, 0.3)
})

test_that("Binary Y with discrete X: CMA/CID variance matches DeLong", {
  skip_if_no_pROC()
  
  # Test with discrete predictors (common in practice)
  test_cases <- list(
    list(name = "5-level predictor", X_levels = 5, seed = 4004),
    list(name = "10-level predictor", X_levels = 10, seed = 4005),
    list(name = "20-level predictor", X_levels = 20, seed = 4006)
  )
  
  n <- 250
  
  for (tc in test_cases) {
    set.seed(tc$seed)
    X <- sample(1:tc$X_levels, n, replace = TRUE)
    Y <- rbinom(n, 1, 0.6)
    
    result_cma <- acor.test(X, Y, method = "cma")
    
    roc_obj <- pROC::roc(Y, X, direction = '<', quiet = TRUE)
    delong_var <- pROC::var(roc_obj)
    
    # Compare estimates
    expect_equal(unname(result_cma$estimate), as.numeric(pROC::auc(roc_obj)), 
                 tolerance = 1e-4,
                 info = sprintf("%s: CMA should equal AUC", tc$name))
    
    # Compare variance
    our_var_scaled <- result_cma$variance / n
    
    expect_equal(our_var_scaled, delong_var, 
                 tolerance =0.05,
                 info = sprintf("%s: CMA variance should match DeLong", tc$name))
    
  }
})


# ============================================================================
# Test D: Compare multivariate (m=2) to pROC DeLong covariance
# ============================================================================

test_that("Binary Y m=2: Covariance structure matches DeLong", {
  skip_if_no_pROC()
  
  set.seed(5001)
  n <- 250
  
  # Create correlated predictors
  X1 <- rnorm(n)
  X2 <- 0.5 * X1 + rnorm(n, sd = sqrt(1 - 0.5^2))
  Y <- rbinom(n, 1, 0.6)
  
  X <- cbind(X1, X2)
  result <- acor.test(X, Y, method = "cma")
  
  # pROC
  roc1 <- pROC::roc(Y, X1, direction = '<', quiet = TRUE)
  roc2 <- pROC::roc(Y, X2, direction = '<', quiet = TRUE)
  
  # Get DeLong covariance
  # pROC::cov computes covariance between two ROC curves
  delong_cov <- pROC::cov(roc1, roc2)
  delong_var1 <- pROC::var(roc1)
  delong_var2 <- pROC::var(roc2)
  
  # Our covariance matrix (scaled)
  our_cov_scaled <- result$variance / n
  
  # Compare diagonal elements
  expect_equal(our_cov_scaled[1, 1], delong_var1, 
               tolerance = 0.05,
               info = "Variance of first predictor should match DeLong")
  
  expect_equal(our_cov_scaled[2, 2], delong_var2, 
               tolerance = 0.05,
               info = "Variance of second predictor should match DeLong")
  
  # Compare off-diagonal (covariance)
  expect_equal(our_cov_scaled[1, 2], delong_cov, 
               tolerance = 0.05 ,  # Relative + absolute tolerance
               info = "Covariance should match DeLong")
})

test_that("Binary Y m=2 with discrete predictors: Covariance matches DeLong", {
  skip_if_no_pROC()
  
  set.seed(5002)
  n <- 300
  
  # Discrete predictors (common in practice, e.g., rating scales)
  X1 <- sample(1:10, n, replace = TRUE)
  X2 <- sample(1:5, n, replace = TRUE)
  Y <- rbinom(n, 1, 0.5)
  
  X <- cbind(X1, X2)
  result <- acor.test(X, Y, method = "cma")
  
  # pROC
  roc1 <- pROC::roc(Y, X1, direction = '<', quiet = TRUE)
  roc2 <- pROC::roc(Y, X2, direction = '<', quiet = TRUE)
  
  delong_var1 <- pROC::var(roc1)
  delong_var2 <- pROC::var(roc2)
  delong_cov <- pROC::cov(roc1, roc2)
  
  # Our covariance matrix (scaled)
  our_cov_scaled <- result$variance / n
  
  # Compare estimates first
  expect_equal(unname(result$estimate[1]), as.numeric(pROC::auc(roc1)), 
               tolerance = 1e-3,
               info = "First predictor AUC should match")
  expect_equal(unname(result$estimate[2]), as.numeric(pROC::auc(roc2)), 
               tolerance = 1e-3,
               info = "Second predictor AUC should match")
  
  # Compare variances (with tolerance for discrete data)
  expect_equal(our_cov_scaled[1, 1], delong_var1, 
               tolerance = 0.05,
               info = "Variance of discrete predictor 1 should match DeLong")
  expect_equal(our_cov_scaled[2, 2], delong_var2, 
               tolerance = 0.05,
               info = "Variance of discrete predictor 2 should match DeLong")
})


# ============================================================================
# Test E: Simulation-based validation of variance estimators
# ============================================================================

test_that("Variance estimators are unbiased (simulation)", {
  set.seed(6001)
  n <- 100
  n_sim <- 500
  
  estimates <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    X <- rnorm(n)
    Y <- rnorm(n)
    estimates[i] <- acor.test(X, Y, method = "agc")$estimate
  }
  
  # Empirical variance
  emp_var <- var(estimates)
  
  # Get one theoretical variance estimate
  X <- rnorm(n)
  Y <- rnorm(n)
  theory_var <- acor.test(X, Y, method = "agc")$variance_ind / n
  
  # They should be in the same ballpark
  ratio <- emp_var / theory_var
  
  expect_gt(ratio, 0.5)
  expect_lt(ratio, 2.0)
})

test_that("Independence variance is correct under true independence", {
  set.seed(6002)
  n <- 150
  n_sim <- 300
  
  # Collect p-values under independence
  p_values_ind <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    X <- rnorm(n)
    Y <- rnorm(n)  # Independent
    p_values_ind[i] <- acor.test(X, Y, method = "agc")$p.value_ind
  }
  
  # P-values should be uniform under null
  # Check that Type I error rate is close to nominal
  alpha <- 0.05
  rejection_rate <- mean(p_values_ind < alpha)
  
  # Should be close to 0.05 (allow some Monte Carlo error)
  expect_gt(rejection_rate, 0.02)
  expect_lt(rejection_rate, 0.10)
  
  # KS test for uniformity
  ks_test <- ks.test(p_values_ind, "punif")
  expect_gt(ks_test$p.value, 0.01)
})


# ============================================================================
# Test F: Edge cases and robustness
# ============================================================================

test_that("Variance estimation works with many ties", {
  set.seed(7001)
  n <- 100
  
  # Highly discrete data
  X <- sample(1:5, n, replace = TRUE)
  Y <- sample(1:3, n, replace = TRUE)
  
  result_agc <- acor.test(X, Y, method = "agc")
  result_akc <- acor.test(X, Y, method = "akc")
  
  # Should not error and should give valid results
  expect_true(!is.na(result_agc$variance))
  expect_true(!is.na(result_agc$variance_ind))
  expect_true(result_agc$variance > 0)
  expect_true(result_agc$variance_ind > 0)
  
  expect_true(!is.na(result_akc$variance))
  expect_true(!is.na(result_akc$variance_ind))
})

test_that("Variance estimation works across discrete configurations", {
  # Test various discrete data configurations matching test_akc_fast.R patterns
  test_cases <- list(
    list(name = "Binary Y", n = 100, seed = 7002,
         gen_X = function(n) rnorm(n),
         gen_Y = function(n) rbinom(n, 1, 0.6)),
    list(name = "Discrete Y (5 levels)", n = 100, seed = 7003,
         gen_X = function(n) rnorm(n),
         gen_Y = function(n) sample(1:5, n, replace = TRUE)),
    list(name = "Discrete X (10 levels)", n = 100, seed = 7004,
         gen_X = function(n) sample(1:10, n, replace = TRUE),
         gen_Y = function(n) rnorm(n)),
    list(name = "Discrete X and Y", n = 100, seed = 7005,
         gen_X = function(n) sample(1:5, n, replace = TRUE),
         gen_Y = function(n) sample(1:3, n, replace = TRUE)),
    list(name = "Large discrete", n = 200, seed = 7006,
         gen_X = function(n) sample(1:20, n, replace = TRUE),
         gen_Y = function(n) sample(1:10, n, replace = TRUE)),
    list(name = "Ties in X only", n = 100, seed = 7007,
         gen_X = function(n) { x <- rnorm(n); x[seq(1, n, by = 5)] <- x[1]; x },
         gen_Y = function(n) rnorm(n) + rnorm(n, sd = 0.0001)),
    list(name = "Many ties in Y", n = 100, seed = 7008,
         gen_X = function(n) rnorm(n),
         gen_Y = function(n) sample(1:3, n, replace = TRUE))
  )
  
  for (tc in test_cases) {
    set.seed(tc$seed)
    X <- tc$gen_X(tc$n)
    Y <- tc$gen_Y(tc$n)
    
    result_agc <- acor.test(X, Y, method = "agc")
    result_akc <- acor.test(X, Y, method = "akc")
    result_cma <- acor.test(X, Y, method = "cma")
    result_cid <- acor.test(X, Y, method = "cid")
    
    # All methods should produce valid results
    for (result in list(result_agc, result_akc, result_cma, result_cid)) {
      expect_true(!is.na(result$estimate), 
                  info = sprintf("%s: estimate should be valid", tc$name))
      expect_true(!is.na(result$variance) && result$variance > 0, 
                  info = sprintf("%s: variance should be positive", tc$name))
      expect_true(!is.na(result$variance_ind) && result$variance_ind > 0, 
                  info = sprintf("%s: variance_ind should be positive", tc$name))
      expect_true(!is.na(result$p.value) && result$p.value >= 0 && result$p.value <= 1, 
                  info = sprintf("%s: p.value should be in [0,1]", tc$name))
      expect_true(!is.na(result$p.value_ind) && result$p.value_ind >= 0 && result$p.value_ind <= 1, 
                  info = sprintf("%s: p.value_ind should be in [0,1]", tc$name))
    }
  }
})

test_that("Variance estimation is stable across sample sizes", {
  set.seed(7009)
  
  sample_sizes <- c(50, 100, 200, 500)
  
  for (n in sample_sizes) {
    X <- rnorm(n)
    Y <- rnorm(n)
    
    result <- acor.test(X, Y, method = "agc")
    
    # Variance / n should decrease approximately as 1/n
    # So variance itself should be roughly constant (asymptotic variance)
    expect_true(!is.na(result$variance), info = paste("n =", n))
    expect_true(result$variance > 0, info = paste("n =", n))
    expect_true(result$variance < 10, info = paste("n =", n))  # Reasonable upper bound
  }
})

test_that("Variance estimation is stable for discrete data across sample sizes", {
  set.seed(7010)
  
  sample_sizes <- c(50, 100, 200, 500)
  
  for (n in sample_sizes) {
    # Discrete data
    X <- sample(1:10, n, replace = TRUE)
    Y <- sample(1:5, n, replace = TRUE)
    
    result_agc <- acor.test(X, Y, method = "agc")
    result_akc <- acor.test(X, Y, method = "akc")
    
    # Should produce valid results
    expect_true(!is.na(result_agc$variance) && result_agc$variance > 0, 
                info = paste("AGC n =", n))
    expect_true(!is.na(result_akc$variance) && result_akc$variance > 0, 
                info = paste("AKC n =", n))
  }
})


# ============================================================================
# Test G: Chi-square test validation for m >= 2
# ============================================================================

test_that("m=2: Chi-square statistic equals z-squared from difference test", {
  # For m=2, the chi-square test should be equivalent to squaring the z-statistic
  # from the difference test: chi_sq = z^2
  
  set.seed(8001)
  n <- 200
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  Y <- rnorm(n)
  
  X <- cbind(X1, X2)
  result <- acor.test(X, Y, method = "agc")
  
  # Manually compute z-statistic for difference (old approach)
  diff <- result$estimate[1] - result$estimate[2]
  var_diff <- result$variance[1, 1] + result$variance[2, 2] - 2 * result$variance[1, 2]
  se_diff <- sqrt(var_diff / n)
  z_stat_manual <- diff / se_diff
  
  # Chi-square should equal z^2
  expect_equal(result$statistic, z_stat_manual^2, 
               tolerance = 1e-10,
               info = "Chi-square statistic should equal z^2 for m=2")
  
  # Degrees of freedom should be 1 for m=2
  expect_equal(result$df, 1,
               info = "Degrees of freedom should be 1 for m=2")
  
  # P-values should match (chi-sq with df=1 vs two-sided z-test)
  p_from_z <- 2 * (1 - pnorm(abs(z_stat_manual)))
  expect_equal(result$p.value, p_from_z, 
               tolerance = 1e-10,
               info = "Chi-square p-value should match two-sided z-test p-value for m=2")
})

test_that("m=2: Chi-square test with independence variance equals z-squared", {
  set.seed(8002)
  n <- 200
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  Y <- rnorm(n)
  
  X <- cbind(X1, X2)
  result <- acor.test(X, Y, method = "agc")
  
  # Manually compute z-statistic for difference using independence variance
  diff <- result$estimate[1] - result$estimate[2]
  var_diff_ind <- result$variance_ind[1, 1] + result$variance_ind[2, 2] - 2 * result$variance_ind[1, 2]
  se_diff_ind <- sqrt(var_diff_ind / n)
  z_stat_ind_manual <- diff / se_diff_ind
  
  # Chi-square_ind should equal z^2
  expect_equal(result$statistic_ind, z_stat_ind_manual^2, 
               tolerance = 1e-10,
               info = "Chi-square_ind statistic should equal z^2 for m=2")
  
  # P-values should match
  p_from_z_ind <- 2 * (1 - pnorm(abs(z_stat_ind_manual)))
  expect_equal(result$p.value_ind, p_from_z_ind, 
               tolerance = 1e-10,
               info = "Chi-square_ind p-value should match two-sided z-test p-value for m=2")
})

test_that("m=2: Chi-square equivalence holds for all methods", {
  set.seed(8003)
  n <- 150
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  Y <- rnorm(n)
  
  X <- cbind(X1, X2)
  
  for (method in c("akc", "agc", "cid", "cma")) {
    result <- acor.test(X, Y, method = method)
    
    # Compute z-statistic manually
    diff <- result$estimate[1] - result$estimate[2]
    var_diff <- result$variance[1, 1] + result$variance[2, 2] - 2 * result$variance[1, 2]
    z_stat_manual <- diff / sqrt(var_diff / n)
    
    expect_equal(result$statistic, z_stat_manual^2, 
                 tolerance = 1e-10,
                 info = sprintf("%s: Chi-square should equal z^2", toupper(method)))
  }
})

test_that("m=2 binary Y: Chi-square test matches pROC DeLong test", {
  skip_if_no_pROC()
  
  set.seed(8004)
  n <- 200
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  Y <- rbinom(n, 1, 0.6)
  
  # Add signal to X1
  X1 <- X1 + 0.5 * Y
  
  X <- cbind(X1, X2)
  result <- acor.test(X, Y, method = "cma")
  
  # pROC DeLong test
  roc1 <- pROC::roc(Y, X1, direction = "<", quiet = TRUE)
  roc2 <- pROC::roc(Y, X2, direction = "<", quiet = TRUE)
  delong_test <- pROC::roc.test(roc1, roc2, method = "delong")
  
  # DeLong uses z-statistic, so z^2 should match our chi-square
  expect_equal(result$statistic, as.numeric(delong_test$statistic)^2,
               tolerance = 0.1,  # Allow some difference due to variance estimation
               info = "Chi-square should approximately equal DeLong z^2")
  
  # P-values should be similar
  expect_equal(result$p.value, delong_test$p.value, 
               tolerance = 0.05,
               info = "P-value should be similar to DeLong p-value")
})


# ============================================================================
# Test H: Tests for m > 2 (three or more predictors)
# ============================================================================

test_that("m=3: Chi-square test runs without error", {
  set.seed(9001)
  n <- 200
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  X3 <- rnorm(n)
  Y <- rnorm(n)
  
  X <- cbind(X1, X2, X3)
  
  # Should run without error for all methods
  for (method in c("akc", "agc", "cid", "cma")) {
    result <- acor.test(X, Y, method = method)
    
    expect_true(!is.na(result$statistic),
                info = sprintf("%s m=3: statistic should be valid", toupper(method)))
    expect_true(!is.na(result$p.value) && result$p.value >= 0 && result$p.value <= 1,
                info = sprintf("%s m=3: p.value should be in [0,1]", toupper(method)))
    expect_true(result$df >= 1,
                info = sprintf("%s m=3: df should be >= 1", toupper(method)))
  }
})

test_that("m=3: Degrees of freedom is correct", {
  # For m=3 predictors, we have 3 pairwise comparisons
  # The contrast matrix L has 3 rows, but rank should be 2 (m-1)
  # because the constraints are linearly dependent
  
  set.seed(9002)
  n <- 200
  X <- matrix(rnorm(n * 3), ncol = 3)
  Y <- rnorm(n)
  
  result <- acor.test(X, Y, method = "agc")
  
  # df should be m-1 = 2 (not 3) because rank of contrast matrix is m-1
  expect_equal(result$df, 2,
               info = "Degrees of freedom for m=3 should be 2 (m-1)")
})

test_that("m=4: Degrees of freedom is m-1", {
  set.seed(9003)
  n <- 300
  X <- matrix(rnorm(n * 4), ncol = 4)
  Y <- rnorm(n)
  
  result <- acor.test(X, Y, method = "agc")
  
  # df should be m-1 = 3
  expect_equal(result$df, 3,
               info = "Degrees of freedom for m=4 should be 3 (m-1)")
})

test_that("m=3: Results table has correct structure", {
  set.seed(9004)
  n <- 200
  X <- matrix(rnorm(n * 3), ncol = 3)
  Y <- rnorm(n)
  
  result <- acor.test(X, Y, method = "agc")
  
  # Results table should have 3 rows (one per predictor)
  expect_equal(nrow(result$results), 3,
               info = "Results table should have 3 rows for m=3")
  
  # Check column names
  expected_cols <- c("predictor", "estimate", "statistic", "statistic_ind", 
                     "p.value", "p.value_ind", "CI_lower", "CI_upper")
  expect_equal(names(result$results), expected_cols,
               info = "Results table should have correct column names")
  
  # Predictor names should be X1, X2, X3
  expect_equal(result$results$predictor, c("X1", "X2", "X3"),
               info = "Predictor names should be X1, X2, X3")
})

test_that("m=3: Individual estimates match marginal tests", {
  set.seed(9005)
  n <- 200
  X <- matrix(rnorm(n * 3), ncol = 3)
  Y <- rnorm(n)
  
  result_multi <- acor.test(X, Y, method = "agc")
  
  # Individual estimates should match single-predictor results
  for (i in 1:3) {
    result_single <- acor.test(X[, i], Y, method = "agc")
    
    expect_equal(result_multi$results$estimate[i], result_single$estimate,
                 tolerance = 1e-10,
                 info = sprintf("X%d estimate should match single-predictor result", i))
  }
})

test_that("m=3: Pairwise differences are correctly computed", {
  set.seed(9006)
  n <- 200
  X <- matrix(rnorm(n * 3), ncol = 3)
  Y <- rnorm(n)
  
  result <- acor.test(X, Y, method = "agc")
  
  # Check pairwise differences
  # Order should be: (1-2), (1-3), (2-3)
  est <- result$estimate
  expected_diffs <- c(est[1] - est[2], est[1] - est[3], est[2] - est[3])
  
  expect_equal(result$pairwise_differences, expected_diffs,
               tolerance = 1e-10,
               info = "Pairwise differences should be correctly computed")
})

test_that("m=3: Contrast matrix has correct structure", {
  set.seed(9007)
  n <- 200
  X <- matrix(rnorm(n * 3), ncol = 3)
  Y <- rnorm(n)
  
  result <- acor.test(X, Y, method = "agc")
  
  # Contrast matrix should be 3x3 for m=3 (3 pairwise comparisons)
  expect_equal(dim(result$contrast_matrix), c(3, 3),
               info = "Contrast matrix should be 3x3 for m=3")
  
  # Each row should sum to 0 (contrast)
  row_sums <- rowSums(result$contrast_matrix)
  expect_equal(row_sums, rep(0, 3),
               info = "Each row of contrast matrix should sum to 0")
})

test_that("m=3 binary Y: Test produces valid results", {
  set.seed(9008)
  n <- 250
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  X3 <- rnorm(n)
  Y <- rbinom(n, 1, 0.6)
  
  # Add signal to X1
  X1 <- X1 + 0.5 * Y
  
  X <- cbind(X1, X2, X3)
  result <- acor.test(X, Y, method = "cma")
  
  # X1 should have highest estimate (has signal)
  expect_true(result$results$estimate[1] > result$results$estimate[2],
              info = "X1 (with signal) should have higher estimate than X2")
  expect_true(result$results$estimate[1] > result$results$estimate[3],
              info = "X1 (with signal) should have higher estimate than X3")
  
  # Global test should likely be significant
  expect_true(result$p.value < 0.1,
              info = "Global test should detect difference when one predictor has signal")
})

test_that("m=3: Type I error is controlled under null", {
  # Under H0 (all correlations equal), p-values should be uniform
  set.seed(9009)
  n <- 150
  n_sim <- 100
  
  p_values <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    # All independent - all correlations should be ~0
    X <- matrix(rnorm(n * 3), ncol = 3)
    Y <- rnorm(n)
    
    result <- acor.test(X, Y, method = "agc")
    p_values[i] <- result$p.value
  }
  
  # Type I error at alpha = 0.05 should be close to 0.05
  rejection_rate <- mean(p_values < 0.05)
  
  expect_true(rejection_rate > 0.01,
              info = "Type I error should be > 1%")
  expect_true(rejection_rate < 0.15,
              info = "Type I error should be < 15%")
})

test_that("m=5: Chi-square test works with many predictors", {
  set.seed(9010)
  n <- 300
  X <- matrix(rnorm(n * 5), ncol = 5)
  Y <- rnorm(n)
  
  result <- acor.test(X, Y, method = "agc")
  
  # Should produce valid results
  
  expect_true(!is.na(result$statistic),
              info = "m=5: statistic should be valid")
  expect_equal(result$df, 4,
               info = "m=5: df should be 4 (m-1)")
  expect_equal(nrow(result$results), 5,
               info = "m=5: results table should have 5 rows")
  expect_equal(length(result$pairwise_differences), 10,
               info = "m=5: should have 10 pairwise differences")
})

