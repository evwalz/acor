# ============================================================================
# Tests for HAC (time series) variance paths and variance convergence
# ============================================================================
#
# Coverage gaps addressed:
#   1. Univariate HAC via acor.test() public API (IID = FALSE)
#   2. Multivariate HAC via acor.test() (IID = FALSE, m >= 2)
#   3. HAC version consistency (original vs v2, original vs binary)
#   4. Multivariate HAC version consistency
#   5. Variance / n convergence as n grows
#   6. Binary kernel dispatch: _binary path matches general path
#
# ============================================================================

library(testthat)
library(acor)


# ============================================================================
# Section 1: Univariate HAC through acor.test() public API
# ============================================================================

test_that("acor.test with IID=FALSE runs without error (AKC, continuous)", {
  set.seed(10001)
  n <- 200
  X <- rnorm(n)
  Y <- rnorm(n)
  
  result <- acor.test(X, Y, method = "akc", IID = FALSE)
  
  expect_true(!is.na(result$estimate))
  expect_true(!is.na(result$variance) && result$variance > 0)
  expect_true(!is.na(result$variance_ind) && result$variance_ind > 0)
  expect_true(!is.na(result$p.value) && result$p.value >= 0 && result$p.value <= 1)
  expect_true(!is.na(result$p.value_ind) && result$p.value_ind >= 0 && result$p.value_ind <= 1)
  expect_false(result$IID)
})

test_that("acor.test with IID=FALSE runs without error (AGC, continuous)", {
  set.seed(10002)
  n <- 200
  X <- rnorm(n)
  Y <- rnorm(n)
  
  result <- acor.test(X, Y, method = "agc", IID = FALSE)
  
  expect_true(!is.na(result$estimate))
  expect_true(!is.na(result$variance) && result$variance > 0)
  expect_true(!is.na(result$variance_ind) && result$variance_ind > 0)
  expect_false(result$IID)
})

test_that("acor.test with IID=FALSE runs without error (CID, continuous)", {
  set.seed(10003)
  n <- 200
  X <- rnorm(n)
  Y <- rnorm(n)
  
  result <- acor.test(X, Y, method = "cid", IID = FALSE)
  
  expect_true(!is.na(result$estimate))
  expect_true(!is.na(result$variance) && result$variance > 0)
  expect_false(result$IID)
})

test_that("acor.test with IID=FALSE runs without error (CMA, continuous)", {
  set.seed(10004)
  n <- 200
  X <- rnorm(n)
  Y <- rnorm(n)
  
  result <- acor.test(X, Y, method = "cma", IID = FALSE)
  
  expect_true(!is.na(result$estimate))
  expect_true(!is.na(result$variance) && result$variance > 0)
  expect_false(result$IID)
})

test_that("HAC estimate equals IID estimate (point estimate unchanged)", {
  set.seed(10005)
  n <- 200
  X <- rnorm(n)
  Y <- rnorm(n)
  
  for (method in c("akc", "agc", "cid", "cma")) {
    r_iid <- acor.test(X, Y, method = method, IID = TRUE)
    r_hac <- acor.test(X, Y, method = method, IID = FALSE)
    
    expect_equal(r_hac$estimate, r_iid$estimate, tolerance = 1e-10,
                 info = sprintf("%s: HAC and IID should have same point estimate", method))
  }
})

test_that("HAC variance differs from IID variance for autocorrelated data", {
  set.seed(10006)
  n <- 300
  
  # Generate autocorrelated data (AR(1) with rho = 0.8)
  e_x <- arima.sim(model = list(ar = 0.8), n = n)
  e_y <- arima.sim(model = list(ar = 0.8), n = n)
  X <- as.numeric(e_x)
  Y <- as.numeric(e_y)
  
  for (method in c("akc", "agc")) {
    r_iid <- acor.test(X, Y, method = method, IID = TRUE)
    r_hac <- acor.test(X, Y, method = method, IID = FALSE)
    
    # HAC variance should generally be larger than IID for autocorrelated data
    # (not guaranteed in every sample, but the point is they should differ)
    expect_false(isTRUE(all.equal(r_hac$variance, r_iid$variance)),
                 info = sprintf("%s: HAC and IID variances should differ for autocorrelated data",
                                method))
  }
})

test_that("HAC works with binary Y", {
  set.seed(10007)
  n <- 200
  X <- rnorm(n)
  Y <- rbinom(n, 1, 0.6)
  
  for (method in c("akc", "agc", "cid", "cma")) {
    result <- acor.test(X, Y, method = method, IID = FALSE)
    
    expect_true(!is.na(result$estimate),
                info = sprintf("%s binary HAC: estimate valid", method))
    expect_true(!is.na(result$variance) && result$variance > 0,
                info = sprintf("%s binary HAC: variance positive", method))
    expect_true(!is.na(result$p.value) && result$p.value >= 0 && result$p.value <= 1,
                info = sprintf("%s binary HAC: p-value in [0,1]", method))
  }
})

test_that("HAC works with discrete data", {
  set.seed(10008)
  n <- 200
  
  test_cases <- list(
    list(name = "Discrete X (10), continuous Y",
         X = sample(1:10, n, replace = TRUE), Y = rnorm(n)),
    list(name = "Continuous X, discrete Y (5)",
         X = rnorm(n), Y = sample(1:5, n, replace = TRUE)),
    list(name = "Discrete X (5), discrete Y (3)",
         X = sample(1:5, n, replace = TRUE), Y = sample(1:3, n, replace = TRUE))
  )
  
  for (tc in test_cases) {
    for (method in c("akc", "agc")) {
      result <- acor.test(tc$X, tc$Y, method = method, IID = FALSE)
      
      expect_true(!is.na(result$variance) && result$variance > 0,
                  info = sprintf("%s %s: HAC variance positive", method, tc$name))
      expect_true(!is.na(result$p.value),
                  info = sprintf("%s %s: HAC p-value valid", method, tc$name))
    }
  }
})

test_that("HAC works with negative correlation", {
  set.seed(10009)
  n <- 200
  X <- rnorm(n)
  Y <- -0.7 * X + rnorm(n, sd = 0.5)
  
  for (method in c("akc", "agc")) {
    result <- acor.test(X, Y, method = method, IID = FALSE)
    
    expect_true(result$estimate < 0,
                info = sprintf("%s: estimate should be negative for negative association", method))
    expect_true(result$variance > 0,
                info = sprintf("%s: HAC variance should be positive even for negative correlation", method))
  }
})


# ============================================================================
# Section 2: Multivariate HAC through acor.test()
# ============================================================================

test_that("m=2: acor.test with IID=FALSE runs for all methods", {
  set.seed(20001)
  n <- 200
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  Y <- rnorm(n)
  X <- cbind(X1, X2)
  
  for (method in c("akc", "agc", "cid", "cma")) {
    result <- acor.test(X, Y, method = method, IID = FALSE)
    
    expect_true(!is.na(result$statistic),
                info = sprintf("%s m=2 HAC: statistic valid", method))
    expect_true(!is.na(result$p.value) && result$p.value >= 0 && result$p.value <= 1,
                info = sprintf("%s m=2 HAC: p-value in [0,1]", method))
    expect_equal(result$df, 1,
                 info = sprintf("%s m=2 HAC: df should be 1", method))
    expect_false(result$IID)
  }
})

test_that("m=3: acor.test with IID=FALSE runs correctly", {
  set.seed(20002)
  n <- 200
  X <- matrix(rnorm(n * 3), ncol = 3)
  Y <- rnorm(n)
  
  for (method in c("akc", "agc")) {
    result <- acor.test(X, Y, method = method, IID = FALSE)
    
    expect_true(!is.na(result$statistic),
                info = sprintf("%s m=3 HAC: statistic valid", method))
    expect_equal(result$df, 2,
                 info = sprintf("%s m=3 HAC: df should be 2", method))
    expect_equal(nrow(result$results), 3,
                 info = sprintf("%s m=3 HAC: results table has 3 rows", method))
  }
})

test_that("m=2 HAC: estimates match IID estimates", {
  set.seed(20003)
  n <- 200
  X <- cbind(rnorm(n), rnorm(n))
  Y <- rnorm(n)
  
  for (method in c("akc", "agc")) {
    r_iid <- acor.test(X, Y, method = method, IID = TRUE)
    r_hac <- acor.test(X, Y, method = method, IID = FALSE)
    
    expect_equal(r_hac$estimate, r_iid$estimate, tolerance = 1e-10,
                 info = sprintf("%s m=2: HAC and IID point estimates should match", method))
  }
})

test_that("m=2 HAC: covariance matrix is positive definite", {
  set.seed(20004)
  n <- 200
  X <- cbind(rnorm(n), rnorm(n))
  Y <- rnorm(n)
  
  for (method in c("akc", "agc")) {
    result <- acor.test(X, Y, method = method, IID = FALSE)
    
    # Variance is a 2x2 matrix for m=2
    eigenvalues <- eigen(result$variance, only.values = TRUE)$values
    expect_true(all(eigenvalues > 0),
                info = sprintf("%s m=2 HAC: covariance matrix should be positive definite", method))
  }
})

test_that("m=2 HAC binary Y: valid results", {
  set.seed(20005)
  n <- 250
  X1 <- rnorm(n) + 0.5 * rbinom(n, 1, 0.6)
  X2 <- rnorm(n)
  Y <- rbinom(n, 1, 0.6)
  X <- cbind(X1, X2)
  
  for (method in c("akc", "agc", "cid", "cma")) {
    result <- acor.test(X, Y, method = method, IID = FALSE)
    
    expect_true(!is.na(result$statistic),
                info = sprintf("%s m=2 binary HAC: statistic valid", method))
    expect_true(!is.na(result$p.value),
                info = sprintf("%s m=2 binary HAC: p-value valid", method))
  }
})

test_that("m=2 HAC: chi-square equals z-squared", {
  # Same structural test as in test_acor_test0.R but with IID=FALSE
  set.seed(20006)
  n <- 200
  X <- cbind(rnorm(n), rnorm(n))
  Y <- rnorm(n)
  
  for (method in c("akc", "agc")) {
    result <- acor.test(X, Y, method = method, IID = FALSE)
    
    diff <- result$estimate[1] - result$estimate[2]
    var_diff <- result$variance[1, 1] + result$variance[2, 2] - 2 * result$variance[1, 2]
    z_stat_manual <- diff / sqrt(var_diff / n)
    
    expect_equal(result$statistic, z_stat_manual^2, tolerance = 1e-10,
                 info = sprintf("%s m=2 HAC: chi-square should equal z^2", method))
  }
})


# ============================================================================
# Section 3: HAC version consistency (internal Sigma functions)
# ============================================================================
# These test that original, v2, and binary HAC kernels give identical results.
# Uses internal functions via acor:::.

test_that("AKC univariate HAC: all versions agree (continuous)", {
  set.seed(30001)
  n <- 200
  X <- rnorm(n)
  Y <- rnorm(n)
  
  r_orig <- acor:::compute_akc_variance_auto(X, Y, IID = FALSE, version = "original")
  r_v2   <- acor:::compute_akc_variance_auto(X, Y, IID = FALSE, version = "v2")
  
  expect_equal(r_v2$akc,     r_orig$akc,     tolerance = 1e-10)
  expect_equal(r_v2$var,     r_orig$var,     tolerance = 1e-10)
  expect_equal(r_v2$var_ind, r_orig$var_ind, tolerance = 1e-10)
})

test_that("AKC univariate HAC: all versions agree (discrete)", {
  set.seed(30002)
  n <- 200
  X <- sample(1:10, n, replace = TRUE)
  Y <- sample(1:5, n, replace = TRUE)
  
  r_orig <- acor:::compute_akc_variance_auto(X, Y, IID = FALSE, version = "original")
  r_v1   <- acor:::compute_akc_variance_auto(X, Y, IID = FALSE, version = "v1")
  r_v2   <- acor:::compute_akc_variance_auto(X, Y, IID = FALSE, version = "v2")
  
  expect_equal(r_v1$var,     r_orig$var,     tolerance = 1e-10)
  expect_equal(r_v2$var,     r_orig$var,     tolerance = 1e-10)
  expect_equal(r_v1$var_ind, r_orig$var_ind, tolerance = 1e-10)
  expect_equal(r_v2$var_ind, r_orig$var_ind, tolerance = 1e-10)
})

test_that("AGC univariate HAC: original vs v2 agree (continuous)", {
  set.seed(30003)
  n <- 200
  X <- rnorm(n)
  Y <- rnorm(n)
  y_rank <- rank(Y, ties.method = "average")
  x_rank <- rank(X, ties.method = "average")
  
  r_orig <- acor:::compute_agc_variance_auto(y_rank, x_rank, IID = FALSE, version = "original")
  r_v2   <- acor:::compute_agc_variance_auto(y_rank, x_rank, IID = FALSE, version = "v2")
  
  expect_equal(r_v2$agc,     r_orig$agc,     tolerance = 1e-10)
  expect_equal(r_v2$var,     r_orig$var,     tolerance = 1e-10)
  expect_equal(r_v2$var_ind, r_orig$var_ind, tolerance = 1e-10)
})

test_that("AGC univariate HAC: original vs v2 agree (discrete)", {
  set.seed(30004)
  n <- 200
  X <- sample(1:10, n, replace = TRUE)
  Y <- sample(1:5, n, replace = TRUE)
  y_rank <- rank(Y, ties.method = "average")
  x_rank <- rank(X, ties.method = "average")
  
  r_orig <- acor:::compute_agc_variance_auto(y_rank, x_rank, IID = FALSE, version = "original")
  r_v2   <- acor:::compute_agc_variance_auto(y_rank, x_rank, IID = FALSE, version = "v2")
  
  expect_equal(r_v2$var,     r_orig$var,     tolerance = 1e-10)
  expect_equal(r_v2$var_ind, r_orig$var_ind, tolerance = 1e-10)
})


# ============================================================================
# Section 4: Multivariate HAC version consistency
# ============================================================================

test_that("AKC multivariate HAC m=3: all versions agree (continuous)", {
  set.seed(40001)
  n <- 200
  X <- matrix(rnorm(n * 3), ncol = 3)
  Y <- rnorm(n)
  
  r_orig <- acor:::compute_akc_multivariate_variance_auto(X, Y, IID = FALSE, version = "original")
  r_v2   <- acor:::compute_akc_multivariate_variance_auto(X, Y, IID = FALSE, version = "v2")
  
  expect_equal(r_v2$akc_vector, r_orig$akc_vector, tolerance = 1e-10)
  expect_equal(r_v2$Sigma,      r_orig$Sigma,      tolerance = 1e-10)
  expect_equal(r_v2$Sigma_ind,  r_orig$Sigma_ind,  tolerance = 1e-10)
})

test_that("AKC multivariate HAC m=3: all versions agree (discrete)", {
  set.seed(40002)
  n <- 200
  X <- matrix(sample(1:10, n * 3, replace = TRUE), ncol = 3)
  Y <- sample(1:5, n, replace = TRUE)
  
  r_orig <- acor:::compute_akc_multivariate_variance_auto(X, Y, IID = FALSE, version = "original")
  r_v1   <- acor:::compute_akc_multivariate_variance_auto(X, Y, IID = FALSE, version = "v1")
  r_v2   <- acor:::compute_akc_multivariate_variance_auto(X, Y, IID = FALSE, version = "v2")
  
  expect_equal(r_v1$Sigma,     r_orig$Sigma,     tolerance = 1e-10)
  expect_equal(r_v2$Sigma,     r_orig$Sigma,     tolerance = 1e-10)
  expect_equal(r_v1$Sigma_ind, r_orig$Sigma_ind, tolerance = 1e-10)
  expect_equal(r_v2$Sigma_ind, r_orig$Sigma_ind, tolerance = 1e-10)
})

test_that("AGC multivariate HAC m=3: original vs v2 agree (continuous)", {
  set.seed(40003)
  n <- 200
  X <- matrix(rnorm(n * 3), ncol = 3)
  Y <- rnorm(n)
  y_rank <- rank(Y, ties.method = "average")
  xarray_ranks <- matrix(0, nrow = 3, ncol = n)
  for (j in 1:3) xarray_ranks[j, ] <- rank(X[, j], ties.method = "average")
  
  r_orig <- acor:::compute_agc_multivariate_variance_auto(y_rank, xarray_ranks,
                                                          IID = FALSE, version = "original")
  r_v2   <- acor:::compute_agc_multivariate_variance_auto(y_rank, xarray_ranks,
                                                          IID = FALSE, version = "v2")
  
  expect_equal(r_v2$agc_vector, r_orig$agc_vector, tolerance = 1e-10)
  expect_equal(r_v2$Sigma,      r_orig$Sigma,      tolerance = 1e-10)
  expect_equal(r_v2$Sigma_ind,  r_orig$Sigma_ind,  tolerance = 1e-10)
})

test_that("AGC multivariate HAC m=3: original vs v2 agree (binary Y)", {
  set.seed(40004)
  n <- 200
  X <- matrix(rnorm(n * 3), ncol = 3)
  Y <- rbinom(n, 1, 0.5)
  y_rank <- rank(Y, ties.method = "average")
  xarray_ranks <- matrix(0, nrow = 3, ncol = n)
  for (j in 1:3) xarray_ranks[j, ] <- rank(X[, j], ties.method = "average")
  
  r_orig <- acor:::compute_agc_multivariate_variance_auto(y_rank, xarray_ranks,
                                                          IID = FALSE, version = "original")
  r_v2   <- acor:::compute_agc_multivariate_variance_auto(y_rank, xarray_ranks,
                                                          IID = FALSE, version = "v2")
  
  expect_equal(r_v2$Sigma,     r_orig$Sigma,     tolerance = 1e-10)
  expect_equal(r_v2$Sigma_ind, r_orig$Sigma_ind, tolerance = 1e-10)
})

test_that("AKC multivariate HAC m=5: original vs v2 agree", {
  set.seed(40005)
  n <- 200
  X <- matrix(rnorm(n * 5), ncol = 5)
  Y <- rnorm(n)
  
  r_orig <- acor:::compute_akc_multivariate_variance_auto(X, Y, IID = FALSE, version = "original")
  r_v2   <- acor:::compute_akc_multivariate_variance_auto(X, Y, IID = FALSE, version = "v2")
  
  expect_equal(r_v2$Sigma,     r_orig$Sigma,     tolerance = 1e-10)
  expect_equal(r_v2$Sigma_ind, r_orig$Sigma_ind, tolerance = 1e-10)
})


# ============================================================================
# Section 5: Binary kernel dispatch — binary path matches general path
# ============================================================================

test_that("AKC: binary kernel matches general kernel for binary Y (IID)", {
  set.seed(50001)
  n <- 200
  X <- rnorm(n)
  Y <- rbinom(n, 1, 0.6)
  
  r_orig   <- acor:::compute_akc_variance_auto(X, Y, IID = TRUE, version = "original")
  r_binary <- acor:::compute_akc_variance_auto(X, Y, IID = TRUE, version = "binary")
  
  expect_equal(r_binary$akc,     r_orig$akc,     tolerance = 1e-10)
  expect_equal(r_binary$var,     r_orig$var,     tolerance = 1e-10)
  expect_equal(r_binary$var_ind, r_orig$var_ind, tolerance = 1e-10)
})

test_that("AKC: binary kernel matches general kernel for binary Y (HAC)", {
  set.seed(50002)
  n <- 200
  X <- rnorm(n)
  Y <- rbinom(n, 1, 0.5)
  
  r_orig   <- acor:::compute_akc_variance_auto(X, Y, IID = FALSE, version = "original")
  r_binary <- acor:::compute_akc_variance_auto(X, Y, IID = FALSE, version = "binary")
  
  expect_equal(r_binary$akc,     r_orig$akc,     tolerance = 1e-10)
  expect_equal(r_binary$var,     r_orig$var,     tolerance = 1e-10)
  expect_equal(r_binary$var_ind, r_orig$var_ind, tolerance = 1e-10)
})

test_that("AGC: binary kernel matches general kernel for binary Y (IID)", {
  set.seed(50003)
  n <- 200
  X <- rnorm(n)
  Y <- rbinom(n, 1, 0.6)
  y_rank <- rank(Y, ties.method = "average")
  x_rank <- rank(X, ties.method = "average")
  
  r_orig   <- acor:::compute_agc_variance_auto(y_rank, x_rank, IID = TRUE, version = "original")
  r_binary <- acor:::compute_agc_variance_auto(y_rank, x_rank, IID = TRUE, version = "binary")
  
  expect_equal(r_binary$agc,     r_orig$agc,     tolerance = 1e-10)
  expect_equal(r_binary$var,     r_orig$var,     tolerance = 1e-10)
  expect_equal(r_binary$var_ind, r_orig$var_ind, tolerance = 1e-10)
})

test_that("AGC: binary kernel matches general kernel for binary Y (HAC)", {
  set.seed(50004)
  n <- 200
  X <- rnorm(n)
  Y <- rbinom(n, 1, 0.5)
  y_rank <- rank(Y, ties.method = "average")
  x_rank <- rank(X, ties.method = "average")
  
  r_orig   <- acor:::compute_agc_variance_auto(y_rank, x_rank, IID = FALSE, version = "original")
  r_binary <- acor:::compute_agc_variance_auto(y_rank, x_rank, IID = FALSE, version = "binary")
  
  expect_equal(r_binary$agc,     r_orig$agc,     tolerance = 1e-10)
  expect_equal(r_binary$var,     r_orig$var,     tolerance = 1e-10)
  expect_equal(r_binary$var_ind, r_orig$var_ind, tolerance = 1e-10)
})

test_that("AKC multivariate: binary kernel matches general for binary Y (IID)", {
  set.seed(50005)
  n <- 200
  X <- matrix(rnorm(n * 3), ncol = 3)
  Y <- rbinom(n, 1, 0.6)
  
  r_orig   <- acor:::compute_akc_multivariate_variance_auto(X, Y, IID = TRUE, version = "original")
  r_binary <- acor:::compute_akc_multivariate_variance_auto(X, Y, IID = TRUE, version = "binary")
  
  expect_equal(r_binary$Sigma,     r_orig$Sigma,     tolerance = 1e-10)
  expect_equal(r_binary$Sigma_ind, r_orig$Sigma_ind, tolerance = 1e-10)
})

test_that("AKC multivariate: binary kernel matches general for binary Y (HAC)", {
  set.seed(50006)
  n <- 200
  X <- matrix(rnorm(n * 3), ncol = 3)
  Y <- rbinom(n, 1, 0.5)
  
  r_orig   <- acor:::compute_akc_multivariate_variance_auto(X, Y, IID = FALSE, version = "original")
  r_binary <- acor:::compute_akc_multivariate_variance_auto(X, Y, IID = FALSE, version = "binary")
  
  expect_equal(r_binary$Sigma,     r_orig$Sigma,     tolerance = 1e-10)
  expect_equal(r_binary$Sigma_ind, r_orig$Sigma_ind, tolerance = 1e-10)
})


# ============================================================================
# Section 6: Variance / n convergence — SE shrinks with n
# ============================================================================

test_that("IID variance: SE decreases as n grows (AGC, continuous)", {
  se_values <- numeric(4)
  sample_sizes <- c(100, 200, 500, 1000)
  
  for (i in seq_along(sample_sizes)) {
    set.seed(60000 + i)
    n <- sample_sizes[i]
    X <- rnorm(n)
    Y <- 0.5 * X + rnorm(n, sd = 0.8)
    
    result <- acor.test(X, Y, method = "agc", IID = TRUE)
    se_values[i] <- sqrt(result$variance / n)
  }
  
  # SE should be strictly decreasing as n grows
  for (i in 2:length(se_values)) {
    expect_true(se_values[i] < se_values[i - 1],
                info = sprintf("SE at n=%d should be less than SE at n=%d",
                               sample_sizes[i], sample_sizes[i - 1]))
  }
})

test_that("IID variance: SE decreases as n grows (AKC, continuous)", {
  se_values <- numeric(4)
  sample_sizes <- c(100, 200, 500, 1000)
  
  for (i in seq_along(sample_sizes)) {
    set.seed(60100 + i)
    n <- sample_sizes[i]
    X <- rnorm(n)
    Y <- 0.5 * X + rnorm(n, sd = 0.8)
    
    result <- acor.test(X, Y, method = "akc", IID = TRUE)
    se_values[i] <- sqrt(result$variance / n)
  }
  
  for (i in 2:length(se_values)) {
    expect_true(se_values[i] < se_values[i - 1],
                info = sprintf("AKC SE at n=%d should be less than SE at n=%d",
                               sample_sizes[i], sample_sizes[i - 1]))
  }
})

test_that("IID variance: asymptotic variance stabilises as n grows", {
  # The asymptotic variance (not divided by n) should converge to a constant
  var_values <- numeric(4)
  sample_sizes <- c(200, 500, 1000, 2000)
  
  for (i in seq_along(sample_sizes)) {
    set.seed(60200 + i)
    n <- sample_sizes[i]
    X <- rnorm(n)
    Y <- 0.5 * X + rnorm(n, sd = 0.8)
    
    result <- acor.test(X, Y, method = "agc", IID = TRUE)
    var_values[i] <- result$variance
  }
  
  # All asymptotic variances should be in a reasonable range of each other
  # (not a strict test, but checks they don't diverge or collapse)
  cv <- sd(var_values) / mean(var_values)
  expect_true(cv < 0.5,
              info = "Coefficient of variation of asymptotic variance should be small across n")
})

test_that("HAC variance: SE decreases as n grows (AGC)", {
  se_values <- numeric(3)
  sample_sizes <- c(200, 500, 1000)
  
  for (i in seq_along(sample_sizes)) {
    set.seed(60300 + i)
    n <- sample_sizes[i]
    X <- rnorm(n)
    Y <- 0.5 * X + rnorm(n, sd = 0.8)
    
    result <- acor.test(X, Y, method = "agc", IID = FALSE)
    se_values[i] <- sqrt(result$variance / n)
  }
  
  for (i in 2:length(se_values)) {
    expect_true(se_values[i] < se_values[i - 1],
                info = sprintf("HAC SE at n=%d should be less than SE at n=%d",
                               sample_sizes[i], sample_sizes[i - 1]))
  }
})

test_that("Independence variance: SE decreases as n grows", {
  se_values <- numeric(4)
  sample_sizes <- c(100, 200, 500, 1000)
  
  for (i in seq_along(sample_sizes)) {
    set.seed(60400 + i)
    n <- sample_sizes[i]
    X <- rnorm(n)
    Y <- rnorm(n)  # independent
    
    result <- acor.test(X, Y, method = "agc", IID = TRUE)
    se_values[i] <- sqrt(result$variance_ind / n)
  }
  
  for (i in 2:length(se_values)) {
    expect_true(se_values[i] < se_values[i - 1],
                info = sprintf("Independence SE at n=%d should be less than SE at n=%d",
                               sample_sizes[i], sample_sizes[i - 1]))
  }
})

test_that("Multivariate variance: SE decreases as n grows (m=2)", {
  se_values <- numeric(3)
  sample_sizes <- c(200, 500, 1000)
  
  for (i in seq_along(sample_sizes)) {
    set.seed(60500 + i)
    n <- sample_sizes[i]
    X <- cbind(rnorm(n), rnorm(n))
    Y <- 0.3 * X[, 1] + rnorm(n)
    
    result <- acor.test(X, Y, method = "agc", IID = TRUE)
    # Check SE of first predictor
    se_values[i] <- sqrt(result$variance[1, 1] / n)
  }
  
  for (i in 2:length(se_values)) {
    expect_true(se_values[i] < se_values[i - 1],
                info = sprintf("Multivariate SE at n=%d should be less than SE at n=%d",
                               sample_sizes[i], sample_sizes[i - 1]))
  }
})

test_that("IID variance: empirical variance matches theoretical (simulation)", {
  # More focused version: run many replicates, check that the average
  # estimated SE is close to the empirical SE of the estimates
  set.seed(60600)
  n <- 200
  n_sim <- 300
  
  estimates <- numeric(n_sim)
  avg_se <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    X <- rnorm(n)
    Y <- 0.3 * X + rnorm(n, sd = 0.9)
    
    result <- acor.test(X, Y, method = "agc", IID = TRUE)
    estimates[i] <- result$estimate
    avg_se[i] <- sqrt(result$variance / n)
  }
  
  empirical_se <- sd(estimates)
  mean_estimated_se <- mean(avg_se)
  
  # The ratio should be close to 1
  ratio <- mean_estimated_se / empirical_se
  expect_true(ratio > 0.7, info = "Estimated SE should not drastically underestimate empirical SE")
  expect_true(ratio < 1.5, info = "Estimated SE should not drastically overestimate empirical SE")
})