# ============================================================================
# Tests for edge cases, input validation, Fisher CI, boundary estimates,
# print methods, acor() structure, and CID/CMA transformations
# ============================================================================

library(testthat)
library(acor)

# ============================================================================
# Section 1: Input validation — acor()
# ============================================================================

test_that("acor() rejects Y with only 1 unique value", {
  expect_error(acor(rnorm(10), rep(5, 10), method = "akc"),
               "at least 2 distinct")
})

test_that("acor() rejects mismatched lengths", {
  expect_error(acor(rnorm(10), rnorm(15), method = "akc"),
               "same number of observations")
})

test_that("acor() rejects NA values in X", {
  x <- rnorm(10)
  x[3] <- NA
  expect_error(acor(x, rnorm(10), method = "akc"), "NA")
})

test_that("acor() rejects NA values in Y", {
  y <- rnorm(10)
  y[5] <- NA
  expect_error(acor(rnorm(10), y, method = "akc"), "NA")
})

test_that("acor() rejects non-numeric X", {
  expect_error(acor(letters[1:10], rnorm(10), method = "akc"))
})

test_that("acor() rejects non-numeric Y", {
  expect_error(acor(rnorm(10), letters[1:10], method = "akc"),
               "numeric vector")
})

test_that("acor() rejects multi-column Y data frames", {
  y <- data.frame(a = 1:5, b = 5:1)
  expect_error(acor(1:5, y, method = "akc"), "multi-column data frame")
})

test_that("acor() works with n = 2 (minimum valid)", {
  result <- acor(c(1, 2), c(3, 4), method = "akc")
  expect_true(is.finite(result$estimate))
})

test_that("acor() works with n = 3", {
  result <- acor(c(1, 2, 3), c(1, 3, 2), method = "agc")
  expect_true(is.finite(result$estimate))
})

# ============================================================================
# Section 2: Input validation — acor.test()
# ============================================================================

test_that("acor.test() rejects Y with only 1 unique value", {
  expect_error(acor.test(rnorm(10), rep(5, 10), method = "akc"),
               "at least 2 distinct")
})

test_that("acor.test() rejects mismatched lengths", {
  expect_error(acor.test(rnorm(10), rnorm(15), method = "akc"),
               "same number of observations")
})

test_that("acor.test() rejects NA in data", {
  x <- rnorm(20)
  y <- rnorm(20)
  x[7] <- NA
  expect_error(acor.test(x, y, method = "akc"), "NA")
})

test_that("acor.test() rejects multi-column Y data frames", {
  y <- data.frame(a = 1:5, b = 5:1)
  expect_error(acor.test(1:5, y, method = "akc"), "multi-column data frame")
})

test_that("acor.test() requires at least 3 observations", {
  expect_error(acor.test(c(1, 2), c(1, 2), method = "akc"),
               "At least 3 observations")
})

test_that("acor.test() rejects invalid conf.level", {
  expect_error(acor.test(rnorm(20), rnorm(20), method = "akc", conf.level = 2),
               "conf.level")
})

test_that("acor.test() rejects multivariate pearson inference", {
  X <- cbind(rnorm(20), rnorm(20))
  y <- rnorm(20)
  expect_error(acor.test(X, y, method = "pearson"),
               "Multivariate Pearson inference is not yet implemented")
})

test_that("acor.test() rejects multivariate rho_b inference", {
  X <- cbind(rnorm(20), rnorm(20))
  y <- rnorm(20)
  expect_error(acor.test(X, y, method = "rho_b"),
               "Multivariate rho_b inference is not yet implemented")
})

test_that("acor.test() rejects multivariate tau_b inference", {
  X <- cbind(rnorm(20), rnorm(20))
  y <- rnorm(20)
  expect_error(acor.test(X, y, method = "tau_b"),
               "Multivariate tau_b inference is not yet implemented")
})

test_that("acor.test() rejects multivariate gamma inference", {
  X <- cbind(rnorm(20), rnorm(20))
  y <- rnorm(20)
  expect_error(acor.test(X, y, method = "gamma"),
               "Multivariate gamma inference is not yet implemented")
})

test_that("acor.test() rejects HAC inference with fewer than 3 observations", {
  expect_error(acor.test(c(1, 2), c(1, 2), method = "akc", IID = FALSE),
               "At least 3 observations")
  expect_error(acor.test(cbind(c(1, 2), c(2, 1)), c(1, 2), method = "rho_a", IID = FALSE),
               "At least 3 observations")
})


# ============================================================================
# Section 3: Fisher transformation
# ============================================================================

test_that("Fisher CI differs from standard CI (m = 1, AKC)", {
  set.seed(5001)
  x <- rnorm(200)
  y <- 0.5 * x + rnorm(200, sd = 0.5)
  
  res_std <- acor.test(x, y, method = "akc", fisher = FALSE)
  res_fish <- acor.test(x, y, method = "akc", fisher = TRUE)
  
  # Estimates should be identical
  
  expect_equal(unname(res_fish$estimate), unname(res_std$estimate))
  
  # CIs should differ
  expect_false(isTRUE(all.equal(res_fish$conf.int, res_std$conf.int)),
               info = "Fisher CI should differ from standard CI")
  
  # Fisher CI should be within [-1, 1] for AKC
  expect_true(res_fish$conf.int[1] >= -1)
  expect_true(res_fish$conf.int[2] <= 1)
})

test_that("Fisher CI differs from standard CI (m = 1, CMA)", {
  set.seed(5002)
  x <- rnorm(200)
  y <- 0.5 * x + rnorm(200, sd = 0.5)
  
  res_std <- acor.test(x, y, method = "cma", fisher = FALSE)
  res_fish <- acor.test(x, y, method = "cma", fisher = TRUE)
  
  expect_equal(unname(res_fish$estimate), unname(res_std$estimate))
  expect_false(isTRUE(all.equal(res_fish$conf.int, res_std$conf.int)))
  
  # Fisher CI should be within [0, 1] for CMA
  expect_true(res_fish$conf.int[1] >= 0)
  expect_true(res_fish$conf.int[2] <= 1)
})

test_that("Fisher CI works for multivariate (m >= 2)", {
  set.seed(5003)
  y <- rnorm(150)
  X <- cbind(y + rnorm(150, sd = 0.3), rnorm(150))
  
  res_fish <- acor.test(X, y, method = "agc", fisher = TRUE)
  res_std <- acor.test(X, y, method = "agc", fisher = FALSE)
  
  # Should have results table with CIs
  expect_true("CI_lower" %in% names(res_fish$results))
  expect_true("CI_upper" %in% names(res_fish$results))
  
  # Fisher CIs should differ from standard
  expect_false(isTRUE(all.equal(res_fish$results$CI_lower,
                                res_std$results$CI_lower)))
})

# ============================================================================
# Section 4: Boundary estimates (atanh clamp)
# ============================================================================

test_that("Perfect concordance (AKC near +1) does not produce Inf", {
  x <- 1:20
  y <- 1:20
  result <- acor.test(x, y, method = "akc", fisher = TRUE)
  expect_true(all(is.finite(result$conf.int)))
  expect_true(is.finite(result$statistic))
})

test_that("Perfect discordance (AKC near -1) does not produce Inf", {
  x <- 1:20
  y <- 20:1
  result <- acor.test(x, y, method = "akc", fisher = TRUE)
  expect_true(all(is.finite(result$conf.int)))
  expect_true(is.finite(result$statistic))
})

test_that("Perfect concordance with CID (near 1) does not produce Inf", {
  x <- 1:20
  y <- 1:20
  result <- acor.test(x, y, method = "cid", fisher = TRUE)
  expect_true(all(is.finite(result$conf.int)))
})

test_that("Perfect concordance with CMA (near 1) does not produce Inf", {
  x <- 1:20
  y <- 1:20
  result <- acor.test(x, y, method = "cma", fisher = TRUE)
  expect_true(all(is.finite(result$conf.int)))
})

test_that("Multivariate with near-perfect correlation and Fisher does not error", {
  x1 <- 1:30
  x2 <- rnorm(30)
  y <- 1:30
  X <- cbind(x1, x2)
  result <- acor.test(X, y, method = "akc", fisher = TRUE)
  expect_true(all(is.finite(result$results$CI_lower)))
  expect_true(all(is.finite(result$results$CI_upper)))
})

# ============================================================================
# Section 5: print methods
# ============================================================================

test_that("print.acor_htest works for m = 1 without error", {
  set.seed(6001)
  res <- acor.test(rnorm(50), rnorm(50), method = "akc")
  expect_output(print(res), "AKC test")
  expect_output(print(res), "estimate")
  expect_output(print(res), "statistic")
})

test_that("print.acor_htest shows the correct one-sided alternative", {
  x <- 1:10
  y <- 10:1
  res <- acor.test(x, y, method = "akc", alternative = "less")
  expect_output(print(res), "is less than")
})

test_that("print.acor_htest works for m >= 2 without error", {
  set.seed(6002)
  y <- rnorm(50)
  X <- cbind(rnorm(50), rnorm(50))
  res <- acor.test(X, y, method = "agc")
  expect_output(print(res), "AGC test")
  expect_output(print(res), "Chi-squared")
  expect_output(print(res), "Individual predictors")
})

test_that("print.htest fallback works for m = 1", {
  set.seed(6003)
  res <- acor.test(rnorm(50), rnorm(50), method = "akc")
  # Should not error when calling print.htest directly
  expect_output(stats:::print.htest(res))
})

test_that("m = 1 object inherits from htest", {
  set.seed(6004)
  res <- acor.test(rnorm(50), rnorm(50), method = "akc")
  expect_s3_class(res, "htest")
  expect_s3_class(res, "acor_htest")
})

test_that("m >= 2 object is acor_htest but not htest", {
  set.seed(6005)
  y <- rnorm(50)
  X <- cbind(rnorm(50), rnorm(50))
  res <- acor.test(X, y, method = "akc")
  expect_s3_class(res, "acor_htest")
  expect_false(inherits(res, "htest"))
})

test_that("identical predictor columns give zero pairwise statistic and p-value 1", {
  set.seed(6006)
  y <- rnorm(40)
  x <- rnorm(40)
  X <- cbind(x, x)
  res <- acor.test(X, y, method = "akc")
  expect_equal(res$pairwise_results$difference, 0)
  expect_equal(res$pairwise_results$statistic, 0)
  expect_equal(res$pairwise_results$p.value, 1)
})

# ============================================================================
# Section 6: acor() return structure
# ============================================================================

test_that("acor() returns correct structure for all methods (m = 1)", {
  set.seed(7001)
  x <- rnorm(80)
  y <- rnorm(80)
  
  for (method in c("akc", "agc", "cid", "cma")) {
    res <- acor(x, y, method = method)
    expect_true(is.list(res), info = method)
    expect_true("estimate" %in% names(res), info = method)
    expect_true("method" %in% names(res), info = method)
    expect_equal(res$method, method, info = method)
    expect_length(res$estimate, 1)
    expect_true(is.finite(res$estimate), info = method)
  }
})

test_that("acor() returns correct structure for all methods (m = 3)", {
  set.seed(7002)
  X <- matrix(rnorm(240), ncol = 3)
  y <- rnorm(80)
  
  for (method in c("akc", "agc", "cid", "cma")) {
    res <- acor(X, y, method = method)
    expect_true(is.list(res), info = method)
    expect_length(res$estimate, 3)
    expect_true(all(is.finite(res$estimate)), info = method)
  }
})

# ============================================================================
# Section 7: CID/CMA transformation consistency
# ============================================================================

test_that("CID = (AKC + 1) / 2 for acor() (m = 1)", {
  set.seed(8001)
  x <- rnorm(100)
  y <- 0.5 * x + rnorm(100)
  
  akc <- acor(x, y, method = "akc")$estimate
  cid <- acor(x, y, method = "cid")$estimate
  
  expect_equal(cid, (akc + 1) / 2, tolerance = 1e-12)
})

test_that("CMA = (AGC + 1) / 2 for acor() (m = 1)", {
  set.seed(8002)
  x <- rnorm(100)
  y <- 0.5 * x + rnorm(100)
  
  agc <- acor(x, y, method = "agc")$estimate
  cma <- acor(x, y, method = "cma")$estimate
  
  expect_equal(cma, (agc + 1) / 2, tolerance = 1e-12)
})

test_that("CID = (AKC + 1) / 2 for acor() (m = 3)", {
  set.seed(8003)
  X <- matrix(rnorm(300), ncol = 3)
  y <- rnorm(100)
  
  akcs <- acor(X, y, method = "akc")$estimate
  cids <- acor(X, y, method = "cid")$estimate
  
  expect_equal(cids, (akcs + 1) / 2, tolerance = 1e-12)
})

test_that("CMA = (AGC + 1) / 2 for acor() (m = 3)", {
  set.seed(8004)
  X <- matrix(rnorm(300), ncol = 3)
  y <- rnorm(100)
  
  agcs <- acor(X, y, method = "agc")$estimate
  cmas <- acor(X, y, method = "cma")$estimate
  
  expect_equal(cmas, (agcs + 1) / 2, tolerance = 1e-12)
})

test_that("CID/AKC transformation holds for acor.test() estimates (m = 1)", {
  set.seed(8005)
  x <- rnorm(100)
  y <- 0.5 * x + rnorm(100)
  
  res_akc <- acor.test(x, y, method = "akc")
  res_cid <- acor.test(x, y, method = "cid")
  
  expect_equal(unname(res_cid$estimate), (unname(res_akc$estimate) + 1) / 2,
               tolerance = 1e-12)
})

test_that("CMA/AGC transformation holds for acor.test() estimates (m = 1)", {
  set.seed(8006)
  x <- rnorm(100)
  y <- 0.5 * x + rnorm(100)
  
  res_agc <- acor.test(x, y, method = "agc")
  res_cma <- acor.test(x, y, method = "cma")
  
  expect_equal(unname(res_cma$estimate), (unname(res_agc$estimate) + 1) / 2,
               tolerance = 1e-12)
})

test_that("CID/AKC variance transformation: Var(CID) = Var(AKC) / 4", {
  set.seed(8007)
  x <- rnorm(150)
  y <- 0.5 * x + rnorm(150)
  
  res_akc <- acor.test(x, y, method = "akc")
  res_cid <- acor.test(x, y, method = "cid")
  
  expect_equal(res_cid$variance, res_akc$variance / 4, tolerance = 1e-12)
})

test_that("CMA/AGC variance transformation: Var(CMA) = Var(AGC) / 4", {
  set.seed(8008)
  x <- rnorm(150)
  y <- 0.5 * x + rnorm(150)
  
  res_agc <- acor.test(x, y, method = "agc")
  res_cma <- acor.test(x, y, method = "cma")
  
  expect_equal(res_cma$variance, res_agc$variance / 4, tolerance = 1e-12)
})

# ============================================================================
# Section 8: IID = FALSE multivariate through public API
# ============================================================================

test_that("acor.test() with IID = FALSE works for m >= 2 (AKC)", {
  set.seed(9001)
  y <- rnorm(150)
  X <- cbind(y + rnorm(150, sd = 0.5), rnorm(150))
  
  result <- acor.test(X, y, method = "akc", IID = FALSE)
  expect_true(is.finite(result$statistic))
  expect_true(result$p.value >= 0 && result$p.value <= 1)
  expect_length(result$estimate, 2)
  expect_false(result$IID)
})

test_that("acor.test() with IID = FALSE works for m >= 2 (AGC)", {
  set.seed(9002)
  y <- rnorm(150)
  X <- cbind(y + rnorm(150, sd = 0.5), rnorm(150))
  
  result <- acor.test(X, y, method = "agc", IID = FALSE)
  expect_true(is.finite(result$statistic))
  expect_true(result$p.value >= 0 && result$p.value <= 1)
  expect_false(result$IID)
})

test_that("acor.test() with IID = FALSE works for m >= 2 (CID)", {
  set.seed(9003)
  y <- rnorm(150)
  X <- cbind(y + rnorm(150, sd = 0.5), rnorm(150), rnorm(150))
  
  result <- acor.test(X, y, method = "cid", IID = FALSE)
  expect_true(is.finite(result$statistic))
  expect_length(result$estimate, 3)
})

test_that("acor.test() with IID = FALSE works for m >= 2 (CMA)", {
  set.seed(9004)
  y <- rnorm(150)
  X <- cbind(y + rnorm(150, sd = 0.5), rnorm(150), rnorm(150))
  
  result <- acor.test(X, y, method = "cma", IID = FALSE)
  expect_true(is.finite(result$statistic))
  expect_length(result$estimate, 3)
})


# ============================================================================
# Section 9: Infinitesimal jackknife (IJ) variance — validation & smoke tests
# ============================================================================

test_that("acor.test() rejects variance='ij' for unsupported methods", {
  set.seed(1234)
  x <- rnorm(50); y <- rnorm(50)
  for (m in c("tau_a", "tau_b", "rho_a", "rho_b", "gamma", "pearson")) {
    expect_error(acor.test(x, y, method = m, variance = "ij"),
                 "only supported for methods")
  }
})

test_that("acor.test() with variance='ij' and IID=FALSE works for AKC", {
  set.seed(1235)
  x <- rnorm(200); y <- x + rnorm(200)
  
  res <- acor.test(x, y, method = "akc", variance = "ij", IID = FALSE)
  expect_s3_class(res, "acor_htest")
  expect_true(is.finite(res$statistic))
  expect_true(res$variance > 0)
  expect_true(res$variance_ind > 0)
  expect_false(res$IID)
})

test_that("acor.test() with variance='ij' and IID=FALSE works for AGC", {
  set.seed(1235)
  x <- rnorm(200); y <- x + rnorm(200)
  
  res <- acor.test(x, y, method = "agc", variance = "ij", IID = FALSE)
  expect_true(is.finite(res$statistic))
  expect_true(res$variance > 0)
  expect_false(res$IID)
})

test_that("IJ HAC variance >= IJ IID variance for AKC", {
  set.seed(1240)
  n <- 300
  x <- as.numeric(arima.sim(list(ar = 0.5), n = n))
  y <- x + as.numeric(arima.sim(list(ar = 0.5), n = n))
  
  res_iid <- acor.test(x, y, method = "akc", variance = "ij", IID = TRUE)
  res_hac <- acor.test(x, y, method = "akc", variance = "ij", IID = FALSE)
  
  expect_equal(unname(res_iid$estimate), unname(res_hac$estimate), tolerance = 1e-10)
  expect_true(res_hac$variance >= res_iid$variance * 0.5)
})

test_that("IJ HAC works for multivariate AKC", {
  set.seed(1241)
  y <- rnorm(150)
  X <- cbind(y + rnorm(150, sd = 0.5), rnorm(150))
  
  res <- acor.test(X, y, method = "akc", variance = "ij", IID = FALSE)
  expect_length(res$estimate, 2)
  expect_true(is.matrix(res$variance))
  expect_true(all(diag(res$variance) > 0))
})

test_that("acor.test() with variance='ij' works for multivariate AKC", {
  set.seed(1236)
  y <- rnorm(150)
  X <- cbind(y + rnorm(150, sd = 0.5), rnorm(150), y + rnorm(150))
  
  res <- acor.test(X, y, method = "akc", variance = "ij")
  expect_length(res$estimate, 3)
  expect_true(is.matrix(res$variance))
  expect_equal(dim(res$variance), c(3, 3))
  expect_true(all(diag(res$variance) > 0))
  expect_true(is.finite(res$statistic))
  expect_true(res$p.value >= 0 && res$p.value <= 1)
})

test_that("acor.test() with variance='ij' works for multivariate AGC", {
  set.seed(1237)
  y <- rnorm(150)
  X <- cbind(y + rnorm(150, sd = 0.5), rnorm(150))
  
  res <- acor.test(X, y, method = "agc", variance = "ij")
  expect_length(res$estimate, 2)
  expect_true(is.matrix(res$variance))
  expect_equal(dim(res$variance), c(2, 2))
  expect_true(all(diag(res$variance) > 0))
})

test_that("acor.test() with variance='ij' works for multivariate CID", {
  set.seed(1238)
  y <- rnorm(150)
  X <- cbind(y + rnorm(150, sd = 0.5), rnorm(150))
  
  res <- acor.test(X, y, method = "cid", variance = "ij")
  expect_length(res$estimate, 2)
  expect_true(all(res$estimate > 0 & res$estimate < 1))
})

test_that("IJ multivariate diagonal matches univariate IJ variance", {
  set.seed(1239)
  y <- rnorm(200)
  x1 <- y + rnorm(200, sd = 0.5)
  x2 <- rnorm(200)
  
  res_mv <- acor.test(cbind(x1, x2), y, method = "akc", variance = "ij")
  res_u1 <- acor.test(x1, y, method = "akc", variance = "ij")
  res_u2 <- acor.test(x2, y, method = "akc", variance = "ij")
  
  expect_equal(res_mv$variance[1, 1], res_u1$variance, tolerance = 1e-10)
  expect_equal(res_mv$variance[2, 2], res_u2$variance, tolerance = 1e-10)
})

test_that("acor.test() with variance='ij' works for AKC", {
  set.seed(2001)
  x <- rnorm(200); y <- x + rnorm(200)
  
  res <- acor.test(x, y, method = "akc", variance = "ij")
  expect_s3_class(res, "acor_htest")
  expect_true(is.finite(res$statistic))
  expect_true(res$p.value >= 0 && res$p.value <= 1)
  expect_true(res$variance > 0)
  expect_true(res$variance_ind > 0)
  expect_true(!is.null(res$statistic_ind))
  expect_true(!is.null(res$p.value_ind))
  expect_equal(res$variance_method, "ij")
})

test_that("acor.test() with variance='ij' works for AGC", {
  set.seed(2002)
  x <- rnorm(200); y <- x + rnorm(200)
  
  res <- acor.test(x, y, method = "agc", variance = "ij")
  expect_s3_class(res, "acor_htest")
  expect_true(is.finite(res$statistic))
  expect_true(res$p.value >= 0 && res$p.value <= 1)
  expect_true(res$variance > 0)
  expect_true(res$variance_ind > 0)
  expect_equal(res$variance_method, "ij")
})

test_that("acor.test() with variance='ij' works for CID and CMA", {
  set.seed(2003)
  x <- rnorm(200); y <- x + rnorm(200)
  
  res_cid <- acor.test(x, y, method = "cid", variance = "ij")
  res_cma <- acor.test(x, y, method = "cma", variance = "ij")
  
  expect_true(res_cid$estimate > 0 && res_cid$estimate < 1)
  expect_true(res_cma$estimate > 0 && res_cma$estimate < 1)
  expect_true(res_cid$variance_ind > 0)
  expect_true(res_cma$variance_ind > 0)
})

test_that("IJ and delta give same point estimate for AKC", {
  set.seed(2004)
  x <- rnorm(200); y <- x + rnorm(200)
  
  res_delta <- acor.test(x, y, method = "akc", variance = "delta")
  res_ij    <- acor.test(x, y, method = "akc", variance = "ij")
  
  expect_equal(unname(res_ij$estimate), unname(res_delta$estimate),
               tolerance = 1e-10)
})

test_that("IJ and delta give same point estimate for AGC", {
  set.seed(2005)
  x <- rnorm(200); y <- x + rnorm(200)
  
  res_delta <- acor.test(x, y, method = "agc", variance = "delta")
  res_ij    <- acor.test(x, y, method = "agc", variance = "ij")
  
  expect_equal(unname(res_ij$estimate), unname(res_delta$estimate),
               tolerance = 1e-10)
})

test_that("IJ and delta give identical independence variance for AKC", {
  set.seed(2006)
  x <- rnorm(200); y <- x + rnorm(200)
  
  res_delta <- acor.test(x, y, method = "akc", variance = "delta")
  res_ij    <- acor.test(x, y, method = "akc", variance = "ij")
  
  expect_equal(res_ij$variance_ind, res_delta$variance_ind, tolerance = 1e-12)
})

test_that("IJ and delta give identical independence variance for AGC", {
  set.seed(2007)
  x <- rnorm(200); y <- x + rnorm(200)
  
  res_delta <- acor.test(x, y, method = "agc", variance = "delta")
  res_ij    <- acor.test(x, y, method = "agc", variance = "ij")
  
  expect_equal(res_ij$variance_ind, res_delta$variance_ind, tolerance = 1e-12)
})