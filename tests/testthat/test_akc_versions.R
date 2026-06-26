# ============================================================================
# Tests: AKC Kernel Consistency — Auto-dispatch vs Original Reference
# ============================================================================
#
# Verify that auto-dispatched optimized kernels produce numerically identical
# estimates and variances to the O(n^2) original reference kernels in
# helper-original-kernels.R.
#
# ============================================================================

library(testthat)
library(acor)

TOL <- 1e-10

run_uni <- function(X, Y, IID = TRUE) {
  acor:::compute_akc_variance_auto(X, Y, IID = IID)
}

run_mv <- function(X, Y, IID = TRUE) {
  acor:::compute_akc_multivariate_variance_auto(X, Y, IID = IID)
}


# ============================================================================
# Section 1: Univariate IID
# ============================================================================

test_that("Uni IID | continuous X & Y: auto-dispatch matches original", {
  set.seed(101); n <- 300
  X <- rnorm(n); Y <- rnorm(n)

  r_orig <- Sigma_akc(X, Y)
  r_auto <- run_uni(X, Y)

  expect_equal(r_auto$akc,     r_orig$akc,     tolerance = TOL)
  expect_equal(r_auto$var,     r_orig$var,     tolerance = TOL)
  expect_equal(r_auto$var_ind, r_orig$var_ind, tolerance = TOL)
})

test_that("Uni IID | binary Y: auto-dispatch matches original", {
  set.seed(102); n <- 300
  X <- rnorm(n); Y <- rbinom(n, 1, 0.5)

  r_orig <- Sigma_akc(X, Y)
  r_auto <- run_uni(X, Y)

  expect_equal(r_auto$akc,     r_orig$akc,     tolerance = TOL)
  expect_equal(r_auto$var,     r_orig$var,     tolerance = TOL)
  expect_equal(r_auto$var_ind, r_orig$var_ind, tolerance = TOL)
})

test_that("Uni IID | discrete / heavy ties: auto-dispatch matches original", {
  set.seed(104); n <- 300
  X <- sample(1:5, n, replace = TRUE)
  Y <- sample(1:3, n, replace = TRUE)

  r_orig <- Sigma_akc(X, Y)
  r_auto <- run_uni(X, Y)

  expect_equal(r_auto$akc,     r_orig$akc,     tolerance = TOL)
  expect_equal(r_auto$var,     r_orig$var,     tolerance = TOL)
  expect_equal(r_auto$var_ind, r_orig$var_ind, tolerance = TOL)
})


# ============================================================================
# Section 2: Univariate HAC
# ============================================================================

test_that("Uni HAC | continuous X & Y: auto-dispatch matches original", {
  set.seed(201); n <- 300
  X <- rnorm(n); Y <- rnorm(n)

  r_orig <- Sigma_akc_ts(X, Y)
  r_auto <- run_uni(X, Y, IID = FALSE)

  expect_equal(r_auto$akc,     r_orig$akc,     tolerance = TOL)
  expect_equal(r_auto$var,     r_orig$var,     tolerance = TOL)
  expect_equal(r_auto$var_ind, r_orig$var_ind, tolerance = TOL)
})


# ============================================================================
# Section 3: Multivariate IID (m = 3)
# ============================================================================

test_that("MV IID m=3 | continuous X & Y: auto-dispatch matches original", {
  set.seed(301); n <- 300
  X <- matrix(rnorm(n * 3), ncol = 3); Y <- rnorm(n)

  r_orig <- Sigma_akc_multivariate(X, Y)
  r_auto <- run_mv(X, Y)

  expect_equal(r_auto$akc_vector, r_orig$akc_vector, tolerance = TOL)
  expect_equal(r_auto$Sigma,      r_orig$Sigma,      tolerance = TOL)
  expect_equal(r_auto$Sigma_ind,  r_orig$Sigma_ind,  tolerance = TOL)
})


# ============================================================================
# Section 4: Multivariate HAC (m = 3)
# ============================================================================

test_that("MV HAC m=3 | continuous X & Y: auto-dispatch matches original", {
  set.seed(401); n <- 300
  X <- matrix(rnorm(n * 3), ncol = 3); Y <- rnorm(n)

  r_orig <- Sigma_akc_multivariate_ts(X, Y)
  r_auto <- run_mv(X, Y, IID = FALSE)

  expect_equal(r_auto$akc_vector, r_orig$akc_vector, tolerance = TOL)
  expect_equal(r_auto$Sigma,      r_orig$Sigma,      tolerance = TOL)
  expect_equal(r_auto$Sigma_ind,  r_orig$Sigma_ind,  tolerance = TOL)
})
