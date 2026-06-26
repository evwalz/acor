# ============================================================================
# Tests: AGC Version Consistency (original vs v2)
# ============================================================================
#
# Verify that "original" and "v2" (Fenwick tree) produce numerically
# identical AGC estimates and variances across representative scenarios.
#
# ============================================================================

library(testthat)
library(acor)

TOL <- 1e-10
TOL_PVAL <- 0.01

run_agc_uni <- function(y_rank, x_rank, version, IID = TRUE) {
  if (IID) {
    if (version == "original") {
      Sigma_agc(y_rank, x_rank)
    } else {
      acor:::Sigma_agc_v2(y_rank, x_rank)
    }
  } else {
    if (version == "original") {
      Sigma_agc_ts(y_rank, x_rank)
    } else {
      acor:::Sigma_agc_ts_v2(y_rank, x_rank)
    }
  }
}

run_agc_mv <- function(y_rank, xarray_ranks, version, IID = TRUE) {
  if (IID) {
    if (version == "original") {
      Sigma_agc_multivariate(y_rank, xarray_ranks)
    } else {
      acor:::Sigma_agc_multivariate_v2(y_rank, xarray_ranks)
    }
  } else {
    if (version == "original") {
      Sigma_agc_multivariate_ts(y_rank, xarray_ranks)
    } else {
      acor:::Sigma_agc_multivariate_ts_v2(y_rank, xarray_ranks)
    }
  }
}

make_ranks <- function(X, Y) {
  y_rank <- rank(Y, ties.method = "average")
  if (is.vector(X)) {
    x_rank <- rank(X, ties.method = "average")
    return(list(y_rank = y_rank, x_rank = x_rank))
  }
  m <- ncol(X)
  n <- nrow(X)
  xarray_ranks <- matrix(0, nrow = m, ncol = n)
  for (i in 1:m) {
    xarray_ranks[i, ] <- rank(X[, i], ties.method = "average")
  }
  list(y_rank = y_rank, xarray_ranks = xarray_ranks)
}


# ============================================================================
# Section 1: Univariate IID
# ============================================================================

test_that("AGC Uni IID | continuous X & Y: original vs v2 agree", {
  set.seed(101); n <- 300
  X <- rnorm(n); Y <- rnorm(n)
  rk <- make_ranks(X, Y)

  r_orig <- run_agc_uni(rk$y_rank, rk$x_rank, "original")
  r_v2   <- run_agc_uni(rk$y_rank, rk$x_rank, "v2")

  expect_equal(r_v2$agc,     r_orig$agc,     tolerance = TOL)
  expect_equal(r_v2$var,     r_orig$var,     tolerance = TOL)
  expect_equal(r_v2$var_ind, r_orig$var_ind, tolerance = TOL)
})

test_that("AGC Uni IID | binary Y: original vs v2 agree", {
  set.seed(102); n <- 300
  X <- rnorm(n); Y <- rbinom(n, 1, 0.5)
  rk <- make_ranks(X, Y)

  r_orig <- run_agc_uni(rk$y_rank, rk$x_rank, "original")
  r_v2   <- run_agc_uni(rk$y_rank, rk$x_rank, "v2")

  expect_equal(r_v2$agc,     r_orig$agc,     tolerance = TOL)
  expect_equal(r_v2$var,     r_orig$var,     tolerance = TOL)
  expect_equal(r_v2$var_ind, r_orig$var_ind, tolerance = TOL)
})

test_that("AGC Uni IID | discrete X & Y: original vs v2 agree", {
  set.seed(104); n <- 300
  X <- sample(1:5, n, replace = TRUE); Y <- sample(1:3, n, replace = TRUE)
  rk <- make_ranks(X, Y)

  r_orig <- run_agc_uni(rk$y_rank, rk$x_rank, "original")
  r_v2   <- run_agc_uni(rk$y_rank, rk$x_rank, "v2")

  expect_equal(r_v2$agc,     r_orig$agc,     tolerance = TOL)
  expect_equal(r_v2$var,     r_orig$var,     tolerance = TOL)
  expect_equal(r_v2$var_ind, r_orig$var_ind, tolerance = TOL)
})

test_that("AGC Uni IID | heavy ties in X: original vs v2 agree", {
  set.seed(105); n <- 300
  X <- rnorm(n); X[seq(1, n, by = 4)] <- X[1]; Y <- rnorm(n)
  rk <- make_ranks(X, Y)

  r_orig <- run_agc_uni(rk$y_rank, rk$x_rank, "original")
  r_v2   <- run_agc_uni(rk$y_rank, rk$x_rank, "v2")

  expect_equal(r_v2$agc,     r_orig$agc,     tolerance = TOL)
  expect_equal(r_v2$var,     r_orig$var,     tolerance = TOL)
  expect_equal(r_v2$var_ind, r_orig$var_ind, tolerance = TOL)
})


# ============================================================================
# Section 2: Univariate HAC
# ============================================================================

test_that("AGC Uni HAC | continuous X & Y: original vs v2 agree", {
  set.seed(201); n <- 300
  X <- rnorm(n); Y <- rnorm(n)
  rk <- make_ranks(X, Y)

  r_orig <- run_agc_uni(rk$y_rank, rk$x_rank, "original", IID = FALSE)
  r_v2   <- run_agc_uni(rk$y_rank, rk$x_rank, "v2",      IID = FALSE)

  expect_equal(r_v2$agc,     r_orig$agc,     tolerance = TOL)
  expect_equal(r_v2$var,     r_orig$var,     tolerance = TOL)
  expect_equal(r_v2$var_ind, r_orig$var_ind, tolerance = TOL)
})

test_that("AGC Uni HAC | binary Y: original vs v2 agree", {
  set.seed(202); n <- 300
  X <- rnorm(n); Y <- rbinom(n, 1, 0.6)
  rk <- make_ranks(X, Y)

  r_orig <- run_agc_uni(rk$y_rank, rk$x_rank, "original", IID = FALSE)
  r_v2   <- run_agc_uni(rk$y_rank, rk$x_rank, "v2",      IID = FALSE)

  expect_equal(r_v2$agc,     r_orig$agc,     tolerance = TOL)
  expect_equal(r_v2$var,     r_orig$var,     tolerance = TOL)
  expect_equal(r_v2$var_ind, r_orig$var_ind, tolerance = TOL)
})


# ============================================================================
# Section 3: Multivariate IID (m = 3)
# ============================================================================

test_that("AGC MV IID m=3 | continuous X & Y: original vs v2 agree", {
  set.seed(301); n <- 300
  X <- matrix(rnorm(n * 3), ncol = 3); Y <- rnorm(n)
  rk <- make_ranks(X, Y)

  r_orig <- run_agc_mv(rk$y_rank, rk$xarray_ranks, "original")
  r_v2   <- run_agc_mv(rk$y_rank, rk$xarray_ranks, "v2")

  expect_equal(r_v2$agc_vector, r_orig$agc_vector, tolerance = TOL)
  expect_equal(r_v2$Sigma,      r_orig$Sigma,      tolerance = TOL)
  expect_equal(r_v2$Sigma_ind,  r_orig$Sigma_ind,  tolerance = TOL)
})


# ============================================================================
# Section 4: Multivariate HAC (m = 3)
# ============================================================================

test_that("AGC MV HAC m=3 | continuous X & Y: original vs v2 agree", {
  set.seed(401); n <- 300
  X <- matrix(rnorm(n * 3), ncol = 3); Y <- rnorm(n)
  rk <- make_ranks(X, Y)

  r_orig <- run_agc_mv(rk$y_rank, rk$xarray_ranks, "original", IID = FALSE)
  r_v2   <- run_agc_mv(rk$y_rank, rk$xarray_ranks, "v2",      IID = FALSE)

  expect_equal(r_v2$agc_vector, r_orig$agc_vector, tolerance = TOL)
  expect_equal(r_v2$Sigma,      r_orig$Sigma,      tolerance = TOL)
  expect_equal(r_v2$Sigma_ind,  r_orig$Sigma_ind,  tolerance = TOL)
})


# ============================================================================
# Section 5: Variance validation against known references
# ============================================================================

test_that("v2: var_ind = 1 for continuous X & Y | n=300", {
  set.seed(811); n <- 300
  X <- rnorm(n); Y <- rnorm(n)
  rk <- make_ranks(X, Y)

  r_v2 <- run_agc_uni(rk$y_rank, rk$x_rank, "v2")
  expect_equal(r_v2$var_ind, 1, tolerance = 1e-12)
})

test_that("v2: ind p-value matches Spearman | n=500, two.sided", {
  set.seed(821); n <- 500
  X <- rnorm(n); Y <- rnorm(n)

  agc_result <- acor.test(X, Y, method = "agc", alternative = "two.sided")
  sp_result  <- cor.test(X, Y, method = "spearman", exact = FALSE)

  expect_equal(agc_result$p.value_ind, sp_result$p.value, tolerance = TOL_PVAL)
})


# ============================================================================
# Section 6: Kernel / Sigma correctness (v2 and binary vs original)
# ============================================================================

test_that("kernel_agc_v2_cpp matches original kernel", {
  test_cases <- list(
    list(name = "Continuous, no ties", n = 150, seed = 9001,
         gen_X = function(n) rnorm(n) + rnorm(n, sd = 0.0001),
         gen_Y = function(n, X) 0.5 * X + rnorm(n, sd = 0.8) + rnorm(n, sd = 0.0001)),
    list(name = "Discrete X and Y", n = 150, seed = 9004,
         gen_X = function(n) sample(1:8, n, replace = TRUE),
         gen_Y = function(n, X) sample(1:4, n, replace = TRUE)),
    list(name = "Many ties (3x3)", n = 200, seed = 9005,
         gen_X = function(n) sample(1:3, n, replace = TRUE),
         gen_Y = function(n, X) sample(1:3, n, replace = TRUE))
  )

  for (tc in test_cases) {
    set.seed(tc$seed)
    X <- tc$gen_X(tc$n)
    Y <- tc$gen_Y(tc$n, X)

    y_rank <- rank(Y, ties.method = "average")
    x_rank <- rank(X, ties.method = "average")

    result <- acor:::comp_rho_agc(y_rank, x_rank)
    rho <- result$rho

    kp_original <- compute_kp_original(y_rank, x_rank, rho)
    kp_v2 <- acor:::kernel_agc_v2_cpp(x_rank, y_rank, rho)

    expect_equal(kp_v2, kp_original, tolerance = 1e-10,
                 info = sprintf("%s: v2 kernel should match original", tc$name))
  }
})

test_that("kernel_agc_binary matches original kernel", {
  test_cases <- list(
    list(name = "Binary Y, continuous X", n = 150, seed = 9101,
         gen_X = function(n) rnorm(n),
         gen_Y = function(n, X) rbinom(n, 1, plogis(0.5 * X))),
    list(name = "Binary Y, discrete X", n = 150, seed = 9102,
         gen_X = function(n) sample(1:10, n, replace = TRUE),
         gen_Y = function(n, X) rbinom(n, 1, 0.4))
  )

  for (tc in test_cases) {
    set.seed(tc$seed)
    X <- tc$gen_X(tc$n)
    Y <- tc$gen_Y(tc$n, X)

    y_rank <- rank(Y, ties.method = "average")
    x_rank <- rank(X, ties.method = "average")

    result <- acor:::comp_rho_agc(y_rank, x_rank)
    rho <- result$rho

    kp_original <- compute_kp_original(y_rank, x_rank, rho)
    kp_binary <- acor:::kernel_agc_binary(x_rank, y_rank, rho)

    expect_equal(kp_binary, kp_original, tolerance = 1e-10,
                 info = sprintf("%s: binary kernel should match original", tc$name))
  }
})

test_that("Sigma_agc_v2 and Sigma_agc_binary match Sigma_agc (univariate IID)", {
  test_cases <- list(
    list(name = "Continuous", n = 150, seed = 9201, binary_y = FALSE),
    list(name = "Binary Y", n = 150, seed = 9203, binary_y = TRUE)
  )

  for (tc in test_cases) {
    set.seed(tc$seed)
    X <- rnorm(tc$n)
    Y <- if (tc$binary_y) rbinom(tc$n, 1, plogis(0.5 * X)) else 0.5 * X + rnorm(tc$n, sd = 0.8)

    y_rank <- rank(Y, ties.method = "average")
    x_rank <- rank(X, ties.method = "average")

    res_orig <- Sigma_agc(y_rank, x_rank)
    res_v2 <- acor:::Sigma_agc_v2(y_rank, x_rank)

    expect_equal(res_v2$agc, res_orig$agc, tolerance = 1e-10,
                 info = sprintf("%s: agc v2 should match original", tc$name))
    expect_equal(res_v2$var, res_orig$var, tolerance = 1e-10,
                 info = sprintf("%s: var v2 should match original", tc$name))

    if (tc$binary_y) {
      res_bin <- acor:::Sigma_agc_binary(y_rank, x_rank)
      expect_equal(res_bin$agc, res_orig$agc, tolerance = 1e-10,
                   info = sprintf("%s: agc binary should match original", tc$name))
      expect_equal(res_bin$var, res_orig$var, tolerance = 1e-10,
                   info = sprintf("%s: var binary should match original", tc$name))
    }
  }
})

test_that("Sigma_agc_ts_v2 matches Sigma_agc_ts (HAC)", {
  set.seed(9301); n <- 150
  X <- rnorm(n); Y <- 0.5 * X + rnorm(n, sd = 0.8)

  y_rank <- rank(Y, ties.method = "average")
  x_rank <- rank(X, ties.method = "average")

  res_orig <- Sigma_agc_ts(y_rank, x_rank)
  res_v2 <- acor:::Sigma_agc_ts_v2(y_rank, x_rank)

  expect_equal(res_v2$agc, res_orig$agc, tolerance = 1e-10)
  expect_equal(res_v2$var, res_orig$var, tolerance = 1e-10)
})

test_that("Sigma_agc_multivariate_v2 matches original (multivariate IID)", {
  set.seed(9402); n <- 100; m <- 3
  X <- matrix(rnorm(n * m), nrow = n, ncol = m)
  Y <- 0.5 * rowMeans(X) + rnorm(n, sd = 0.8)

  y_rank <- rank(Y, ties.method = "average")
  xarray_ranks <- matrix(0, nrow = m, ncol = n)
  for (j in 1:m) {
    xarray_ranks[j, ] <- rank(X[, j], ties.method = "average")
  }

  res_orig <- Sigma_agc_multivariate(y_rank, xarray_ranks)
  res_v2 <- acor:::Sigma_agc_multivariate_v2(y_rank, xarray_ranks)

  expect_equal(res_v2$agc_vector, res_orig$agc_vector, tolerance = 1e-10)
  expect_equal(res_v2$Sigma, res_orig$Sigma, tolerance = 1e-10)
})

test_that("Sigma_agc_multivariate_ts_v2 matches original (multivariate HAC)", {
  set.seed(9501); n <- 100; m <- 2
  X <- matrix(rnorm(n * m), nrow = n, ncol = m)
  Y <- 0.5 * rowMeans(X) + rnorm(n, sd = 0.8)

  y_rank <- rank(Y, ties.method = "average")
  xarray_ranks <- matrix(0, nrow = m, ncol = n)
  for (j in 1:m) {
    xarray_ranks[j, ] <- rank(X[, j], ties.method = "average")
  }

  res_orig <- Sigma_agc_multivariate_ts(y_rank, xarray_ranks)
  res_v2 <- acor:::Sigma_agc_multivariate_ts_v2(y_rank, xarray_ranks)

  expect_equal(res_v2$agc_vector, res_orig$agc_vector, tolerance = 1e-10)
  expect_equal(res_v2$Sigma, res_orig$Sigma, tolerance = 1e-10)
})
