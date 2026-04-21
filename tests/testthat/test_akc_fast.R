# ============================================================================
# Tests for Optimized AKC Kernel Computations
# ============================================================================
#
# Tests correctness and performance of:
# 1. V1: O(R×M + n) using unique values
# 2. V2: O(n log n) using Fenwick trees
#
# Compared against original O(n²) implementations
#
# Runtime benchmarks have been moved to benchmarks/run_benchmarks.R
#
# ============================================================================

library(testthat)
library(acor) 

# ============================================================================
# Helper: Original (loop-based) vectorized kernel functions for comparison
# ============================================================================

# Vectorized original K_tau (loop-based)
K_tau_vec_original <- function(X, Y, tau_XY) {
  n <- length(X)
  result <- numeric(n)
  for (i in 1:n) {
    result[i] <- K_tau(X[i], Y[i], X, Y, tau_XY)
  }
  result
}

# Vectorized original K_p (loop-based)
K_p_vec_original <- function(Y, tau_y) {
  n <- length(Y)
  result <- numeric(n)
  for (i in 1:n) {
    result[i] <- K_p(Y[i], Y, tau_y)
  }
  result
}

# Vectorized original F_bar (loop-based)
F_bar_vec_original <- function(X) {
  sapply(X, function(x) F_bar(x, X))
}

# Vectorized original H_bar (loop-based)
H_bar_vec_original <- function(X, Y) {
  n <- length(X)
  sapply(1:n, function(i) H_bar(X[i], Y[i], X, Y))
}

# ============================================================================
# TEST 1: Correctness of F_bar_vec_v1
# ============================================================================

test_that("F_bar_vec_v1 produces identical results to original", {

  test_cases <- list(
    list(name = "Continuous, no ties", n = 100, seed = 1001,
         gen_X = function(n) rnorm(n) + rnorm(n, sd = 0.0001)),
    list(name = "With ties", n = 100, seed = 1002,
         gen_X = function(n) { x <- rnorm(n); x[seq(1, n, by = 5)] <- x[1]; x }),
    list(name = "Discrete (10 levels)", n = 100, seed = 1003,
         gen_X = function(n) sample(1:10, n, replace = TRUE)),
    list(name = "Binary", n = 100, seed = 1004,
         gen_X = function(n) rbinom(n, 1, 0.6))
  )

  for (tc in test_cases) {
    set.seed(tc$seed)
    X <- tc$gen_X(tc$n)

    F_orig <- F_bar_vec_original(X)
    F_v1 <- acor:::F_bar_vec(X)

    expect_equal(F_v1, F_orig, tolerance = 1e-10,
                 info = sprintf("%s: F_bar_vec_v1 should match original", tc$name))
  }
})

# ============================================================================
# TEST 2: Correctness of K_p_vec_v1 and K_p_vec_v2
# ============================================================================

test_that("K_p_vec_v1 and K_p_vec_v2 produce identical results to original", {

  test_cases <- list(
    list(name = "Continuous Y", n = 100, seed = 2001,
         gen_Y = function(n) rnorm(n) + rnorm(n, sd = 0.0001)),
    list(name = "Binary Y", n = 100, seed = 2002,
         gen_Y = function(n) rbinom(n, 1, 0.6)),
    list(name = "Discrete Y (5 levels)", n = 100, seed = 2003,
         gen_Y = function(n) sample(1:5, n, replace = TRUE)),
    list(name = "Many ties in Y", n = 100, seed = 2004,
         gen_Y = function(n) sample(1:3, n, replace = TRUE))
  )

  for (tc in test_cases) {
    set.seed(tc$seed)
    Y <- tc$gen_Y(tc$n)

    # Compute tau_Y using the original function
    tau_Y_result <- acor:::tau_Y_func(Y)
    tau_Y <- tau_Y_result$expectation

    K_p_orig <- K_p_vec_original(Y, tau_Y)
    K_p_v1 <- acor:::K_p_vec(Y, tau_Y)
    K_p_v2 <- acor:::K_p_vec(Y, tau_Y)

    expect_equal(K_p_v1, K_p_orig, tolerance = 1e-10,
                 info = sprintf("%s: K_p_vec_v1 should match original", tc$name))
    expect_equal(K_p_v2, K_p_orig, tolerance = 1e-10,
                 info = sprintf("%s: K_p_vec_v2 should match original", tc$name))
  }
})

# ============================================================================
# TEST 3: Correctness of H_bar_vec_v1 and H_bar_vec_v2
# ============================================================================

test_that("H_bar_vec_v1 and H_bar_vec_v2_cpp produce identical results to original", {

  test_cases <- list(
    list(name = "Continuous, no ties", n = 100, seed = 3001,
         gen_X = function(n) rnorm(n) + rnorm(n, sd = 0.0001),
         gen_Y = function(n) rnorm(n) + rnorm(n, sd = 0.0001)),
    list(name = "Binary Y", n = 100, seed = 3002,
         gen_X = function(n) rnorm(n),
         gen_Y = function(n) rbinom(n, 1, 0.6)),
    list(name = "Discrete X and Y", n = 100, seed = 3003,
         gen_X = function(n) sample(1:5, n, replace = TRUE),
         gen_Y = function(n) sample(1:3, n, replace = TRUE)),
    list(name = "Ties in X only", n = 100, seed = 3004,
         gen_X = function(n) { x <- rnorm(n); x[seq(1, n, by = 5)] <- x[1]; x },
         gen_Y = function(n) rnorm(n) + rnorm(n, sd = 0.0001)),
    list(name = "Large discrete", n = 200, seed = 3005,
         gen_X = function(n) sample(1:20, n, replace = TRUE),
         gen_Y = function(n) sample(1:10, n, replace = TRUE))
  )

  for (tc in test_cases) {
    set.seed(tc$seed)
    X <- tc$gen_X(tc$n)
    Y <- tc$gen_Y(tc$n)

    H_orig <- H_bar_vec_original(X, Y)
    H_v1 <- acor:::H_bar_vec_v1(X, Y)
    H_v2 <- acor:::H_bar_vec_v2_cpp(X, Y)

    expect_equal(H_v1, H_orig, tolerance = 1e-10,
                 info = sprintf("%s: H_bar_vec_v1 should match original", tc$name))
    expect_equal(H_v2, H_orig, tolerance = 1e-10,
                 info = sprintf("%s: H_bar_vec_v2_cpp should match original", tc$name))
  }
})

# ============================================================================
# TEST 4: Correctness of K_tau_vec_v1 and K_tau_vec_v2
# ============================================================================

test_that("K_tau_vec_v1 and K_tau_vec_v2 produce identical results to original", {

  test_cases <- list(
    list(name = "Continuous, no ties", n = 100, seed = 4001,
         gen_X = function(n) rnorm(n) + rnorm(n, sd = 0.0001),
         gen_Y = function(n, X) 0.5 * X + rnorm(n, sd = 0.8) + rnorm(n, sd = 0.0001)),
    list(name = "Binary Y", n = 100, seed = 4002,
         gen_X = function(n) rnorm(n),
         gen_Y = function(n, X) rbinom(n, 1, 0.6)),
    list(name = "Discrete X and Y", n = 100, seed = 4003,
         gen_X = function(n) sample(1:10, n, replace = TRUE),
         gen_Y = function(n, X) sample(1:5, n, replace = TRUE)),
    list(name = "Ties in X", n = 100, seed = 4004,
         gen_X = function(n) { x <- rnorm(n); x[seq(1, n, by = 5)] <- x[1]; x },
         gen_Y = function(n, X) rnorm(n)),
    list(name = "Larger sample, continuous", n = 300, seed = 4005,
         gen_X = function(n) rnorm(n) + rnorm(n, sd = 0.0001),
         gen_Y = function(n, X) 0.5 * X + rnorm(n, sd = 0.8) + rnorm(n, sd = 0.0001))
  )

  for (tc in test_cases) {
    set.seed(tc$seed)
    X <- tc$gen_X(tc$n)
    Y <- tc$gen_Y(tc$n, X)

    akc_result <- kendall_tau_sign(X, Y)
    tau_XY <- akc_result$expectation

    K_tau_orig <- K_tau_vec_original(X, Y, tau_XY)
    K_tau_v1 <- acor:::K_tau_vec_v1(X, Y, tau_XY)
    K_tau_v2 <- acor:::K_tau_vec_v2(X, Y, tau_XY)

    expect_equal(K_tau_v1, K_tau_orig, tolerance = 1e-10,
                 info = sprintf("%s: K_tau_vec_v1 should match original", tc$name))
    expect_equal(K_tau_v2, K_tau_orig, tolerance = 1e-10,
                 info = sprintf("%s: K_tau_vec_v2 should match original", tc$name))
  }
})

# ============================================================================
# TEST 5: Correctness of Sigma_akc_multivariate_v1 and Sigma_akc_multivariate_v2
# ============================================================================

test_that("Sigma_akc_multivariate_v1 and _v2 produce identical results to original", {

  test_cases <- list(
    list(name = "Continuous, single predictor", n = 100, m = 1, seed = 5001, binary_y = FALSE),
    list(name = "Continuous, multiple predictors", n = 100, m = 3, seed = 5002, binary_y = FALSE),
    list(name = "Binary Y, single predictor", n = 100, m = 1, seed = 5003, binary_y = TRUE),
    list(name = "Binary Y, multiple predictors", n = 100, m = 3, seed = 5004, binary_y = TRUE),
    list(name = "Larger sample", n = 200, m = 2, seed = 5005, binary_y = FALSE)
  )

  for (tc in test_cases) {
    set.seed(tc$seed)
    X <- matrix(rnorm(tc$n * tc$m), nrow = tc$n, ncol = tc$m)

    if (tc$binary_y) {
      Y <- rbinom(tc$n, 1, 0.6)
    } else {
      Y <- 0.5 * rowMeans(X) + rnorm(tc$n, sd = 0.8)
    }

    res_orig <- Sigma_akc_multivariate(X, Y)
    res_v1 <- acor:::Sigma_akc_multivariate_v1(X, Y)
    res_v2 <- acor:::Sigma_akc_multivariate_v2(X, Y)

    expect_equal(res_v1$akc_vector, res_orig$akc_vector, tolerance = 1e-10,
                 info = sprintf("%s: akc_vector v1 should match original", tc$name))
    expect_equal(res_v1$Sigma, res_orig$Sigma, tolerance = 1e-10,
                 info = sprintf("%s: Sigma v1 should match original", tc$name))

    expect_equal(res_v2$akc_vector, res_orig$akc_vector, tolerance = 1e-10,
                 info = sprintf("%s: akc_vector v2 should match original", tc$name))
    expect_equal(res_v2$Sigma, res_orig$Sigma, tolerance = 1e-10,
                 info = sprintf("%s: Sigma v2 should match original", tc$name))
  }
})

# ============================================================================
# TEST 6: Correctness of Sigma_akc_multivariate_ts_v1 and _ts_v2 (HAC)
# ============================================================================

test_that("Sigma_akc_multivariate_ts_v1 and _ts_v2 produce identical results to original", {

  test_cases <- list(
    list(name = "Continuous, single predictor", n = 100, m = 1, seed = 6001, binary_y = FALSE),
    list(name = "Continuous, multiple predictors", n = 100, m = 3, seed = 6002, binary_y = FALSE),
    list(name = "Binary Y, single predictor", n = 100, m = 1, seed = 6003, binary_y = TRUE),
    list(name = "Binary Y, multiple predictors", n = 100, m = 3, seed = 6004, binary_y = TRUE),
    list(name = "Larger sample", n = 200, m = 2, seed = 6005, binary_y = FALSE)
  )

  for (tc in test_cases) {
    set.seed(tc$seed)
    X <- matrix(rnorm(tc$n * tc$m), nrow = tc$n, ncol = tc$m)

    if (tc$binary_y) {
      Y <- rbinom(tc$n, 1, 0.6)
    } else {
      Y <- 0.5 * rowMeans(X) + rnorm(tc$n, sd = 0.8)
    }

    res_orig <- Sigma_akc_multivariate_ts(X, Y)
    res_v1 <- Sigma_akc_multivariate_ts_v1(X, Y)
    res_v2 <- Sigma_akc_multivariate_ts_v2(X, Y)

    expect_equal(res_v1$akc_vector, res_orig$akc_vector, tolerance = 1e-10,
                 info = sprintf("%s: akc_vector v1 should match original", tc$name))
    expect_equal(res_v1$Sigma, res_orig$Sigma, tolerance = 1e-10,
                 info = sprintf("%s: Sigma v1 should match original", tc$name))

    expect_equal(res_v2$akc_vector, res_orig$akc_vector, tolerance = 1e-10,
                 info = sprintf("%s: akc_vector v2 should match original", tc$name))
    expect_equal(res_v2$Sigma, res_orig$Sigma, tolerance = 1e-10,
                 info = sprintf("%s: Sigma v2 should match original", tc$name))
  }
})

# tau_Y_func uses pair_tie_proportion_cpp so Y ties match kendall_tau_sign_cpp.

test_that("tau_Y_func p_tie_y matches pair_tie_proportion_cpp (Kendall-consistent)", {
  set.seed(4001)
  for (i in 1:25) {
    Y <- rnorm(80 + i)
    ty <- acor:::tau_Y_func(Y)
    p_cpp <- pair_tie_proportion_cpp(as.numeric(Y))
    expect_equal(ty$p_tie_y, p_cpp, tolerance = 1e-15)
    expect_equal(ty$expectation, 1 - p_cpp, tolerance = 1e-15)
  }
})

test_that("tau_Y_func matches table(Y) for exact-duplicate numeric Y", {
  Y <- rep(c(0, 1.5, pi), c(4, 10, 6))
  n <- length(Y)
  np <- n * (n - 1) / 2
  f <- as.numeric(table(Y))
  p_tab <- sum(f * (f - 1) / 2) / np
  ty <- acor:::tau_Y_func(Y)
  expect_equal(ty$p_tie_y, p_tab, tolerance = 1e-15)
  expect_equal(ty$p_tie_y, pair_tie_proportion_cpp(as.numeric(Y)), tolerance = 1e-15)
})
