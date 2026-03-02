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
# ============================================================================

library(testthat)
library(acor) 

# ============================================================================
# Helper: Timing function
# ============================================================================

time_execution <- function(expr, replications = 5, envir = parent.frame()) {
  times <- numeric(replications)
  for (i in 1:replications) {
    start <- Sys.time()
    eval(expr, envir = envir)
    end <- Sys.time()
    times[i] <- as.numeric(difftime(end, start, units = "secs"))
  }
  list(
    mean = mean(times),
    median = median(times),
    sd = sd(times),
    min = min(times),
    max = max(times)
  )
}

# ============================================================================
# Helper: Original (loop-based) vectorized kernel functions for comparison
# ============================================================================

# Vectorized original K_tau (loop-based)
K_tau_vec_original <- function(X, Y, tau_XY) {
  n <- length(X)
  result <- numeric(n)
  for (i in 1:n) {
    result[i] <- acor:::K_tau(X[i], Y[i], X, Y, tau_XY)
  }
  result
}

# Vectorized original K_p (loop-based)
K_p_vec_original <- function(Y, tau_y) {
  n <- length(Y)
  result <- numeric(n)
  for (i in 1:n) {
    result[i] <- acor:::K_p(Y[i], Y, tau_y)
  }
  result
}

# Vectorized original F_bar (loop-based)
F_bar_vec_original <- function(X) {
  sapply(X, function(x) acor:::F_bar(x, X))
}

# Vectorized original H_bar (loop-based)
H_bar_vec_original <- function(X, Y) {
  n <- length(X)
  sapply(1:n, function(i) acor:::H_bar(X[i], Y[i], X, Y))
}

Sigma_akc_v1_wrapper <- function(X, Y) {
  if (acor:::is_binary(Y)) {
    tau_Y_result <- acor:::tau_Y_func_binary(Y)
    akc_result <- acor:::kendall_tau_sign_binary(X, Y)
  } else {
    tau_Y_result <- acor:::tau_Y_func(Y)
    akc_result <- acor:::kendall_tau_sign(X, Y)
  }
  tau_Y <- tau_Y_result$expectation
  p_Y <- tau_Y_result$p_tie_y
  tau_XY <- akc_result$expectation
  akc <- akc_result$tau

  var <- acor:::Sigma_akc_v1(X, Y, tau_XY, tau_Y, p_Y)
  list(akc = akc, var = var)
}

Sigma_akc_v2_wrapper <- function(X, Y) {
  if (acor:::is_binary(Y)) {
    tau_Y_result <- acor:::tau_Y_func_binary(Y)
    akc_result <- acor:::kendall_tau_sign_binary(X, Y)
  } else {
    tau_Y_result <- acor:::tau_Y_func(Y)
    akc_result <- acor:::kendall_tau_sign(X, Y)
  }
  tau_Y <- tau_Y_result$expectation
  p_Y <- tau_Y_result$p_tie_y
  tau_XY <- akc_result$expectation
  akc <- akc_result$tau

  var <- acor:::Sigma_akc_v2(X, Y, tau_XY, tau_Y, p_Y)
  list(akc = akc, var = var)
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
    F_v1 <- acor:::F_bar_vec_v1(X)

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
    K_p_v1 <- acor:::K_p_vec_v1(Y, tau_Y)
    K_p_v2 <- acor:::K_p_vec_v2(Y, tau_Y)

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

    # Use kendall_tau_sign to get tau_XY
    akc_result <- acor:::kendall_tau_sign(X, Y)
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

    res_orig <- acor:::Sigma_akc_multivariate(X, Y)
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

# ============================================================================
# TEST 7: Runtime benchmark - K_tau kernel computation
# ============================================================================

test_that("Runtime benchmark: K_tau - Original vs V1 vs V2", {
  skip_on_cran()
  cat("\n")
  cat(paste(rep("=", 80), collapse = ""))
  cat("\nRUNTIME BENCHMARK: K_tau Implementations\n")
  cat(paste(rep("=", 80), collapse = ""))
  cat("\n\n")

  n_replications <- 5
  sample_sizes <- c(500, 1000, 5000, 8000)

  # -------------------------------------------------------------------------
  # Scenario 1: Continuous data (no ties)
  # -------------------------------------------------------------------------
  cat("SCENARIO 1: Continuous data (no ties)\n")
  cat(sprintf("%8s %12s %12s %12s %12s %12s\n",
              "n", "Original", "V1", "V2", "V1 speedup", "V2 speedup"))
  cat(paste(rep("-", 80), collapse = ""))
  cat("\n")

  for (n in sample_sizes) {
    set.seed(7000 + n)
    X <- rnorm(n) + rnorm(n, sd = 0.0001)
    Y <- 0.5 * X + rnorm(n, sd = 0.8) + rnorm(n, sd = 0.0001)

    akc_result <- acor:::kendall_tau_sign(X, Y)
    tau_XY <- akc_result$expectation

    time_orig <- time_execution(quote(K_tau_vec_original(X, Y, tau_XY)),
                                replications = n_replications)$median
    time_v1 <- time_execution(quote(acor:::K_tau_vec_v1(X, Y, tau_XY)),
                              replications = n_replications)$median
    time_v2 <- time_execution(quote(acor:::K_tau_vec_v2(X, Y, tau_XY)),
                              replications = n_replications)$median

    speedup_v1 <- if (time_v1 > 1e-9) time_orig / time_v1 else NA
    speedup_v2 <- if (time_v2 > 1e-9) time_orig / time_v2 else NA

    cat(sprintf("%8d %12.4f %12.4f %12.4f %12.1fx %12.1fx\n",
                n, time_orig, time_v1, time_v2,
                ifelse(is.na(speedup_v1), 0, speedup_v1),
                ifelse(is.na(speedup_v2), 0, speedup_v2)))
  }

  cat("\n")

  # -------------------------------------------------------------------------
  # Scenario 2: Binary Y
  # -------------------------------------------------------------------------
  cat("SCENARIO 2: Binary Y\n")
  cat(sprintf("%8s %12s %12s %12s %12s %12s\n",
              "n", "Original", "V1", "V2", "V1 speedup", "V2 speedup"))
  cat(paste(rep("-", 80), collapse = ""))
  cat("\n")

  for (n in sample_sizes) {
    set.seed(8000 + n)
    X <- rnorm(n)
    Y <- rbinom(n, 1, 0.6)

    akc_result <- acor:::kendall_tau_sign_binary(X, Y)
    tau_XY <- akc_result$expectation

    time_orig <- time_execution(quote(K_tau_vec_original(X, Y, tau_XY)),
                                replications = n_replications)$median
    time_v1 <- time_execution(quote(acor:::K_tau_vec_v1(X, Y, tau_XY)),
                              replications = n_replications)$median
    time_v2 <- time_execution(quote(acor:::K_tau_vec_v2(X, Y, tau_XY)),
                              replications = n_replications)$median

    speedup_v1 <- if (time_v1 > 1e-9) time_orig / time_v1 else NA
    speedup_v2 <- if (time_v2 > 1e-9) time_orig / time_v2 else NA

    cat(sprintf("%8d %12.4f %12.4f %12.4f %12.1fx %12.1fx\n",
                n, time_orig, time_v1, time_v2,
                ifelse(is.na(speedup_v1), 0, speedup_v1),
                ifelse(is.na(speedup_v2), 0, speedup_v2)))
  }

  cat("\n")

  # -------------------------------------------------------------------------
  # Scenario 3: Discrete X and Y (10 levels each)
  # -------------------------------------------------------------------------
  cat("SCENARIO 3: Discrete X and Y (10 levels each)\n")
  cat(sprintf("%8s %12s %12s %12s %12s %12s\n",
              "n", "Original", "V1", "V2", "V1 speedup", "V2 speedup"))
  cat(paste(rep("-", 80), collapse = ""))
  cat("\n")

  for (n in sample_sizes) {
    set.seed(9000 + n)
    X <- sample(1:10, n, replace = TRUE)
    Y <- sample(1:10, n, replace = TRUE)

    akc_result <- acor:::kendall_tau_sign(X, Y)
    tau_XY <- akc_result$expectation

    time_orig <- time_execution(quote(K_tau_vec_original(X, Y, tau_XY)),
                                replications = n_replications)$median
    time_v1 <- time_execution(quote(acor:::K_tau_vec_v1(X, Y, tau_XY)),
                              replications = n_replications)$median
    time_v2 <- time_execution(quote(acor:::K_tau_vec_v2(X, Y, tau_XY)),
                              replications = n_replications)$median

    speedup_v1 <- if (time_v1 > 1e-9) time_orig / time_v1 else NA
    speedup_v2 <- if (time_v2 > 1e-9) time_orig / time_v2 else NA

    cat(sprintf("%8d %12.4f %12.4f %12.4f %12.1fx %12.1fx\n",
                n, time_orig, time_v1, time_v2,
                ifelse(is.na(speedup_v1), 0, speedup_v1),
                ifelse(is.na(speedup_v2), 0, speedup_v2)))
  }

  cat("\n")
  expect_true(TRUE, info = "Benchmark completed")
})

# ============================================================================
# TEST 8: Runtime benchmark - Full multivariate IID variance
# ============================================================================

test_that("Runtime benchmark: Sigma_akc_multivariate - Original vs V1 vs V2", {
  skip_on_cran()
  cat("\n")
  cat(paste(rep("=", 80), collapse = ""))
  cat("\nRUNTIME BENCHMARK: Sigma_akc_multivariate (IID)\n")
  cat(paste(rep("=", 80), collapse = ""))
  cat("\n\n")

  n_replications <- 5
  sample_sizes <- c(500, 1000, 5000)
  m <- 3  # number of predictors

  # -------------------------------------------------------------------------
  # Continuous data
  # -------------------------------------------------------------------------
  cat("Continuous data (m = 3 predictors):\n")
  cat(sprintf("%8s %12s %12s %12s %12s %12s\n",
              "n", "Original", "V1", "V2", "V1 speedup", "V2 speedup"))
  cat(paste(rep("-", 80), collapse = ""))
  cat("\n")

  for (n in sample_sizes) {
    set.seed(10000 + n)
    X <- matrix(rnorm(n * m), nrow = n, ncol = m)
    Y <- 0.5 * rowMeans(X) + rnorm(n, sd = 0.8)

    time_orig <- time_execution(quote(acor:::Sigma_akc_multivariate(X, Y)),
                                replications = n_replications)$median
    time_v1 <- time_execution(quote(acor:::Sigma_akc_multivariate_v1(X, Y)),
                              replications = n_replications)$median
    time_v2 <- time_execution(quote(acor:::Sigma_akc_multivariate_v2(X, Y)),
                              replications = n_replications)$median

    speedup_v1 <- if (time_v1 > 1e-9) time_orig / time_v1 else NA
    speedup_v2 <- if (time_v2 > 1e-9) time_orig / time_v2 else NA

    cat(sprintf("%8d %12.4f %12.4f %12.4f %12.1fx %12.1fx\n",
                n, time_orig, time_v1, time_v2,
                ifelse(is.na(speedup_v1), 0, speedup_v1),
                ifelse(is.na(speedup_v2), 0, speedup_v2)))
  }

  cat("\n")

  # -------------------------------------------------------------------------
  # Binary Y
  # -------------------------------------------------------------------------
  cat("Binary Y (m = 3 predictors):\n")
  cat(sprintf("%8s %12s %12s %12s %12s %12s\n",
              "n", "Original", "V1", "V2", "V1 speedup", "V2 speedup"))
  cat(paste(rep("-", 80), collapse = ""))
  cat("\n")

  for (n in sample_sizes) {
    set.seed(11000 + n)
    X <- matrix(rnorm(n * m), nrow = n, ncol = m)
    Y <- rbinom(n, 1, 0.6)

    time_orig <- time_execution(quote(acor:::Sigma_akc_multivariate(X, Y)),
                                replications = n_replications)$median
    time_v1 <- time_execution(quote(acor:::Sigma_akc_multivariate_v1(X, Y)),
                              replications = n_replications)$median
    time_v2 <- time_execution(quote(acor:::Sigma_akc_multivariate_v2(X, Y)),
                              replications = n_replications)$median

    speedup_v1 <- if (time_v1 > 1e-9) time_orig / time_v1 else NA
    speedup_v2 <- if (time_v2 > 1e-9) time_orig / time_v2 else NA

    cat(sprintf("%8d %12.4f %12.4f %12.4f %12.1fx %12.1fx\n",
                n, time_orig, time_v1, time_v2,
                ifelse(is.na(speedup_v1), 0, speedup_v1),
                ifelse(is.na(speedup_v2), 0, speedup_v2)))
  }

  cat("\n")
  expect_true(TRUE, info = "Benchmark completed")
})

# ============================================================================
# TEST 9: Runtime benchmark - Full multivariate HAC variance
# ============================================================================

test_that("Runtime benchmark: Sigma_akc_multivariate_ts (HAC) - Original vs V1 vs V2", {
  skip_on_cran()
  cat("\n")
  cat(paste(rep("=", 80), collapse = ""))
  cat("\nRUNTIME BENCHMARK: Sigma_akc_multivariate_ts (HAC)\n")
  cat(paste(rep("=", 80), collapse = ""))
  cat("\n\n")

  n_replications <- 5
  sample_sizes <- c(500, 1000, 5000)#, 8000)
  m <- 3  # number of predictors

  # -------------------------------------------------------------------------
  # Continuous data
  # -------------------------------------------------------------------------
  cat("Continuous data (m = 3 predictors):\n")
  cat(sprintf("%8s %12s %12s %12s %12s %12s\n",
              "n", "Original", "V1", "V2", "V1 speedup", "V2 speedup"))
  cat(paste(rep("-", 80), collapse = ""))
  cat("\n")

  for (n in sample_sizes) {
    set.seed(12000 + n)
    X <- matrix(rnorm(n * m), nrow = n, ncol = m)
    Y <- 0.5 * rowMeans(X) + rnorm(n, sd = 0.8)

    time_orig <- time_execution(quote(Sigma_akc_multivariate_ts(X, Y)),
                                replications = n_replications)$median
    time_v1 <- time_execution(quote(Sigma_akc_multivariate_ts_v1(X, Y)),
                              replications = n_replications)$median
    time_v2 <- time_execution(quote(Sigma_akc_multivariate_ts_v2(X, Y)),
                              replications = n_replications)$median

    speedup_v1 <- if (time_v1 > 1e-9) time_orig / time_v1 else NA
    speedup_v2 <- if (time_v2 > 1e-9) time_orig / time_v2 else NA

    cat(sprintf("%8d %12.4f %12.4f %12.4f %12.1fx %12.1fx\n",
                n, time_orig, time_v1, time_v2,
                ifelse(is.na(speedup_v1), 0, speedup_v1),
                ifelse(is.na(speedup_v2), 0, speedup_v2)))
  }

  cat("\n")

  # -------------------------------------------------------------------------
  # Binary Y
  # -------------------------------------------------------------------------
  cat("Binary Y (m = 3 predictors):\n")
  cat(sprintf("%8s %12s %12s %12s %12s %12s\n",
              "n", "Original", "V1", "V2", "V1 speedup", "V2 speedup"))
  cat(paste(rep("-", 80), collapse = ""))
  cat("\n")

  for (n in sample_sizes) {
    set.seed(13000 + n)
    X <- matrix(rnorm(n * m), nrow = n, ncol = m)
    Y <- rbinom(n, 1, 0.6)

    time_orig <- time_execution(quote(Sigma_akc_multivariate_ts(X, Y)),
                                replications = n_replications)$median
    time_v1 <- time_execution(quote(Sigma_akc_multivariate_ts_v1(X, Y)),
                              replications = n_replications)$median
    time_v2 <- time_execution(quote(Sigma_akc_multivariate_ts_v2(X, Y)),
                              replications = n_replications)$median

    speedup_v1 <- if (time_v1 > 1e-9) time_orig / time_v1 else NA
    speedup_v2 <- if (time_v2 > 1e-9) time_orig / time_v2 else NA

    cat(sprintf("%8d %12.4f %12.4f %12.4f %12.1fx %12.1fx\n",
                n, time_orig, time_v1, time_v2,
                ifelse(is.na(speedup_v1), 0, speedup_v1),
                ifelse(is.na(speedup_v2), 0, speedup_v2)))
  }

  cat("\n")
  expect_true(TRUE, info = "Benchmark completed")
})

# ============================================================================
# TEST 10: Computational complexity analysis
# ============================================================================

test_that("V1 vs V2 comparison: Binary Y with large n", {
  skip_on_cran()
  cat("\n")
  cat(paste(rep("=", 80), collapse = ""))
  cat("\nV1 vs V2 COMPARISON: Binary Y (n = 100000)\n")
  cat(paste(rep("=", 80), collapse = ""))
  cat("\n\n")

  n <- 100000
  n_replications <- 3

  set.seed(20001)
  X <- rnorm(n)  # continuous X, no ties
  Y <- rbinom(n, 1, 0.6)  # binary Y

  akc_result <- acor:::kendall_tau_sign_binary(X, Y)
  tau_XY <- akc_result$expectation

  cat(sprintf("n = %d, unique(Y) = %d, unique(X) = %d\n\n", n, length(unique(Y)), length(unique(X))))

  # Time V1
  time_v1 <- time_execution(quote(acor:::K_tau_vec_v1(X, Y, tau_XY)),
                            replications = n_replications)$median

  # Time V2
  time_v2 <- time_execution(quote(acor:::K_tau_vec_v2(X, Y, tau_XY)),
                            replications = n_replications)$median

  cat(sprintf("K_tau_vec_v1: %.4f seconds\n", time_v1))
  cat(sprintf("K_tau_vec_v2: %.4f seconds\n", time_v2))

  if (time_v1 < time_v2) {
    cat(sprintf("\n--> V1 is %.1fx faster than V2 for binary Y\n", time_v2 / time_v1))
  } else {
    cat(sprintf("\n--> V2 is %.1fx faster than V1 for binary Y\n", time_v1 / time_v2))
  }

  cat("\n")
  expect_true(TRUE, info = "Binary Y comparison completed")
})

test_that("V1 vs V2 comparison: Y with 3 distinct values and large n", {
  skip_on_cran()
  cat("\n")
  cat(paste(rep("=", 80), collapse = ""))
  cat("\nV1 vs V2 COMPARISON: Y with 3 distinct values (n = 100000)\n")
  cat(paste(rep("=", 80), collapse = ""))
  cat("\n\n")

  n <- 100000
  n_replications <- 3

  set.seed(20002)
  X <- rnorm(n)  # continuous X, no ties
  Y <- sample(1:3, n, replace = TRUE)  # Y with 3 distinct values

  akc_result <- acor:::kendall_tau_sign(X, Y)
  tau_XY <- akc_result$expectation

  cat(sprintf("n = %d, unique(Y) = %d, unique(X) = %d\n\n", n, length(unique(Y)), length(unique(X))))

  # Time V1
  time_v1 <- time_execution(quote(acor:::K_tau_vec_v1(X, Y, tau_XY)),
                            replications = n_replications)$median

  # Time V2
  time_v2 <- time_execution(quote(acor:::K_tau_vec_v2(X, Y, tau_XY)),
                            replications = n_replications)$median

  cat(sprintf("K_tau_vec_v1: %.4f seconds\n", time_v1))
  cat(sprintf("K_tau_vec_v2: %.4f seconds\n", time_v2))

  if (time_v1 < time_v2) {
    cat(sprintf("\n--> V1 is %.1fx faster than V2 for Y with 3 values\n", time_v2 / time_v1))
  } else {
    cat(sprintf("\n--> V2 is %.1fx faster than V1 for Y with 3 values\n", time_v1 / time_v2))
  }

  cat("\n")
  expect_true(TRUE, info = "Y with 3 values comparison completed")
})



test_that("Computational complexity analysis", {
  skip_on_cran()
  cat("\n")
  cat(paste(rep("=", 80), collapse = ""))
  cat("\nCOMPUTATIONAL COMPLEXITY ANALYSIS\n")
  cat(paste(rep("=", 80), collapse = ""))
  cat("\n\n")

  sample_sizes <- c(500, 1000, 2000, 5000, 8000)
  n_replications <- 3

  results <- data.frame(
    n = integer(),
    method = character(),
    time = numeric(),
    stringsAsFactors = FALSE
  )

  cat("Timing K_tau_vec for different sample sizes (continuous data):\n\n")

  for (n in sample_sizes) {
    set.seed(14000 + n)
    X <- rnorm(n) + rnorm(n, sd = 0.0001)
    Y <- rnorm(n) + rnorm(n, sd = 0.0001)

    akc_result <- acor:::kendall_tau_sign(X, Y)
    tau_XY <- akc_result$expectation

    time_orig <- time_execution(quote(K_tau_vec_original(X, Y, tau_XY)),
                                replications = n_replications)$median
    time_v1 <- time_execution(quote(acor:::K_tau_vec_v1(X, Y, tau_XY)),
                              replications = n_replications)$median
    time_v2 <- time_execution(quote(acor:::K_tau_vec_v2(X, Y, tau_XY)),
                              replications = n_replications)$median

    results <- rbind(results,
                     data.frame(n = n, method = "Original", time = time_orig),
                     data.frame(n = n, method = "V1", time = time_v1),
                     data.frame(n = n, method = "V2", time = time_v2))

    cat(sprintf("n = %5d: Original = %.4fs, V1 = %.4fs, V2 = %.4fs\n",
                n, time_orig, time_v1, time_v2))
  }

  cat("\nEstimated complexity exponents (time ~ n^power):\n\n")

  for (method in c("Original", "V1", "V2")) {
    subset_data <- results[results$method == method, ]
    valid_idx <- subset_data$time > 0 & is.finite(subset_data$time)
    subset_data <- subset_data[valid_idx, ]

    if (nrow(subset_data) >= 3) {
      log_model <- lm(log(time) ~ log(n), data = subset_data)
      power <- coef(log_model)[2]
      r_squared <- summary(log_model)$r.squared
      cat(sprintf("  %s: O(n^%.2f), R² = %.4f\n", method, power, r_squared))
    } else {
      cat(sprintf("  %s: Insufficient data for complexity estimation\n", method))
    }
  }

  cat("\nExpected complexities:\n")
  cat("  Original: O(n^2) - nested loops\n")
  cat("  V1: O(R×M + n) - for continuous data R≈M≈n so ~O(n^2), but faster constants\n")
  cat("  V2: O(n log n) - Fenwick tree\n")
  cat("\n")

  expect_true(TRUE, info = "Complexity analysis completed")
})


# ============================================================================
# TEST 11: Univariate Sigma_akc runtime benchmark
# ============================================================================

test_that("Runtime benchmark: Sigma_akc vs Sigma_akc_v1 vs Sigma_akc_v2", {
  skip_on_cran()
  cat("\n")
  cat(paste(rep("=", 80), collapse = ""))
  cat("\nRUNTIME BENCHMARK: Sigma_akc (Univariate)\n")
  cat(paste(rep("=", 80), collapse = ""))
  cat("\n\n")

  n_replications <- 3

  results <- data.frame(
    scenario = character(),
    n = integer(),
    M = integer(),
    M_over_n = numeric(),
    time_orig = numeric(),
    time_v1 = numeric(),
    time_v2 = numeric(),
    winner = character(),
    stringsAsFactors = FALSE
  )

  # -------------------------------------------------------------------------
  # Binary Y
  # -------------------------------------------------------------------------
  cat("Binary Y:\n")
  cat(sprintf("%8s %12s %12s %12s %12s %12s\n",
              "n", "Original", "V1", "V2", "V1 speedup", "V2 speedup"))
  cat(paste(rep("-", 80), collapse = ""))
  cat("\n")

  for (n in c(1000, 2000, 5000, 10000)) {
    set.seed(50000 + n)
    X <- rnorm(n)
    Y <- rbinom(n, 1, 0.6)

    M <- length(unique(Y))

    # Only run original for small n (it's slow)
    if (n <= 1000) {
      time_orig <- time_execution(quote(acor:::Sigma_akc(X, Y)),
                                  replications = n_replications)$median
    } else {
      time_orig <- NA
    }

    time_v1 <- time_execution(quote(Sigma_akc_v1_wrapper(X, Y)),
                              replications = n_replications)$median
    time_v2 <- time_execution(quote(Sigma_akc_v2_wrapper(X, Y)),
                              replications = n_replications)$median

    winner <- if (time_v1 < time_v2) "v1" else "v2"

    speedup_v1 <- if (!is.na(time_orig) && time_v1 > 1e-9) time_orig / time_v1 else NA
    speedup_v2 <- if (!is.na(time_orig) && time_v2 > 1e-9) time_orig / time_v2 else NA

    cat(sprintf("%8d %12.4f %12.4f %12.4f %12s %12s\n",
                n,
                ifelse(is.na(time_orig), NA, time_orig),
                time_v1, time_v2,
                ifelse(is.na(speedup_v1), "-", sprintf("%.1fx", speedup_v1)),
                ifelse(is.na(speedup_v2), "-", sprintf("%.1fx", speedup_v2))))

    results <- rbind(results, data.frame(
      scenario = "binary", n = n, M = M, M_over_n = M/n,
      time_orig = ifelse(is.na(time_orig), NA, time_orig),
      time_v1 = time_v1, time_v2 = time_v2, winner = winner
    ))
  }

  # -------------------------------------------------------------------------
  # Discrete Y with various M values
  # -------------------------------------------------------------------------
  cat("\nDiscrete Y (various M levels):\n")
  cat(sprintf("%8s %8s %12s %12s %12s %12s %12s\n",
              "n", "M", "Original", "V1", "V2", "V1 speedup", "V2 speedup"))
  cat(paste(rep("-", 90), collapse = ""))
  cat("\n")

  M_values <- c(50, 100, 250, 500)

  for (n in c(1000, 2000, 5000, 10000)) {
    for (M_target in M_values) {
      if (M_target >= n) next  # Skip if M >= n

      set.seed(51000 + n + M_target)
      X <- rnorm(n)
      Y <- sample(1:M_target, n, replace = TRUE)

      M <- length(unique(Y))

      # Only run original for small n
      if (n <= 1000) {
        time_orig <- time_execution(quote(acor:::Sigma_akc(X, Y)),
                                    replications = n_replications)$median
      } else {
        time_orig <- NA
      }

      time_v1 <- time_execution(quote(Sigma_akc_v1_wrapper(X, Y)),
                                replications = n_replications)$median
      time_v2 <- time_execution(quote(Sigma_akc_v2_wrapper(X, Y)),
                                replications = n_replications)$median

      winner <- if (time_v1 < time_v2) "v1" else "v2"

      speedup_v1 <- if (!is.na(time_orig) && time_v1 > 1e-9) time_orig / time_v1 else NA
      speedup_v2 <- if (!is.na(time_orig) && time_v2 > 1e-9) time_orig / time_v2 else NA

      cat(sprintf("%8d %8d %12.4f %12.4f %12.4f %12s %12s\n",
                  n, M,
                  ifelse(is.na(time_orig), NA, time_orig),
                  time_v1, time_v2,
                  ifelse(is.na(speedup_v1), "-", sprintf("%.1fx", speedup_v1)),
                  ifelse(is.na(speedup_v2), "-", sprintf("%.1fx", speedup_v2))))

      results <- rbind(results, data.frame(
        scenario = "discrete", n = n, M = M, M_over_n = M/n,
        time_orig = ifelse(is.na(time_orig), NA, time_orig),
        time_v1 = time_v1, time_v2 = time_v2, winner = winner
      ))
    }
  }

  # -------------------------------------------------------------------------
  # Continuous Y
  # -------------------------------------------------------------------------
  cat("\nContinuous Y:\n")
  cat(sprintf("%8s %12s %12s %12s %12s %12s\n",
              "n", "Original", "V1", "V2", "V1 speedup", "V2 speedup"))
  cat(paste(rep("-", 80), collapse = ""))
  cat("\n")

  for (n in c(1000, 2000, 5000)) {
    set.seed(52000 + n)
    X <- rnorm(n)
    Y <- rnorm(n)

    M <- length(unique(Y))

    # Only run original for small n
    if (n <= 1000) {
      time_orig <- time_execution(quote(acor:::Sigma_akc(X, Y)),
                                  replications = n_replications)$median
    } else {
      time_orig <- NA
    }

    time_v1 <- time_execution(quote(Sigma_akc_v1_wrapper(X, Y)),
                              replications = n_replications)$median
    time_v2 <- time_execution(quote(Sigma_akc_v2_wrapper(X, Y)),
                              replications = n_replications)$median

    winner <- if (time_v1 < time_v2) "v1" else "v2"

    speedup_v1 <- if (!is.na(time_orig) && time_v1 > 1e-9) time_orig / time_v1 else NA
    speedup_v2 <- if (!is.na(time_orig) && time_v2 > 1e-9) time_orig / time_v2 else NA

    cat(sprintf("%8d %12.4f %12.4f %12.4f %12s %12s\n",
                n,
                ifelse(is.na(time_orig), NA, time_orig),
                time_v1, time_v2,
                ifelse(is.na(speedup_v1), "-", sprintf("%.1fx", speedup_v1)),
                ifelse(is.na(speedup_v2), "-", sprintf("%.1fx", speedup_v2))))

    results <- rbind(results, data.frame(
      scenario = "continuous", n = n, M = M, M_over_n = M/n,
      time_orig = ifelse(is.na(time_orig), NA, time_orig),
      time_v1 = time_v1, time_v2 = time_v2, winner = winner
    ))
  }

  # -------------------------------------------------------------------------
  # Summary
  # -------------------------------------------------------------------------
  cat("\n")
  cat(paste(rep("-", 80), collapse = ""))
  cat("\nSummary (V1 vs V2):\n")
  cat(sprintf("  V1 wins: %d/%d (%.1f%%)\n",
              sum(results$winner == "v1"), nrow(results),
              100 * mean(results$winner == "v1")))
  cat(sprintf("  V2 wins: %d/%d (%.1f%%)\n",
              sum(results$winner == "v2"), nrow(results),
              100 * mean(results$winner == "v2")))

  cat("\nBy scenario:\n")
  for (sc in unique(results$scenario)) {
    subset <- results[results$scenario == sc, ]
    cat(sprintf("  %s: V1 wins %d/%d\n", sc,
                sum(subset$winner == "v1"), nrow(subset)))
  }

  cat("\n")
  expect_true(TRUE, info = "Univariate Sigma_akc benchmark completed")
})


# ============================================================================
# TEST 12: Multivariate Sigma_akc_multivariate runtime benchmark
# ============================================================================

test_that("Runtime benchmark: Sigma_akc_multivariate vs _v1 vs _v2", {
  skip_on_cran()
  cat("\n")
  cat(paste(rep("=", 80), collapse = ""))
  cat("\nRUNTIME BENCHMARK: Sigma_akc_multivariate\n")
  cat(paste(rep("=", 80), collapse = ""))
  cat("\n\n")

  n_replications <- 3
  m <- 3  # number of predictors

  results <- data.frame(
    scenario = character(),
    n = integer(),
    M = integer(),
    M_over_n = numeric(),
    time_orig = numeric(),
    time_v1 = numeric(),
    time_v2 = numeric(),
    winner = character(),
    stringsAsFactors = FALSE
  )

  # -------------------------------------------------------------------------
  # Binary Y
  # -------------------------------------------------------------------------
  cat("Binary Y (m=3 predictors):\n")
  cat(sprintf("%8s %12s %12s %12s %12s %12s\n",
              "n", "Original", "V1", "V2", "V1 speedup", "V2 speedup"))
  cat(paste(rep("-", 80), collapse = ""))
  cat("\n")

  for (n in c(1000, 2000, 5000, 10000)) {
    set.seed(60000 + n)
    X <- matrix(rnorm(n * m), nrow = n, ncol = m)
    Y <- rbinom(n, 1, 0.6)

    M <- length(unique(Y))

    # Only run original for small n
    if (n <= 500) {
      time_orig <- time_execution(quote(acor:::Sigma_akc_multivariate(X, Y)),
                                  replications = n_replications)$median
    } else {
      time_orig <- NA
    }

    time_v1 <- time_execution(quote(acor:::Sigma_akc_multivariate_v1(X, Y)),
                              replications = n_replications)$median
    time_v2 <- time_execution(quote(acor:::Sigma_akc_multivariate_v2(X, Y)),
                              replications = n_replications)$median

    winner <- if (time_v1 < time_v2) "v1" else "v2"

    speedup_v1 <- if (!is.na(time_orig) && time_v1 > 1e-9) time_orig / time_v1 else NA
    speedup_v2 <- if (!is.na(time_orig) && time_v2 > 1e-9) time_orig / time_v2 else NA

    cat(sprintf("%8d %12.4f %12.4f %12.4f %12s %12s\n",
                n,
                ifelse(is.na(time_orig), NA, time_orig),
                time_v1, time_v2,
                ifelse(is.na(speedup_v1), "-", sprintf("%.1fx", speedup_v1)),
                ifelse(is.na(speedup_v2), "-", sprintf("%.1fx", speedup_v2))))

    results <- rbind(results, data.frame(
      scenario = "binary", n = n, M = M, M_over_n = M/n,
      time_orig = ifelse(is.na(time_orig), NA, time_orig),
      time_v1 = time_v1, time_v2 = time_v2, winner = winner
    ))
  }

  # -------------------------------------------------------------------------
  # Discrete Y with various M values
  # -------------------------------------------------------------------------
  cat("\nDiscrete Y (m=3 predictors, various M levels):\n")
  cat(sprintf("%8s %8s %12s %12s %12s %12s %12s\n",
              "n", "M", "Original", "V1", "V2", "V1 speedup", "V2 speedup"))
  cat(paste(rep("-", 90), collapse = ""))
  cat("\n")

  M_values <- c(50, 100, 250, 500)

  for (n in c(1000, 2000, 5000, 10000)) {
    for (M_target in M_values) {
      if (M_target >= n) next

      set.seed(61000 + n + M_target)
      X <- matrix(rnorm(n * m), nrow = n, ncol = m)
      Y <- sample(1:M_target, n, replace = TRUE)

      M <- length(unique(Y))

      # Only run original for small n
      if (n <= 500) {
        time_orig <- time_execution(quote(acor:::Sigma_akc_multivariate(X, Y)),
                                    replications = n_replications)$median
      } else {
        time_orig <- NA
      }

      time_v1 <- time_execution(quote(acor:::Sigma_akc_multivariate_v1(X, Y)),
                                replications = n_replications)$median
      time_v2 <- time_execution(quote(acor:::Sigma_akc_multivariate_v2(X, Y)),
                                replications = n_replications)$median

      winner <- if (time_v1 < time_v2) "v1" else "v2"

      speedup_v1 <- if (!is.na(time_orig) && time_v1 > 1e-9) time_orig / time_v1 else NA
      speedup_v2 <- if (!is.na(time_orig) && time_v2 > 1e-9) time_orig / time_v2 else NA

      cat(sprintf("%8d %8d %12.4f %12.4f %12.4f %12s %12s\n",
                  n, M,
                  ifelse(is.na(time_orig), NA, time_orig),
                  time_v1, time_v2,
                  ifelse(is.na(speedup_v1), "-", sprintf("%.1fx", speedup_v1)),
                  ifelse(is.na(speedup_v2), "-", sprintf("%.1fx", speedup_v2))))

      results <- rbind(results, data.frame(
        scenario = "discrete", n = n, M = M, M_over_n = M/n,
        time_orig = ifelse(is.na(time_orig), NA, time_orig),
        time_v1 = time_v1, time_v2 = time_v2, winner = winner
      ))
    }
  }

  # -------------------------------------------------------------------------
  # Continuous Y
  # -------------------------------------------------------------------------
  cat("\nContinuous Y (m=3 predictors):\n")
  cat(sprintf("%8s %12s %12s %12s %12s %12s\n",
              "n", "Original", "V1", "V2", "V1 speedup", "V2 speedup"))
  cat(paste(rep("-", 80), collapse = ""))
  cat("\n")

  for (n in c(1000, 2000, 5000)) {
    set.seed(62000 + n)
    X <- matrix(rnorm(n * m), nrow = n, ncol = m)
    Y <- rnorm(n)

    M <- length(unique(Y))

    # Only run original for small n
    if (n <= 500) {
      time_orig <- time_execution(quote(acor:::Sigma_akc_multivariate(X, Y)),
                                  replications = n_replications)$median
    } else {
      time_orig <- NA
    }

    time_v1 <- time_execution(quote(acor:::Sigma_akc_multivariate_v1(X, Y)),
                              replications = n_replications)$median
    time_v2 <- time_execution(quote(acor:::Sigma_akc_multivariate_v2(X, Y)),
                              replications = n_replications)$median

    winner <- if (time_v1 < time_v2) "v1" else "v2"

    speedup_v1 <- if (!is.na(time_orig) && time_v1 > 1e-9) time_orig / time_v1 else NA
    speedup_v2 <- if (!is.na(time_orig) && time_v2 > 1e-9) time_orig / time_v2 else NA

    cat(sprintf("%8d %12.4f %12.4f %12.4f %12s %12s\n",
                n,
                ifelse(is.na(time_orig), NA, time_orig),
                time_v1, time_v2,
                ifelse(is.na(speedup_v1), "-", sprintf("%.1fx", speedup_v1)),
                ifelse(is.na(speedup_v2), "-", sprintf("%.1fx", speedup_v2))))

    results <- rbind(results, data.frame(
      scenario = "continuous", n = n, M = M, M_over_n = M/n,
      time_orig = ifelse(is.na(time_orig), NA, time_orig),
      time_v1 = time_v1, time_v2 = time_v2, winner = winner
    ))
  }

  # -------------------------------------------------------------------------
  # Summary
  # -------------------------------------------------------------------------
  cat("\n")
  cat(paste(rep("-", 80), collapse = ""))
  cat("\nSummary (V1 vs V2):\n")
  cat(sprintf("  V1 wins: %d/%d (%.1f%%)\n",
              sum(results$winner == "v1"), nrow(results),
              100 * mean(results$winner == "v1")))
  cat(sprintf("  V2 wins: %d/%d (%.1f%%)\n",
              sum(results$winner == "v2"), nrow(results),
              100 * mean(results$winner == "v2")))

  cat("\nBy scenario:\n")
  for (sc in unique(results$scenario)) {
    subset <- results[results$scenario == sc, ]
    cat(sprintf("  %s: V1 wins %d/%d\n", sc,
                sum(subset$winner == "v1"), nrow(subset)))
  }

  cat("\n")
  expect_true(TRUE, info = "Multivariate Sigma_akc benchmark completed")
})


# ============================================================================
# TEST 13: Threshold analysis - find optimal M/n cutoff
# ============================================================================

test_that("Threshold analysis: optimal M/n cutoff for V1 vs V2", {
  skip_on_cran()
  cat("\n")
  cat(paste(rep("=", 80), collapse = ""))
  cat("\nTHRESHOLD ANALYSIS: Finding optimal M/n cutoff\n")
  cat(paste(rep("=", 80), collapse = ""))
  cat("\n\n")

  n_replications <- 3

  # Collect results across many M/n ratios
  results <- data.frame(
    n = integer(),
    M = integer(),
    M_over_n = numeric(),
    time_v1 = numeric(),
    time_v2 = numeric(),
    winner = character(),
    stringsAsFactors = FALSE
  )

  cat("Collecting timing data across M/n ratios...\n\n")
  cat(sprintf("%8s %8s %10s %12s %12s %8s\n",
              "n", "M", "M/n", "V1 time", "V2 time", "Winner"))
  cat(paste(rep("-", 70), collapse = ""))
  cat("\n")

  # Test various M/n ratios
  for (n in c(2000, 5000, 10000)) {
    # M values giving different M/n ratios
    M_values <- c(2, 5, 10, 20, 50, 100, 200, 500, round(n * 0.5), round(n * 0.8))
    M_values <- unique(M_values[M_values < n])
    M_values <- sort(M_values)

    for (M_target in M_values) {
      set.seed(70000 + n + M_target)
      X <- rnorm(n)
      Y <- sample(1:M_target, n, replace = TRUE)

      M <- length(unique(Y))

      time_v1 <- time_execution(quote(Sigma_akc_v1_wrapper(X, Y)),
                                replications = n_replications)$median
      time_v2 <- time_execution(quote(Sigma_akc_v2_wrapper(X, Y)),
                                replications = n_replications)$median

      winner <- if (time_v1 < time_v2) "v1" else "v2"

      cat(sprintf("%8d %8d %10.4f %12.4f %12.4f %8s\n",
                  n, M, M/n, time_v1, time_v2, winner))

      results <- rbind(results, data.frame(
        n = n, M = M, M_over_n = M/n,
        time_v1 = time_v1, time_v2 = time_v2, winner = winner
      ))
    }
  }

  # Also add continuous case
  cat("\nContinuous Y cases:\n")
  for (n in c(2000, 5000)) {
    set.seed(71000 + n)
    X <- rnorm(n)
    Y <- rnorm(n)

    M <- length(unique(Y))

    time_v1 <- time_execution(quote(Sigma_akc_v1_wrapper(X, Y)),
                              replications = n_replications)$median
    time_v2 <- time_execution(quote(Sigma_akc_v2_wrapper(X, Y)),
                              replications = n_replications)$median

    winner <- if (time_v1 < time_v2) "v1" else "v2"

    cat(sprintf("%8d %8d %10.4f %12.4f %12.4f %8s\n",
                n, M, M/n, time_v1, time_v2, winner))

    results <- rbind(results, data.frame(
      n = n, M = M, M_over_n = M/n,
      time_v1 = time_v1, time_v2 = time_v2, winner = winner
    ))
  }

  # -------------------------------------------------------------------------
  # Test different thresholds
  # -------------------------------------------------------------------------
  cat("\n")
  cat(paste(rep("-", 70), collapse = ""))
  cat("\nTesting M/n thresholds (M/n < threshold -> use V1):\n\n")
  cat(sprintf("%12s  %10s  %12s  %12s\n",
              "Threshold", "Accuracy", "Mean Regret", "Max Regret"))
  cat(paste(rep("-", 50), collapse = ""))
  cat("\n")

  for (thresh in c(0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7)) {
    predicted <- ifelse(results$M_over_n < thresh, "v1", "v2")
    accuracy <- mean(predicted == results$winner)

    regret <- sapply(1:nrow(results), function(i) {
      chosen <- if (predicted[i] == "v1") results$time_v1[i] else results$time_v2[i]
      optimal <- min(results$time_v1[i], results$time_v2[i])
      (chosen - optimal) / optimal * 100
    })

    cat(sprintf("%12.2f  %9.1f%%  %11.2f%%  %11.2f%%\n",
                thresh, accuracy * 100, mean(regret), max(regret)))
  }

  # -------------------------------------------------------------------------
  # Find crossover point
  # -------------------------------------------------------------------------
  cat("\n")
  cat(paste(rep("-", 70), collapse = ""))
  cat("\nCrossover analysis (where does V2 start winning?):\n\n")

  # Group by M/n bins
  results$M_bin <- cut(results$M_over_n,
                       breaks = c(0, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0),
                       labels = c("0-0.01", "0.01-0.05", "0.05-0.1",
                                  "0.1-0.2", "0.2-0.5", "0.5-1.0"))

  cat(sprintf("%12s  %10s  %10s  %15s\n",
              "M/n range", "V1 wins", "V2 wins", "V1 win rate"))
  cat(paste(rep("-", 50), collapse = ""))
  cat("\n")

  for (bin in levels(results$M_bin)) {
    subset <- results[results$M_bin == bin, ]
    if (nrow(subset) > 0) {
      v1_wins <- sum(subset$winner == "v1")
      v2_wins <- sum(subset$winner == "v2")
      cat(sprintf("%12s  %10d  %10d  %14.1f%%\n",
                  bin, v1_wins, v2_wins, 100 * v1_wins / nrow(subset)))
    }
  }

  cat("\n")
  expect_true(TRUE, info = "Threshold analysis completed")
})