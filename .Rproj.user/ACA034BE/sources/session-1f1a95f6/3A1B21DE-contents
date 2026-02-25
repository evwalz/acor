# ============================================================================
# Tests: AGC Version Consistency (original vs v2) and Runtime Benchmarking
# ============================================================================
#
# Purpose:
#   1. Verify that "original" (kernel_ties_optim2) and "v2" (kernel_agc_v2)
#      produce numerically identical AGC estimates and variances across all
#      data scenarios.
#   2. Print runtimes to validate that v2 (Fenwick tree, O(n log n)) is
#      faster than original for large n continuous data.
#
# The AGC Sigma functions are called directly to force each version.
#
# ============================================================================

library(testthat)
library(acor)  # Replace with your actual package name

# If agc_functions_v2.R is not yet part of the package, source it:
# source("agc_functions_v2.R")

# Tolerance for numerical equivalence across versions
TOL <- 1e-10

# ---- Shorthand wrappers to force a specific version -----------------------

# Univariate IID
run_agc_uni <- function(y_rank, x_rank, version, IID = TRUE) {
  if (IID) {
    if (version == "original") {
      Sigma_agc(y_rank, x_rank)
    } else {  # v2
      Sigma_agc_v2(y_rank, x_rank)
    }
  } else {
    if (version == "original") {
      Sigma_agc_ts(y_rank, x_rank)
    } else {  # v2
      Sigma_agc_ts_v2(y_rank, x_rank)
    }
  }
}

# Multivariate IID / HAC
run_agc_mv <- function(y_rank, xarray_ranks, version, IID = TRUE) {
  if (IID) {
    if (version == "original") {
      Sigma_agc_multivariate(y_rank, xarray_ranks)
    } else {  # v2
      Sigma_agc_multivariate_v2(y_rank, xarray_ranks)
    }
  } else {
    if (version == "original") {
      Sigma_agc_multivariate_ts(y_rank, xarray_ranks)
    } else {  # v2
      Sigma_agc_multivariate_ts_v2(y_rank, xarray_ranks)
    }
  }
}

# Helper: generate ranked inputs (AGC functions operate on ranks)
make_ranks <- function(X, Y) {
  y_rank <- rank(Y, ties.method = "average")
  if (is.vector(X)) {
    x_rank <- rank(X, ties.method = "average")
    return(list(y_rank = y_rank, x_rank = x_rank))
  } else {
    m <- ncol(X)
    n <- nrow(X)
    xarray_ranks <- matrix(0, nrow = m, ncol = n)
    for (i in 1:m) {
      xarray_ranks[i, ] <- rank(X[, i], ties.method = "average")
    }
    return(list(y_rank = y_rank, xarray_ranks = xarray_ranks))
  }
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

test_that("AGC Uni IID | discrete X (10 levels), continuous Y: original vs v2 agree", {
  set.seed(103); n <- 300
  X <- sample(1:10, n, replace = TRUE); Y <- rnorm(n)
  rk <- make_ranks(X, Y)
  
  r_orig <- run_agc_uni(rk$y_rank, rk$x_rank, "original")
  r_v2   <- run_agc_uni(rk$y_rank, rk$x_rank, "v2")
  
  expect_equal(r_v2$agc,     r_orig$agc,     tolerance = TOL)
  expect_equal(r_v2$var,     r_orig$var,     tolerance = TOL)
  expect_equal(r_v2$var_ind, r_orig$var_ind, tolerance = TOL)
})

test_that("AGC Uni IID | discrete X (5 levels) & discrete Y (3 levels): original vs v2 agree", {
  set.seed(104); n <- 300
  X <- sample(1:5, n, replace = TRUE); Y <- sample(1:3, n, replace = TRUE)
  rk <- make_ranks(X, Y)
  
  r_orig <- run_agc_uni(rk$y_rank, rk$x_rank, "original")
  r_v2   <- run_agc_uni(rk$y_rank, rk$x_rank, "v2")
  
  expect_equal(r_v2$agc,     r_orig$agc,     tolerance = TOL)
  expect_equal(r_v2$var,     r_orig$var,     tolerance = TOL)
  expect_equal(r_v2$var_ind, r_orig$var_ind, tolerance = TOL)
})

test_that("AGC Uni IID | heavy ties in X (~25%): original vs v2 agree", {
  set.seed(105); n <- 300
  X <- rnorm(n); X[seq(1, n, by = 4)] <- X[1]; Y <- rnorm(n)
  rk <- make_ranks(X, Y)
  
  r_orig <- run_agc_uni(rk$y_rank, rk$x_rank, "original")
  r_v2   <- run_agc_uni(rk$y_rank, rk$x_rank, "v2")
  
  expect_equal(r_v2$agc,     r_orig$agc,     tolerance = TOL)
  expect_equal(r_v2$var,     r_orig$var,     tolerance = TOL)
  expect_equal(r_v2$var_ind, r_orig$var_ind, tolerance = TOL)
})

test_that("AGC Uni IID | correlated X & Y: original vs v2 agree", {
  set.seed(106); n <- 300
  X <- rnorm(n); Y <- 0.7 * X + 0.3 * rnorm(n)
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

test_that("AGC Uni HAC | discrete X (10) & Y (5): original vs v2 agree", {
  set.seed(203); n <- 300
  X <- sample(1:10, n, replace = TRUE); Y <- sample(1:5, n, replace = TRUE)
  rk <- make_ranks(X, Y)
  
  r_orig <- run_agc_uni(rk$y_rank, rk$x_rank, "original", IID = FALSE)
  r_v2   <- run_agc_uni(rk$y_rank, rk$x_rank, "v2",      IID = FALSE)
  
  expect_equal(r_v2$agc,     r_orig$agc,     tolerance = TOL)
  expect_equal(r_v2$var,     r_orig$var,     tolerance = TOL)
  expect_equal(r_v2$var_ind, r_orig$var_ind, tolerance = TOL)
})

test_that("AGC Uni HAC | correlated X & Y: original vs v2 agree", {
  set.seed(204); n <- 300
  X <- rnorm(n); Y <- 0.7 * X + 0.3 * rnorm(n)
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

test_that("AGC MV IID m=3 | binary Y: original vs v2 agree", {
  set.seed(302); n <- 300
  X <- matrix(rnorm(n * 3), ncol = 3); Y <- rbinom(n, 1, 0.5)
  rk <- make_ranks(X, Y)
  
  r_orig <- run_agc_mv(rk$y_rank, rk$xarray_ranks, "original")
  r_v2   <- run_agc_mv(rk$y_rank, rk$xarray_ranks, "v2")
  
  expect_equal(r_v2$agc_vector, r_orig$agc_vector, tolerance = TOL)
  expect_equal(r_v2$Sigma,      r_orig$Sigma,      tolerance = TOL)
  expect_equal(r_v2$Sigma_ind,  r_orig$Sigma_ind,  tolerance = TOL)
})

test_that("AGC MV IID m=3 | discrete X (10) & Y (5): original vs v2 agree", {
  set.seed(303); n <- 300
  X <- matrix(sample(1:10, n * 3, replace = TRUE), ncol = 3)
  Y <- sample(1:5, n, replace = TRUE)
  rk <- make_ranks(X, Y)
  
  r_orig <- run_agc_mv(rk$y_rank, rk$xarray_ranks, "original")
  r_v2   <- run_agc_mv(rk$y_rank, rk$xarray_ranks, "v2")
  
  expect_equal(r_v2$agc_vector, r_orig$agc_vector, tolerance = TOL)
  expect_equal(r_v2$Sigma,      r_orig$Sigma,      tolerance = TOL)
  expect_equal(r_v2$Sigma_ind,  r_orig$Sigma_ind,  tolerance = TOL)
})

test_that("AGC MV IID m=3 | correlated X & Y: original vs v2 agree", {
  set.seed(304); n <- 300
  X <- matrix(rnorm(n * 3), ncol = 3)
  Y <- 0.5 * X[, 1] + 0.3 * X[, 2] + 0.2 * rnorm(n)
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

test_that("AGC MV HAC m=3 | binary Y: original vs v2 agree", {
  set.seed(402); n <- 300
  X <- matrix(rnorm(n * 3), ncol = 3); Y <- rbinom(n, 1, 0.6)
  rk <- make_ranks(X, Y)
  
  r_orig <- run_agc_mv(rk$y_rank, rk$xarray_ranks, "original", IID = FALSE)
  r_v2   <- run_agc_mv(rk$y_rank, rk$xarray_ranks, "v2",      IID = FALSE)
  
  expect_equal(r_v2$agc_vector, r_orig$agc_vector, tolerance = TOL)
  expect_equal(r_v2$Sigma,      r_orig$Sigma,      tolerance = TOL)
  expect_equal(r_v2$Sigma_ind,  r_orig$Sigma_ind,  tolerance = TOL)
})

test_that("AGC MV HAC m=3 | discrete X (10) & Y (5): original vs v2 agree", {
  set.seed(403); n <- 300
  X <- matrix(sample(1:10, n * 3, replace = TRUE), ncol = 3)
  Y <- sample(1:5, n, replace = TRUE)
  rk <- make_ranks(X, Y)
  
  r_orig <- run_agc_mv(rk$y_rank, rk$xarray_ranks, "original", IID = FALSE)
  r_v2   <- run_agc_mv(rk$y_rank, rk$xarray_ranks, "v2",      IID = FALSE)
  
  expect_equal(r_v2$agc_vector, r_orig$agc_vector, tolerance = TOL)
  expect_equal(r_v2$Sigma,      r_orig$Sigma,      tolerance = TOL)
  expect_equal(r_v2$Sigma_ind,  r_orig$Sigma_ind,  tolerance = TOL)
})


# ============================================================================
# Section 5: Multivariate with m = 5
# ============================================================================

test_that("AGC MV IID m=5 | continuous X & Y: original vs v2 agree", {
  set.seed(501); n <- 200
  X <- matrix(rnorm(n * 5), ncol = 5); Y <- rnorm(n)
  rk <- make_ranks(X, Y)
  
  r_orig <- run_agc_mv(rk$y_rank, rk$xarray_ranks, "original")
  r_v2   <- run_agc_mv(rk$y_rank, rk$xarray_ranks, "v2")
  
  expect_equal(r_v2$agc_vector, r_orig$agc_vector, tolerance = TOL)
  expect_equal(r_v2$Sigma,      r_orig$Sigma,      tolerance = TOL)
  expect_equal(r_v2$Sigma_ind,  r_orig$Sigma_ind,  tolerance = TOL)
})

test_that("AGC MV HAC m=5 | continuous X & Y: original vs v2 agree", {
  set.seed(502); n <- 200
  X <- matrix(rnorm(n * 5), ncol = 5); Y <- rnorm(n)
  rk <- make_ranks(X, Y)
  
  r_orig <- run_agc_mv(rk$y_rank, rk$xarray_ranks, "original", IID = FALSE)
  r_v2   <- run_agc_mv(rk$y_rank, rk$xarray_ranks, "v2",      IID = FALSE)
  
  expect_equal(r_v2$agc_vector, r_orig$agc_vector, tolerance = TOL)
  expect_equal(r_v2$Sigma,      r_orig$Sigma,      tolerance = TOL)
  expect_equal(r_v2$Sigma_ind,  r_orig$Sigma_ind,  tolerance = TOL)
})


# ============================================================================
# Section 6: Runtime benchmarks
# ============================================================================
# Not testthat assertions — prints a formatted table of median runtimes.
# Run interactively or via Rscript.
#
# Expected pattern:
#   - "original" (kernel_ties_optim2) builds R×M sign matrices → O(R*M*n)
#     which is O(n²) for continuous data but fast for discrete
#   - "v2" (kernel_agc_v2 with Fenwick tree) → O(n log n)
#     faster for large n continuous, comparable or slower for discrete
# ============================================================================

run_agc_benchmarks <- function(reps = 5) {
  
  time_med <- function(expr_fn, reps) {
    median(replicate(reps, system.time(expr_fn())["elapsed"]))
  }
  
  # ---- Univariate scenarios --------------------------------------------------
  scenarios <- list(
    # --- Continuous data: watch original degrade vs v2 ---
    list(label = "n=200,   continuous",        seed = 601, n = 200,
         gen_X = function(n) rnorm(n),
         gen_Y = function(n) rnorm(n)),
    
    list(label = "n=500,   continuous",        seed = 602, n = 500,
         gen_X = function(n) rnorm(n),
         gen_Y = function(n) rnorm(n)),
    
    list(label = "n=2000,  continuous",        seed = 603, n = 2000,
         gen_X = function(n) rnorm(n),
         gen_Y = function(n) rnorm(n)),
    
    list(label = "n=5000,  continuous",        seed = 604, n = 5000,
         gen_X = function(n) rnorm(n),
         gen_Y = function(n) rnorm(n)),
    
    # --- Binary Y ---
    list(label = "n=500,   binary Y",          seed = 606, n = 500,
         gen_X = function(n) rnorm(n),
         gen_Y = function(n) rbinom(n, 1, 0.5)),
    
    list(label = "n=5000,  binary Y",          seed = 607, n = 5000,
         gen_X = function(n) rnorm(n),
         gen_Y = function(n) rbinom(n, 1, 0.5)),
    
    # --- Discrete X ---
    list(label = "n=500,   discrete X (10 lv)",     seed = 608, n = 500,
         gen_X = function(n) sample(1:10, n, replace = TRUE),
         gen_Y = function(n) rnorm(n)),
    
    list(label = "n=500,   discrete X+Y (10x5)",    seed = 609, n = 500,
         gen_X = function(n) sample(1:10, n, replace = TRUE),
         gen_Y = function(n) sample(1:5, n, replace = TRUE)),
    
    list(label = "n=5000,  discrete X+Y (10x5)",    seed = 610, n = 5000,
         gen_X = function(n) sample(1:10, n, replace = TRUE),
         gen_Y = function(n) sample(1:5, n, replace = TRUE)),
    
    list(label = "n=10000, discrete X+Y (10x5)",    seed = 611, n = 10000,
         gen_X = function(n) sample(1:10, n, replace = TRUE),
         gen_Y = function(n) sample(1:5, n, replace = TRUE)),
    
    # --- Correlated ---
    list(label = "n=2000,  correlated (r~0.7)",     seed = 612, n = 2000,
         gen_X = function(n) rnorm(n),
         gen_Y = function(n) { x <- rnorm(n); 0.7 * x + 0.3 * rnorm(n) })
  )
  
  # ---- Print header ----------------------------------------------------------
  cat("\n")
  cat(strrep("=", 68), "\n")
  cat("  AGC Runtime Benchmarks: original vs v2  (IID, univariate)\n")
  cat(sprintf("  Median of %d replications — times in seconds\n", reps))
  cat(strrep("=", 68), "\n")
  cat(sprintf("  %-36s  %10s  %10s\n", "Scenario", "original", "v2"))
  cat(strrep("-", 68), "\n")
  
  for (sc in scenarios) {
    set.seed(sc$seed)
    X <- sc$gen_X(sc$n)
    Y <- sc$gen_Y(sc$n)
    rk <- make_ranks(X, Y)
    
    t_orig <- time_med(function() run_agc_uni(rk$y_rank, rk$x_rank, "original"), reps)
    t_v2   <- time_med(function() run_agc_uni(rk$y_rank, rk$x_rank, "v2"),      reps)
    
    # Mark faster version
    if (t_orig <= t_v2) {
      orig_str <- sprintf("%9.4f*", t_orig)
      v2_str   <- sprintf("%10.4f", t_v2)
    } else {
      orig_str <- sprintf("%10.4f", t_orig)
      v2_str   <- sprintf("%9.4f*", t_v2)
    }
    
    cat(sprintf("  %-36s  %s  %s\n", sc$label, orig_str, v2_str))
  }
  
  cat(strrep("-", 68), "\n")
  cat("  * = faster version\n\n")
  
  
  # ---- Multivariate scenarios ------------------------------------------------
  mv_scenarios <- list(
    list(label = "n=300,  m=3, continuous",   seed = 701, n = 300, m = 3,
         gen_X = function(n, m) matrix(rnorm(n * m), ncol = m),
         gen_Y = function(n) rnorm(n)),
    
    list(label = "n=1000, m=3, continuous",   seed = 702, n = 1000, m = 3,
         gen_X = function(n, m) matrix(rnorm(n * m), ncol = m),
         gen_Y = function(n) rnorm(n)),
    
    list(label = "n=5000, m=3, continuous",   seed = 703, n = 5000, m = 3,
         gen_X = function(n, m) matrix(rnorm(n * m), ncol = m),
         gen_Y = function(n) rnorm(n)),
    
    list(label = "n=300,  m=3, binary Y",     seed = 704, n = 300, m = 3,
         gen_X = function(n, m) matrix(rnorm(n * m), ncol = m),
         gen_Y = function(n) rbinom(n, 1, 0.5)),
    
    list(label = "n=1000, m=3, binary Y",     seed = 705, n = 1000, m = 3,
         gen_X = function(n, m) matrix(rnorm(n * m), ncol = m),
         gen_Y = function(n) rbinom(n, 1, 0.5)),
    
    list(label = "n=300,  m=5, continuous",   seed = 706, n = 300, m = 5,
         gen_X = function(n, m) matrix(rnorm(n * m), ncol = m),
         gen_Y = function(n) rnorm(n)),
    
    list(label = "n=1000, m=5, continuous",   seed = 707, n = 1000, m = 5,
         gen_X = function(n, m) matrix(rnorm(n * m), ncol = m),
         gen_Y = function(n) rnorm(n))
  )
  
  cat(strrep("=", 68), "\n")
  cat("  AGC Runtime Benchmarks: original vs v2  (IID, multivariate)\n")
  cat(sprintf("  Median of %d replications — times in seconds\n", reps))
  cat(strrep("=", 68), "\n")
  cat(sprintf("  %-36s  %10s  %10s\n", "Scenario", "original", "v2"))
  cat(strrep("-", 68), "\n")
  
  for (sc in mv_scenarios) {
    set.seed(sc$seed)
    X <- sc$gen_X(sc$n, sc$m)
    Y <- sc$gen_Y(sc$n)
    rk <- make_ranks(X, Y)
    
    t_orig <- time_med(function() run_agc_mv(rk$y_rank, rk$xarray_ranks, "original"), reps)
    t_v2   <- time_med(function() run_agc_mv(rk$y_rank, rk$xarray_ranks, "v2"),      reps)
    
    if (t_orig <= t_v2) {
      orig_str <- sprintf("%9.4f*", t_orig)
      v2_str   <- sprintf("%10.4f", t_v2)
    } else {
      orig_str <- sprintf("%10.4f", t_orig)
      v2_str   <- sprintf("%9.4f*", t_v2)
    }
    
    cat(sprintf("  %-36s  %s  %s\n", sc$label, orig_str, v2_str))
  }
  
  cat(strrep("-", 68), "\n")
  cat("  * = faster version\n\n")
  
  invisible(NULL)
}

# Run benchmarks when script is sourced directly (not during testthat)
if (!interactive() || identical(Sys.getenv("RUN_BENCHMARKS"), "1")) {
  run_agc_benchmarks()
}


# # ============================================================================
# # Section 7: Variance Validation — v2 against known references and original
# # ============================================================================
# #
# # These tests validate the VARIANCE estimator from the v2 kernel against:
# #   (a) Known analytical values (var_ind = 1 for continuous data)
# #   (b) cor.test(method = "spearman") p-values, which implicitly validate
# #       the independence variance (since p-values are derived from it)
# #   (c) The original kernel at large n (the regime v2 is designed for),
# #       with relaxed tolerance to separate numerical noise from logic bugs
# #
# # This catches bugs in the v2 kernel that version-vs-version tests at
# # small n (Sections 1–5) might miss.
# # ============================================================================
# 
# # P-values use different asymptotics (normal vs t/AS89), so allow tolerance
# TOL_PVAL <- 0.01
# 
# 
# # ---- 7a: Independence variance = 1 for continuous data ----
# # For continuous X & Y (no ties): zeta_3X = zeta_3Y = 1/n^2
# # so var_ind = (1 - 1/n^2) / (1 - 1/n^2) = 1 exactly.
# 
# test_that("v2: var_ind = 1 for continuous X & Y | n=300", {
#   set.seed(811); n <- 300
#   X <- rnorm(n); Y <- rnorm(n)
#   rk <- make_ranks(X, Y)
#   
#   r_v2 <- run_agc_uni(rk$y_rank, rk$x_rank, "v2")
#   expect_equal(r_v2$var_ind, 1, tolerance = 1e-12)
# })
# 
# test_that("v2: var_ind = 1 for continuous X & Y | n=5000", {
#   set.seed(812); n <- 5000
#   X <- rnorm(n); Y <- rnorm(n)
#   rk <- make_ranks(X, Y)
#   
#   r_v2 <- run_agc_uni(rk$y_rank, rk$x_rank, "v2")
#   expect_equal(r_v2$var_ind, 1, tolerance = 1e-12)
# })
# 
# test_that("v2: var_ind = 1 for continuous X & Y | n=10000", {
#   set.seed(813); n <- 10000
#   X <- rnorm(n); Y <- rnorm(n)
#   rk <- make_ranks(X, Y)
#   
#   r_v2 <- run_agc_uni(rk$y_rank, rk$x_rank, "v2")
#   expect_equal(r_v2$var_ind, 1, tolerance = 1e-12)
# })
# 
# 
# # ---- 7b: P-value from var_ind matches cor.test(spearman) ----
# # Under independence with continuous data, the p-value derived from var_ind
# # in acor.test should closely match the Spearman asymptotic p-value.
# # This validates the full variance → se → z-stat → p-value chain.
# 
# test_that("v2: ind p-value matches Spearman | n=500, two.sided", {
#   set.seed(821); n <- 500
#   X <- rnorm(n); Y <- rnorm(n)
#   
#   agc_result <- acor.test(X, Y, method = "agc", alternative = "two.sided")
#   sp_result  <- cor.test(X, Y, method = "spearman", exact = FALSE)
#   
#   expect_equal(agc_result$p.value_ind, sp_result$p.value, tolerance = TOL_PVAL)
# })
# 
# test_that("v2: ind p-value matches Spearman | n=2000, two.sided", {
#   set.seed(822); n <- 2000
#   X <- rnorm(n); Y <- rnorm(n)
#   
#   agc_result <- acor.test(X, Y, method = "agc", alternative = "two.sided")
#   sp_result  <- cor.test(X, Y, method = "spearman", exact = FALSE)
#   
#   expect_equal(agc_result$p.value_ind, sp_result$p.value, tolerance = TOL_PVAL)
# })
# 
# test_that("v2: ind p-value matches Spearman | n=5000, two.sided", {
#   set.seed(823); n <- 5000
#   X <- rnorm(n); Y <- rnorm(n)
#   
#   agc_result <- acor.test(X, Y, method = "agc", alternative = "two.sided")
#   sp_result  <- cor.test(X, Y, method = "spearman", exact = FALSE)
#   
#   expect_equal(agc_result$p.value_ind, sp_result$p.value, tolerance = TOL_PVAL)
# })
# 
# test_that("v2: ind p-value matches Spearman | n=1000, one-sided greater", {
#   set.seed(831); n <- 1000
#   X <- rnorm(n); Y <- rnorm(n)
#   
#   agc_result <- acor.test(X, Y, method = "agc", alternative = "greater")
#   sp_result  <- cor.test(X, Y, method = "spearman", alternative = "greater",
#                          exact = FALSE)
#   
#   expect_equal(agc_result$p.value_ind, sp_result$p.value, tolerance = TOL_PVAL)
# })
# 
# test_that("v2: ind p-value matches Spearman | n=1000, one-sided less", {
#   set.seed(832); n <- 1000
#   X <- rnorm(n); Y <- rnorm(n)
#   
#   agc_result <- acor.test(X, Y, method = "agc", alternative = "less")
#   sp_result  <- cor.test(X, Y, method = "spearman", alternative = "less",
#                          exact = FALSE)
#   
#   expect_equal(agc_result$p.value_ind, sp_result$p.value, tolerance = TOL_PVAL)
# })
# 
# 
# # ---- 7c: Full variance (var) — original vs v2 at large n ----
# # The full asymptotic variance (not just var_ind) is the main output of the
# # kernel.  At large n with continuous data, the original builds n×n sign
# # matrices while v2 uses Fenwick trees — different numerics, same answer.
# # Relaxed tolerance (1e-6) for var to allow floating-point accumulation
# # differences; tight tolerance for var_ind (analytically determined).
# 
# test_that("v2 var matches original | n=2000, continuous, IID", {
#   set.seed(841); n <- 2000
#   X <- rnorm(n); Y <- rnorm(n)
#   rk <- make_ranks(X, Y)
#   
#   r_orig <- run_agc_uni(rk$y_rank, rk$x_rank, "original")
#   r_v2   <- run_agc_uni(rk$y_rank, rk$x_rank, "v2")
#   
#   expect_equal(r_v2$var,     r_orig$var,     tolerance = 1e-6)
#   expect_equal(r_v2$var_ind, r_orig$var_ind, tolerance = 1e-10)
# })
# 
# test_that("v2 var matches original | n=5000, continuous, IID", {
#   set.seed(842); n <- 5000
#   X <- rnorm(n); Y <- rnorm(n)
#   rk <- make_ranks(X, Y)
#   
#   r_orig <- run_agc_uni(rk$y_rank, rk$x_rank, "original")
#   r_v2   <- run_agc_uni(rk$y_rank, rk$x_rank, "v2")
#   
#   expect_equal(r_v2$var,     r_orig$var,     tolerance = 1e-6)
#   expect_equal(r_v2$var_ind, r_orig$var_ind, tolerance = 1e-10)
# })
# 
# test_that("v2 var matches original | n=10000, continuous, IID", {
#   set.seed(843); n <- 10000
#   X <- rnorm(n); Y <- rnorm(n)
#   rk <- make_ranks(X, Y)
#   
#   r_orig <- run_agc_uni(rk$y_rank, rk$x_rank, "original")
#   r_v2   <- run_agc_uni(rk$y_rank, rk$x_rank, "v2")
#   
#   expect_equal(r_v2$var,     r_orig$var,     tolerance = 1e-6)
#   expect_equal(r_v2$var_ind, r_orig$var_ind, tolerance = 1e-10)
# })
# 
# test_that("v2 var matches original | n=5000, correlated, IID", {
#   set.seed(844); n <- 5000
#   X <- rnorm(n); Y <- 0.7 * X + 0.3 * rnorm(n)
#   rk <- make_ranks(X, Y)
#   
#   r_orig <- run_agc_uni(rk$y_rank, rk$x_rank, "original")
#   r_v2   <- run_agc_uni(rk$y_rank, rk$x_rank, "v2")
#   
#   expect_equal(r_v2$var,     r_orig$var,     tolerance = 1e-6)
#   expect_equal(r_v2$var_ind, r_orig$var_ind, tolerance = 1e-10)
# })
# 
# test_that("v2 var matches original | n=5000, continuous, HAC", {
#   set.seed(845); n <- 5000
#   X <- rnorm(n); Y <- rnorm(n)
#   rk <- make_ranks(X, Y)
#   
#   r_orig <- run_agc_uni(rk$y_rank, rk$x_rank, "original", IID = FALSE)
#   r_v2   <- run_agc_uni(rk$y_rank, rk$x_rank, "v2",      IID = FALSE)
#   
#   expect_equal(r_v2$var,     r_orig$var,     tolerance = 1e-6)
#   expect_equal(r_v2$var_ind, r_orig$var_ind, tolerance = 1e-6)
# })
# 
# 
# # ---- 7d: Multivariate Sigma — original vs v2 at large n ----
# 
# test_that("v2 Sigma matches original | MV m=3, n=2000, continuous, IID", {
#   set.seed(851); n <- 2000
#   X <- matrix(rnorm(n * 3), ncol = 3); Y <- rnorm(n)
#   rk <- make_ranks(X, Y)
#   
#   r_orig <- run_agc_mv(rk$y_rank, rk$xarray_ranks, "original")
#   r_v2   <- run_agc_mv(rk$y_rank, rk$xarray_ranks, "v2")
#   
#   expect_equal(r_v2$Sigma,     r_orig$Sigma,     tolerance = 1e-6)
#   expect_equal(r_v2$Sigma_ind, r_orig$Sigma_ind, tolerance = 1e-10)
# })
# 
# test_that("v2 Sigma matches original | MV m=3, n=5000, continuous, IID", {
#   set.seed(852); n <- 5000
#   X <- matrix(rnorm(n * 3), ncol = 3); Y <- rnorm(n)
#   rk <- make_ranks(X, Y)
#   
#   r_orig <- run_agc_mv(rk$y_rank, rk$xarray_ranks, "original")
#   r_v2   <- run_agc_mv(rk$y_rank, rk$xarray_ranks, "v2")
#   
#   expect_equal(r_v2$Sigma,     r_orig$Sigma,     tolerance = 1e-6)
#   expect_equal(r_v2$Sigma_ind, r_orig$Sigma_ind, tolerance = 1e-10)
# })
# 
# test_that("v2 Sigma matches original | MV m=3, n=5000, continuous, HAC", {
#   set.seed(853); n <- 5000
#   X <- matrix(rnorm(n * 3), ncol = 3); Y <- rnorm(n)
#   rk <- make_ranks(X, Y)
#   
#   r_orig <- run_agc_mv(rk$y_rank, rk$xarray_ranks, "original", IID = FALSE)
#   r_v2   <- run_agc_mv(rk$y_rank, rk$xarray_ranks, "v2",      IID = FALSE)
#   
#   expect_equal(r_v2$Sigma,     r_orig$Sigma,     tolerance = 1e-6)
#   expect_equal(r_v2$Sigma_ind, r_orig$Sigma_ind, tolerance = 1e-6)
# })
# 
# test_that("v2 Sigma matches original | MV m=5, n=2000, continuous, IID", {
#   set.seed(854); n <- 2000
#   X <- matrix(rnorm(n * 5), ncol = 5); Y <- rnorm(n)
#   rk <- make_ranks(X, Y)
#   
#   r_orig <- run_agc_mv(rk$y_rank, rk$xarray_ranks, "original")
#   r_v2   <- run_agc_mv(rk$y_rank, rk$xarray_ranks, "v2")
#   
#   expect_equal(r_v2$Sigma,     r_orig$Sigma,     tolerance = 1e-6)
#   expect_equal(r_v2$Sigma_ind, r_orig$Sigma_ind, tolerance = 1e-10)
# })
# 
# 
