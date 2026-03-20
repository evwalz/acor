# ============================================================================
# Tests: AKC Kernel Consistency — Auto-dispatch vs Original Reference
# ============================================================================
#
# Purpose:
#   1. Verify that the auto-dispatched optimized kernels (v1/v2, selected
#      internally by select_kernel_version()) produce numerically identical
#      estimates and variances to the O(n^2) original reference kernels
#      defined in helper-original-kernels.R.
#   2. Print runtimes to validate that the optimized code is faster than
#      the original reference.
#
# The original (point-wise O(n^2)) kernels in helper-original-kernels.R
# serve as the ground-truth reference. The auto-dispatched functions
# are called without a version argument.
#
# ============================================================================

library(testthat)
library(acor) 

# Tolerance for numerical equivalence across versions
TOL <- 1e-10

# Shorthand wrappers for auto-dispatched functions
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

test_that("Uni IID | discrete X (10 levels), continuous Y: auto-dispatch matches original", {
  
  set.seed(103); n <- 300
  X <- sample(1:10, n, replace = TRUE); Y <- rnorm(n)
  
  r_orig <- Sigma_akc(X, Y)
  r_auto <- run_uni(X, Y)
  
  expect_equal(r_auto$akc,     r_orig$akc,     tolerance = TOL)
  expect_equal(r_auto$var,     r_orig$var,     tolerance = TOL)
  expect_equal(r_auto$var_ind, r_orig$var_ind, tolerance = TOL)
})

test_that("Uni IID | discrete X (5 levels) & discrete Y (3 levels): auto-dispatch matches original", {
  
  set.seed(104); n <- 300
  X <- sample(1:5, n, replace = TRUE); Y <- sample(1:3, n, replace = TRUE)
  
  r_orig <- Sigma_akc(X, Y)
  r_auto <- run_uni(X, Y)
  
  expect_equal(r_auto$akc,     r_orig$akc,     tolerance = TOL)
  expect_equal(r_auto$var,     r_orig$var,     tolerance = TOL)
  expect_equal(r_auto$var_ind, r_orig$var_ind, tolerance = TOL)
})

test_that("Uni IID | heavy ties in X (~25%): auto-dispatch matches original", {
  
  set.seed(105); n <- 300
  X <- rnorm(n); X[seq(1, n, by = 4)] <- X[1]; Y <- rnorm(n)
  
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

test_that("Uni HAC | binary Y: auto-dispatch matches original", {
  set.seed(202); n <- 300
  X <- rnorm(n); Y <- rbinom(n, 1, 0.6)
  
  r_orig <- Sigma_akc_ts(X, Y)
  r_auto <- run_uni(X, Y, IID = FALSE)
  
  expect_equal(r_auto$akc,     r_orig$akc,     tolerance = TOL)
  expect_equal(r_auto$var,     r_orig$var,     tolerance = TOL)
  expect_equal(r_auto$var_ind, r_orig$var_ind, tolerance = TOL)
})

test_that("Uni HAC | discrete X (10) & Y (5): auto-dispatch matches original", {
  set.seed(203); n <- 300
  X <- sample(1:10, n, replace = TRUE); Y <- sample(1:5, n, replace = TRUE)
  
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

test_that("MV IID m=3 | binary Y: auto-dispatch matches original", {
  set.seed(302); n <- 300
  X <- matrix(rnorm(n * 3), ncol = 3); Y <- rbinom(n, 1, 0.5)
  
  r_orig <- Sigma_akc_multivariate(X, Y)
  r_auto <- run_mv(X, Y)
  
  expect_equal(r_auto$akc_vector, r_orig$akc_vector, tolerance = TOL)
  expect_equal(r_auto$Sigma,      r_orig$Sigma,      tolerance = TOL)
  expect_equal(r_auto$Sigma_ind,  r_orig$Sigma_ind,  tolerance = TOL)
})

test_that("MV IID m=3 | discrete X (10) & Y (5): auto-dispatch matches original", {
  set.seed(303); n <- 300
  X <- matrix(sample(1:10, n * 3, replace = TRUE), ncol = 3)
  Y <- sample(1:5, n, replace = TRUE)
  
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

test_that("MV HAC m=3 | binary Y: auto-dispatch matches original", {
  set.seed(402); n <- 300
  X <- matrix(rnorm(n * 3), ncol = 3); Y <- rbinom(n, 1, 0.6)
  
  r_orig <- Sigma_akc_multivariate_ts(X, Y)
  r_auto <- run_mv(X, Y, IID = FALSE)
  
  expect_equal(r_auto$akc_vector, r_orig$akc_vector, tolerance = TOL)
  expect_equal(r_auto$Sigma,      r_orig$Sigma,      tolerance = TOL)
  expect_equal(r_auto$Sigma_ind,  r_orig$Sigma_ind,  tolerance = TOL)
})

test_that("MV HAC m=3 | discrete X (10) & Y (5): auto-dispatch matches original", {
  set.seed(403); n <- 300
  X <- matrix(sample(1:10, n * 3, replace = TRUE), ncol = 3)
  Y <- sample(1:5, n, replace = TRUE)
  
  r_orig <- Sigma_akc_multivariate_ts(X, Y)
  r_auto <- run_mv(X, Y, IID = FALSE)
  
  expect_equal(r_auto$akc_vector, r_orig$akc_vector, tolerance = TOL)
  expect_equal(r_auto$Sigma,      r_orig$Sigma,      tolerance = TOL)
  expect_equal(r_auto$Sigma_ind,  r_orig$Sigma_ind,  tolerance = TOL)
})


# ============================================================================
# Section 5: Runtime benchmarks
# ============================================================================
# Not testthat assertions — prints a formatted table of median runtimes.
# Run interactively to compare auto-dispatched (optimized) vs original.
#
# Expected pattern:
#   - Auto-dispatched should be faster or equal for all scenarios
#   - Advantage grows with n (original is O(n^2))
# ============================================================================

run_benchmarks <- function(reps = 5) {
  
  time_med <- function(expr_fn, reps) {
    median(replicate(reps, system.time(expr_fn())["elapsed"]))
  }
  
  # ---- Scenario definitions ------------------------------------------------
  scenarios <- list(
    list(label = "n=200,   continuous",     seed = 501, n = 200,
         gen_X = function(n) rnorm(n),
         gen_Y = function(n) rnorm(n)),
    
    list(label = "n=500,   continuous",     seed = 502, n = 500,
         gen_X = function(n) rnorm(n),
         gen_Y = function(n) rnorm(n)),
    
    list(label = "n=2000,  continuous",     seed = 503, n = 2000,
         gen_X = function(n) rnorm(n),
         gen_Y = function(n) rnorm(n)),
    
    list(label = "n=5000,  continuous",     seed = 504, n = 5000,
         gen_X = function(n) rnorm(n),
         gen_Y = function(n) rnorm(n)),
    
    list(label = "n=200,   binary Y",       seed = 506, n = 200,
         gen_X = function(n) rnorm(n),
         gen_Y = function(n) rbinom(n, 1, 0.5)),
    
    list(label = "n=500,   binary Y",       seed = 507, n = 500,
         gen_X = function(n) rnorm(n),
         gen_Y = function(n) rbinom(n, 1, 0.5)),
    
    list(label = "n=5000,  binary Y",       seed = 508, n = 5000,
         gen_X = function(n) rnorm(n),
         gen_Y = function(n) rbinom(n, 1, 0.5)),
    
    list(label = "n=500,   discrete X (10 lv)", seed = 509, n = 500,
         gen_X = function(n) sample(1:10, n, replace = TRUE),
         gen_Y = function(n) rnorm(n)),
    
    list(label = "n=500,   discrete X+Y (10x5)", seed = 510, n = 500,
         gen_X = function(n) sample(1:10, n, replace = TRUE),
         gen_Y = function(n) sample(1:5, n, replace = TRUE)),
    
    list(label = "n=5000,  discrete X+Y (10x5)", seed = 511, n = 5000,
         gen_X = function(n) sample(1:10, n, replace = TRUE),
         gen_Y = function(n) sample(1:5, n, replace = TRUE)),
    
    list(label = "n=10000, discrete X+Y (10x5)", seed = 511, n = 10000,
         gen_X = function(n) sample(1:10, n, replace = TRUE),
         gen_Y = function(n) sample(1:5, n, replace = TRUE))
  )
  
  # ---- Print header --------------------------------------------------------
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("  Runtime Benchmarks: original (O(n^2)) vs auto-dispatch (IID, univariate)\n")
  cat(sprintf("  Median of %d replications — times in seconds\n", reps))
  cat(strrep("=", 70), "\n")
  cat(sprintf("  %-32s  %-10s  %8s  %8s\n",
              "Scenario", "selected", "original", "auto"))
  cat(strrep("-", 70), "\n")
  
  for (sc in scenarios) {
    set.seed(sc$seed)
    X <- sc$gen_X(sc$n)
    Y <- sc$gen_Y(sc$n)
    
    selected <- acor:::select_kernel_version(Y, X)
    
    t_orig <- time_med(function() Sigma_akc(X, Y), reps)
    t_auto <- time_med(function() run_uni(X, Y), reps)
    
    cat(sprintf("  %-32s  %-10s  %8.4f  %8.4f\n",
                sc$label, selected, t_orig, t_auto))
  }
  
  cat(strrep("-", 70), "\n\n")
  
  
  # ---- Multivariate section ------------------------------------------------
  mv_scenarios <- list(
    list(label = "n=300,  m=3, continuous",  seed = 601, n = 300, m = 3,
         gen_X = function(n, m) matrix(rnorm(n * m), ncol = m),
         gen_Y = function(n) rnorm(n)),
    
    list(label = "n=1000, m=3, continuous",  seed = 602, n = 1000, m = 3,
         gen_X = function(n, m) matrix(rnorm(n * m), ncol = m),
         gen_Y = function(n) rnorm(n)),
    
    list(label = "n=5000, m=3, continuous",  seed = 603, n = 5000, m = 3,
         gen_X = function(n, m) matrix(rnorm(n * m), ncol = m),
         gen_Y = function(n) rnorm(n)),
    
    list(label = "n=300,  m=3, binary Y",    seed = 604, n = 300, m = 3,
         gen_X = function(n, m) matrix(rnorm(n * m), ncol = m),
         gen_Y = function(n) rbinom(n, 1, 0.5)),
    
    list(label = "n=1000, m=3, binary Y",    seed = 605, n = 1000, m = 3,
         gen_X = function(n, m) matrix(rnorm(n * m), ncol = m),
         gen_Y = function(n) rbinom(n, 1, 0.5)),
    
    list(label = "n=300,  m=5, continuous",  seed = 606, n = 300, m = 5,
         gen_X = function(n, m) matrix(rnorm(n * m), ncol = m),
         gen_Y = function(n) rnorm(n))
  )
  
  cat(strrep("=", 70), "\n")
  cat("  Runtime Benchmarks: original vs auto-dispatch (IID, multivariate)\n")
  cat(sprintf("  Median of %d replications — times in seconds\n", reps))
  cat(strrep("=", 70), "\n")
  cat(sprintf("  %-32s  %-10s  %8s  %8s\n",
              "Scenario", "selected", "original", "auto"))
  cat(strrep("-", 70), "\n")
  
  for (sc in mv_scenarios) {
    set.seed(sc$seed)
    X <- sc$gen_X(sc$n, sc$m)
    Y <- sc$gen_Y(sc$n)
    
    selected <- acor:::select_kernel_version(Y, X)
    
    t_orig <- time_med(function() Sigma_akc_multivariate(X, Y), reps)
    t_auto <- time_med(function() run_mv(X, Y), reps)
    
    cat(sprintf("  %-32s  %-10s  %8.4f  %8.4f\n",
                sc$label, selected, t_orig, t_auto))
  }
  
  cat(strrep("-", 70), "\n\n")
  
  invisible(NULL)
}

# Run benchmarks when script is sourced directly (not during testthat)

#if (interactive() && identical(Sys.getenv("RUN_BENCHMARKS"), "1")) {
#run_benchmarks()
#}


# ============================================================================
# Section 6: acor.test() end-to-end benchmarks vs reference implementations
# ============================================================================
# Compares wall-clock time of acor.test() (AKC and AGC) against:
#   - pROC::roc() + pROC::var()     for binary Y   (DeLong)
#   - survival::concordance()       for continuous/discrete Y
#
# The AKC kernel version selected by select_kernel_version() is printed
# alongside each scenario so you can see which code path is exercised.
#
# Requirements: pROC package (binary Y rows are skipped if not available)
# ============================================================================

run_acor_benchmarks <- function(reps = 5) {
  
  has_pROC <- requireNamespace("pROC", quietly = TRUE)
  
  time_med <- function(expr_fn, reps) {
    median(replicate(reps, system.time(expr_fn())["elapsed"]))
  }
  
  # ---- Scenario definitions ------------------------------------------------
  scenarios <- list(
    
    list(label = "n=200,  binary Y,  cont X",
         seed = 701, n = 200,  Y_type = "binary",
         X_levels = NA),
    
    list(label = "n=500,  binary Y,  cont X",
         seed = 702, n = 500,  Y_type = "binary",
         X_levels = NA),
    
    list(label = "n=2000, binary Y,  cont X",
         seed = 703, n = 2000, Y_type = "binary",
         X_levels = NA),
    
    list(label = "n=5000, binary Y,  cont X",
         seed = 704, n = 5000, Y_type = "binary",
         X_levels = NA),
    
    list(label = "n=500,  binary Y,  disc X (10 lv)",
         seed = 705, n = 500,  Y_type = "binary",
         X_levels = 10),
    
    list(label = "n=2000, binary Y,  disc X (10 lv)",
         seed = 706, n = 2000, Y_type = "binary",
         X_levels = 10),
    
    list(label = "n=200,  cont Y,    cont X",
         seed = 711, n = 200,  Y_type = "continuous",
         X_levels = NA),
    
    list(label = "n=500,  cont Y,    cont X",
         seed = 712, n = 500,  Y_type = "continuous",
         X_levels = NA),
    
    list(label = "n=2000, cont Y,    cont X",
         seed = 713, n = 2000, Y_type = "continuous",
         X_levels = NA),
    
    list(label = "n=5000, cont Y,    cont X",
         seed = 714, n = 5000, Y_type = "continuous",
         X_levels = NA),
    
    list(label = "n=500,  cont Y,    disc X (10 lv)",
         seed = 721, n = 500,  Y_type = "continuous",
         X_levels = 10),
    
    list(label = "n=500,  cont Y,    disc X (5 lv)",
         seed = 722, n = 500,  Y_type = "continuous",
         X_levels = 5),
    
    list(label = "n=2000, cont Y,    disc X (10 lv)",
         seed = 723, n = 2000, Y_type = "continuous",
         X_levels = 10),
    
    list(label = "n=500,  disc Y (5 lv),  cont X",
         seed = 731, n = 500,  Y_type = 5,
         X_levels = NA),
    
    list(label = "n=500,  disc Y (5 lv),  disc X (10 lv)",
         seed = 732, n = 500,  Y_type = 5,
         X_levels = 10),
    
    list(label = "n=2000, disc Y (5 lv),  disc X (10 lv)",
         seed = 733, n = 2000, Y_type = 5,
         X_levels = 10)
  )
  
  # ---- Header --------------------------------------------------------------
  w <- 90
  cat("\n")
  cat(strrep("=", w), "\n")
  cat("  acor.test() end-to-end benchmarks vs reference implementations\n")
  cat(sprintf("  Median of %d replications — times in seconds\n", reps))
  cat("  AKC version = kernel version selected by select_kernel_version()\n")
  cat("  Reference: DeLong (pROC) for binary Y | survival::concordance() for continuous/discrete Y\n")
  cat(strrep("=", w), "\n")
  cat(sprintf("  %-36s  %-10s  %8s  %8s  %8s  %8s\n",
              "Scenario", "AKC ver", "AKC", "AGC", "reference", "ref method"))
  cat(strrep("-", w), "\n")
  
  for (sc in scenarios) {
    set.seed(sc$seed)
    n <- sc$n
    
    X <- if (is.na(sc$X_levels)) rnorm(n) else sample(seq_len(sc$X_levels), n, replace = TRUE)
    
    if (sc$Y_type == "binary") {
      Y          <- rbinom(n, 1, 0.5)
      ref_method <- "DeLong"
      can_ref    <- has_pROC
    } else if (sc$Y_type == "continuous") {
      Y          <- rnorm(n)
      ref_method <- "concordance"
      can_ref    <- TRUE
    } else {
      Y          <- sample(seq_len(sc$Y_type), n, replace = TRUE)
      ref_method <- "concordance"
      can_ref    <- TRUE
    }
    
    akc_ver <- acor:::select_kernel_version(Y, X)
    
    t_akc <- time_med(function() acor.test(X, Y, method = "akc"), reps)
    t_agc <- time_med(function() acor.test(X, Y, method = "agc"), reps)
    
    if (!can_ref) {
      t_ref <- NA_real_
    } else if (ref_method == "DeLong") {
      t_ref <- time_med(function() {
        roc_obj <- pROC::roc(Y, X, direction = "<", quiet = TRUE)
        pROC::var(roc_obj)
      }, reps)
    } else {
      t_ref <- time_med(function() survival::concordance(Y ~ X), reps)
    }
    
    ref_str <- if (is.na(t_ref)) "     N/A" else sprintf("%8.4f", t_ref)
    
    cat(sprintf("  %-36s  %-10s  %8.4f  %8.4f  %s  %s\n",
                sc$label, akc_ver,
                t_akc, t_agc,
                ref_str, ref_method))
  }
  
  cat(strrep("-", w), "\n")
  if (!has_pROC) {
    cat("  NOTE: pROC not installed — DeLong reference rows show N/A\n")
  }
  cat("\n")
  
  invisible(NULL)
}

# Run both benchmark sets when script is sourced directly
#if (interactive() && identical(Sys.getenv("RUN_BENCHMARKS"), "1")) {
#run_acor_benchmarks()
#}
