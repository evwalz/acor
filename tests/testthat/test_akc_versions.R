# ============================================================================
# Tests: Version Consistency and Runtime Benchmarking
# ============================================================================
#
# Purpose:
#   1. Verify that "original", "v1", and "v2" produce numerically identical
#      estimates and variances across all data scenarios.
#   2. Print runtimes to validate that the version selection heuristic makes
#      sense:
#        - v1  should be fastest for highly discrete / many-ties data
#        - v2  should be fastest for large n continuous data
#        - original is the O(n^2) baseline, expected fastest only for small n
#
# The internal dispatcher functions are called directly so that the version
# can be forced, bypassing select_kernel_version().
#
# ============================================================================

library(testthat)
library(acor)  # Replace with your actual package name

# Tolerance for numerical equivalence across versions
TOL <- 1e-10

# Shorthand wrappers to force a specific version
run_uni <- function(X, Y, version, IID = TRUE) {
  acor:::compute_akc_variance_auto(X, Y, IID = IID, version = version)
}

run_mv <- function(X, Y, version, IID = TRUE) {
  acor:::compute_akc_multivariate_variance_auto(X, Y, IID = IID, version = version)
}


# ============================================================================
# Section 1: Univariate IID
# ============================================================================

test_that("Uni IID | continuous X & Y: all versions agree", {

  set.seed(101); n <- 300
  X <- rnorm(n); Y <- rnorm(n)
  
  r_orig <- run_uni(X, Y, "original")
  r_v1   <- run_uni(X, Y, "v1")
  r_v2   <- run_uni(X, Y, "v2")
  
  expect_equal(r_v1$akc,     r_orig$akc,     tolerance = TOL)
  expect_equal(r_v2$akc,     r_orig$akc,     tolerance = TOL)
  expect_equal(r_v1$var,     r_orig$var,     tolerance = TOL)
  expect_equal(r_v2$var,     r_orig$var,     tolerance = TOL)
  expect_equal(r_v1$var_ind, r_orig$var_ind, tolerance = TOL)
  expect_equal(r_v2$var_ind, r_orig$var_ind, tolerance = TOL)
})

test_that("Uni IID | binary Y: all versions agree", {

  set.seed(102); n <- 300
  X <- rnorm(n); Y <- rbinom(n, 1, 0.5)
  
  r_orig <- run_uni(X, Y, "original")
  r_v1   <- run_uni(X, Y, "v1")
  r_v2   <- run_uni(X, Y, "v2")
  
  expect_equal(r_v1$akc,     r_orig$akc,     tolerance = TOL)
  expect_equal(r_v2$akc,     r_orig$akc,     tolerance = TOL)
  expect_equal(r_v1$var,     r_orig$var,     tolerance = TOL)
  expect_equal(r_v2$var,     r_orig$var,     tolerance = TOL)
  expect_equal(r_v1$var_ind, r_orig$var_ind, tolerance = TOL)
  expect_equal(r_v2$var_ind, r_orig$var_ind, tolerance = TOL)
})

test_that("Uni IID | discrete X (10 levels), continuous Y: all versions agree", {

  set.seed(103); n <- 300
  X <- sample(1:10, n, replace = TRUE); Y <- rnorm(n)
  
  r_orig <- run_uni(X, Y, "original")
  r_v1   <- run_uni(X, Y, "v1")
  r_v2   <- run_uni(X, Y, "v2")
  
  expect_equal(r_v1$akc,     r_orig$akc,     tolerance = TOL)
  expect_equal(r_v2$akc,     r_orig$akc,     tolerance = TOL)
  expect_equal(r_v1$var,     r_orig$var,     tolerance = TOL)
  expect_equal(r_v2$var,     r_orig$var,     tolerance = TOL)
  expect_equal(r_v1$var_ind, r_orig$var_ind, tolerance = TOL)
  expect_equal(r_v2$var_ind, r_orig$var_ind, tolerance = TOL)
})

test_that("Uni IID | discrete X (5 levels) & discrete Y (3 levels): all versions agree", {

  set.seed(104); n <- 300
  X <- sample(1:5, n, replace = TRUE); Y <- sample(1:3, n, replace = TRUE)
  
  r_orig <- run_uni(X, Y, "original")
  r_v1   <- run_uni(X, Y, "v1")
  r_v2   <- run_uni(X, Y, "v2")
  
  expect_equal(r_v1$akc,     r_orig$akc,     tolerance = TOL)
  expect_equal(r_v2$akc,     r_orig$akc,     tolerance = TOL)
  expect_equal(r_v1$var,     r_orig$var,     tolerance = TOL)
  expect_equal(r_v2$var,     r_orig$var,     tolerance = TOL)
  expect_equal(r_v1$var_ind, r_orig$var_ind, tolerance = TOL)
  expect_equal(r_v2$var_ind, r_orig$var_ind, tolerance = TOL)
})

test_that("Uni IID | heavy ties in X (~25%): all versions agree", {

  set.seed(105); n <- 300
  X <- rnorm(n); X[seq(1, n, by = 4)] <- X[1]; Y <- rnorm(n)
  
  r_orig <- run_uni(X, Y, "original")
  r_v1   <- run_uni(X, Y, "v1")
  r_v2   <- run_uni(X, Y, "v2")
  
  expect_equal(r_v1$akc,     r_orig$akc,     tolerance = TOL)
  expect_equal(r_v2$akc,     r_orig$akc,     tolerance = TOL)
  expect_equal(r_v1$var,     r_orig$var,     tolerance = TOL)
  expect_equal(r_v2$var,     r_orig$var,     tolerance = TOL)
  expect_equal(r_v1$var_ind, r_orig$var_ind, tolerance = TOL)
  expect_equal(r_v2$var_ind, r_orig$var_ind, tolerance = TOL)
})


# ============================================================================
# Section 2: Univariate HAC
# ============================================================================

test_that("Uni HAC | continuous X & Y: all versions agree", {

  set.seed(201); n <- 300
  X <- rnorm(n); Y <- rnorm(n)
  
  r_orig <- run_uni(X, Y, "original", IID = FALSE)
  r_v1   <- run_uni(X, Y, "v1",       IID = FALSE)
  r_v2   <- run_uni(X, Y, "v2",       IID = FALSE)
  
  expect_equal(r_v1$akc,     r_orig$akc,     tolerance = TOL)
  expect_equal(r_v2$akc,     r_orig$akc,     tolerance = TOL)
  expect_equal(r_v1$var,     r_orig$var,     tolerance = TOL)
  expect_equal(r_v2$var,     r_orig$var,     tolerance = TOL)
  expect_equal(r_v1$var_ind, r_orig$var_ind, tolerance = TOL)
  expect_equal(r_v2$var_ind, r_orig$var_ind, tolerance = TOL)
})

test_that("Uni HAC | binary Y: all versions agree", {
  set.seed(202); n <- 300
  X <- rnorm(n); Y <- rbinom(n, 1, 0.6)
  
  r_orig <- run_uni(X, Y, "original", IID = FALSE)
  r_v1   <- run_uni(X, Y, "v1",       IID = FALSE)
  r_v2   <- run_uni(X, Y, "v2",       IID = FALSE)
  
  expect_equal(r_v1$akc,     r_orig$akc,     tolerance = TOL)
  expect_equal(r_v2$akc,     r_orig$akc,     tolerance = TOL)
  expect_equal(r_v1$var,     r_orig$var,     tolerance = TOL)
  expect_equal(r_v2$var,     r_orig$var,     tolerance = TOL)
  expect_equal(r_v1$var_ind, r_orig$var_ind, tolerance = TOL)
  expect_equal(r_v2$var_ind, r_orig$var_ind, tolerance = TOL)
})

test_that("Uni HAC | discrete X (10) & Y (5): all versions agree", {
  skip_on_cran()
  set.seed(203); n <- 300
  X <- sample(1:10, n, replace = TRUE); Y <- sample(1:5, n, replace = TRUE)
  
  r_orig <- run_uni(X, Y, "original", IID = FALSE)
  r_v1   <- run_uni(X, Y, "v1",       IID = FALSE)
  r_v2   <- run_uni(X, Y, "v2",       IID = FALSE)
  
  expect_equal(r_v1$akc,     r_orig$akc,     tolerance = TOL)
  expect_equal(r_v2$akc,     r_orig$akc,     tolerance = TOL)
  expect_equal(r_v1$var,     r_orig$var,     tolerance = TOL)
  expect_equal(r_v2$var,     r_orig$var,     tolerance = TOL)
  expect_equal(r_v1$var_ind, r_orig$var_ind, tolerance = TOL)
  expect_equal(r_v2$var_ind, r_orig$var_ind, tolerance = TOL)
})


# ============================================================================
# Section 3: Multivariate IID (m = 3)
# ============================================================================

test_that("MV IID m=3 | continuous X & Y: all versions agree", {
  set.seed(301); n <- 300
  X <- matrix(rnorm(n * 3), ncol = 3); Y <- rnorm(n)
  
  r_orig <- run_mv(X, Y, "original")
  r_v1   <- run_mv(X, Y, "v1")
  r_v2   <- run_mv(X, Y, "v2")
  
  expect_equal(r_v1$akc_vector, r_orig$akc_vector, tolerance = TOL)
  expect_equal(r_v2$akc_vector, r_orig$akc_vector, tolerance = TOL)
  expect_equal(r_v1$Sigma,      r_orig$Sigma,      tolerance = TOL)
  expect_equal(r_v2$Sigma,      r_orig$Sigma,      tolerance = TOL)
  expect_equal(r_v1$Sigma_ind,  r_orig$Sigma_ind,  tolerance = TOL)
  expect_equal(r_v2$Sigma_ind,  r_orig$Sigma_ind,  tolerance = TOL)
})

test_that("MV IID m=3 | binary Y: all versions agree", {
  set.seed(302); n <- 300
  X <- matrix(rnorm(n * 3), ncol = 3); Y <- rbinom(n, 1, 0.5)
  
  r_orig <- run_mv(X, Y, "original")
  r_v1   <- run_mv(X, Y, "v1")
  r_v2   <- run_mv(X, Y, "v2")
  
  expect_equal(r_v1$akc_vector, r_orig$akc_vector, tolerance = TOL)
  expect_equal(r_v2$akc_vector, r_orig$akc_vector, tolerance = TOL)
  expect_equal(r_v1$Sigma,      r_orig$Sigma,      tolerance = TOL)
  expect_equal(r_v2$Sigma,      r_orig$Sigma,      tolerance = TOL)
  expect_equal(r_v1$Sigma_ind,  r_orig$Sigma_ind,  tolerance = TOL)
  expect_equal(r_v2$Sigma_ind,  r_orig$Sigma_ind,  tolerance = TOL)
})

test_that("MV IID m=3 | discrete X (10) & Y (5): all versions agree", {
  set.seed(303); n <- 300
  X <- matrix(sample(1:10, n * 3, replace = TRUE), ncol = 3)
  Y <- sample(1:5, n, replace = TRUE)
  
  r_orig <- run_mv(X, Y, "original")
  r_v1   <- run_mv(X, Y, "v1")
  r_v2   <- run_mv(X, Y, "v2")
  
  expect_equal(r_v1$akc_vector, r_orig$akc_vector, tolerance = TOL)
  expect_equal(r_v2$akc_vector, r_orig$akc_vector, tolerance = TOL)
  expect_equal(r_v1$Sigma,      r_orig$Sigma,      tolerance = TOL)
  expect_equal(r_v2$Sigma,      r_orig$Sigma,      tolerance = TOL)
  expect_equal(r_v1$Sigma_ind,  r_orig$Sigma_ind,  tolerance = TOL)
  expect_equal(r_v2$Sigma_ind,  r_orig$Sigma_ind,  tolerance = TOL)
})


# ============================================================================
# Section 4: Multivariate HAC (m = 3)
# ============================================================================

test_that("MV HAC m=3 | continuous X & Y: all versions agree", {
  set.seed(401); n <- 300
  X <- matrix(rnorm(n * 3), ncol = 3); Y <- rnorm(n)
  
  r_orig <- run_mv(X, Y, "original", IID = FALSE)
  r_v1   <- run_mv(X, Y, "v1",       IID = FALSE)
  r_v2   <- run_mv(X, Y, "v2",       IID = FALSE)
  
  expect_equal(r_v1$akc_vector, r_orig$akc_vector, tolerance = TOL)
  expect_equal(r_v2$akc_vector, r_orig$akc_vector, tolerance = TOL)
  expect_equal(r_v1$Sigma,      r_orig$Sigma,      tolerance = TOL)
  expect_equal(r_v2$Sigma,      r_orig$Sigma,      tolerance = TOL)
  expect_equal(r_v1$Sigma_ind,  r_orig$Sigma_ind,  tolerance = TOL)
  expect_equal(r_v2$Sigma_ind,  r_orig$Sigma_ind,  tolerance = TOL)
})

test_that("MV HAC m=3 | binary Y: all versions agree", {
  skip_on_cran()
  set.seed(402); n <- 300
  X <- matrix(rnorm(n * 3), ncol = 3); Y <- rbinom(n, 1, 0.6)
  
  r_orig <- run_mv(X, Y, "original", IID = FALSE)
  r_v1   <- run_mv(X, Y, "v1",       IID = FALSE)
  r_v2   <- run_mv(X, Y, "v2",       IID = FALSE)
  
  expect_equal(r_v1$akc_vector, r_orig$akc_vector, tolerance = TOL)
  expect_equal(r_v2$akc_vector, r_orig$akc_vector, tolerance = TOL)
  expect_equal(r_v1$Sigma,      r_orig$Sigma,      tolerance = TOL)
  expect_equal(r_v2$Sigma,      r_orig$Sigma,      tolerance = TOL)
  expect_equal(r_v1$Sigma_ind,  r_orig$Sigma_ind,  tolerance = TOL)
  expect_equal(r_v2$Sigma_ind,  r_orig$Sigma_ind,  tolerance = TOL)
})

test_that("MV HAC m=3 | discrete X (10) & Y (5): all versions agree", {
  skip_on_cran()
  set.seed(403); n <- 300
  X <- matrix(sample(1:10, n * 3, replace = TRUE), ncol = 3)
  Y <- sample(1:5, n, replace = TRUE)
  
  r_orig <- run_mv(X, Y, "original", IID = FALSE)
  r_v1   <- run_mv(X, Y, "v1",       IID = FALSE)
  r_v2   <- run_mv(X, Y, "v2",       IID = FALSE)
  
  expect_equal(r_v1$akc_vector, r_orig$akc_vector, tolerance = TOL)
  expect_equal(r_v2$akc_vector, r_orig$akc_vector, tolerance = TOL)
  expect_equal(r_v1$Sigma,      r_orig$Sigma,      tolerance = TOL)
  expect_equal(r_v2$Sigma,      r_orig$Sigma,      tolerance = TOL)
  expect_equal(r_v1$Sigma_ind,  r_orig$Sigma_ind,  tolerance = TOL)
  expect_equal(r_v2$Sigma_ind,  r_orig$Sigma_ind,  tolerance = TOL)
})


# ============================================================================
# Section 5: Runtime benchmarks
# ============================================================================
# Not testthat assertions — prints a formatted table of median runtimes.
# Run interactively or via Rscript to validate version selection heuristic.
#
# Expected pattern:
#   - "original" fastest at small n continuous (but O(n^2) kills it for large n)
#   - "v1" fastest whenever data is discrete or binary Y
#   - "v2" fastest at large n continuous (n > ~4000)
# ============================================================================

run_benchmarks <- function(reps = 5) {
  
  time_med <- function(expr_fn, reps) {
    median(replicate(reps, system.time(expr_fn())["elapsed"]))
  }
  
  # ---- Scenario definitions ------------------------------------------------
  scenarios <- list(
    # --- Continuous data: watches original degrade vs v2 improving with n ---
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
    
    #list(label = "n=10000, continuous",     seed = 505, n = 10000,
    #     gen_X = function(n) rnorm(n),
    #     gen_Y = function(n) rnorm(n)),
    
    # --- Binary Y: v1 should dominate across all n -------------------------
    list(label = "n=200,   binary Y",       seed = 506, n = 200,
         gen_X = function(n) rnorm(n),
         gen_Y = function(n) rbinom(n, 1, 0.5)),
    
    list(label = "n=500,   binary Y",       seed = 507, n = 500,
         gen_X = function(n) rnorm(n),
         gen_Y = function(n) rbinom(n, 1, 0.5)),
    
    list(label = "n=5000,  binary Y",       seed = 508, n = 5000,
         gen_X = function(n) rnorm(n),
         gen_Y = function(n) rbinom(n, 1, 0.5)),
    
    # --- Discrete X: v1 advantage from tie structure -----------------------
    list(label = "n=500,   discrete X (10 lv)", seed = 509, n = 500,
         gen_X = function(n) sample(1:10, n, replace = TRUE),
         gen_Y = function(n) rnorm(n)),
    
    list(label = "n=500,   discrete X+Y (10x5)", seed = 510, n = 500,
         gen_X = function(n) sample(1:10, n, replace = TRUE),
         gen_Y = function(n) sample(1:5, n, replace = TRUE)),
    
    list(label = "n=5000,  discrete X+Y (10x5)", seed = 511, n = 5000,
         gen_X = function(n) sample(1:10, n, replace = TRUE),
         gen_Y = function(n) sample(1:5, n, replace = TRUE)),
    
    list(label = "n=10000,  discrete X+Y (10x5)", seed = 511, n = 10000,
         gen_X = function(n) sample(1:10, n, replace = TRUE),
         gen_Y = function(n) sample(1:5, n, replace = TRUE))
  )
  
  # ---- Print header --------------------------------------------------------
  cat("\n")
  cat(strrep("=", 78), "\n")
  cat("  Runtime Benchmarks: original vs v1 vs v2  (IID, univariate)\n")
  cat(sprintf("  Median of %d replications — times in seconds\n", reps))
  cat(strrep("=", 78), "\n")
  cat(sprintf("  %-32s  %-10s  %8s  %8s  %8s\n",
              "Scenario", "selected", "original", "v1", "v2"))
  cat(strrep("-", 78), "\n")
  
  for (sc in scenarios) {
    set.seed(sc$seed)
    X <- sc$gen_X(sc$n)
    Y <- sc$gen_Y(sc$n)
    
    selected <- acor:::select_kernel_version(Y, X)
    
    t_orig <- time_med(function() run_uni(X, Y, "original"), reps)
    t_v1   <- time_med(function() run_uni(X, Y, "v1"),       reps)
    t_v2   <- time_med(function() run_uni(X, Y, "v2"),       reps)
    
    # Mark the selected version with an asterisk
    orig_str <- if (selected == "original") sprintf("%7.4f*", t_orig) else sprintf("%8.4f",  t_orig)
    v1_str   <- if (selected == "v1")       sprintf("%7.4f*", t_v1)   else sprintf("%8.4f",  t_v1)
    v2_str   <- if (selected == "v2")       sprintf("%7.4f*", t_v2)   else sprintf("%8.4f",  t_v2)
    
    cat(sprintf("  %-32s  %-10s  %s  %s  %s\n",
                sc$label, selected, orig_str, v1_str, v2_str))
  }
  
  cat(strrep("-", 78), "\n")
  cat("  * = version that would be selected by select_kernel_version()\n\n")
  
  
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
  
  cat(strrep("=", 78), "\n")
  cat("  Runtime Benchmarks: original vs v1 vs v2  (IID, multivariate)\n")
  cat(sprintf("  Median of %d replications — times in seconds\n", reps))
  cat(strrep("=", 78), "\n")
  cat(sprintf("  %-32s  %-10s  %8s  %8s  %8s\n",
              "Scenario", "selected", "original", "v1", "v2"))
  cat(strrep("-", 78), "\n")
  
  for (sc in mv_scenarios) {
    set.seed(sc$seed)
    X <- sc$gen_X(sc$n, sc$m)
    Y <- sc$gen_Y(sc$n)
    
    selected <- acor:::select_kernel_version(Y, X)
    
    t_orig <- time_med(function() run_mv(X, Y, "original"), reps)
    t_v1   <- time_med(function() run_mv(X, Y, "v1"),       reps)
    t_v2   <- time_med(function() run_mv(X, Y, "v2"),       reps)
    
    orig_str <- if (selected == "original") sprintf("%7.4f*", t_orig) else sprintf("%8.4f",  t_orig)
    v1_str   <- if (selected == "v1")       sprintf("%7.4f*", t_v1)   else sprintf("%8.4f",  t_v1)
    v2_str   <- if (selected == "v2")       sprintf("%7.4f*", t_v2)   else sprintf("%8.4f",  t_v2)
    
    cat(sprintf("  %-32s  %-10s  %s  %s  %s\n",
                sc$label, selected, orig_str, v1_str, v2_str))
  }
  
  cat(strrep("-", 78), "\n")
  cat("  * = version that would be selected by select_kernel_version()\n\n")
  
  invisible(NULL)
}

# Run benchmarks when script is sourced directly (not during testthat)
if (interactive() && identical(Sys.getenv("RUN_BENCHMARKS"), "1")) {
  run_benchmarks()
}


# ============================================================================
# Section 6: acor.test() end-to-end benchmarks vs reference implementations
# ============================================================================
# Compares wall-clock time of acor.test() (AKC and AGC) against:
#   - pROC::roc() + pROC::var()     for binary Y   (DeLong)
#   - cor.test(method="spearman")   for continuous Y
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
  # Discretisation levels for X (NA = continuous rnorm)
  # Y type: "binary", "continuous", or integer = number of discrete levels
  scenarios <- list(
    
    # --- Binary Y: acor.test vs DeLong -------------------------------------
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
    
    # --- Continuous Y: acor.test vs cor.test(spearman) --------------------
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
    
    # --- Discrete X: acor.test vs cor.test(spearman) ----------------------
    list(label = "n=500,  cont Y,    disc X (10 lv)",
         seed = 721, n = 500,  Y_type = "continuous",
         X_levels = 10),
    
    list(label = "n=500,  cont Y,    disc X (5 lv)",
         seed = 722, n = 500,  Y_type = "continuous",
         X_levels = 5),
    
    list(label = "n=2000, cont Y,    disc X (10 lv)",
         seed = 723, n = 2000, Y_type = "continuous",
         X_levels = 10),
    
    # --- Discrete Y: acor.test (no direct reference) ----------------------
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
    
    # Generate X
    X <- if (is.na(sc$X_levels)) rnorm(n) else sample(seq_len(sc$X_levels), n, replace = TRUE)
    
    # Generate Y and determine reference method
    if (sc$Y_type == "binary") {
      Y          <- rbinom(n, 1, 0.5)
      ref_method <- "DeLong"
      can_ref    <- has_pROC
    } else if (sc$Y_type == "continuous") {
      Y          <- rnorm(n)
      ref_method <- "concordance"
      can_ref    <- TRUE
    } else {
      # discrete Y
      Y          <- sample(seq_len(sc$Y_type), n, replace = TRUE)
      ref_method <- "concordance"
      can_ref    <- TRUE
    }
    
    # AKC version that will be selected
    akc_ver <- acor:::select_kernel_version(Y, X)
    
    # Time acor.test
    t_akc <- time_med(function() acor.test(X, Y, method = "akc"), reps)
    t_agc <- time_med(function() acor.test(X, Y, method = "agc"), reps)
    
    # Time reference
    if (!can_ref) {
      t_ref <- NA_real_
    } else if (ref_method == "DeLong") {
      t_ref <- time_med(function() {
        roc_obj <- pROC::roc(Y, X, direction = "<", quiet = TRUE)
        pROC::var(roc_obj)
      }, reps)
    } else {
      # survival::concordance handles continuous, discrete, and binary Y
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
if (interactive() && identical(Sys.getenv("RUN_BENCHMARKS"), "1")) {
  run_acor_benchmarks()
}