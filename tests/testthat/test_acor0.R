# Tests for correlation measure equivalences

library(testthat)
library(acor) 

# Helper function to compute Harrell's C-index
compute_harrell_c <- function(X, Y) {
  n <- length(Y)
  concordant <- 0
  discordant <- 0
  tied_x <- 0
  tied_y <- 0
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (Y[i] != Y[j]) {
        # Y values differ - this is a usable pair
        if (X[i] == X[j]) {
          # Tie in predictions
          tied_x <- tied_x + 1
        } else if ((Y[i] > Y[j] && X[i] > X[j]) || (Y[i] < Y[j] && X[i] < X[j])) {
          concordant <- concordant + 1
        } else {
          discordant <- discordant + 1
        }
      } else {
        tied_y <- tied_y + 1
      }
    }
  }
  
  # Standard formula: add 0.5 for each tie in X
  c_index <- (concordant + 0.5 * tied_x) / (concordant + discordant + tied_x)
  
  return(c_index)
}

# Helper function to compute AUC (equivalent to C-index for binary Y)
compute_auc <- function(X, Y) {
  # For binary Y, AUC is the C-index
  compute_harrell_c(X, Y)
}

# ============================================================================
# Test 1: For binary Y, check that CMA = CID = AUC
# ============================================================================

test_that("Binary Y: CMA equals CID", {
  set.seed(100)
  n <- 100
  X <- rnorm(n)
  Y <- rbinom(n, 1, 0.6)
  
  # Get CMA
  res_cma <- acor(X, Y, method = "cma")
  cma <- unname(res_cma$estimate)
  
  # Get CID
  res_cid <- acor(X, Y, method = "cid")
  cid <- unname(res_cid$estimate)
  
  # They should be equal for binary Y
  expect_equal(cma, cid, tolerance = 1e-10,
               info = "CMA should equal CID for binary Y")
})

test_that("Binary Y: CMA equals AUC", {
  set.seed(101)
  n <- 100
  X <- rnorm(n)
  Y <- rbinom(n, 1, 0.6)
  
  # Get CMA
  res_cma <- acor(X, Y, method = "cma")
  cma <- unname(res_cma$estimate)
  
  # Compute AUC
  auc <- compute_auc(X, Y)
  
  # They should be equal for binary Y
  expect_equal(cma, auc, tolerance = 1e-10,
               info = "CMA should equal AUC for binary Y")
})

test_that("Binary Y: CID equals AUC", {
  set.seed(102)
  n <- 100
  X <- rnorm(n)
  Y <- rbinom(n, 1, 0.7)
  
  # Get CID
  res_cid <- acor(X, Y, method = "cid")
  cid <- unname(res_cid$estimate)
  
  # Compute AUC
  auc <- compute_auc(X, Y)
  
  # They should be equal for binary Y
  expect_equal(cid, auc, tolerance = 1e-10,
               info = "CID should equal AUC for binary Y")
})

test_that("Binary Y with ties in X: CMA = CID = AUC", {
  set.seed(103)
  n <- 100
  X <- rnorm(n)
  X[seq(1, n, by = 5)] <- X[1]  # Create ties
  Y <- rbinom(n, 1, 0.5)
  
  res_cma <- acor(X, Y, method = "cma")
  res_cid <- acor(X, Y, method = "cid")
  auc <- compute_auc(X, Y)
  
  expect_equal(unname(res_cma$estimate), unname(res_cid$estimate), tolerance = 1e-10)
  expect_equal(unname(res_cma$estimate), auc, tolerance = 1e-10)
  expect_equal(unname(res_cid$estimate), auc, tolerance = 1e-10)
})

# ============================================================================
# Test 2: For binary Y, check that variance estimators are equal for AKC and AGC
# ============================================================================

test_that("Binary Y: AKC and AGC have equal variances", {
  set.seed(200)
  n <- 100
  X <- rnorm(n)
  Y <- rbinom(n, 1, 0.6)
  
  # Get AKC variance
  res_akc <- acor(X, Y, method = "akc")
  var_akc <- unname(res_akc$variance)
  
  # Get AGC variance
  res_agc <- acor(X, Y, method = "agc")
  var_agc <- unname(res_agc$variance)
  
  # They should be equal for binary Y
  tolerance = 0.03
  expect_equal(var_akc, var_agc, tolerance = tolerance,
               info = "AKC and AGC variances should be equal for binary Y")
})

test_that("Binary Y with ties: AKC and AGC have equal variances", {
  set.seed(201)
  n <- 100
  X <- rnorm(n)
  X[seq(1, n, by = 7)] <- X[2]  # Create ties
  Y <- rbinom(n, 1, 0.5)
  
  res_akc <- acor(X, Y, method = "akc")
  res_agc <- acor(X, Y, method = "agc")

  tolerance = 0.03
  
  expect_equal(unname(res_akc$variance), unname(res_agc$variance), tolerance = tolerance,
               info = "AKC and AGC variances should be equal for binary Y even with ties in X")
})

test_that("Binary Y multiple predictors: AKC and AGC have equal covariance matrices", {
  set.seed(202)
  n <- 100
  X <- matrix(rnorm(n * 3), ncol = 3)
  Y <- rbinom(n, 1, 0.6)
  
  res_akc <- acor(X, Y, method = "akc")
  res_agc <- acor(X, Y, method = "agc")
  
  tolerance = 0.03
  # Covariance matrices should be equal
  expect_equal(unname(res_akc$variance), unname(res_agc$variance), tolerance = tolerance,
               info = "AKC and AGC covariance matrices should be equal for binary Y")
})

# ============================================================================
# Test 3: No ties in X and Y, check AKC = Kendall's tau and AGC = Spearman's rho
# ============================================================================

test_that("No ties: AKC equals Kendall's tau", {
  set.seed(300)
  n <- 100
  # Add small noise to ensure no ties
  X <- rnorm(n) + rnorm(n, sd = 0.0001)
  Y <- 0.5 * X + rnorm(n, sd = 0.8) + rnorm(n, sd = 0.0001)
  
  # Verify no ties
  expect_equal(length(unique(X)), n, info = "X should have no ties")
  expect_equal(length(unique(Y)), n, info = "Y should have no ties")
  
  # Get AKC
  res_akc <- acor(X, Y, method = "akc")
  akc <- unname(res_akc$estimate)
  
  # Get Kendall's tau
  kendall <- cor(X, Y, method = "kendall")
  
  # They should be equal when there are no ties
  expect_equal(akc, kendall, tolerance = 1e-10,
               info = "AKC should equal Kendall's tau when there are no ties")
})

test_that("No ties: AGC equals Spearman's rho", {
  set.seed(301)
  n <- 100
  # Add small noise to ensure no ties
  X <- rnorm(n) + rnorm(n, sd = 0.0001)
  Y <- 0.5 * X + rnorm(n, sd = 0.8) + rnorm(n, sd = 0.0001)
  
  # Verify no ties
  expect_equal(length(unique(X)), n, info = "X should have no ties")
  expect_equal(length(unique(Y)), n, info = "Y should have no ties")
  
  # Get AGC
  res_agc <- acor(X, Y, method = "agc")
  agc <- unname(res_agc$estimate)
  
  # Get Spearman's rho
  spearman <- cor(X, Y, method = "spearman")
  
  # They should be equal when there are no ties
  expect_equal(agc, spearman, tolerance = 1e-10,
               info = "AGC should equal Spearman's rho when there are no ties")
})

test_that("No ties multiple tests: AKC = Kendall and AGC = Spearman", {
  # Test with different sample sizes and correlations
  test_cases <- list(
    list(n = 50, rho = 0.3, seed = 302),
    list(n = 100, rho = 0.7, seed = 303),
    list(n = 150, rho = -0.5, seed = 304)
  )
  
  for (tc in test_cases) {
    set.seed(tc$seed)
    X <- rnorm(tc$n) + rnorm(tc$n, sd = 0.0001)
    Y <- tc$rho * X + sqrt(1 - tc$rho^2) * rnorm(tc$n) + rnorm(tc$n, sd = 0.0001)
    
    # Verify no ties
    expect_equal(length(unique(X)), tc$n)
    expect_equal(length(unique(Y)), tc$n)
    
    # Test AKC = Kendall
    akc <- unname(acor(X, Y, method = "akc")$estimate)
    kendall <- cor(X, Y, method = "kendall")
    expect_equal(akc, kendall, tolerance = 1e-10,
                 info = sprintf("n=%d, rho=%.1f: AKC = Kendall", tc$n, tc$rho))
    
    # Test AGC = Spearman
    agc <- unname(acor(X, Y, method = "agc")$estimate)
    spearman <- cor(X, Y, method = "spearman")
    expect_equal(agc, spearman, tolerance = 1e-10,
                 info = sprintf("n=%d, rho=%.1f: AGC = Spearman", tc$n, tc$rho))
  }
})

# ============================================================================
# Test 4: Binary Y, check CID = Harrell's C-index (both with and without ties)
# ============================================================================

test_that("Binary Y, no ties: CID equals Harrell's C-index", {
  set.seed(400)
  n <- 100
  X <- rnorm(n) + rnorm(n, sd = 0.0001)
  Y <- rbinom(n, 1, 0.6)
  
  # Verify no ties in X
  expect_equal(length(unique(X)), n, info = "X should have no ties")
  
  # Get CID
  res_cid <- acor(X, Y, method = "cid")
  cid <- unname(res_cid$estimate)
  
  # Compute Harrell's C-index
  c_index <- compute_harrell_c(X, Y)
  
  # They should be equal
  expect_equal(cid, c_index, tolerance = 1e-10,
               info = "CID should equal Harrell's C-index for binary Y with no ties")
})

test_that("Binary Y, with ties in X: CID equals Harrell's C-index", {
  set.seed(401)
  n <- 100
  X <- rnorm(n)
  X[seq(1, n, by = 8)] <- X[1]  # Create ties in X
  Y <- rbinom(n, 1, 0.6)
  
  # Verify ties exist
  expect_lt(length(unique(X)), n)
  
  # Get CID
  res_cid <- acor(X, Y, method = "cid")
  cid <- unname(res_cid$estimate)
  
  # Compute Harrell's C-index
  c_index <- compute_harrell_c(X, Y)
  
  # They should be equal even with ties
  expect_equal(cid, c_index, tolerance = 1e-10,
               info = "CID should equal Harrell's C-index for binary Y even with ties in X")
})

test_that("Binary Y, with ties in Y (shouldn't happen but test edge case)", {
  set.seed(402)
  n <- 100
  X <- rnorm(n)
  Y <- rbinom(n, 1, 0.6)
  
  # Y is binary so has "ties" by definition
  expect_equal(length(unique(Y)), 2, info = "Binary Y has exactly 2 unique values")
  
  # Get CID
  res_cid <- acor(X, Y, method = "cid")
  cid <- unname(res_cid$estimate)
  
  # Compute Harrell's C-index
  c_index <- compute_harrell_c(X, Y)
  
  # They should be equal
  expect_equal(cid, c_index, tolerance = 1e-10,
               info = "CID should equal Harrell's C-index for binary Y")
})

test_that("Binary Y, multiple scenarios: CID = Harrell's C", {
  test_cases <- list(
    list(name = "No ties", create_ties = FALSE, seed = 403),
    list(name = "Few ties", create_ties = TRUE, tie_freq = 10, seed = 404),
    list(name = "Many ties", create_ties = TRUE, tie_freq = 5, seed = 405),
    list(name = "Extreme ties", create_ties = TRUE, tie_freq = 3, seed = 406)
  )
  
  for (tc in test_cases) {
    set.seed(tc$seed)
    n <- 100
    X <- rnorm(n)
    
    if (tc$create_ties) {
      X[seq(1, n, by = tc$tie_freq)] <- X[1]
    } else {
      X <- X + rnorm(n, sd = 0.0001)
    }
    
    Y <- rbinom(n, 1, 0.6)
    
    cid <- unname(acor(X, Y, method = "cid")$estimate)
    c_index <- compute_harrell_c(X, Y)
    
    expect_equal(cid, c_index, tolerance = 1e-10,
                 info = sprintf("%s: CID should equal Harrell's C-index", tc$name))
  }
})

# ============================================================================
# Additional cross-checks for consistency
# ============================================================================

test_that("Binary Y: All equivalent measures agree", {
  set.seed(500)
  n <- 100
  X <- rnorm(n)
  Y <- rbinom(n, 1, 0.6)
  
  # Get all measures
  cma <- unname(acor(X, Y, method = "cma")$estimate)
  cid <- unname(acor(X, Y, method = "cid")$estimate)
  auc <- compute_auc(X, Y)
  c_index <- compute_harrell_c(X, Y)
  
  # All should be equal
  expect_equal(cma, cid, tolerance = 1e-10)
  expect_equal(cma, auc, tolerance = 1e-10)
  expect_equal(cma, c_index, tolerance = 1e-10)
  expect_equal(cid, auc, tolerance = 1e-10)
  expect_equal(cid, c_index, tolerance = 1e-10)
  expect_equal(auc, c_index, tolerance = 1e-10)
})

test_that("Binary Y: AGC = 2*CMA - 1 and AKC = 2*CID - 1", {
  set.seed(501)
  n <- 100
  X <- rnorm(n)
  Y <- rbinom(n, 1, 0.6)
  
  cma <- unname(acor(X, Y, method = "cma")$estimate)
  cid <- unname(acor(X, Y, method = "cid")$estimate)
  agc <- unname(acor(X, Y, method = "agc")$estimate)
  akc <- unname(acor(X, Y, method = "akc")$estimate)
  
  # Check transformations
  expect_equal(agc, 2 * cma - 1, tolerance = 1e-10)
  expect_equal(akc, 2 * cid - 1, tolerance = 1e-10)
  
  # Since CMA = CID for binary Y, AGC should equal AKC
  expect_equal(agc, akc, tolerance = 1e-10)
})

test_that("Consistency across different binary Y probabilities", {
  probs <- c(0.3, 0.5, 0.7, 0.9)
  
  for (i in seq_along(probs)) {
    set.seed(500 + i)
    n <- 100
    X <- rnorm(n)
    Y <- rbinom(n, 1, probs[i])
    
    # Skip if Y is constant (all 0s or all 1s)
    if (length(unique(Y)) < 2) next
    
    cma <- unname(acor(X, Y, method = "cma")$estimate)
    cid <- unname(acor(X, Y, method = "cid")$estimate)
    auc <- compute_auc(X, Y)
    
    expect_equal(cma, cid, tolerance = 1e-10,
                 info = sprintf("p = %.1f: CMA = CID", probs[i]))
    expect_equal(cma, auc, tolerance = 1e-10,
                 info = sprintf("p = %.1f: CMA = AUC", probs[i]))
  }
})



# test_that("Binary Y: AKC and AGC variance difference decreases with n", {
#   # Test multiple sample sizes
#   sample_sizes <- c(100, 500)
#   rel_diffs <- numeric(length(sample_sizes))
  
#   for (i in seq_along(sample_sizes)) {
#     set.seed(200 + i)
#     n <- sample_sizes[i]
#     X <- rnorm(n)
#     Y <- rbinom(n, 1, 0.6)
    
#     res_akc <- acor(X, Y, method = "akc")
#     res_agc <- acor(X, Y, method = "agc")
    
#     var_akc <- unname(res_akc$variance)
#     var_agc <- unname(res_agc$variance)
    
#     # Compute relative difference
#     rel_diffs[i] <- abs(var_akc - var_agc) / var_agc
    
#     cat(sprintf("n = %5d: AKC var = %.6f, AGC var = %.6f, Rel diff = %.6f\n",
#                 n, var_akc, var_agc, rel_diffs[i]))
#   }
  
#   # Check that relative difference is decreasing
#   expect_true(all(diff(rel_diffs) < 0), 
#               info = "Relative difference should decrease as n increases")
  
#   # For largest n, difference should be small
#   expect_lt(rel_diffs[length(rel_diffs)], 0.01,
#             info = "For large n, relative difference should be < 1%")
# })