test_that("kendall_tau_a matches cor(method='kendall') when no ties", {
  set.seed(42)
  n <- 100
  X <- rnorm(n)
  Y <- rnorm(n)
  
  expected <- cor(X, Y, method = "kendall")
  result <- acor:::kendall_tau_a(X, Y)
  
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("kendall_tau_a is correct with ties", {
  set.seed(123)
  n <- 50
  X <- sample(1:10, n, replace = TRUE)
  Y <- sample(1:10, n, replace = TRUE)
  
  # Brute force tau-a
  pairs <- combn(n, 2)
  sgn <- sign(X[pairs[2,]] - X[pairs[1,]]) * sign(Y[pairs[2,]] - Y[pairs[1,]])
  expected <- mean(sgn)
  
  result <- acor:::kendall_tau_a(X, Y)
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("kendall_tau_b matches cor(method='kendall') with ties", {
  set.seed(99)
  n <- 80
  X <- sample(1:5, n, replace = TRUE)
  Y <- sample(1:5, n, replace = TRUE)
  
  expected <- cor(X, Y, method = "kendall")
  result <- acor:::kendall_tau_b(X, Y)
  
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("kendall_tau_b matches cor(method='kendall') without ties", {
  set.seed(7)
  n <- 100
  X <- rnorm(n)
  Y <- rnorm(n)
  
  expected <- cor(X, Y, method = "kendall")
  result <- acor:::kendall_tau_b(X, Y)
  
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("kendall_tau_a and kendall_tau_b agree when no ties", {
  set.seed(55)
  n <- 100
  X <- rnorm(n)
  Y <- rnorm(n)
  
  expect_equal(acor:::kendall_tau_a(X, Y), acor:::kendall_tau_b(X, Y), tolerance = 1e-10)
})

test_that("goodman_kruskal_gamma is correct via brute force", {
  set.seed(77)
  n <- 50
  X <- sample(1:10, n, replace = TRUE)
  Y <- sample(1:10, n, replace = TRUE)
  
  # Brute force
  pairs <- combn(n, 2)
  dx <- X[pairs[2,]] - X[pairs[1,]]
  dy <- Y[pairs[2,]] - Y[pairs[1,]]
  C <- sum(dx * dy > 0)
  D <- sum(dx * dy < 0)
  expected <- (C - D) / (C + D)
  
  result <- acor:::goodman_kruskal_gamma(X, Y)
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("goodman_kruskal_gamma equals tau-a when no ties", {
  set.seed(33)
  n <- 100
  X <- rnorm(n)
  Y <- rnorm(n)
  
  # No ties means C + D = choose(n,2), so gamma = tau-a
  expect_equal(acor:::goodman_kruskal_gamma(X, Y), acor:::kendall_tau_a(X, Y), tolerance = 1e-10)
})

test_that("all measures handle perfect concordance", {
  X <- 1:20
  Y <- 1:20
  
  expect_equal(acor:::kendall_tau_a(X, Y), 1.0)
  expect_equal(acor:::kendall_tau_b(X, Y), 1.0)
  expect_equal(acor:::goodman_kruskal_gamma(X, Y), 1.0)
})

test_that("all measures handle perfect discordance", {
  X <- 1:20
  Y <- 20:1
  
  expect_equal(acor:::kendall_tau_a(X, Y), -1.0)
  expect_equal(acor:::kendall_tau_b(X, Y), -1.0)
  expect_equal(acor:::goodman_kruskal_gamma(X, Y), -1.0)
})

test_that("duplicate export attribute on kendall_tau_a does not cause issues", {
  # kendall_tau_a has // [[Rcpp::export]] twice in the source file
  # Just verify it's callable
  expect_true(is.numeric(acor:::kendall_tau_a(1:10, 10:1)))
})

test_that("comp_spearman_rho_b matches cor(method='spearman') without ties", {
  set.seed(42)
  n <- 100
  X <- rnorm(n)
  Y <- rnorm(n)
  
  x_rank <- rank(X)
  y_rank <- rank(Y)
  
  expected <- cor(X, Y, method = "spearman")
  expect_equal(acor:::comp_spearman_rho_b(y_rank, x_rank), expected, tolerance = 1e-10)
})

test_that("comp_spearman_rho_b matches cor(method='spearman') with ties", {
  set.seed(123)
  n <- 80
  X <- sample(1:5, n, replace = TRUE)
  Y <- sample(1:5, n, replace = TRUE)
  
  x_rank <- rank(X)
  y_rank <- rank(Y)
  
  expected <- cor(X, Y, method = "spearman")
  expect_equal(acor:::comp_spearman_rho_b(y_rank, x_rank), expected, tolerance = 1e-10)
})

test_that("comp_spearman_rho_a matches cor(method='spearman') without ties", {
  set.seed(55)
  n <- 100
  X <- rnorm(n)
  Y <- rnorm(n)
  
  x_rank <- rank(X)
  y_rank <- rank(Y)
  
  expected <- ((n^3-n)/n^3)*cor(X, Y, method = "spearman")
  expect_equal(acor:::comp_spearman_rho_a(y_rank, x_rank), expected, tolerance = 1e-10)
})

test_that("comp_spearman_rho_a and comp_spearman_rho_b differ with ties", {
  set.seed(77)
  n <- 100
  X <- sample(1:5, n, replace = TRUE)
  Y <- sample(1:5, n, replace = TRUE)
  
  x_rank <- rank(X)
  y_rank <- rank(Y)
  
  expect_false(isTRUE(all.equal(
    acor:::comp_spearman_rho_a(y_rank, x_rank),
    acor:::comp_spearman_rho_b(y_rank, x_rank),
    tolerance = 1e-10
  )))
})

test_that("comp_spearman_rho_a and comp_spearman_rho_b agree without ties", {
  set.seed(77)
  n <- 100
  X <- rnorm(n)
  Y <- rnorm(n)
  
  x_rank <- rank(X)
  y_rank <- rank(Y)
  
  rho_a <- acor:::comp_spearman_rho_a(y_rank, x_rank)
  rho_b <- acor:::comp_spearman_rho_b(y_rank, x_rank)
  expect_equal(rho_b / rho_a, n^2 / (n^2 - 1), tolerance = 1e-10)
})

test_that("comp_spearman_rho_a differs from cor(method='spearman') with ties", {
  set.seed(99)
  n <- 80
  X <- sample(1:5, n, replace = TRUE)
  Y <- sample(1:5, n, replace = TRUE)
  
  x_rank <- rank(X)
  y_rank <- rank(Y)
  
  rho_a <- acor:::comp_spearman_rho_a(y_rank, x_rank)
  rho_cor <- ((n^3-n)/n^3)*cor(X, Y, method = "spearman")
  
  # rho_a is uncorrected, so it should NOT match cor() when there are ties
  expect_false(isTRUE(all.equal(rho_a, rho_cor, tolerance = 1e-10)))
})

test_that("comp_pearson_rho matches cor(method='pearson') without ties", {
  set.seed(42)
  n <- 100
  X <- rnorm(n)
  Y <- rnorm(n)
  
  expected <- cor(X, Y, method = "pearson")
  expect_equal(acor:::comp_pearson_rho(X, Y), expected, tolerance = 1e-10)
})

test_that("comp_pearson_rho matches cor(method='pearson') with ties", {
  set.seed(123)
  n <- 80
  X <- sample(1:5, n, replace = TRUE)
  Y <- sample(1:5, n, replace = TRUE)
  
  expected <- cor(X, Y, method = "pearson")
  expect_equal(acor:::comp_pearson_rho(X, Y), expected, tolerance = 1e-10)
})

test_that("comp_pearson_rho handles perfect correlation", {
  X <- 1:20
  Y <- 1:20
  expect_equal(acor:::comp_pearson_rho(X, Y), 1.0)
})

test_that("comp_pearson_rho handles perfect negative correlation", {
  X <- 1:20
  Y <- 20:1
  expect_equal(acor:::comp_pearson_rho(X, Y), -1.0)
})


test_that("kendall_tau_b_mod matches brute force triple-tie version", {
  set.seed(42)
  n <- 50
  X <- sample(1:5, n, replace = TRUE)
  Y <- sample(1:5, n, replace = TRUE)
  
  # Brute force O(n^3) triple tie counting
  ties_x <- 0
  ties_y <- 0
  for (i in 3:n) {
    for (j in 2:(i - 1)) {
      for (k in 1:(j - 1)) {
        ties_x <- ties_x + ifelse(X[i] == X[j] & X[j] == X[k], 1, 0)
        ties_y <- ties_y + ifelse(Y[i] == Y[j] & Y[j] == Y[k], 1, 0)
      }
    }
  }
  ties_x <- ties_x / choose(n, 3)
  ties_y <- ties_y / choose(n, 3)
  
  tau_a <- kendall_tau_a(X, Y)
  expected <- tau_a / sqrt((1 - ties_x) * (1 - ties_y))
  
  expect_equal(kendall_tau_b_mod(X, Y), expected, tolerance = 1e-10)
})

test_that("tau-a IID variance matches brute force (no ties)", {
  set.seed(42)
  n <- 100
  X <- rnorm(n)
  Y <- rnorm(n)
  
  result <- acor:::compute_tau_a_variance(X, Y, IID = TRUE)
  
  G_XY <- Vectorize(function(x_val, y_val) (mean(X <= x_val & Y <= y_val) + mean(X <= x_val & Y < y_val) + mean(X < x_val & Y <= y_val) + mean(X < x_val & Y < y_val)) / 4)
  G_X <- Vectorize(function(x_val) (mean(X < x_val) + mean(X <= x_val)) / 2)
  G_Y <- Vectorize(function(y_val) (mean(Y < y_val) + mean(Y <= y_val)) / 2)
  
  tau <- acor:::kendall_tau_a(X, Y)
  var_expected <- 4 * mean((4 * G_XY(X, Y) - 2 * (G_X(X) + G_Y(Y)) + 1 - tau)^2)
  X_TieProb3 <- sum((table(X) / n)^3)
  Y_TieProb3 <- sum((table(Y) / n)^3)
  var_ind_expected <- 4 / 9 * (1 - X_TieProb3) * (1 - Y_TieProb3)
  
  expect_equal(result$tau_a, tau, tolerance = 1e-10)
  expect_equal(result$var, var_expected, tolerance = 1e-10)
  expect_equal(result$var_ind, var_ind_expected, tolerance = 1e-10)
})

test_that("tau-a IID variance matches brute force (with ties)", {
  set.seed(123)
  n <- 50
  X <- sample(1:5, n, replace = TRUE)
  Y <- sample(1:5, n, replace = TRUE)
  
  result <- acor:::compute_tau_a_variance(X, Y, IID = TRUE)
  
  G_XY <- Vectorize(function(x_val, y_val) (mean(X <= x_val & Y <= y_val) + mean(X <= x_val & Y < y_val) + mean(X < x_val & Y <= y_val) + mean(X < x_val & Y < y_val)) / 4)
  G_X <- Vectorize(function(x_val) (mean(X < x_val) + mean(X <= x_val)) / 2)
  G_Y <- Vectorize(function(y_val) (mean(Y < y_val) + mean(Y <= y_val)) / 2)
  
  tau <- acor:::kendall_tau_a(X, Y)
  var_expected <- 4 * mean((4 * G_XY(X, Y) - 2 * (G_X(X) + G_Y(Y)) + 1 - tau)^2)
  
  X_TieProb3 <- sum((table(X) / n)^3)
  Y_TieProb3 <- sum((table(Y) / n)^3)
  var_ind_expected <- 4 / 9 * (1 - X_TieProb3) * (1 - Y_TieProb3)
  
  expect_equal(result$tau_a, tau, tolerance = 1e-10)
  expect_equal(result$var, var_expected, tolerance = 1e-10)
  expect_equal(result$var_ind, var_ind_expected, tolerance = 1e-10)
})

test_that("tau-a HAC variance matches brute force (no ties)", {
  set.seed(42)
  n <- 100
  X <- rnorm(n)
  Y <- rnorm(n)
  
  result <- acor:::compute_tau_a_variance(X, Y, IID = FALSE)
  
  tau <- acor:::kendall_tau_a(X, Y)
  var_expected <- Tau_LRV(X, Y, tau)
  var_ind_expected <- Tau_ind_LRV(X, Y)
  
  expect_equal(result$tau_a, tau, tolerance = 1e-10)
  expect_equal(result$var, var_expected, tolerance = 1e-10)
  expect_equal(result$var_ind, var_ind_expected, tolerance = 1e-10)
})

test_that("tau-a HAC variance matches brute force (with ties)", {
  set.seed(123)
  n <- 50
  X <- sample(1:5, n, replace = TRUE)
  Y <- sample(1:5, n, replace = TRUE)
  
  result <- acor:::compute_tau_a_variance(X, Y, IID = FALSE)
  
  tau <- acor:::kendall_tau_a(X, Y)
  var_expected <- Tau_LRV(X, Y, tau)
  var_ind_expected <- Tau_ind_LRV(X, Y)
  
  expect_equal(result$tau_a, tau, tolerance = 1e-10)
  expect_equal(result$var, var_expected, tolerance = 1e-10)
  expect_equal(result$var_ind, var_ind_expected, tolerance = 1e-10)
})


test_that("multivariate tau-a IID covariance matches brute force (no ties)", {
  set.seed(42)
  n <- 100
  X <- matrix(rnorm(200), ncol = 2)
  Y <- rnorm(n)
  
  result <- acor:::compute_tau_a_multivariate_variance(X, Y, IID = TRUE)
  
  G_Y <- Vectorize(function(y_val) (mean(Y < y_val) + mean(Y <= y_val)) / 2)
  g_Y <- G_Y(Y)
  
  K_tau_mat <- matrix(0, nrow = n, ncol = 2)
  for (k in 1:2) {
    Xk <- X[, k]
    G_XY <- Vectorize(function(x_val, y_val) (mean(Xk <= x_val & Y <= y_val) + mean(Xk <= x_val & Y < y_val) + mean(Xk < x_val & Y <= y_val) + mean(Xk < x_val & Y < y_val)) / 4)
    G_X <- Vectorize(function(x_val) (mean(Xk < x_val) + mean(Xk <= x_val)) / 2)
    tau_k <- acor:::kendall_tau_a(Xk, Y)
    K_tau_mat[, k] <- 4 * G_XY(Xk, Y) - 2 * (G_X(Xk) + g_Y) + 1 - tau_k
  }
  
  Sigma_expected <- 4 * (t(K_tau_mat) %*% K_tau_mat) / n
  
  expect_equal(result$tau_a_vector, c(acor:::kendall_tau_a(X[,1], Y), acor:::kendall_tau_a(X[,2], Y)), tolerance = 1e-10)
  expect_equal(result$Sigma, Sigma_expected, tolerance = 1e-10)
})

test_that("multivariate tau-a IID covariance matches brute force (with ties)", {
  set.seed(123)
  n <- 50
  X <- cbind(sample(1:5, n, replace = TRUE), sample(1:5, n, replace = TRUE))
  Y <- sample(1:5, n, replace = TRUE)
  
  result <- acor:::compute_tau_a_multivariate_variance(X, Y, IID = TRUE)
  
  G_Y <- Vectorize(function(y_val) (mean(Y < y_val) + mean(Y <= y_val)) / 2)
  g_Y <- G_Y(Y)
  
  K_tau_mat <- matrix(0, nrow = n, ncol = 2)
  for (k in 1:2) {
    Xk <- X[, k]
    G_XY <- Vectorize(function(x_val, y_val) (mean(Xk <= x_val & Y <= y_val) + mean(Xk <= x_val & Y < y_val) + mean(Xk < x_val & Y <= y_val) + mean(Xk < x_val & Y < y_val)) / 4)
    G_X <- Vectorize(function(x_val) (mean(Xk < x_val) + mean(Xk <= x_val)) / 2)
    tau_k <- acor:::kendall_tau_a(Xk, Y)
    K_tau_mat[, k] <- 4 * G_XY(Xk, Y) - 2 * (G_X(Xk) + g_Y) + 1 - tau_k
  }
  
  Sigma_expected <- 4 * (t(K_tau_mat) %*% K_tau_mat) / n
  
  expect_equal(result$Sigma, Sigma_expected, tolerance = 1e-10)
})

test_that("multivariate tau-a HAC covariance matches brute force (no ties)", {
  set.seed(42)
  n <- 100
  X <- matrix(rnorm(200), ncol = 2)
  Y <- rnorm(n)
  
  result <- acor:::compute_tau_a_multivariate_variance(X, Y, IID = FALSE)
  
  G_Y <- Vectorize(function(y_val) (mean(Y < y_val) + mean(Y <= y_val)) / 2)
  g_Y <- G_Y(Y)
  
  K_tau_mat <- matrix(0, nrow = n, ncol = 2)
  for (k in 1:2) {
    Xk <- X[, k]
    G_XY <- Vectorize(function(x_val, y_val) (mean(Xk <= x_val & Y <= y_val) + mean(Xk <= x_val & Y < y_val) + mean(Xk < x_val & Y <= y_val) + mean(Xk < x_val & Y < y_val)) / 4)
    G_X <- Vectorize(function(x_val) (mean(Xk < x_val) + mean(Xk <= x_val)) / 2)
    tau_k <- acor:::kendall_tau_a(Xk, Y)
    K_tau_mat[, k] <- 4 * G_XY(Xk, Y) - 2 * (G_X(Xk) + g_Y) + 1 - tau_k
  }
  
  # IID component
  Sigma_iid <- 4 * (t(K_tau_mat) %*% K_tau_mat) / n
  
  # HAC correction
  b <- floor(2 * n^(1/3))
  Sigma_hac <- matrix(0, nrow = 2, ncol = 2)
  for (h in 1:b) {
    omega <- 1 - h / (b + 1)
    K_lag <- K_tau_mat[1:(n - h), , drop = FALSE]
    K_lead <- K_tau_mat[(h + 1):n, , drop = FALSE]
    autocov_h <- (t(K_lag) %*% K_lead + t(K_lead) %*% K_lag) / n
    Sigma_hac <- Sigma_hac + omega * autocov_h
  }
  
  Sigma_expected <- Sigma_iid + 4 * Sigma_hac
  
  expect_equal(result$Sigma, Sigma_expected, tolerance = 1e-10)
})

test_that("multivariate tau-a HAC covariance matches brute force (with ties)", {
  set.seed(123)
  n <- 50
  X <- cbind(sample(1:5, n, replace = TRUE), sample(1:5, n, replace = TRUE))
  Y <- sample(1:5, n, replace = TRUE)
  
  result <- acor:::compute_tau_a_multivariate_variance(X, Y, IID = FALSE)
  
  G_Y <- Vectorize(function(y_val) (mean(Y < y_val) + mean(Y <= y_val)) / 2)
  g_Y <- G_Y(Y)
  
  K_tau_mat <- matrix(0, nrow = n, ncol = 2)
  for (k in 1:2) {
    Xk <- X[, k]
    G_XY <- Vectorize(function(x_val, y_val) (mean(Xk <= x_val & Y <= y_val) + mean(Xk <= x_val & Y < y_val) + mean(Xk < x_val & Y <= y_val) + mean(Xk < x_val & Y < y_val)) / 4)
    G_X <- Vectorize(function(x_val) (mean(Xk < x_val) + mean(Xk <= x_val)) / 2)
    tau_k <- acor:::kendall_tau_a(Xk, Y)
    K_tau_mat[, k] <- 4 * G_XY(Xk, Y) - 2 * (G_X(Xk) + g_Y) + 1 - tau_k
  }
  
  Sigma_iid <- 4 * (t(K_tau_mat) %*% K_tau_mat) / n
  
  b <- floor(2 * n^(1/3))
  Sigma_hac <- matrix(0, nrow = 2, ncol = 2)
  for (h in 1:b) {
    omega <- 1 - h / (b + 1)
    K_lag <- K_tau_mat[1:(n - h), , drop = FALSE]
    K_lead <- K_tau_mat[(h + 1):n, , drop = FALSE]
    autocov_h <- (t(K_lag) %*% K_lead + t(K_lead) %*% K_lag) / n
    Sigma_hac <- Sigma_hac + omega * autocov_h
  }
  
  Sigma_expected <- Sigma_iid + 4 * Sigma_hac
  
  expect_equal(result$Sigma, Sigma_expected, tolerance = 1e-10)
})

### SPearman rho ####
# === Compare rho_a to AGC in no-ties case ===

test_that("rho_a variance approximately equals AGC variance when no ties (IID)", {
  set.seed(42)
  n <- 100
  X <- rnorm(n)
  Y <- rnorm(n)
  
  x_rank <- rank(X, ties.method = "average")
  y_rank <- rank(Y, ties.method = "average")
  
  result_rho <- acor:::compute_rho_a_variance(x_rank, y_rank, IID = TRUE)
  result_agc <- acor:::compute_agc_variance_auto(y_rank, x_rank, IID = TRUE, version = "v2")
  
  expect_equal(result_rho$rho_a, result_agc$agc, tolerance = 1e-3)
  expect_equal(result_rho$var, result_agc$var, tolerance = 1e-3)
  expect_equal(result_rho$var_ind, result_agc$var_ind, tolerance = 1e-3)
})

test_that("rho_a variance approximately equals AGC variance when no ties (HAC)", {
  set.seed(42)
  n <- 100
  X <- rnorm(n)
  Y <- rnorm(n)
  
  x_rank <- rank(X, ties.method = "average")
  y_rank <- rank(Y, ties.method = "average")
  
  result_rho <- acor:::compute_rho_a_variance(x_rank, y_rank, IID = FALSE)
  result_agc <- acor:::compute_agc_variance_auto(y_rank, x_rank, IID = FALSE, version = "v2")
  
  expect_equal(result_rho$rho_a, result_agc$agc, tolerance = 1e-3)
  expect_equal(result_rho$var, result_agc$var, tolerance = 1e-3)
  expect_equal(result_rho$var_ind, result_agc$var_ind, tolerance = 1e-3)
})

# === Compare rho_a to brute force (IID) ===

test_that("rho_a IID variance matches brute force (no ties)", {
  set.seed(42)
  n <- 100
  X <- rnorm(n)
  Y <- rnorm(n)
  
  x_rank <- rank(X, ties.method = "average")
  y_rank <- rank(Y, ties.method = "average")
  
  result <- acor:::compute_rho_a_variance(x_rank, y_rank, IID = TRUE)
  
  G_XY <- Vectorize(function(x_val, y_val) (mean(X <= x_val & Y <= y_val) + mean(X <= x_val & Y < y_val) + mean(X < x_val & Y <= y_val) + mean(X < x_val & Y < y_val)) / 4)
  G_X <- Vectorize(function(x_val) (mean(X < x_val) + mean(X <= x_val)) / 2)
  G_Y <- Vectorize(function(y_val) (mean(Y < y_val) + mean(Y <= y_val)) / 2)
  g_x <- Vectorize(function(x_val) mean(G_XY(x_val, Y)))
  g_y <- Vectorize(function(y_val) mean(G_XY(X, y_val)))
  
  rho <- 12 * (n - 1) / n^3 * stats::cov(X, Y, method = "spearman")
  G_XX <- G_X(X)
  G_YY <- G_Y(Y)
  var_expected <- 9 * mean((4 * (g_x(X) + g_y(Y) + G_XX * G_YY - G_XX - G_YY) + 1 - rho)^2)
  var_ind_expected <- 1
  
  expect_equal(result$rho_a, rho, tolerance = 1e-10)
  expect_equal(result$var, var_expected, tolerance = 1e-10)
  expect_equal(result$var_ind, var_ind_expected, tolerance = 1e-3)
})

test_that("rho_a IID variance matches brute force (with ties)", {
  set.seed(123)
  n <- 50
  X <- sample(1:5, n, replace = TRUE)
  Y <- sample(1:5, n, replace = TRUE)
  
  x_rank <- rank(X, ties.method = "average")
  y_rank <- rank(Y, ties.method = "average")
  
  result <- acor:::compute_rho_a_variance(x_rank, y_rank, IID = TRUE)
  
  G_XY <- Vectorize(function(x_val, y_val) (mean(X <= x_val & Y <= y_val) + mean(X <= x_val & Y < y_val) + mean(X < x_val & Y <= y_val) + mean(X < x_val & Y < y_val)) / 4)
  G_X <- Vectorize(function(x_val) (mean(X < x_val) + mean(X <= x_val)) / 2)
  G_Y <- Vectorize(function(y_val) (mean(Y < y_val) + mean(Y <= y_val)) / 2)
  g_x <- Vectorize(function(x_val) mean(G_XY(x_val, Y)))
  g_y <- Vectorize(function(y_val) mean(G_XY(X, y_val)))
  
  rho <- 12 * (n - 1) / n^3 * stats::cov(X, Y, method = "spearman")
  G_XX <- G_X(X)
  G_YY <- G_Y(Y)
  var_expected <- 9 * mean((4 * (g_x(X) + g_y(Y) + G_XX * G_YY - G_XX - G_YY) + 1 - rho)^2)
  
  X_TieProb3 <- sum((table(X) / n)^3)
  Y_TieProb3 <- sum((table(Y) / n)^3)
  var_ind_expected <- (1 - X_TieProb3) * (1 - Y_TieProb3)
  
  expect_equal(result$rho_a, rho, tolerance = 1e-10)
  expect_equal(result$var, var_expected, tolerance = 1e-10)
  expect_equal(result$var_ind, var_ind_expected, tolerance = 1e-10)
})

# === Compare rho_a to brute force (HAC) ===

test_that("rho_a HAC variance matches brute force (no ties)", {
  set.seed(42)
  n <- 100
  X <- rnorm(n)
  Y <- rnorm(n)
  
  x_rank <- rank(X, ties.method = "average")
  y_rank <- rank(Y, ties.method = "average")
  
  result <- acor:::compute_rho_a_variance(x_rank, y_rank, IID = FALSE)
  
  rho <- 12 * (n - 1) / n^3 * stats::cov(X, Y, method = "spearman")
  var_expected <- SRho_LRV(X, Y, rho)
  var_ind_expected <- 9 / 4 * Tau_ind_LRV(X, Y)
  
  expect_equal(result$rho_a, rho, tolerance = 1e-10)
  expect_equal(result$var, var_expected, tolerance = 1e-10)
  expect_equal(result$var_ind, var_ind_expected, tolerance = 1e-10)
})

test_that("rho_a HAC variance matches brute force (with ties)", {
  set.seed(123)
  n <- 50
  X <- sample(1:5, n, replace = TRUE)
  Y <- sample(1:5, n, replace = TRUE)
  
  x_rank <- rank(X, ties.method = "average")
  y_rank <- rank(Y, ties.method = "average")
  
  result <- acor:::compute_rho_a_variance(x_rank, y_rank, IID = FALSE)
  
  rho <- 12 * (n - 1) / n^3 * stats::cov(X, Y, method = "spearman")
  var_expected <- SRho_LRV(X, Y, rho)
  var_ind_expected <- 9 / 4 * Tau_ind_LRV(X, Y)
  
  expect_equal(result$rho_a, rho, tolerance = 1e-10)
  expect_equal(result$var, var_expected, tolerance = 1e-10)
  expect_equal(result$var_ind, var_ind_expected, tolerance = 1e-10)
})

# === acor.test integration ===

test_that("acor.test works for rho_a single predictor", {
  set.seed(42)
  n <- 100
  X <- rnorm(n)
  Y <- rnorm(n)
  
  result <- acor:::acor.test(X, Y, method = "rho_a")
  
  expect_s3_class(result, "acor_htest")
  expect_equal(unname(result$estimate), acor(X, Y, method = "rho_a")$estimate, tolerance = 1e-10)
  expect_true(is.numeric(result$p.value))
  expect_true(result$p.value >= 0 && result$p.value <= 1)
})

test_that("acor.test works for rho_a multiple predictors", {
  set.seed(42)
  n <- 100
  X <- matrix(rnorm(200), ncol = 2)
  Y <- rnorm(n)
  
  result <- acor:::acor.test(X, Y, method = "rho_a")
  
  expect_s3_class(result, "acor_htest")
  expect_length(result$estimate, 2)
  expect_true(is.numeric(result$p.value))
})
