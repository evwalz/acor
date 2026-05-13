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

test_that("compute_pearson_variance matches cor() and returns positive variances", {
  set.seed(314)
  n <- 200
  X <- rnorm(n)
  Y <- 0.4 * X + rnorm(n, sd = 0.7)
  
  result <- acor:::compute_pearson_variance(X, Y, IID = TRUE)
  
  expect_equal(result$estimate, cor(X, Y, method = "pearson"), tolerance = 1e-10)
  expect_true(result$var > 0)
  expect_equal(result$var_ind, 1)
})

test_that("pearson HAC variance matches legacy delta-method implementation", {
  set.seed(316)
  n <- 200
  X <- cumsum(rnorm(n))
  Y <- 0.35 * X + cumsum(rnorm(n))

  result <- acor:::compute_pearson_variance(X, Y, IID = FALSE)
  r <- acor:::comp_pearson_rho(X, Y)

  Rho_LRV <- function(X, Y) {
    n <- length(X)
    b <- floor(2 * n^(1 / 3))
    h <- 1:(n - 1)
    w <- pmax(1 - abs(h) / (b + 1), 0)

    mean_x <- mean(X)
    mean_y <- mean(Y)
    sigma_xy <- (n - 1) / n * stats::cov(X, Y)
    var_x <- (n - 1) / n * stats::var(X)
    var_y <- (n - 1) / n * stats::var(Y)

    x_autoc <- stats::acf(X, plot = FALSE, type = "covariance",
                          demean = TRUE, lag.max = n - 1)$acf
    y_autoc <- stats::acf(Y, plot = FALSE, type = "covariance",
                          demean = TRUE, lag.max = n - 1)$acf
    x2_autoc <- stats::acf(X^2, plot = FALSE, type = "covariance",
                           demean = TRUE, lag.max = n - 1)$acf
    y2_autoc <- stats::acf(Y^2, plot = FALSE, type = "covariance",
                           demean = TRUE, lag.max = n - 1)$acf
    xy_autoc <- stats::acf(X * Y, plot = FALSE, type = "covariance",
                           demean = TRUE, lag.max = n - 1)$acf

    xy_crossc <- stats::ccf(X, Y, plot = FALSE, type = "covariance",
                            demean = TRUE, lag.max = n - 1)$acf
    xx2_crossc <- stats::ccf(X, X^2, plot = FALSE, type = "covariance",
                             demean = TRUE, lag.max = n - 1)$acf
    yx2_crossc <- stats::ccf(Y, X^2, plot = FALSE, type = "covariance",
                             demean = TRUE, lag.max = n - 1)$acf
    xy2_crossc <- stats::ccf(X, Y^2, plot = FALSE, type = "covariance",
                             demean = TRUE, lag.max = n - 1)$acf
    yy2_crossc <- stats::ccf(Y, Y^2, plot = FALSE, type = "covariance",
                             demean = TRUE, lag.max = n - 1)$acf
    x2y2_crossc <- stats::ccf(X^2, Y^2, plot = FALSE, type = "covariance",
                              demean = TRUE, lag.max = n - 1)$acf
    xxy_crossc <- stats::ccf(X, X * Y, plot = FALSE, type = "covariance",
                             demean = TRUE, lag.max = n - 1)$acf
    yxy_crossc <- stats::ccf(Y, X * Y, plot = FALSE, type = "covariance",
                             demean = TRUE, lag.max = n - 1)$acf
    x2xy_crossc <- stats::ccf(X^2, X * Y, plot = FALSE, type = "covariance",
                              demean = TRUE, lag.max = n - 1)$acf
    y2xy_crossc <- stats::ccf(Y^2, X * Y, plot = FALSE, type = "covariance",
                              demean = TRUE, lag.max = n - 1)$acf

    x_lrv <- sum(x_autoc[1], 2 * (w * x_autoc[-1]))
    y_lrv <- sum(y_autoc[1], 2 * (w * y_autoc[-1]))
    x2_lrv <- sum(x2_autoc[1], 2 * (w * x2_autoc[-1]))
    y2_lrv <- sum(y2_autoc[1], 2 * (w * y2_autoc[-1]))
    xy_lrv <- sum(xy_autoc[1], 2 * (w * xy_autoc[-1]))
    xy_lrv_c <- sum(c(sort(w), 1, w) * xy_crossc)
    xx2_lrv_c <- sum(c(sort(w), 1, w) * xx2_crossc)
    yx2_lrv_c <- sum(c(sort(w), 1, w) * yx2_crossc)
    xy2_lrv_c <- sum(c(sort(w), 1, w) * xy2_crossc)
    yy2_lrv_c <- sum(c(sort(w), 1, w) * yy2_crossc)
    x2y2_lrv_c <- sum(c(sort(w), 1, w) * x2y2_crossc)
    xxy_lrv_c <- sum(c(sort(w), 1, w) * xxy_crossc)
    yxy_lrv_c <- sum(c(sort(w), 1, w) * yxy_crossc)
    x2xy_lrv_c <- sum(c(sort(w), 1, w) * x2xy_crossc)
    y2xy_lrv_c <- sum(c(sort(w), 1, w) * y2xy_crossc)

    Sigma <- diag(c(x_lrv, y_lrv, x2_lrv, y2_lrv, xy_lrv))
    Sigma[upper.tri(Sigma)] <- c(
      xy_lrv_c, xx2_lrv_c, yx2_lrv_c, xy2_lrv_c, yy2_lrv_c,
      x2y2_lrv_c, xxy_lrv_c, yxy_lrv_c, x2xy_lrv_c, y2xy_lrv_c
    )
    Sigma[lower.tri(Sigma)] <- t(Sigma)[lower.tri(Sigma)]

    A <- matrix(c(
      -2 * mean_x, -mean_y,
      0, 0,
      -mean_x, -2 * mean_y,
      1, 0, 0,
      0, 0, 1,
      0, 1, 0
    ), nrow = 3)
    B <- c(
      -0.5 * sigma_xy / sqrt(var_y) * sqrt(var_x)^(-3),
      1 / sqrt(var_x * var_y),
      -0.5 * sigma_xy / sqrt(var_x) * sqrt(var_y)^(-3)
    )

    as.numeric(B %*% A %*% Sigma %*% t(A) %*% B)
  }

  expect_equal(result$estimate, r, tolerance = 1e-10)
  expect_equal(result$var, Rho_LRV(X, Y), tolerance = 1e-10)
})

test_that("acor.test works for pearson single predictor", {
  set.seed(315)
  n <- 500
  X <- rnorm(n)
  Y <- 0.3 * X + rnorm(n, sd = 0.8)
  
  result <- acor:::acor.test(X, Y, method = "pearson", variance = "plugin")
  ref <- cor.test(X, Y, method = "pearson")
  
  expect_equal(unname(result$estimate), unname(ref$estimate), tolerance = 1e-10)
  expect_equal(result$p.value_ind, ref$p.value, tolerance = 0.03)
  expect_true(is.finite(result$conf.int[1]))
  expect_true(is.finite(result$conf.int[2]))
})


test_that("multivariate pearson IID diagonal matches univariate", {
  set.seed(400)
  n <- 200
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  Y <- 0.4 * X1 - 0.2 * X2 + rnorm(n)
  X_mat <- cbind(X1, X2)

  res_multi <- acor:::compute_pearson_multivariate_variance(X_mat, Y, IID = TRUE)
  res_uni1 <- acor:::compute_pearson_variance(X1, Y, IID = TRUE)
  res_uni2 <- acor:::compute_pearson_variance(X2, Y, IID = TRUE)

  expect_equal(res_multi$estimate_vector[1], res_uni1$estimate, tolerance = 1e-12)
  expect_equal(res_multi$estimate_vector[2], res_uni2$estimate, tolerance = 1e-12)
  expect_equal(res_multi$Sigma[1, 1], res_uni1$var, tolerance = 1e-12)
  expect_equal(res_multi$Sigma[2, 2], res_uni2$var, tolerance = 1e-12)
  expect_equal(res_multi$Sigma[1, 2], res_multi$Sigma[2, 1], tolerance = 1e-14)
})


test_that("multivariate pearson IID independence covariance is corr(X)", {
  set.seed(401)
  n <- 300
  X1 <- rnorm(n)
  X2 <- 0.5 * X1 + rnorm(n)
  X3 <- rnorm(n)
  X_mat <- cbind(X1, X2, X3)
  Y <- rnorm(n)

  Sigma_ind <- acor:::ind_covariance_pearson_iid(X_mat)

  expect_equal(diag(Sigma_ind), rep(1, 3), tolerance = 1e-12)

  X_std <- scale(X_mat, center = TRUE, scale = FALSE)
  X_std <- X_std / rep(sqrt(colMeans(X_std^2)), each = n)
  expected_corr <- (t(X_std) %*% X_std) / n
  expect_equal(Sigma_ind, unname(as.matrix(expected_corr)), tolerance = 1e-12)
})


test_that("multivariate pearson HAC diagonal matches univariate HAC", {
  set.seed(402)
  n <- 200
  X1 <- cumsum(rnorm(n))
  X2 <- rnorm(n)
  Y <- 0.3 * X1 + rnorm(n)
  X_mat <- cbind(X1, X2)

  res_multi <- acor:::compute_pearson_multivariate_variance(X_mat, Y, IID = FALSE)
  res_uni1 <- acor:::compute_pearson_variance(X1, Y, IID = FALSE)
  res_uni2 <- acor:::compute_pearson_variance(X2, Y, IID = FALSE)

  expect_equal(res_multi$estimate_vector[1], res_uni1$estimate, tolerance = 1e-12)
  expect_equal(res_multi$estimate_vector[2], res_uni2$estimate, tolerance = 1e-12)
  expect_equal(res_multi$Sigma[1, 1], res_uni1$var, tolerance = 1e-10)
  expect_equal(res_multi$Sigma[2, 2], res_uni2$var, tolerance = 1e-10)
})


test_that("acor.test works for multivariate pearson IID", {
  set.seed(403)
  n <- 200
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  X3 <- rnorm(n)
  Y <- 0.3 * X1 + 0.2 * X2 + rnorm(n)
  X_mat <- cbind(X1, X2, X3)

  result <- acor.test(X_mat, Y, method = "pearson", IID = TRUE, variance = "plugin")

  expect_length(result$estimate, 3)
  expect_equal(nrow(result$pairwise_results), 3)
  expect_true(all(is.finite(result$p.value_ind)))
  expect_true(is.matrix(result$variance))
  expect_equal(dim(result$variance), c(3, 3))
  expect_true(all(diag(result$variance) > 0))
})


test_that("acor.test multivariate pearson HAC returns valid output", {
  set.seed(404)
  n <- 200
  X1 <- as.numeric(arima.sim(model = list(ar = 0.3), n = n))
  X2 <- rnorm(n)
  Y <- 0.3 * X1 + rnorm(n)
  X_mat <- cbind(X1, X2)

  result <- acor.test(X_mat, Y, method = "pearson", IID = FALSE, variance = "plugin")

  expect_length(result$estimate, 2)
  expect_true(all(diag(result$variance) > 0))
  expect_true(all(is.finite(result$p.value_ind)))
})


Rhob_ind_LRV_ref <- function(X, Y) {
  n <- length(X)
  b <- floor(2 * n^(1 / 3))
  h <- 1:(n - 1)
  w <- pmax(1 - abs(h) / (b + 1), 0)

  x_grade <- (rank(X) - 0.5) / n - 0.5
  y_grade <- (rank(Y) - 0.5) / n - 0.5
  x_autoc <- n / (n - 1) *
    stats::acf(x_grade, plot = FALSE, type = "covariance",
               demean = FALSE, lag.max = n - 1)$acf / stats::var(x_grade)
  y_autoc <- n / (n - 1) *
    stats::acf(y_grade, plot = FALSE, type = "covariance",
               demean = FALSE, lag.max = n - 1)$acf / stats::var(y_grade)

  sum(x_autoc[1] * y_autoc[1], 2 * (w * x_autoc[-1] * y_autoc[-1]))
}

Rhob_LRV_ref <- function(X, Y, spearman, spearman_X, spearman_Y) {
  n <- length(X)
  b <- floor(2 * n^(1 / 3))
  h <- 1:(n - 1)
  w <- pmax(1 - abs(h) / (b + 1), 0)

  G_XY <- Vectorize(function(x_val, y_val) {
    (mean(X <= x_val & Y <= y_val) +
       mean(X <= x_val & Y < y_val) +
       mean(X < x_val & Y <= y_val) +
       mean(X < x_val & Y < y_val)) / 4
  })
  G_X <- Vectorize(function(x_val) (mean(X < x_val) + mean(X <= x_val)) / 2)
  G_Y <- Vectorize(function(y_val) (mean(Y < y_val) + mean(Y <= y_val)) / 2)
  g_x <- Vectorize(function(x_val) mean(G_XY(x_val, Y)))
  g_y <- Vectorize(function(y_val) mean(G_XY(X, y_val)))
  x_eq <- Vectorize(function(x_val) mean(X == x_val))
  y_eq <- Vectorize(function(y_val) mean(Y == y_val))

  G_XX <- G_X(X)
  G_YY <- G_Y(Y)
  k_XY_rho <- 4 * (g_x(X) + g_y(Y) + G_XX * G_YY - G_XX - G_YY) + 1 - spearman
  k_X_rho <- 1 - x_eq(X)^2 - spearman_X
  k_Y_rho <- 1 - y_eq(Y)^2 - spearman_Y

  k_XY_rho_autoc <- stats::acf(k_XY_rho, plot = FALSE, type = "covariance",
                               demean = FALSE, lag.max = n - 1)$acf
  k_X_rho_autoc <- stats::acf(k_X_rho, plot = FALSE, type = "covariance",
                              demean = FALSE, lag.max = n - 1)$acf
  k_Y_rho_autoc <- stats::acf(k_Y_rho, plot = FALSE, type = "covariance",
                              demean = FALSE, lag.max = n - 1)$acf
  k_X_rho_crossc <- stats::ccf(k_XY_rho, k_X_rho, plot = FALSE,
                               type = "covariance", lag.max = n - 1)$acf
  k_Y_rho_crossc <- stats::ccf(k_XY_rho, k_Y_rho, plot = FALSE,
                               type = "covariance", lag.max = n - 1)$acf
  k_XY_rho_crossc <- stats::ccf(k_X_rho, k_Y_rho, plot = FALSE,
                                type = "covariance", lag.max = n - 1)$acf

  sigma_rho_sq <- 9 * sum(k_XY_rho_autoc[1], 2 * (w * k_XY_rho_autoc[-1]))
  sigma_rhoX_sq <- 9 * sum(k_X_rho_autoc[1], 2 * (w * k_X_rho_autoc[-1]))
  sigma_rhoY_sq <- 9 * sum(k_Y_rho_autoc[1], 2 * (w * k_Y_rho_autoc[-1]))
  sigma_rhorhoX <- 9 * sum(c(sort(w), 1, w) * k_X_rho_crossc)
  sigma_rhorhoY <- 9 * sum(c(sort(w), 1, w) * k_Y_rho_crossc)
  sigma_rhoXrhoY <- 9 * sum(c(sort(w), 1, w) * k_XY_rho_crossc)

  (sigma_rho_sq -
     spearman * (sigma_rhorhoX / spearman_X + sigma_rhorhoY / spearman_Y) +
     spearman^2 / 4 * (
       sigma_rhoX_sq / spearman_X^2 +
         sigma_rhoY_sq / spearman_Y^2 +
         (2 * sigma_rhoXrhoY) / spearman_Y / spearman_X
     )) / (spearman_X * spearman_Y)
}

test_that("compute_rho_b_variance matches cor() and legacy IID formula", {
  set.seed(317)
  n <- 150
  X <- sample(1:8, n, replace = TRUE)
  Y <- sample(1:8, n, replace = TRUE)

  result <- acor:::compute_rho_b_variance(X, Y, IID = TRUE)
  rho <- 12 * (n - 1) / n^3 * stats::cov(X, Y, method = "spearman")
  rho_x <- 12 * (n - 1) / n^3 * stats::cov(X, X, method = "spearman")
  rho_y <- 12 * (n - 1) / n^3 * stats::cov(Y, Y, method = "spearman")

  G_XY <- Vectorize(function(x_val, y_val) {
    (mean(X <= x_val & Y <= y_val) +
       mean(X <= x_val & Y < y_val) +
       mean(X < x_val & Y <= y_val) +
       mean(X < x_val & Y < y_val)) / 4
  })
  G_X <- Vectorize(function(x_val) (mean(X < x_val) + mean(X <= x_val)) / 2)
  G_Y <- Vectorize(function(y_val) (mean(Y < y_val) + mean(Y <= y_val)) / 2)
  g_x <- Vectorize(function(x_val) mean(G_XY(x_val, Y)))
  g_y <- Vectorize(function(y_val) mean(G_XY(X, y_val)))
  x_eq <- Vectorize(function(x_val) mean(X == x_val))
  y_eq <- Vectorize(function(y_val) mean(Y == y_val))

  g_xX <- g_x(X)
  g_yY <- g_y(Y)
  G_XX <- G_X(X)
  G_YY <- G_Y(Y)
  x_eqX <- x_eq(X)
  y_eqY <- y_eq(Y)
  var_rho <- 9 * mean((4 * (g_xX + g_yY + G_XX * G_YY - G_XX - G_YY) + 1 - rho)^2)
  var_rho_x <- 9 * mean((1 - x_eqX^2 - rho_x)^2)
  var_rho_y <- 9 * mean((1 - y_eqY^2 - rho_y)^2)
  var_rhorho_x <- 9 * mean((4 * (g_xX + g_yY + G_XX * G_YY - G_XX - G_YY) + 1 - rho) *
                             (1 - x_eqX^2 - rho_x))
  var_rhorho_y <- 9 * mean((4 * (g_xX + g_yY + G_XX * G_YY - G_XX - G_YY) + 1 - rho) *
                             (1 - y_eqY^2 - rho_y))
  var_rho_xrho_y <- 9 * mean((1 - x_eqX^2 - rho_x) * (1 - y_eqY^2 - rho_y))
  var_expected <- (var_rho - rho * (var_rhorho_x / rho_x + var_rhorho_y / rho_y) +
                     0.25 * rho^2 * (var_rho_x / rho_x^2 + var_rho_y / rho_y^2 +
                                       2 * var_rho_xrho_y / (rho_x * rho_y))) /
    (rho_x * rho_y)

  expect_equal(result$rho_b, cor(X, Y, method = "spearman"), tolerance = 1e-10)
  expect_equal(result$var, var_expected, tolerance = 1e-10)
  expect_equal(result$var_ind, 1)
})

test_that("compute_rho_b_variance matches legacy HAC formulas", {
  set.seed(318)
  n <- 150
  X <- cumsum(sample(-1:2, n, replace = TRUE))
  Y <- cumsum(sample(-1:2, n, replace = TRUE))

  result <- acor:::compute_rho_b_variance(X, Y, IID = FALSE)
  rho <- 12 * (n - 1) / n^3 * stats::cov(X, Y, method = "spearman")
  rho_x <- 12 * (n - 1) / n^3 * stats::cov(X, X, method = "spearman")
  rho_y <- 12 * (n - 1) / n^3 * stats::cov(Y, Y, method = "spearman")

  expect_equal(result$rho_b, cor(X, Y, method = "spearman"), tolerance = 1e-10)
  expect_equal(result$var, Rhob_LRV_ref(X, Y, rho, rho_x, rho_y), tolerance = 1e-10)
  expect_equal(result$var_ind, Rhob_ind_LRV_ref(X, Y), tolerance = 1e-10)
})

test_that("acor.test works for rho_b single predictor", {
  set.seed(319)
  n <- 250
  X <- sample(1:10, n, replace = TRUE)
  Y <- sample(1:10, n, replace = TRUE)

  result <- acor:::acor.test(X, Y, method = "rho_b")

  expect_s3_class(result, "acor_htest")
  expect_equal(unname(result$estimate), cor(X, Y, method = "spearman"), tolerance = 1e-10)
  expect_true(is.numeric(result$p.value))
  expect_true(result$p.value >= 0 && result$p.value <= 1)
  expect_true(is.finite(result$conf.int[1]))
  expect_true(is.finite(result$conf.int[2]))
})


test_that("multivariate rho_b m=1 matches univariate variance", {
  set.seed(801)
  n <- 120
  X <- matrix(sample(1:8, n, replace = TRUE), ncol = 1)
  Y <- sample(1:6, n, replace = TRUE)
  uni <- acor:::compute_rho_b_variance(X[, 1], Y, IID = TRUE)
  multi <- acor:::compute_rho_b_multivariate_variance(X, Y, IID = TRUE)
  expect_equal(multi$rho_b_vector, uni$rho_b, tolerance = 1e-12)
  expect_equal(as.numeric(multi$Sigma), uni$var, tolerance = 1e-10)
  expect_equal(as.numeric(multi$Sigma_ind), uni$var_ind, tolerance = 1e-10)
})


test_that("multivariate rho_b IID diagonal matches univariate per predictor", {
  set.seed(802)
  n <- 160
  X1 <- sample(1:8, n, replace = TRUE)
  X2 <- sample(1:7, n, replace = TRUE)
  Y <- sample(1:6, n, replace = TRUE)
  X_mat <- cbind(X1, X2)
  res_multi <- acor:::compute_rho_b_multivariate_variance(X_mat, Y, IID = TRUE)
  res_u1 <- acor:::compute_rho_b_variance(X1, Y, IID = TRUE)
  res_u2 <- acor:::compute_rho_b_variance(X2, Y, IID = TRUE)
  expect_equal(diag(res_multi$Sigma)[1], res_u1$var, tolerance = 1e-9)
  expect_equal(diag(res_multi$Sigma)[2], res_u2$var, tolerance = 1e-9)
  expect_equal(diag(res_multi$Sigma_ind)[1], res_u1$var_ind, tolerance = 1e-10)
  expect_equal(diag(res_multi$Sigma_ind)[2], res_u2$var_ind, tolerance = 1e-10)
})


test_that("multivariate rho_b m=1 HAC close to univariate (HAC estimator tol)", {
  set.seed(803)
  n <- 100
  X <- matrix(sample(1:8, n, replace = TRUE), ncol = 1)
  Y <- sample(1:6, n, replace = TRUE)
  uni <- acor:::compute_rho_b_variance(X[, 1], Y, IID = FALSE)
  multi <- acor:::compute_rho_b_multivariate_variance(X, Y, IID = FALSE)
  expect_equal(multi$rho_b_vector, uni$rho_b, tolerance = 1e-12)
  expect_equal(as.numeric(multi$Sigma), uni$var, tolerance = 0.003)
  expect_equal(as.numeric(multi$Sigma_ind), uni$var_ind, tolerance = 1e-9)
})


test_that("multivariate rho_b HAC diagonal close to univariate per predictor (tol)", {
  set.seed(804)
  n <- 140
  X1 <- sample(1:8, n, replace = TRUE)
  X2 <- sample(1:7, n, replace = TRUE)
  Y <- sample(1:6, n, replace = TRUE)
  X_mat <- cbind(X1, X2)
  res_multi <- acor:::compute_rho_b_multivariate_variance(X_mat, Y, IID = FALSE)
  res_u1 <- acor:::compute_rho_b_variance(X1, Y, IID = FALSE)
  res_u2 <- acor:::compute_rho_b_variance(X2, Y, IID = FALSE)
  expect_equal(diag(res_multi$Sigma)[1], res_u1$var, tolerance = 0.003)
  expect_equal(diag(res_multi$Sigma)[2], res_u2$var, tolerance = 0.003)
  expect_equal(diag(res_multi$Sigma_ind)[1], res_u1$var_ind, tolerance = 1e-9)
  expect_equal(diag(res_multi$Sigma_ind)[2], res_u2$var_ind, tolerance = 1e-9)
})


test_that("acor.test works for rho_b multivariate", {
  set.seed(805)
  n <- 200
  X <- cbind(sample(1:9, n, replace = TRUE), sample(1:7, n, replace = TRUE))
  Y <- sample(1:6, n, replace = TRUE)
  result <- acor:::acor.test(X, Y, method = "rho_b")
  expect_s3_class(result, "acor_htest")
  expect_length(result$estimate, 2)
  expect_true(is.numeric(result$p.value))
  expect_true(result$p.value >= 0 && result$p.value <= 1)
})


test_that("acor.test works for rho_b multivariate HAC", {
  set.seed(806)
  n <- 200
  X <- cbind(sample(1:9, n, replace = TRUE), sample(1:7, n, replace = TRUE))
  Y <- sample(1:6, n, replace = TRUE)
  result <- acor:::acor.test(X, Y, method = "rho_b", IID = FALSE)
  expect_s3_class(result, "acor_htest")
  expect_length(result$estimate, 2)
  expect_false(result$IID)
  expect_equal(dim(result$variance), c(2, 2))
  expect_true(all(diag(result$variance) > 0))
})


test_that("kendall_tau_b equals tau_xy / sqrt(tau_self_X * tau_self_Y)", {
  set.seed(700)
  n <- 90
  X <- sample(1:7, n, replace = TRUE)
  Y <- sample(1:5, n, replace = TRUE)
  tau <- acor:::compute_kendall(X, Y)$expectation
  sx <- acor:::tau_b_tie_preamble(X)$tau_self
  sy <- acor:::tau_b_tie_preamble(Y)$tau_self
  expect_equal(acor:::kendall_tau_b(X, Y), tau / sqrt(sx * sy), tolerance = 1e-12)
})


test_that("compute_tau_b_variance IID matches grad-prime Sigma grad (sqrt tau_b)", {
  set.seed(701)
  n <- 120
  X <- sample(1:6, n, replace = TRUE)
  Y <- sample(1:5, n, replace = TRUE)
  tau <- acor:::compute_kendall(X, Y)$expectation
  pre_x <- acor:::tau_b_tie_preamble(X)
  pre_y <- acor:::tau_b_tie_preamble(Y)
  k_tau <- acor:::K_tau_vec_v2(X, Y, tau)
  grad <- acor:::tau_b_gradient(tau, pre_x$tau_self, pre_y$tau_self)
  S <- acor:::tau_b_iid_covariance(k_tau, pre_x$k_self, pre_y$k_self)
  v_manual <- as.numeric(t(grad) %*% S %*% grad)
  v_pkg <- acor:::compute_tau_b_variance(X, Y, IID = TRUE)$var
  expect_equal(v_pkg, v_manual, tolerance = 1e-12)
})


test_that("multivariate tau_b m=1 matches univariate variance", {
  set.seed(702)
  n <- 100
  X <- matrix(sample(1:6, n, replace = TRUE), ncol = 1)
  Y <- sample(1:5, n, replace = TRUE)
  uni <- acor:::compute_tau_b_variance(X[, 1], Y, IID = TRUE)
  multi <- acor:::compute_tau_b_multivariate_variance(X, Y, IID = TRUE)
  expect_equal(multi$tau_b_vector, uni$tau_b, tolerance = 1e-12)
  expect_equal(as.numeric(multi$Sigma), uni$var, tolerance = 1e-12)
  expect_equal(as.numeric(multi$Sigma_ind), uni$var_ind, tolerance = 1e-12)
})


test_that("multivariate tau_b IID diagonal matches univariate per predictor", {
  set.seed(703)
  n <- 150
  X1 <- sample(1:6, n, replace = TRUE)
  X2 <- sample(1:7, n, replace = TRUE)
  Y <- sample(1:5, n, replace = TRUE)
  X_mat <- cbind(X1, X2)
  res_multi <- acor:::compute_tau_b_multivariate_variance(X_mat, Y, IID = TRUE)
  res_u1 <- acor:::compute_tau_b_variance(X1, Y, IID = TRUE)
  res_u2 <- acor:::compute_tau_b_variance(X2, Y, IID = TRUE)
  expect_equal(diag(res_multi$Sigma)[1], res_u1$var, tolerance = 1e-10)
  expect_equal(diag(res_multi$Sigma)[2], res_u2$var, tolerance = 1e-10)
  expect_equal(diag(res_multi$Sigma_ind)[1], res_u1$var_ind, tolerance = 1e-10)
  expect_equal(diag(res_multi$Sigma_ind)[2], res_u2$var_ind, tolerance = 1e-10)
})


test_that("multivariate tau_b m=1 HAC matches univariate (within HAC estimator tol)", {
  set.seed(705)
  n <- 100
  X <- matrix(sample(1:6, n, replace = TRUE), ncol = 1)
  Y <- sample(1:5, n, replace = TRUE)
  uni <- acor:::compute_tau_b_variance(X[, 1], Y, IID = FALSE)
  multi <- acor:::compute_tau_b_multivariate_variance(X, Y, IID = FALSE)
  expect_equal(multi$tau_b_vector, uni$tau_b, tolerance = 1e-12)
  # Univariate uses acf/ccf LRV; stacked-Psi HAC uses lag cross-products — same
  # estimand, small numerical mismatch (see package stress checks ~2e-3 abs).
  expect_equal(as.numeric(multi$Sigma), uni$var, tolerance = 0.003)
  expect_equal(as.numeric(multi$Sigma_ind), uni$var_ind, tolerance = 1e-10)
})


test_that("multivariate tau_b HAC diagonal matches univariate per predictor (tol)", {
  set.seed(706)
  n <- 150
  X1 <- sample(1:6, n, replace = TRUE)
  X2 <- sample(1:7, n, replace = TRUE)
  Y <- sample(1:5, n, replace = TRUE)
  X_mat <- cbind(X1, X2)
  res_multi <- acor:::compute_tau_b_multivariate_variance(X_mat, Y, IID = FALSE)
  res_u1 <- acor:::compute_tau_b_variance(X1, Y, IID = FALSE)
  res_u2 <- acor:::compute_tau_b_variance(X2, Y, IID = FALSE)
  expect_equal(diag(res_multi$Sigma)[1], res_u1$var, tolerance = 0.003)
  expect_equal(diag(res_multi$Sigma)[2], res_u2$var, tolerance = 0.003)
  expect_equal(diag(res_multi$Sigma_ind)[1], res_u1$var_ind, tolerance = 1e-10)
  expect_equal(diag(res_multi$Sigma_ind)[2], res_u2$var_ind, tolerance = 1e-10)
})


test_that("acor.test works for tau_b multivariate", {
  set.seed(704)
  n <- 200
  X <- cbind(sample(1:8, n, replace = TRUE), sample(1:6, n, replace = TRUE))
  Y <- sample(1:5, n, replace = TRUE)
  result <- acor:::acor.test(X, Y, method = "tau_b")
  expect_s3_class(result, "acor_htest")
  expect_length(result$estimate, 2)
  expect_true(is.numeric(result$p.value))
  expect_true(result$p.value >= 0 && result$p.value <= 1)
})


test_that("acor.test works for tau_b multivariate HAC", {
  set.seed(707)
  n <- 200
  X <- cbind(sample(1:8, n, replace = TRUE), sample(1:6, n, replace = TRUE))
  Y <- sample(1:5, n, replace = TRUE)
  result <- acor:::acor.test(X, Y, method = "tau_b", IID = FALSE)
  expect_s3_class(result, "acor_htest")
  expect_length(result$estimate, 2)
  expect_false(result$IID)
  expect_true(is.matrix(result$variance))
  expect_equal(dim(result$variance), c(2, 2))
  expect_true(all(diag(result$variance) > 0))
  expect_true(is.numeric(result$p.value))
  expect_true(result$p.value >= 0 && result$p.value <= 1)
})


test_that("acor.test works for tau_b single predictor", {
  set.seed(322)
  n <- 250
  X <- sample(1:10, n, replace = TRUE)
  Y <- sample(1:10, n, replace = TRUE)

  result <- acor:::acor.test(X, Y, method = "tau_b")

  expect_s3_class(result, "acor_htest")
  expect_equal(unname(result$estimate), cor(X, Y, method = "kendall"), tolerance = 1e-10)
  expect_true(is.numeric(result$p.value))
  expect_true(result$p.value >= 0 && result$p.value <= 1)
  expect_true(is.finite(result$conf.int[1]))
  expect_true(is.finite(result$conf.int[2]))
})


test_that("acor.test works for gamma single predictor", {
  set.seed(324)
  n <- 250
  X <- sample(1:10, n, replace = TRUE)
  Y <- sample(1:10, n, replace = TRUE)

  result <- acor:::acor.test(X, Y, method = "gamma")

  expect_s3_class(result, "acor_htest")
  expect_equal(unname(result$estimate), acor:::goodman_kruskal_gamma(X, Y), tolerance = 1e-10)
  expect_true(is.numeric(result$p.value))
  expect_true(result$p.value >= 0 && result$p.value <= 1)
  expect_true(is.finite(result$conf.int[1]))
  expect_true(is.finite(result$conf.int[2]))
})


test_that("multivariate gamma IID diagonal matches univariate", {
  set.seed(510)
  n <- 200
  X1 <- sample(1:5, n, replace = TRUE)
  X2 <- sample(1:8, n, replace = TRUE)
  Y <- sample(1:4, n, replace = TRUE)
  X_mat <- cbind(X1, X2)

  res_multi <- acor:::compute_gamma_multivariate_variance(X_mat, Y, IID = TRUE)
  res_uni1 <- acor:::compute_gamma_variance(X1, Y, IID = TRUE)
  res_uni2 <- acor:::compute_gamma_variance(X2, Y, IID = TRUE)

  expect_equal(res_multi$gamma_vector[1], res_uni1$gamma, tolerance = 1e-12)
  expect_equal(res_multi$gamma_vector[2], res_uni2$gamma, tolerance = 1e-12)
  expect_equal(res_multi$Sigma[1, 1], res_uni1$var, tolerance = 1e-12)
  expect_equal(res_multi$Sigma[2, 2], res_uni2$var, tolerance = 1e-12)
  expect_equal(res_multi$Sigma_ind[1, 1], res_uni1$var_ind, tolerance = 1e-12)
  expect_equal(res_multi$Sigma_ind[2, 2], res_uni2$var_ind, tolerance = 1e-12)
  expect_equal(res_multi$Sigma[1, 2], res_multi$Sigma[2, 1], tolerance = 1e-14)
})


test_that("multivariate gamma HAC diagonal matches univariate", {
  set.seed(511)
  n <- 200
  X1 <- sample(1:5, n, replace = TRUE)
  X2 <- sample(1:8, n, replace = TRUE)
  Y <- sample(1:4, n, replace = TRUE)
  X_mat <- cbind(X1, X2)

  res_multi <- acor:::compute_gamma_multivariate_variance(X_mat, Y, IID = FALSE)
  res_uni1 <- acor:::compute_gamma_variance(X1, Y, IID = FALSE)
  res_uni2 <- acor:::compute_gamma_variance(X2, Y, IID = FALSE)

  expect_equal(res_multi$gamma_vector[1], res_uni1$gamma, tolerance = 1e-12)
  expect_equal(res_multi$gamma_vector[2], res_uni2$gamma, tolerance = 1e-12)
  expect_equal(res_multi$Sigma[1, 1], res_uni1$var, tolerance = 1e-4)
  expect_equal(res_multi$Sigma[2, 2], res_uni2$var, tolerance = 1e-4)
})


test_that("acor.test works for multivariate gamma IID", {
  set.seed(512)
  n <- 200
  X1 <- sample(1:5, n, replace = TRUE)
  X2 <- sample(1:8, n, replace = TRUE)
  X3 <- sample(1:3, n, replace = TRUE)
  Y <- sample(1:4, n, replace = TRUE)
  X_mat <- cbind(X1, X2, X3)

  result <- acor.test(X_mat, Y, method = "gamma", IID = TRUE)

  expect_length(result$estimate, 3)
  expect_equal(nrow(result$pairwise_results), 3)
  expect_true(all(is.finite(result$p.value_ind)))
  expect_true(is.matrix(result$variance))
  expect_equal(dim(result$variance), c(3, 3))
  expect_true(all(diag(result$variance) > 0))
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
  result_agc <- acor:::compute_agc_variance_auto(y_rank, x_rank, IID = TRUE)
  
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
  result_agc <- acor:::compute_agc_variance_auto(y_rank, x_rank, IID = FALSE)
  
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
  
  result <- acor:::acor.test(X, Y, method = "rho_a", variance = "plugin")
  
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
  
  result <- acor:::acor.test(X, Y, method = "rho_a", variance = "plugin")
  
  expect_s3_class(result, "acor_htest")
  expect_length(result$estimate, 2)
  expect_true(is.numeric(result$p.value))
})


