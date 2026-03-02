library(testthat)
library(acor)

set.seed(321)

test_that("AGC correlates with Spearman rho (sign/direction)", {
  x <- rnorm(150)
  y <- 0.4 * x + rnorm(150, sd = 0.6)
  agc <- unname(acor(x, y, method = "agc")$estimate)
  rho <- cor(x, y, method = "spearman")
  expect_true(is.finite(agc))
  expect_true(sign(agc) == sign(rho) || abs(agc) < 1e-6)
})

test_that("AKC correlates with Kendall tau (sign/direction)", {
  x <- rnorm(150)
  y <- -0.3 * x + rnorm(150, sd = 0.5)
  akc <- unname(acor(x, y, method = "akc")$estimate)
  kt <- cor(x, y, method = "kendall")
  expect_true(is.finite(akc))
  expect_true(sign(akc) == sign(kt) || abs(akc) < 1e-6)
})

test_that("CMA roughly matches AUC for binary outcomes", {
  skip_if_not_installed("pROC")
  x <- rnorm(200)
  p <- plogis(0.6 * x)
  y <- rbinom(200, 1, p)
  cma <- unname(acor(x, y, method = "cma")$estimate)
  auc <- as.numeric(pROC::roc(y, x, quiet = TRUE)$auc)
  expect_true(abs(cma - auc) < 0.1)
})

test_that("acor.test returns finite result under independence", {
  x <- rnorm(120)
  y <- rnorm(120)
  res <- acor.test(x, y, method = "akc")
  expect_true(is.finite(res$p.value))
  expect_true(res$p.value >= 0 && res$p.value <= 1)
})

test_that("acor.test handles two predictors and returns covariance-backed stat", {
  set.seed(99)
  y <- rnorm(120)
  x1 <- y + rnorm(120, sd = 0.2)
  x2 <- rnorm(120)
  X <- cbind(x1, x2)
  res <- acor.test(X, y, method = "akc")
  expect_length(res$estimate, 2)
  expect_true(is.finite(res$statistic))
  expect_true(res$p.value >= 0 && res$p.value <= 1)
})