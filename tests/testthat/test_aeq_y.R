# Optional aeq-style collapse on numeric Y (acor.test / acor, all methods)

library(testthat)
library(acor)

test_that("bucket_y_aeq_style matches reference aeqSurv output for numeric Y", {
  skip_if_not_installed("survival")
  set.seed(11)
  y <- c(rnorm(30), 1.0, 1.0 + 1e-10, 2.0, pi, pi + 2e-10)
  tol <- sqrt(.Machine$double.eps)
  y_ref <- survival::aeqSurv(survival::Surv(y), tolerance = tol)[, 1]
  y_b <- acor:::bucket_y_aeq_style(y, tol)
  expect_equal(y_b, y_ref, tolerance = 1e-15)
})

test_that("acor.test aeq_y FALSE vs TRUE can change CID estimate on near-ties", {
  set.seed(12)
  x <- rnorm(25)
  y <- c(rnorm(23), 1.0, 1.0 + 1e-10)
  t0 <- acor.test(x, y, method = "cid", aeq_y = FALSE)
  t1 <- acor.test(x, y, method = "cid", aeq_y = TRUE)
  expect_true(is.finite(t0$estimate))
  expect_true(is.finite(t1$estimate))
})

test_that("acor aeq_y uses same Y as manual bucket for CID", {
  skip_if_not_installed("survival")
  set.seed(13)
  x <- rnorm(20)
  y <- c(rnorm(18), 0.5, 0.5 + 1e-10)
  tol <- sqrt(.Machine$double.eps)
  y_b <- acor:::bucket_y_aeq_style(y, tol)
  r0 <- acor(x, y, method = "cid", aeq_y = FALSE)$estimate
  r1 <- acor(x, y, method = "cid", aeq_y = TRUE)$estimate
  r_b <- acor(x, y_b, method = "cid", aeq_y = FALSE)$estimate
  expect_equal(r1, r_b, tolerance = 1e-12)
  expect_true(abs(r0 - r1) > 1e-10)
})

test_that("acor.test pearson can differ when aeq_y merges near-duplicate Y", {
  set.seed(14)
  x <- rnorm(15)
  y <- c(rnorm(13), 1.0, 1.0 + 1e-10)
  p0 <- acor.test(x, y, method = "pearson", aeq_y = FALSE, variance = "plugin")$estimate
  p1 <- acor.test(x, y, method = "pearson", aeq_y = TRUE, variance = "plugin")$estimate
  expect_true(is.finite(p0) && is.finite(p1))
})
