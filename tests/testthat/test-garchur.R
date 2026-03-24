test_that("garchur returns expected structure", {
  set.seed(1)
  y <- cumsum(rnorm(80))
  res <- garchur(y, breaks = 2, model = "ct")
  expect_s3_class(res, "garchur")
  expect_true(is.numeric(res$stat))
  expect_true(length(res$break_dates) == 2)
  expect_true(is.numeric(res$cv5))
  expect_true(res$nobs == 80)
})

test_that("garchur constant model works", {
  set.seed(2)
  y <- cumsum(rnorm(60))
  res <- garchur(y, breaks = 1, model = "c")
  expect_s3_class(res, "garchur")
  expect_equal(length(res$break_dates), 1)
})

test_that("garchur three breaks", {
  set.seed(3)
  y <- cumsum(rnorm(100))
  res <- garchur(y, breaks = 3, model = "ct")
  expect_equal(length(res$break_dates), 3)
})

test_that("garchur GARCH parameters are non-negative", {
  set.seed(4)
  y <- cumsum(rnorm(80))
  res <- garchur(y, breaks = 2)
  expect_true(res$alpha >= 0)
  expect_true(res$beta  >= 0)
  expect_true(res$kappa >= 0)
})

test_that("garchur input validation", {
  expect_error(garchur(rnorm(20)))
  expect_error(garchur(rnorm(50), breaks = 4))
  expect_error(garchur(rnorm(50), trim = 0.5))
})

test_that("print.garchur produces output", {
  set.seed(5)
  y <- cumsum(rnorm(60))
  res <- garchur(y, breaks = 2)
  out <- capture.output(print(res))
  expect_true(length(out) > 0)
})
