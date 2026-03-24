test_that("urstat runs on a random walk with ADF and PP", {
  set.seed(1)
  x <- cumsum(rnorm(50))
  res <- urstat(x, tests = c("ADF", "PP"), strategy = FALSE, stars = FALSE)
  expect_type(res, "list")
  expect_length(res, 1L)
  r <- res[[1]]
  expect_named(r, c("adf", "pp", "decision", "process"), ignore.order = TRUE)
  expect_true(is.character(r$decision))
})

test_that("urstat returns KPSS results", {
  set.seed(2)
  x <- rnorm(50)
  res <- urstat(x, tests = "KPSS", strategy = FALSE, stars = FALSE)
  r <- res[[1]]
  expect_named(r, c("kpss", "decision", "process"), ignore.order = TRUE)
  expect_false(is.na(r$kpss$level$c$stat))
})

test_that("urstat handles list input", {
  set.seed(3)
  xl <- list(a = rnorm(50), b = cumsum(rnorm(50)))
  res <- urstat(xl, tests = "ADF", strategy = FALSE, stars = FALSE)
  expect_length(res, 2L)
  expect_named(res, c("a", "b"))
})

test_that("urstat returns I(0) for white noise", {
  set.seed(4)
  x <- rnorm(80)
  res <- urstat(x, tests = "ADF", strategy = FALSE, stars = FALSE)
  expect_equal(res[[1]]$decision, "I(0)")
})

test_that(".urs_stars returns correct strings", {
  expect_equal(unitrootests:::.urs_stars(0.001), "***")
  expect_equal(unitrootests:::.urs_stars(0.03),  "**")
  expect_equal(unitrootests:::.urs_stars(0.08),  "*")
  expect_equal(unitrootests:::.urs_stars(0.20),  "")
  expect_equal(unitrootests:::.urs_stars(NA),    "")
})
