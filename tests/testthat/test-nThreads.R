test_that("nThreads() takes the option, else min(8, max(2, floor(cores - load)))", {
  withr::local_options(fireSense.nThreads = 4)
  expect_identical(nThreads(), 4L)
  withr::local_options(fireSense.nThreads = 0)
  expect_identical(nThreads(), 1L)
  withr::local_options(fireSense.nThreads = NULL)
  testthat::local_mocked_bindings(.loadAvg1 = function() 0)
  expect_identical(nThreads(cores = 80), 8L)
  testthat::local_mocked_bindings(.loadAvg1 = function() 79)
  expect_identical(nThreads(cores = 80), 2L)
  testthat::local_mocked_bindings(.loadAvg1 = function() 75)
  expect_identical(nThreads(cores = 80), 5L)
  testthat::local_mocked_bindings(.loadAvg1 = function() 0)
  expect_identical(nThreads(cores = 4), 4L)
})

test_that("nThreads() falls back to min(8, cores) when the load average cannot be read", {
  withr::local_options(fireSense.nThreads = NULL)
  testthat::local_mocked_bindings(.loadAvg1 = function() 0)
  expect_identical(nThreads(cores = 80), 8L)
  testthat::local_mocked_bindings(.loadAvg1 = function() NA_real_)
  expect_identical(nThreads(cores = 80), 8L)
  expect_identical(nThreads(cores = 4), 4L)
})

test_that(".loadAvg1() is a number or NA", {
  expect_true(is.numeric(.loadAvg1()) && length(.loadAvg1()) == 1L)
})

test_that("the global xgboost thread count is 1 after load, and setXgbThreads() follows the option", {
  skip_if_not_installed("xgboost")
  withr::defer(xgboost::xgb.set.config(nthread = 1L))
  ## .onLoad ran when the package was loaded for this test
  expect_identical(xgboost::xgb.get.config()$nthread, 1L)
  withr::local_options(fireSense.nThreads = 2)
  setXgbThreads()
  expect_identical(xgboost::xgb.get.config()$nthread, 2L)
  withr::local_options(fireSense.nThreads = NULL)
  setXgbThreads()
  expect_identical(xgboost::xgb.get.config()$nthread, 1L)

  x <- matrix(rnorm(200), ncol = 2)
  b <- xgboost::xgboost(x, rpois(100, 1), objective = "count:poisson", nrounds = 2, nthreads = 1)
  f <- withr::local_tempfile(fileext = ".rds")
  saveRDS(b, f)
  expect_equal(unname(predict(readRDS(f), x)), unname(predict(b, x)), tolerance = 1e-6)
})
