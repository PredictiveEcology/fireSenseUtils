## 2026-09-28: drawing the DEoptim progress figures after every generation took 8.3 s of each 53 s
## generation (16% of the wall time) in a FireSense fit. runDEoptim() now passes `plotEvery` to
## clusters::DEoptimIterative2(); how often figures are drawn must not change the fit's cache key.

callRunDEoptimPlotEvery <- function(cachePath, ...) {
  lower <- stats::setNames(c(0.25, 0.2, 0.1, 0), c("maxAsymptote", "hillSlope1", "inflectionPoint1", "x"))
  suppressMessages(runDEoptim(
    landscape = NULL, annualDTx1000 = NULL, nonAnnualDTx1000 = NULL,
    fireBufferedListDT = NULL, historicalFires = NULL,
    itermax = 5, trace = FALSE, strategy = 2L, cores = c("hostA", "hostB"),
    paths = list(cachePath = cachePath),
    lower = lower, upper = lower + 1, mutuallyExclusive = NULL, formulaToFit = NULL,
    objFunCoresInternal = 1L, covMinMax = NULL, maxFireSpread = 0.3, Nreps = 1L,
    ## the default logPath (and so visualizeDEoptim) is a new tempfile each call, which is in the key
    logPath = file.path(cachePath, "fit.log"), .verbose = FALSE, ...))
}

mockDEoptimPlotEvery <- function(seen, env = parent.frame()) {
  testthat::local_mocked_bindings(
    clusterSetup = function(...) list(itermax = 5, trace = FALSE, strategy = 2L, NP = 40L),
    DEoptimIterative2 = function(fn, lower, upper, control, ...) {
      seen$calls <- seen$calls + 1L
      seen$plotEvery <- list(...)$plotEvery
      list()
    },
    .package = "clusters", .env = env)
  testthat::local_mocked_bindings(termsInDEoptim = function(...) invisible(NULL), .env = env)
}

test_that("plotEvery reaches clusters::DEoptimIterative2(), 25 by default", {
  seen <- new.env(); seen$calls <- 0L
  mockDEoptimPlotEvery(seen)
  withr::local_options(reproducible.useCache = FALSE)
  callRunDEoptimPlotEvery(withr::local_tempdir())
  expect_identical(seen$plotEvery, 25L)
  callRunDEoptimPlotEvery(withr::local_tempdir(), plotEvery = 5L)
  expect_identical(seen$plotEvery, 5L)
})

test_that("changing plotEvery does not change the fit's cache key", {
  seen <- new.env(); seen$calls <- 0L
  mockDEoptimPlotEvery(seen)
  withr::local_options(reproducible.useCache = TRUE)
  cp <- withr::local_tempdir()
  callRunDEoptimPlotEvery(cp, plotEvery = 25L)
  expect_identical(seen$calls, 1L)
  callRunDEoptimPlotEvery(cp, plotEvery = 1L)       # a cache hit: the fit is not run again
  expect_identical(seen$calls, 1L)
  callRunDEoptimPlotEvery(cp, plotEvery = 1L, thresh = 551)  # control: a fit argument does refit
  expect_identical(seen$calls, 2L)
})
