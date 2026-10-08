## 2026-09-28: drawing the DEoptim progress figures after every generation took 8.3 s of each 53 s
## generation (16% of the wall time) in a FireSense fit. runDEoptim() now passes `plotEvery` to
## clusters::DEoptimIterative(); how often figures are drawn must not change the fit's cache key.

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
    clusterSetup = function(...) structure(list(itermax = 5, trace = FALSE, strategy = 2L, NP = 40L),
                                           objsDigest = seen$dataDigest),
    DEoptimIterative = function(fn, lower, upper, control, ...) {
      seen$calls <- seen$calls + 1L
      seen$plotEvery <- list(...)$plotEvery
      list()
    },
    .package = "clusters", .env = env)
  testthat::local_mocked_bindings(termsInDEoptim = function(...) invisible(NULL), .env = env)
}

test_that("plotEvery reaches clusters::DEoptimIterative(), 25 by default", {
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


## 2026-09-29: the data reaches the workers through clusterSetup(objsNeeded), not as an argument of the
## cached call, so the two held-out folds of an ELF (same settings, different years) shared one fit.
test_that("two fits that differ only in the shipped data do not share the fit's cache", {
  seen <- new.env(); seen$calls <- 0L
  mockDEoptimPlotEvery(seen)
  withr::local_options(reproducible.useCache = TRUE)
  cp <- withr::local_tempdir()
  seen$dataDigest <- "digest-of-fold-1"
  callRunDEoptimPlotEvery(cp)
  expect_identical(seen$calls, 1L)
  seen$dataDigest <- "digest-of-fold-2"            # only the shipped data differ
  callRunDEoptimPlotEvery(cp)
  expect_identical(seen$calls, 2L)
  seen$dataDigest <- "digest-of-fold-1"            # the same data again: a cache hit
  callRunDEoptimPlotEvery(cp)
  expect_identical(seen$calls, 2L)
})

## The progress file (one row per generation) goes where the caller says, and is not part of the fit.
test_that("progressFile reaches clusters::DEoptimIterative(), NULL by default", {
  seen <- new.env(); seen$calls <- 0L
  testthat::local_mocked_bindings(
    clusterSetup = function(...) list(itermax = 5, trace = FALSE, strategy = 2L, NP = 40L),
    DEoptimIterative = function(fn, lower, upper, control, ..., progressFile = "unset") {
      seen$progressFile <- progressFile
      list()
    },
    .package = "clusters")
  testthat::local_mocked_bindings(termsInDEoptim = function(...) invisible(NULL))
  withr::local_options(reproducible.useCache = FALSE)
  callRunDEoptimPlotEvery(withr::local_tempdir(), progressFile = "out/DEoptimProgress_a.csv")
  expect_identical(seen$progressFile, "out/DEoptimProgress_a.csv")
  callRunDEoptimPlotEvery(withr::local_tempdir())
  expect_null(seen$progressFile)
})

test_that("changing progressFile does not change the fit's cache key", {
  seen <- new.env(); seen$calls <- 0L
  mockDEoptimPlotEvery(seen)
  withr::local_options(reproducible.useCache = TRUE)
  cp <- withr::local_tempdir()
  callRunDEoptimPlotEvery(cp, progressFile = "a.csv")
  callRunDEoptimPlotEvery(cp, progressFile = "b.csv")  # a cache hit
  expect_identical(seen$calls, 1L)
})
