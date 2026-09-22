## runDEoptim() took no likelihood options, so every fit -- and the re-score of its final population,
## which picks the parameter sets written to the ledger -- ran .objfunSpreadFit() with its defaults
## (weighted = TRUE, i.e. a log(size) weight; sizeLik = "kde"), whatever the caller wanted. The
## cross-validated comparison of 2026-09-21 chose sizeLik = "t" with no weight; the module had no way
## to ask for it.

likArgs <- list(sizeLik = "t", sizeLikDf = 7, weighted = FALSE, adWeight = 12)

callRunDEoptimLik <- function(...) {
  lower <- stats::setNames(c(0.25, 0.2, 0.1, 0), c("maxAsymptote", "hillSlope1", "inflectionPoint1", "x"))
  withr::local_options(reproducible.useCache = FALSE)
  suppressMessages(runDEoptim(
    landscape = NULL, annualDTx1000 = NULL, nonAnnualDTx1000 = NULL,
    fireBufferedListDT = NULL, historicalFires = NULL,
    itermax = 5, trace = FALSE, strategy = 2L, cores = c("hostA", "hostB"),
    paths = list(cachePath = withr::local_tempdir()),
    lower = lower, upper = lower + 1, mutuallyExclusive = NULL, formulaToFit = NULL,
    objFunCoresInternal = 1L, covMinMax = NULL, maxFireSpread = 0.3, Nreps = 1L,
    .verbose = FALSE, ...))
}

test_that("the likelihood options reach the objective during the fit", {
  seen <- new.env()
  testthat::local_mocked_bindings(
    clusterSetup = function(...) list(itermax = 5, trace = FALSE, strategy = 2L, NP = 40L),
    DEoptimIterative2 = function(fn, lower, upper, control, ...) {
      seen$dots <- list(...)
      list()
    },
    .package = "clusters")
  testthat::local_mocked_bindings(termsInDEoptim = function(...) invisible(NULL))

  do.call(callRunDEoptimLik, likArgs)
  expect_identical(seen$dots[names(likArgs)], likArgs)
})

test_that("the same options reach the re-score of the final population", {
  ## The re-score chooses the ledger's parameter sets; scoring them with a different objective from
  ## the one they were fitted under would pick the wrong members.
  seen <- new.env()
  finalPop <- matrix(c(0.26, 1, 1, 0.5), nrow = 1)
  testthat::local_mocked_bindings(
    clusterSetup = function(...) list(itermax = 5, trace = FALSE, strategy = 2L, NP = 40L,
                                      cluster = NULL),
    DEoptimIterative2 = function(fn, lower, upper, control, ...) list(list(member = list(pop = finalPop))),
    .package = "clusters")
  testthat::local_mocked_bindings(
    termsInDEoptim = function(...) invisible(NULL),
    rescorePopulation = function(pop, fn, reps, cl, seed = 1L, fnArgs = list()) {
      seen$fnArgs <- fnArgs
      data.table::data.table(member = 1L, rep = 1L, value = 1)
    })

  do.call(callRunDEoptimLik, likArgs)
  expect_identical(seen$fnArgs[names(likArgs)], likArgs)
})

test_that("runDEoptim's defaults are the objective's defaults", {
  ## so a caller who sets nothing gets exactly what .objfunSpreadFit() would have used on its own
  f <- formals(runDEoptim)[names(likArgs)]
  o <- formals(.objfunSpreadFit)[names(likArgs)]
  expect_identical(f, o)
})
