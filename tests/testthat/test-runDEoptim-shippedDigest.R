## 2026-09-29 gave the DEoptim fit a key on the data shipped to the workers (clusters::shippedObjectsDigest()),
## because the data reach the workers through clusterSetup(objsNeeded), not as arguments of the cached call. The
## re-score, the coefficient profile and the simulated members run the objective on the same shipped data and had
## no such key: a new fit on other data reused the old data's re-score.

shippedDigestRun <- function(seen, cachePath, digest, ...) {
  finalPop <- matrix(c(0.26, 1, 1, 0.5, 0.26, 1, 1, 0.7), nrow = 2, byrow = TRUE)
  testthat::local_mocked_bindings(
    clusterSetup = function(...) structure(list(itermax = 5, trace = FALSE, strategy = 2L, NP = 40L, cluster = NULL),
                                           objsDigest = digest),
    DEoptimIterative = function(...) list(list(member = list(pop = finalPop))),
    .package = "clusters", .env = parent.frame())
  testthat::local_mocked_bindings(
    termsInDEoptim = function(...) invisible(NULL),
    rescorePopulation = function(pop, fn, reps, cl, seed = 1L, fnArgs = list()) {
      seen$rescore <- seen$rescore + 1L
      data.table::data.table(member = 1:2, rep = 1L, value = c(2, 1))
    },
    profileCoefficients = function(best, pop, fn, reps, cl, fnArgs, ...) { seen$profile <- seen$profile + 1L; "profile" },
    simulateFireSizes = function(pop, fn, fnArgs, cl, ...) { seen$sims <- seen$sims + 1L; "sims" },
    .env = parent.frame())
  lower <- stats::setNames(c(0.25, 0.2, 0.1, 0), c("maxAsymptote", "hillSlope1", "inflectionPoint1", "x"))
  suppressMessages(runDEoptim(
    landscape = NULL, annualDTx1000 = NULL, nonAnnualDTx1000 = NULL,
    fireBufferedListDT = NULL, historicalFires = NULL,
    itermax = 5, trace = FALSE, strategy = 2L, cores = c("hostA", "hostB"),
    paths = list(cachePath = cachePath),
    lower = lower, upper = lower + 1, mutuallyExclusive = NULL, formulaToFit = NULL,
    objFunCoresInternal = 1L, covMinMax = NULL, maxFireSpread = 0.3, Nreps = 1L,
    logPath = file.path(cachePath, "fit.log"), .verbose = FALSE, profileReps = 2L, simulateMembers = 1L, ...))
}

test_that("the re-score, profile and simulated members are not shared by fits on other shipped data", {
  seen <- new.env(); seen$rescore <- seen$profile <- seen$sims <- 0L
  withr::local_options(reproducible.useCache = TRUE)
  cp <- withr::local_tempdir()
  shippedDigestRun(seen, cp, "digest-of-data-1")
  expect_identical(c(seen$rescore, seen$profile, seen$sims), c(1L, 1L, 1L))
  shippedDigestRun(seen, cp, "digest-of-data-2")           # only the shipped data differ
  expect_identical(c(seen$rescore, seen$profile, seen$sims), c(2L, 2L, 2L))
  shippedDigestRun(seen, cp, "digest-of-data-1")           # the first data again: cache hits
  expect_identical(c(seen$rescore, seen$profile, seen$sims), c(2L, 2L, 2L))
})

test_that("the whole-landscape tables reach the workers and the re-score only when given", {
  seen <- new.env(); seen$rescore <- seen$profile <- seen$sims <- 0L
  seen$objs <- list()
  withr::local_options(reproducible.useCache = FALSE)
  testthat::local_mocked_bindings(
    clusterSetup = function(..., objsNeeded) { seen$objs <- c(seen$objs, list(unlist(objsNeeded)))
      structure(list(itermax = 5, trace = FALSE, strategy = 2L, NP = 40L, cluster = NULL), objsDigest = "d") },
    DEoptimIterative = function(..., cutMaxTries = "unset") { seen$fitTries <- cutMaxTries
      list(list(member = list(pop = matrix(c(0.26, 1, 1, 0.5), nrow = 1)))) },
    .package = "clusters")
  testthat::local_mocked_bindings(
    termsInDEoptim = function(...) invisible(NULL),
    rescorePopulation = function(pop, fn, reps, cl, seed = 1L, fnArgs = list()) {
      seen$rescoreArgs <- fnArgs; data.table::data.table(member = 1L, rep = 1L, value = 1) })
  lower <- stats::setNames(c(0.25, 0.2, 0.1, 0), c("maxAsymptote", "hillSlope1", "inflectionPoint1", "x"))
  run <- function(...) suppressMessages(runDEoptim(
    landscape = NULL, annualDTx1000 = NULL, nonAnnualDTx1000 = NULL, fireBufferedListDT = NULL, historicalFires = NULL,
    itermax = 5, trace = FALSE, strategy = 2L, cores = c("hostA", "hostB"), paths = list(cachePath = withr::local_tempdir()),
    lower = lower, upper = lower + 1, mutuallyExclusive = NULL, formulaToFit = NULL, objFunCoresInternal = 1L,
    covMinMax = NULL, maxFireSpread = 0.3, Nreps = 1L, .verbose = FALSE, ...))
  run()
  expect_false(any(c("annualDTx1000Whole", "nonAnnualDTx1000Whole") %in% seen$objs[[1]]))
  expect_null(seen$fitTries)                       # NULL (not the default 50) reaches the fit, and is omitted from its key
  expect_false("cutMaxTries" %in% names(seen$rescoreArgs))
  run(annualDTx1000Whole = list(), nonAnnualDTx1000Whole = list(), penaliseRunaways = FALSE, cutMaxTries = 7L)
  expect_true(all(c("annualDTx1000Whole", "nonAnnualDTx1000Whole") %in% seen$objs[[2]]))
  expect_identical(seen$fitTries, 7L)
  expect_identical(seen$rescoreArgs$cutMaxTries, 7L)
  expect_error(run(annualDTx1000Whole = list(), nonAnnualDTx1000Whole = list()), "penaliseRunaways = FALSE")
})
