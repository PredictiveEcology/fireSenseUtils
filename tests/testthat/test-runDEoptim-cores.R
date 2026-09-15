## runDEoptim() sizes the DEoptim cluster and hands its settings to clusters' iterative runner.
## FireSense phase 2 (2026-09-15): the cluster request was hard-coded at 100 workers whatever the
## number of parameters, and NP had to be "around 10 x parameters and EXACTLY the number of
## cores allocated" (Eliot). NP now comes from the built cluster, in clusters::clusterSetup().

callRunDEoptim <- function(seen, npar = 12L, ...) {
  lower <- stats::setNames(rep(0, npar), paste0("p", seq_len(npar)))
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

mockCluster <- function(seen, builtWorkers = 40L) {
  testthat::local_mocked_bindings(
    clusterSetup = function(..., nCoresNeeded, NP) {
      seen$nCoresNeeded <- nCoresNeeded
      ## what clusterSetup() returns once a cluster of `builtWorkers` is running
      list(itermax = 5, trace = FALSE, strategy = 2L, NP = builtWorkers)
    },
    DEoptimIterative2 = function(fn, lower, upper, control, ...) {
      seen$control <- control
      list()
    },
    .package = "clusters", .env = parent.frame())
  testthat::local_mocked_bindings(termsInDEoptim = function(...) invisible(NULL), .env = parent.frame())
}

test_that("runDEoptim asks for about 10 workers per estimated parameter", {
  seen <- new.env()
  mockCluster(seen)
  callRunDEoptim(seen, npar = 12L)
  expect_identical(as.integer(seen$nCoresNeeded), 120L)
})

test_that("the cluster request can be overridden", {
  seen <- new.env()
  mockCluster(seen)
  callRunDEoptim(seen, npar = 12L, nCoresNeeded = 60L)
  expect_identical(as.integer(seen$nCoresNeeded), 60L)
})

test_that("DEoptim settings reach clusterSetup(), and .c is DEoptim's c, not an objective-function argument", {
  ## Eliot, 2026-09-15: "Any user passed args should pass into the DEoptim processes." `.c` was
  ## sent to the objective function (which ignores it), so DEoptim always used its default c.
  seen <- new.env()
  testthat::local_mocked_bindings(
    clusterSetup = function(..., nCoresNeeded, NP, controlArgs = NULL) {
      seen$controlArgs <- controlArgs
      list(itermax = 5, trace = FALSE, strategy = 6L, NP = 40L)
    },
    DEoptimIterative2 = function(fn, lower, upper, control, ...) {
      seen$dots <- names(list(...))
      list()
    },
    .package = "clusters")
  testthat::local_mocked_bindings(termsInDEoptim = function(...) invisible(NULL))
  callRunDEoptim(seen, npar = 12L, .c = 0.9, DEoptimControl = list(CR = 0.7, F = 0.6, p = 0.3))
  expect_equal(seen$controlArgs, list(c = 0.9, CR = 0.7, F = 0.6, p = 0.3))
  expect_false(".c" %in% seen$dots)
  ## a c given in DEoptimControl wins over .c
  callRunDEoptim(seen, npar = 12L, .c = 0.9, DEoptimControl = list(c = 0.2))
  expect_equal(seen$controlArgs$c, 0.2)
})

test_that("DEoptim gets the NP of the cluster that was built and the caller's strategy", {
  seen <- new.env()
  mockCluster(seen, builtWorkers = 57L)
  callRunDEoptim(seen, npar = 12L)
  expect_identical(as.integer(seen$control$NP), 57L)
  expect_identical(as.integer(seen$control$strategy), 2L)
})

test_that("runDEoptim runs where R has no OpenMP", {
  ## CI macOS, 2026-09-15: RhpcBLASctl::omp_get_max_threads() is NA there, and `if (origOmp > 1)`
  ## stopped runDEoptim with "missing value where TRUE/FALSE needed".
  seen <- new.env()
  mockCluster(seen)
  testthat::local_mocked_bindings(omp_get_max_threads = function() NA_integer_)
  expect_no_error(callRunDEoptim(seen, npar = 12L))
  expect_identical(as.integer(seen$nCoresNeeded), 120L)
})
