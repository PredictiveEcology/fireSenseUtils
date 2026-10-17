## Shared by test-runDEoptim-plotEvery.R and test-runDEoptim-convergence.R: a runDEoptim() call with
## placeholder data, and a mock of the cluster set-up and the fit.

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
