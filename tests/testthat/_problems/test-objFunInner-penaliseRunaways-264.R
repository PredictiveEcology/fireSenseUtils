# Extracted from test-objFunInner-penaliseRunaways.R:264

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "fireSenseUtils", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
r60 <- terra::rast(nrows = 60, ncols = 60, xmin = 0, xmax = 60, ymin = 0, ymax = 60)
cellOf60 <- function(row, col) (row - 1L) * 60L + col
runawayFixture <- function(n = 10L) {
  n <- rep_len(n, 5L)
  tops <- c(2L, 2L, 22L, 22L, 42L); lefts <- c(2L, 22L, 2L, 22L, 2L)
  buf <- data.table::rbindlist(lapply(1:5, function(i)
    data.table::data.table(ids = i, pixelID = as.vector(outer(tops[i] + 0:(n[i] - 1L), lefts[i] + 0:(n[i] - 1L),
                                                              cellOf60)), buffer = 0L)))
  centre <- cellOf60(tops + n %/% 2L, lefts + n %/% 2L)
  list(buf = buf, fires = data.table::data.table(cells = centre, size = c(25L, 30L, 20L, 28L, 22L), ids = 1:5))
}
runawayObjective <- function(sp, seed = 1, fx = runawayFixture(), Nreps = 6L, ...) {
  pix <- fx$buf$pixelID
  local_mocked_bindings(
    paramsSeparate = function(...) list(logisticPars = c(0.27, 1, 1), covPars = 1),
    spreadProbFromIntegerCovs = function(...) data.table::data.table(pixelID = pix, cov = 0),
    logisticAll = function(...) seq(sp[1], sp[length(sp)], length.out = length(pix)), # a flat surface is refused as not spread out
    .package = "fireSenseUtils"
  )
  set.seed(seed)
  fireSenseUtils:::objFunInner(
    yr = "2004", annDTx1000 = NULL, par = 1, parsModel = NULL,
    annualFires = fx$fires, nonAnnualDTx1000 = NULL, annualFireBufferedDT = data.table::copy(fx$buf),
    indexNonAnnual = NULL, colsToUse = "cov", covMinMax = NULL,
    mutuallyExclusive = NULL, doAssertions = FALSE, maxFireSpread = 1,
    lowerSpreadProb = 0.05, cells = numeric(3600), lanscape1stQuantileThresh = 1,
    weighted = FALSE, r = r60, Nreps = Nreps,
    doSNLL_FSTest = TRUE, doMADTest = TRUE, doADTest = TRUE, doYearArea = TRUE,
    plot.it = FALSE, verbose = 0, ...
  )
}
runawayEdgeRuns <- function(nRing, touch, ...) {
  pix <- seq_len(nRing + 10L)
  fx <- list(buf = data.table::data.table(ids = 1L, pixelID = pix, buffer = 0L, edge = pix <= nRing),
             fires = data.table::data.table(cells = nRing + 5L, size = 20L, ids = 1L))
  local_mocked_bindings(
    cropToCells = function(r, cells, ...) list(r = r, ncell = terra::ncell(r), toCrop = identity, toFull = identity),
    .package = "fireSenseUtils")
  local_mocked_bindings(
    spreadCpp = function(landscape, loci, ...)
      data.table::data.table(initialLocus = loci, indices = c(loci, seq_len(touch))),
    .package = "SpaDES.tools")
  runawayObjective(c(0.06, 0.08), fx = fx, ...)$runaways
}

# test -------------------------------------------------------------------------
skip_if_not_installed("SpaDES.tools")
expect_equal(runawayEdgeRuns(100L, 0L), 0)
expect_equal(runawayEdgeRuns(100L, 1L), 0)
expect_equal(runawayEdgeRuns(100L, 2L), 0)
expect_equal(runawayEdgeRuns(100L, 3L), 6)
