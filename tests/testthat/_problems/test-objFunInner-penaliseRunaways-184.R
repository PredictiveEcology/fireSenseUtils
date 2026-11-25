# Extracted from test-objFunInner-penaliseRunaways.R:184

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

# test -------------------------------------------------------------------------
seen <- list()
local_mocked_bindings(
    objFunInner = function(runawaySize, runawayBufferMultiple = "absent", censorRunaways, ...) {
      seen[[length(seen) + 1L]] <<- list(size = runawaySize, mult = runawayBufferMultiple)
      list(SNLL_FS = 0)
    },
    .package = "fireSenseUtils"
  )
land <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20, vals = 1)
run <- function(...) {
    seen <<- list()
    fireSenseUtils::.objfunSpreadFit(
      par = c(0.27, 1, 1, 1), landscape = land,
      annualDTx1000 = list(year2004 = data.table::data.table(pixelID = 1:400, cov = 0L)),
      nonAnnualDTx1000 = list(`year2004` = data.table::data.table(pixelID = 1:400)),
      formulaToFit = "~ 0 + cov",
      historicalFires = list(year2004 = data.frame(size = 40L, cells = 1L, ids = 1L)),
      fireBufferedListDT = list(year2004 = data.table::data.table(ids = 1L, pixelID = 1:400)),
      indexNonAnnual = data.table::data.table(ind = 1L, date = "2004"), doAssertions = FALSE, ...
    )
    seen
  }
