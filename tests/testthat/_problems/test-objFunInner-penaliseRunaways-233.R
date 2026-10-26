# Extracted from test-objFunInner-penaliseRunaways.R:233

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
mk <- function(snll, gated = FALSE, thresh = Inf, ...) {
    local_mocked_bindings(
      objFunInner = function(yr, ...) list(SNLL_FS = snll[[yr]], runaways = 0L, nSims = 0L, gated = gated[[yr]]),
      .package = "fireSenseUtils")
    land <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20, vals = 1)
    yrs <- c("2001", "2002", "2003")
    dt <- function(...) data.table::data.table(...)
    fireSenseUtils::.objfunSpreadFit(
      par = c(0.27, 1, 1, 1), landscape = land,
      annualDTx1000 = stats::setNames(lapply(yrs, function(y) dt(pixelID = 1:2, cov = 100L)), yrs),
      nonAnnualDTx1000 = list(`2001_2003` = dt(pixelID = 1:2)),
      formulaToFit = "~ 0 + cov",
      historicalFires = list(`2001` = data.frame(size = c(100, 200), cells = 1:2, ids = 1:2),
                             `2002` = data.frame(size = c(150, 250), cells = 3:4, ids = 1:2),
                             `2003` = data.frame(size = c(5, 6), cells = 5:6, ids = 1:2)),
      fireBufferedListDT = stats::setNames(lapply(yrs, function(y) dt(pixelID = 1:2, ids = 1L)), yrs),
      tests = "snll_fs", Nreps = 1L, doAssertions = FALSE, verbose = 0, thresh = thresh,
      returnTerms = TRUE, ...)
  }
snll <- c("2001" = 100, "2002" = 300, "2003" = 40)
ok <- mk(snll, gated = c("2001" = FALSE, "2002" = FALSE, "2003" = FALSE))
expect_equal(ok[["firstBlockSNLL"]], (100 + 300) / 2)
