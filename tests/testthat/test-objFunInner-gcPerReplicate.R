## objFunInner() collects garbage once per replicate (gc(full = FALSE)) so a fit worker's peak memory stays lower
## (ELF 12.1_2: VmHWM 5.5 GB to 4.7 GB, same time per evaluation). It must not change the objective: a fixed toy
## year, simulated with the real spreadCpp(), gives the same fires and the same objective value.

toyYear <- function(yearSpreadSD = 0, returnSims = FALSE) {
  colsToUse <- spreadDesignCols(formC)
  p <- fireSenseUtils:::splitSpreadPar(c(maxAsymptote = 0.27, clim = 1, fuel = 1.5, youngAge = -1))$par
  set.seed(7)
  fireSenseUtils:::objFunInner(
    yr = "2004", annDTx1000 = data.table::copy(climFix$annualDTx1000[["2004"]]), par = p, parsModel = length(colsToUse),
    annualFires = data.table::as.data.table(climFix$historicalFires[["2004"]]),
    nonAnnualDTx1000 = data.table::copy(climFix$nonAnnualDTx1000),
    annualFireBufferedDT = dtC(ids = 1L, pixelID = 1:nPixC),
    indexNonAnnual = dtC(ind = 1L, date = "2000"), colsToUse = colsToUse, covMinMax = NULL,
    mutuallyExclusive = NULL, doAssertions = FALSE, maxFireSpread = 0.28, lowerSpreadProb = f,
    cells = numeric(nPixC), lanscape1stQuantileThresh = Inf, weighted = TRUE,
    r = terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20), Nreps = 5L,
    doSNLL_FSTest = TRUE, doMADTest = FALSE, doADTest = FALSE, plot.it = FALSE, verbose = 0,
    returnSims = returnSims, yearSpreadSD = yearSpreadSD)
}

test_that("objFunInner's simulated fires and objective are unchanged by the per-replicate gc", {
  sims <- toyYear(returnSims = TRUE)$sims
  expect_identical(NROW(sims), 5L)
  expect_equal(sum(sims$sim), 17)
  obj <- toyYear()
  expect_identical(obj$nSims, 5L)
  expect_false(obj$gated)
  expect_equal(obj$SNLL_FS, 261.22520981670976, tolerance = 1e-12)   # value from before the gc was added
})
