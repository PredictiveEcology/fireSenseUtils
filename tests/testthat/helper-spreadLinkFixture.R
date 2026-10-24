## Shared fixture of the spread-link tests (test-climateOutside.R, test-powerLaw.R): a toy landscape with two
## climate columns, one fuel column and youngAge, and objFunInner() for one year with spreadCpp() mocked.

dtC <- function(...) data.table::data.table(...)
set.seed(21)
nPixC <- 400L
## two climate columns (annual), one fuel column and youngAge (non-annual); 2005 is the hot year
climFix <- list(
  annualDTx1000 = list(`2004` = dtC(pixelID = 1:nPixC, clim = sample(200:500, nPixC, TRUE), clim2 = sample(0:1000, nPixC, TRUE)),
                       `2005` = dtC(pixelID = 1:nPixC, clim = sample(600:1000, nPixC, TRUE), clim2 = sample(0:1000, nPixC, TRUE))),
  nonAnnualDTx1000 = list(`2000` = dtC(pixelID = 1:nPixC, fuel = sample(0:1000, nPixC, TRUE),
                                       youngAge = as.integer(runif(nPixC) < 0.1) * 1000L)),
  historicalFires = list(`2004` = data.frame(size = 50, cells = 1L, ids = 1L),
                         `2005` = data.frame(size = 80, cells = 2L, ids = 1L)))
formC <- "~ 0 + clim + fuel + youngAge"
## a design matrix with named columns, as spreadProbCovsOfYear() gives it
matOf <- function(yr, cols = c("clim", "fuel", "youngAge")) {
  d <- climFix$nonAnnualDTx1000[[1]][climFix$annualDTx1000[[yr]], on = "pixelID"]
  m <- as.matrix(d[, cols, with = FALSE]) / 1000
  colnames(m) <- cols
  m
}
lp3 <- c(0.27, 1, 1) # maxAsymptote, hillSlope1, inflectionPoint1
f <- 0.13

## objFunInner() for one year, the real covariate pipeline and link; spreadCpp() hands back the spreadProb it got
spreadProbSeenC <- function(par, yr, climateRef = NULL, yearSpreadSD = 0, lanscape1stQuantileThresh = Inf,
                            maxFireSpread = 0.28, ...) {
  seen <- new.env()
  seen$reached <- FALSE
  local_mocked_bindings(
    spreadCpp = function(...) {
      a <- list(...)
      seen$sp <- a$spreadProb; seen$reached <- TRUE
      dtC(initialLocus = a$loci, indices = a$loci, id = seq_along(a$loci), active = FALSE)
    }, .package = "SpaDES.tools", .env = parent.frame())
  colsToUse <- spreadDesignCols(formC)
  p <- fireSenseUtils:::splitSpreadPar(par)$par
  out <- fireSenseUtils:::objFunInner(
    yr = yr, annDTx1000 = data.table::copy(climFix$annualDTx1000[[yr]]), par = p, parsModel = length(colsToUse),
    annualFires = data.table::as.data.table(climFix$historicalFires[[yr]]),
    nonAnnualDTx1000 = data.table::copy(climFix$nonAnnualDTx1000),
    annualFireBufferedDT = dtC(ids = 1L, pixelID = 1:nPixC),
    indexNonAnnual = dtC(ind = 1L, date = "2000"), colsToUse = colsToUse, covMinMax = NULL,
    mutuallyExclusive = NULL, doAssertions = FALSE, maxFireSpread = maxFireSpread, lowerSpreadProb = f,
    cells = numeric(nPixC), lanscape1stQuantileThresh = lanscape1stQuantileThresh, weighted = TRUE,
    r = terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20), Nreps = 1L,
    doSNLL_FSTest = TRUE, doMADTest = FALSE, plot.it = FALSE, verbose = 0,
    returnSims = TRUE, climateRef = climateRef, yearSpreadSD = yearSpreadSD, ...)
  list(sp = seen$sp, reached = seen$reached, pSummary = out$pSummary)
}

