## spreadProbGates() applies the objective's own spreadProb gates ("Too burny a landscape", "Not spread
## out enough", median out of range) to a parameter set without simulating a fire. The objective calls the
## same spreadProbGateTest(), so a calibration that screens draws with it screens them as the objective will.

set.seed(1)
nPix <- 400L
dtx <- function(...) data.table::data.table(...)
gateFix <- list(
  annualDTx1000 = list(`2004` = dtx(pixelID = 1:nPix, cov1 = sample(0:1000, nPix, TRUE)),
                       `2005` = dtx(pixelID = 1:nPix, cov1 = sample(0:1000, nPix, TRUE)),
                       `2006` = dtx(pixelID = 1:nPix, cov1 = sample(0:1000, nPix, TRUE))),
  nonAnnualDTx1000 = list(`2000` = dtx(pixelID = 1:nPix, cov2 = sample(0:1000, nPix, TRUE))),
  historicalFires = list(`2004` = data.frame(size = 50, cells = 1L, ids = 1L),
                         `2005` = data.frame(size = 80, cells = 2L, ids = 1L),
                         `2006` = data.frame(size = 5, cells = 3L, ids = 1L)),
  formulaToFit = "~ 0 + cov1 + cov2")
drawPar <- function() c(maxAsymptote = runif(1, 0.2, 0.3),
                        cov1 = runif(1, -50, 50), cov2 = runif(1, -50, 50))

## did the objective refuse the year? objFunInner() reaches spread() only if the gates pass
objectiveRefused <- function(par, yr) {
  reached <- FALSE
  local_mocked_bindings(
    spreadCpp = function(...) {
      reached <<- TRUE
      a <- list(...)
      data.table::data.table(initialLocus = a$loci, indices = a$loci, id = seq_along(a$loci), active = FALSE)
    }, .package = "SpaDES.tools")
  p <- fireSenseUtils:::splitSpreadPar(par)$par
  ps <- fireSenseUtils:::paramsSeparate(p, 2L)
  fireSenseUtils:::objFunInner(
    yr = yr, annDTx1000 = data.table::copy(gateFix$annualDTx1000[[yr]]), par = p, parsModel = 2L,
    annualFires = data.table::as.data.table(gateFix$historicalFires[[yr]]),
    nonAnnualDTx1000 = data.table::copy(gateFix$nonAnnualDTx1000),
    annualFireBufferedDT = dtx(ids = 1L, pixelID = 1:nPix),
    indexNonAnnual = dtx(ind = 1L, date = "2000"), colsToUse = c("cov1", "cov2"), covMinMax = NULL,
    mutuallyExclusive = NULL, doAssertions = FALSE, maxFireSpread = 0.28, lowerSpreadProb = 0.13,
    cells = numeric(nPix), lanscape1stQuantileThresh = 0.265, weighted = TRUE,
    r = terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20), Nreps = 1L,
    doSNLL_FSTest = TRUE, doMADTest = FALSE, doADTest = FALSE, plot.it = FALSE, verbose = 0,
    returnSims = TRUE)
  !reached
}

test_that("the objective's gate decision and spreadProbGates() agree, parameter set by parameter set", {
  pars <- replicate(60, drawPar(), simplify = FALSE)
  g <- do.call(spreadProbGates, c(list(par = pars), gateFix))
  expect_s3_class(g, "data.frame")
  expect_setequal(c("pass", "burny", "notSpread", "medianOK"), names(g))
  expect_true(any(g$pass) && any(!g$pass))   # the fixture exercises both outcomes
  yrs <- c("2004", "2005")                   # the two largest fire years: the first block
  for (i in seq_along(pars)) {
    refused <- vapply(yrs, function(y) objectiveRefused(pars[[i]], y), logical(1))
    expect_identical(g$pass[i], !any(refused), info = paste("draw", i))
  }
  ## a single set gives the same answer, with the years it tested
  one <- do.call(spreadProbGates, c(list(par = pars[[3]]), gateFix))
  expect_identical(one$pass, g$pass[3])
  expect_setequal(one$years, yrs)
})
