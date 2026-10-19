## Climate outside the maxAsymptote ceiling (spreadProbFromCovs(climateRef =)). Under the old link climate is a
## covariate inside the logistic, so where fuel puts a pixel near the ceiling a dry year adds almost nothing
## (2026-10 hindcasts: 59-77% of pixels within 0.01 of the ceiling every year). Under the new link
## p0 = f + (maxAsymptote - f) * plogis(x_fuel) and p = plogis(qlogis(p0) + beta_clim * (climate - climateRef)),
## uncapped; the year effect eps_Y is added on the same logit afterwards. Off (climateRef = NULL) is bit-identical.

test_that("climateRef = NULL is logisticAll(), bit for bit, for both logistic forms", {
  mat <- matOf("2004")
  b <- c(clim = 2, fuel = 3, youngAge = -1)
  expect_identical(spreadProbFromCovs(mat, lp3, b, f), logisticAll(lp3, mat, b, f))
  expect_identical(spreadProbFromCovs(mat, lp3, b, f, climateRef = NULL), logisticAll(lp3, mat, b, f))
  expect_identical(spreadProbFromCovs(mat, lp3, b, f, climateRef = list()), logisticAll(lp3, mat, b, f))
  up <- c(0.27, 1, 1, -0.3)
  expect_identical(spreadProbFromCovs(mat, up, b, f, link = "logistic3pUpper"),
                   logisticAll(up, mat, b, f, link = "logistic3pUpper"))
  expect_false(formals(spreadProbFromCovs)$climateAtRef)
  expect_null(formals(spreadProbFromCovs)$climateRef)
})

test_that("with climateRef, climate acts after the fuel ceiling: p0 at the reference, uncapped above it", {
  mat <- matOf("2005")
  b <- c(clim = 4, fuel = 3, youngAge = -1)
  ref <- list(clim = 0.3)
  p <- spreadProbFromCovs(mat, lp3, b, f, climateRef = ref)
  ## by hand
  p0 <- f + (lp3[1] - f) * plogis(mat[, "fuel"] * 3 + mat[, "youngAge"] * -1)
  expect_equal(as.numeric(p), plogis(qlogis(p0) + 4 * (mat[, "clim"] - 0.3)))
  ## climateAtRef gives p0, the part the ceiling applies to
  expect_equal(as.numeric(spreadProbFromCovs(mat, lp3, b, f, climateRef = ref, climateAtRef = TRUE)), unname(p0))
  expect_true(all(p0 <= lp3[1]))
  ## the hot year: most pixels are above maxAsymptote -- the ceiling no longer binds for climate
  expect_gt(mean(p > lp3[1]), 0.9)
  expect_gt(max(p), lp3[1] + 0.1)
  ## at climate = climateRef the new link is p0 exactly
  matRef <- mat; matRef[, "clim"] <- 0.3
  expect_equal(as.numeric(spreadProbFromCovs(matRef, lp3, b, f, climateRef = ref)), unname(p0))
  ## a higher climate value raises p, monotonically, for every pixel
  matHi <- mat; matHi[, "clim"] <- mat[, "clim"] + 0.1
  expect_true(all(spreadProbFromCovs(matHi, lp3, b, f, climateRef = ref) > p))
  ## the same shape as logisticAll(): one column, one row per pixel
  expect_identical(dim(p), c(nPixC, 1L))
})

test_that("several climate columns, column order, an intercept, and the names are checked", {
  mat <- cbind(`(Intercept)` = 1, matOf("2004", c("fuel", "clim", "clim2", "youngAge")))
  b <- c(`(Intercept)` = -0.5, fuel = 3, clim = 2, clim2 = 1, youngAge = -1)
  ref <- c(clim2 = 0.5, clim = 0.4) # not in the matrix's order
  p <- spreadProbFromCovs(mat, lp3, b, f, climateRef = ref)
  p0 <- f + (lp3[1] - f) * plogis(-0.5 + 3 * mat[, "fuel"] - mat[, "youngAge"])
  expect_equal(as.numeric(p), plogis(qlogis(p0) + 2 * (mat[, "clim"] - 0.4) + 1 * (mat[, "clim2"] - 0.5)))
  expect_error(spreadProbFromCovs(mat, lp3, b, f, climateRef = c(CMD = 0.4)), "must be a column of `mat`")
  expect_error(spreadProbFromCovs(unname(mat), lp3, b, f, climateRef = ref), "must be a column of `mat`")
})

test_that("spreadClimateCols() and spreadClimateRef(): the annual non-indicator covariates and their design-scale means", {
  expect_identical(spreadClimateCols(formC, climFix$annualDTx1000), "clim")
  expect_identical(spreadClimateCols("~ 1 + clim2 + clim + youngAge + fuel", climFix$annualDTx1000), c("clim2", "clim"))
  expect_identical(spreadClimateCols("~ 0 + fuel + youngAge", climFix$annualDTx1000), character(0))
  ## youngAge annual: an indicator, not climate
  ann <- lapply(climFix$annualDTx1000, function(d) dtC(d, youngAge = 0L))
  expect_identical(spreadClimateCols("~ 0 + clim + youngAge", ann), "clim")
  expect_null(spreadClimateRef(climFix$annualDTx1000, climFix$nonAnnualDTx1000, "~ 0 + fuel + youngAge"))
  ## the mean over every pixel-year, rescaled with covMinMax
  mm <- dtC(clim = c(0, 2), clim2 = c(0, 1), fuel = c(0, 1), youngAge = c(0, 1))
  ref <- spreadClimateRef(climFix$annualDTx1000, climFix$nonAnnualDTx1000, formC, covMinMax = mm)
  expect_named(ref, "clim")
  expect_equal(ref$clim, mean(c(climFix$annualDTx1000[[1]]$clim, climFix$annualDTx1000[[2]]$clim) / 1000 / 2))
  ## with an intercept the reference is on the centred scale: the centre itself, so 0
  formI <- "~ 1 + clim + fuel + youngAge"
  centre <- spreadCovCentre(climFix$annualDTx1000, climFix$nonAnnualDTx1000, formI, covMinMax = mm, mutuallyExclusive = NULL)
  refI <- spreadClimateRef(climFix$annualDTx1000, climFix$nonAnnualDTx1000, formI, covMinMax = mm,
                           mutuallyExclusive = NULL, covCentre = centre)
  expect_equal(refI$clim, 0)
  expect_equal(centre$clim, ref$clim)
})

test_that("the objective's spreadProb under the new link is the hand-built one, with the year effect on the same logit", {
  skip_if_not_installed("SpaDES.tools")
  par <- c(maxAsymptote = 0.27, clim = 4, fuel = 3, youngAge = -1)
  ref <- list(clim = 0.3)
  mat <- matOf("2005")
  p0 <- f + (0.27 - f) * plogis(3 * mat[, "fuel"] - mat[, "youngAge"])
  expected <- plogis(qlogis(p0) + 4 * (mat[, "clim"] - 0.3))
  got <- spreadProbSeenC(par, "2005", climateRef = ref)
  expect_true(got$reached)
  expect_equal(got$sp[-2L], unname(expected[-2L]))   # cell 2 is the ignition
  expect_equal(got$sp[2L], 1)
  ## without climateRef: the old link, climate inside the logistic, which the ceiling caps (a smaller climate
  ## coefficient, or the hot year is so uniformly at the ceiling that the objective's gate refuses it)
  parOld <- c(maxAsymptote = 0.27, clim = 1, fuel = 3, youngAge = -1)
  old <- spreadProbSeenC(parOld, "2005")
  pOld <- as.numeric(f + (0.27 - f) * plogis(mat %*% c(1, 3, -1)))   # every pixel, the ignition cell included
  expect_equal(old$sp[-2L], pOld[-2L])
  expect_true(all(old$sp[-2L] <= 0.27))
  expect_gt(mean(got$sp[-2L] > 0.27), 0.9)
  ## the year effect eps_Y is one shift of the same logit, after climate
  set.seed(3)
  withEps <- spreadProbSeenC(par, "2005", climateRef = ref, yearSpreadSD = 1)
  eps <- qlogis(withEps$sp[-2L]) - qlogis(got$sp[-2L])
  expect_lt(diff(range(eps)), 1e-9)
  expect_gt(abs(eps[1]), 1e-6)
  ## the saturation summary describes p0, the part the ceiling applies to
  expect_equal(got$pSummary$atCeiling, sum(p0 >= 0.27 - 0.002))
  expect_equal(old$pSummary$atCeiling, sum(pOld >= 0.27 - 0.002))
})

test_that("the gates look at p0: a hot year whose p is past maxFireSpread is still simulated", {
  skip_if_not_installed("SpaDES.tools")
  par <- c(maxAsymptote = 0.27, clim = 4, fuel = 3, youngAge = -1)
  ref <- list(clim = 0.3)
  mat <- matOf("2005")
  ## under the new link the median p of the hot year is above maxFireSpread ...
  pNew <- plogis(qlogis(f + (0.27 - f) * plogis(3 * mat[, "fuel"] - mat[, "youngAge"])) + 4 * (mat[, "clim"] - 0.3))
  expect_gt(median(pNew), 0.28)
  ## ... and the objective simulates it anyway, because the gate is on p0
  got <- spreadProbSeenC(par, "2005", climateRef = ref, lanscape1stQuantileThresh = 0.265)
  expect_true(got$reached)
  ## spreadProbGates() agrees with the objective, with and without climateRef, year by year
  gatesOf <- function(climateRef) do.call(spreadProbGates, c(list(
    par = par, formulaToFit = formC, climateRef = climateRef, years = c("2004", "2005"), mutuallyExclusive = NULL,
    maxFireSpread = 0.28, lowerSpreadProb = f), climFix[c("annualDTx1000", "nonAnnualDTx1000", "historicalFires")]))
  expect_true(gatesOf(ref)$pass)
  for (yr in c("2004", "2005"))
    expect_true(spreadProbSeenC(par, yr, climateRef = ref, lanscape1stQuantileThresh = 0.265)$reached)
  ## a set that fails on p0 fails for the gates and the objective alike
  burny <- c(maxAsymptote = 0.27, clim = 4, fuel = 30, youngAge = -1)
  gb2 <- do.call(spreadProbGates, c(list(
    par = burny, formulaToFit = formC, climateRef = ref, years = "2005", mutuallyExclusive = NULL,
    maxFireSpread = 0.28, lowerSpreadProb = f), climFix[c("annualDTx1000", "nonAnnualDTx1000", "historicalFires")]))
  expect_identical(gb2$pass, spreadProbSeenC(burny, "2005", climateRef = ref, lanscape1stQuantileThresh = 0.265)$reached)
  expect_false(gb2$pass)
})

test_that(".objfunSpreadFit() with climateRef = NULL is the objective as it was, bit for bit", {
  skip_if_not_installed("SpaDES.tools")
  r <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20, vals = 1)
  fb <- lapply(climFix$historicalFires, function(h) dtC(pixelID = 1:nPixC, buffer = as.integer(1:nPixC <= 30), ids = 1L))
  call <- function(...) fireSenseUtils::.objfunSpreadFit(
    par = c(maxAsymptote = 0.27, clim = 2, fuel = 3, youngAge = -1), landscape = r,
    annualDTx1000 = climFix$annualDTx1000, nonAnnualDTx1000 = climFix$nonAnnualDTx1000,
    formulaToFit = formC, historicalFires = climFix$historicalFires, fireBufferedListDT = fb,
    mutuallyExclusive = NULL, tests = "snll_fs", Nreps = 3L, doAssertions = FALSE, verbose = 0,
    lanscape1stQuantileThresh = Inf, thresh = Inf, ...)
  set.seed(7); v0 <- call()
  set.seed(7); v1 <- call(climateRef = NULL)
  expect_identical(v1, v0)
  expect_true(is.finite(v0) && v0 < 1e6)
  ## and climateRef changes what the fires see: the spread probabilities the simulations summarise
  set.seed(7); s0 <- call(returnSims = TRUE)
  set.seed(7); s2 <- call(climateRef = list(clim = 0.3), returnSims = TRUE)
  expect_false(identical(attr(s2, "spreadProb"), attr(s0, "spreadProb")))
  expect_null(formals(fireSenseUtils::.objfunSpreadFit)$climateRef)
})

test_that("spreadFitValidationData() predicts p with the fit's link and records climateRef", {
  skip_if_not_installed("SpaDES.tools")
  r <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20, vals = 1)
  fb <- lapply(climFix$historicalFires, function(h) dtC(pixelID = 1:nPixC, buffer = as.integer(1:nPixC <= 30), ids = 1L))
  par <- c(maxAsymptote = 0.27, clim = 4, fuel = 3, youngAge = -1)
  ref <- list(clim = 0.3)
  d <- spreadFitValidationData(par, r, climFix$annualDTx1000, climFix$nonAnnualDTx1000, formC,
                               climFix$historicalFires, fb, mutuallyExclusive = NULL, Nreps = 2L,
                               climateRef = ref, seed = 1)
  mat <- matOf("2005")
  p0 <- f + (0.27 - f) * plogis(3 * mat[, "fuel"] - mat[, "youngAge"])
  expected <- plogis(qlogis(p0) + 4 * (mat[, "clim"] - 0.3))
  d5 <- d[d$year == "2005"]
  expect_equal(d5$p, unname(expected[d5$pixelID]))
  expect_identical(attr(d, "spreadFitValidation")$climateRef, ref)
  ## the response curves use the same function (no error, a ggplot)
  expect_s3_class(plotSpreadFitResponse(d), "ggplot")
})

test_that("runDEoptim() hands climateRef to the fit and the re-score, and omits it from the key only when NULL", {
  seen <- new.env()
  lower <- stats::setNames(c(0.25, 0, 0), c("maxAsymptote", "clim", "fuel"))
  testthat::local_mocked_bindings(
    clusterSetup = function(...) list(itermax = 5, trace = FALSE, strategy = 2L, NP = 40L, cluster = NULL),
    DEoptimIterative = function(fn, lower, ...) { seen$fit <- list(...)$climateRef; list(list(member = list(pop = matrix(lower + 0.1, 1)))) },
    .package = "clusters")
  testthat::local_mocked_bindings(
    termsInDEoptim = function(...) invisible(NULL),
    rescorePopulation = function(pop, fn, reps, cl, seed = 1L, fnArgs = list()) {
      seen$rescore <- fnArgs; data.table::data.table(member = 1L, rep = 1L, value = 1) })
  seen$omit <- list()
  testthat::local_mocked_bindings(Cache = function(FUN, ..., omitArgs = NULL, cachePath = NULL, .functionName = NULL,
                                                   .cacheExtra = NULL, useCache = NULL) {
    seen$omit[[length(seen$omit) + 1L]] <- omitArgs
    if (is.function(FUN)) FUN(...) else FUN
  }, .package = "fireSenseUtils")
  run <- function(...) suppressMessages(runDEoptim(
    landscape = NULL, annualDTx1000 = NULL, nonAnnualDTx1000 = NULL, fireBufferedListDT = NULL,
    historicalFires = NULL, itermax = 5, trace = FALSE, strategy = 2L, cores = c("hostA", "hostB"),
    paths = list(cachePath = withr::local_tempdir()), lower = lower, upper = lower + 1,
    mutuallyExclusive = NULL, formulaToFit = "~ 0 + clim + fuel", objFunCoresInternal = 1L, covMinMax = NULL,
    maxFireSpread = 0.3, Nreps = 1L, .verbose = FALSE, ...))
  run(climateRef = list(clim = 0.4))
  expect_identical(seen$fit, list(clim = 0.4))
  expect_identical(seen$rescore$climateRef, list(clim = 0.4))
  expect_false("climateRef" %in% seen$omit[[1]])
  seen$omit <- list()
  run()
  expect_null(seen$fit)
  expect_false("climateRef" %in% names(seen$rescore))
  expect_true("climateRef" %in% seen$omit[[1]])
  expect_null(formals(runDEoptim)$climateRef)
})
