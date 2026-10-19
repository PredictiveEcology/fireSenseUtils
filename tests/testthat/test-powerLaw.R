## The "powerLaw" spread link (spreadLink = "powerLaw" in fireSense_spreadFit), from the percolation calibration:
## p = pc - (pc - floor) * exp(-x / gamma), x the linear predictor, pc (the maxAsymptote slot) fitted, gamma fixed.
## The year effect eps_Y is added to x, not to the logit of p. The default link ("ceiling") is untouched.
## Fixtures (climFix, matOf, spreadProbSeenC, ...) are in helper-spreadLinkFixture.R.

pcT <- 0.27

test_that("powerLawLink() is pc - (pc - floor) * exp(-x / gamma), with p below 0 set to 0", {
  expect_equal(powerLawLink(0, pcT, 1.9, f), f)
  expect_equal(powerLawLink(c(-1, 0.5, 2), 0.25, 2, 0.13), 0.25 - 0.12 * exp(-c(-1, 0.5, 2) / 2))
  expect_equal(powerLawLink(1e3, pcT, 1.9, f), pcT)
  ## log(pc - p) is linear in x with slope -1 / gamma
  x <- seq(-0.5, 3, length.out = 8)
  expect_equal(diff(log(pcT - powerLawLink(x, pcT, 1.9, f))) / diff(x), rep(-1 / 1.9, 7))
  ## below the floor for x < 0, and 0 once that goes negative
  expect_lt(powerLawLink(-0.1, pcT, 1.9, f), f)
  expect_identical(powerLawLink(c(-50, -5), 0.25, 1.9, 0.13), c(0, 0))
  ## a one-column matrix (mat %*% beta) gives a plain vector
  expect_identical(powerLawLink(matrix(c(0, 1), ncol = 1), pcT, 1.9, f), powerLawLink(c(0, 1), pcT, 1.9, f))
  ## a bigger gamma flattens the curve
  expect_gt(powerLawLink(1, pcT, 1, f), powerLawLink(1, pcT, 3, f))
})

test_that("link = 'powerLaw' goes through logisticAll() and spreadProbFromCovs(); the default link is unchanged", {
  mat <- matOf("2004")
  b <- c(clim = 2, fuel = 3, youngAge = -1)
  x <- as.vector(mat %*% b)
  pl <- spreadProbFromCovs(mat, lp3, b, f, link = "powerLaw", powerLawGamma = 1.5)
  expect_equal(as.numeric(pl), powerLawLink(x, lp3[1], 1.5, f))
  expect_equal(as.numeric(logisticAll(lp3, mat, b, f, link = "powerLaw", gamma = 1.5)), as.numeric(pl))
  ## the default gamma is the calibration's
  expect_identical(spreadPowerLawGamma, 1.9)
  expect_identical(spreadProbFromCovs(mat, lp3, b, f, link = "powerLaw"),
                   spreadProbFromCovs(mat, lp3, b, f, link = "powerLaw", powerLawGamma = spreadPowerLawGamma))
  ## "ceiling" output: the logistic, bit for bit, whatever gamma is
  expected <- f + (lp3[1] - f) * plogis(x)
  expect_equal(as.numeric(logisticAll(lp3, mat, b, f)), expected)
  expect_identical(spreadProbFromCovs(mat, lp3, b, f, powerLawGamma = 7), logisticAll(lp3, mat, b, f))
  up <- c(0.27, 1, 1, -0.3)
  expect_identical(spreadProbFromCovs(mat, up, b, f, link = "logistic3pUpper", powerLawGamma = 7),
                   logisticAll(up, mat, b, f, link = "logistic3pUpper"))
  expect_error(spreadProbFromCovs(mat, lp3, b, f, link = "powerLaw", climateRef = list(clim = 0.3)), "climateRef")
  expect_true("powerLaw" %in% spreadLinkOptions)
  expect_identical(spreadLinkOptions[1], "ceiling")
})

test_that("the objective's spreadProb under powerLaw is the hand-built one, with the year effect on x", {
  skip_if_not_installed("SpaDES.tools")
  par <- c(maxAsymptote = pcT, clim = 1, fuel = 1.5, youngAge = -1)
  mat <- matOf("2005")
  x <- as.vector(mat %*% par[-1])
  got <- spreadProbSeenC(par, "2005", link = "powerLaw", powerLawGamma = 1.9)
  expect_true(got$reached)
  expect_equal(got$sp[-2L], powerLawLink(x, pcT, 1.9, f)[-2L])   # cell 2 is the ignition
  expect_equal(got$sp[2L], 1)
  ## year effect: every pixel shifts by the same eps on x ...
  set.seed(3)
  withEps <- spreadProbSeenC(par, "2005", link = "powerLaw", powerLawGamma = 1.9, yearSpreadSD = 0.5)
  keep <- got$sp[-2L] > 0 & withEps$sp[-2L] > 0
  epsX <- function(p) -1.9 * (log(pcT - p[keep]))   # x = -gamma * log(pc - p) + const
  dx <- epsX(withEps$sp[-2L]) - epsX(got$sp[-2L])
  expect_gt(sum(keep), 100)
  expect_lt(diff(range(dx)), 1e-8)
  expect_gt(abs(dx[1]), 1e-6)
  ## ... and so not on the logit of p
  dl <- qlogis(withEps$sp[-2L][keep]) - qlogis(got$sp[-2L][keep])
  expect_gt(diff(range(dl)), 1e-3)
  ## the year effect cannot take p above pc, so no year percolates
  expect_true(all(withEps$sp[-2L] <= pcT))
  ## with the default link the year effect is still on the logit
  parOld <- c(maxAsymptote = pcT, clim = 1, fuel = 1.5, youngAge = -1)
  set.seed(3)
  old0 <- spreadProbSeenC(parOld, "2005")
  set.seed(3)
  old1 <- spreadProbSeenC(parOld, "2005", yearSpreadSD = 0.5)
  de <- qlogis(old1$sp[-2L]) - qlogis(old0$sp[-2L])
  expect_lt(diff(range(de)), 1e-9)
})

test_that("spreadProbGates() and spreadFitValidationData() use the powerLaw link and its gamma", {
  skip_if_not_installed("SpaDES.tools")
  par <- c(maxAsymptote = pcT, clim = 1, fuel = 1.5, youngAge = -1)
  g <- function(gamma, thresh = 0.265) do.call(spreadProbGates, c(list(
    par = par, formulaToFit = formC, link = "powerLaw", powerLawGamma = gamma, years = "2004", mutuallyExclusive = NULL,
    maxFireSpread = 0.28, lowerSpreadProb = f, lanscape1stQuantileThresh = thresh),
    climFix[c("annualDTx1000", "nonAnnualDTx1000", "historicalFires")]))
  expect_identical(g(1.9)$pass, spreadProbSeenC(par, "2004", link = "powerLaw", powerLawGamma = 1.9,
                                                lanscape1stQuantileThresh = 0.265)$reached)
  r <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20, vals = 1)
  fb <- lapply(climFix$historicalFires, function(h) dtC(pixelID = 1:nPixC, buffer = as.integer(1:nPixC <= 30), ids = 1L))
  d <- spreadFitValidationData(par, r, climFix$annualDTx1000, climFix$nonAnnualDTx1000, formC,
                               climFix$historicalFires, fb, mutuallyExclusive = NULL, Nreps = 2L,
                               link = "powerLaw", powerLawGamma = 1.5, seed = 1)
  mat <- matOf("2004")
  d4 <- d[d$year == "2004"]
  expect_equal(d4$p, powerLawLink(as.vector(mat %*% par[-1]), pcT, 1.5, f)[d4$pixelID])
  expect_identical(attr(d, "spreadFitValidation")$powerLawGamma, 1.5)
  expect_s3_class(plotSpreadFitResponse(d), "ggplot")
})

test_that(".objfunSpreadFit() and runDEoptim() carry powerLawGamma; the default link's call and key do not change", {
  skip_if_not_installed("SpaDES.tools")
  expect_identical(formals(fireSenseUtils::.objfunSpreadFit)$powerLawGamma, quote(spreadPowerLawGamma))
  r <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20, vals = 1)
  fb <- lapply(climFix$historicalFires, function(h) dtC(pixelID = 1:nPixC, buffer = as.integer(1:nPixC <= 30), ids = 1L))
  call <- function(...) fireSenseUtils::.objfunSpreadFit(
    par = c(maxAsymptote = pcT, clim = 1, fuel = 1.5, youngAge = -1), landscape = r,
    annualDTx1000 = climFix$annualDTx1000, nonAnnualDTx1000 = climFix$nonAnnualDTx1000,
    formulaToFit = formC, historicalFires = climFix$historicalFires, fireBufferedListDT = fb,
    mutuallyExclusive = NULL, tests = "snll_fs", Nreps = 3L, doAssertions = FALSE, verbose = 0,
    lanscape1stQuantileThresh = Inf, thresh = Inf, ...)
  set.seed(7); v0 <- call()
  set.seed(7); v1 <- call(powerLawGamma = 4)
  expect_identical(v1, v0)                       # gamma is unused without the link
  set.seed(7); s1 <- call(link = "powerLaw", returnSims = TRUE)
  set.seed(7); s2 <- call(link = "powerLaw", powerLawGamma = 4, returnSims = TRUE)
  expect_false(identical(attr(s1, "spreadProb"), attr(s2, "spreadProb")))

  seen <- new.env()
  lower <- stats::setNames(c(0.14, 0, 0), c("maxAsymptote", "clim", "fuel"))
  testthat::local_mocked_bindings(
    clusterSetup = function(...) list(itermax = 5, trace = FALSE, strategy = 2L, NP = 40L, cluster = NULL),
    DEoptimIterative = function(fn, lower, ...) { seen$fit <- list(...); list(list(member = list(pop = matrix(lower + 0.1, 1)))) },
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
  run(link = "powerLaw", powerLawGamma = 2.2)
  expect_identical(seen$fit$powerLawGamma, 2.2)
  expect_identical(seen$rescore$powerLawGamma, 2.2)
  expect_false("powerLawGamma" %in% seen$omit[[1]])
  seen$omit <- list()
  run()
  expect_true("powerLawGamma" %in% seen$omit[[1]])
  expect_false("powerLawGamma" %in% names(seen$rescore))
  expect_identical(powerLawOmitArgs(NULL), "powerLawGamma")
  expect_identical(powerLawOmitArgs("logistic3pUpper"), "powerLawGamma")
  expect_identical(powerLawOmitArgs("powerLaw"), character())
})
