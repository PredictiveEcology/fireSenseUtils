## The "powerLawSmooth" and "richards" spread links: a lower asymptote L (0 <= L <= floor). The "powerLaw" link treats
## the floor as p at x = 0 and lets x < 0 fall to 0, so a cell with a low linear predictor cannot burn.
## Fixtures (climFix, matOf, lp3, f, ...) are in helper-spreadLinkFixture.R; wlObjective() in helper-wholeLandscape.R.

pcS <- 0.27
Lm <- spreadLowerAsymptoteMin   # a lower asymptote is never below it (checkLinkParams())

## spreadCpp() with the links 2 and 3 (SpaDES.tools PR "powerLawSmooth" and "richards"); linkProbCpp() is its test hook
skipUnlessSmoothLinks <- function() {
  skip_if_not_installed("SpaDES.tools")
  skip_if_not(
                exists("linkProbCpp", asNamespace("SpaDES.tools")), "SpaDES.tools without the links powerLawSmooth and richards")
}

test_that("with L = 0.13 both links give 0.13 at x = -50; with L = the minimum both stay above 0; the power law is 0", {
  expect_identical(powerLawLink(-50, pcS, 1.9, f), 0)
  expect_equal(powerLawSmoothLink(-50, pcS, 1.9, f, lower = 0.13), 0.13)
  expect_equal(richardsLink(-50, pcS, lower = 0.13), 0.13)
  expect_gt(powerLawSmoothLink(-50, pcS, 1.9, f, lower = Lm), 0)
  expect_gt(richardsLink(-50, pcS, lower = Lm), 0)
  ## never below L, and it falls monotonically with x
  x <- seq(-60, 60, length.out = 241)
  for (L in c(Lm, 0.05, 0.13)) {
    ps <- powerLawSmoothLink(x, pcS, 1.9, f, lower = L)
    pr <- richardsLink(x, pcS, lower = L)
    expect_true(all(ps >= L) && all(pr >= L))
    expect_true(all(diff(ps) >= 0) && all(diff(pr) >= 0))
    expect_true(all(ps <= pcS) && all(pr <= pcS))
  }
})

test_that("powerLawSmooth is the power law for x >= 0, continuous at 0, and the floor for L = floor", {
  x <- c(0, 1e-12, 0.01, 0.5, 2, 10, 100)
  for (shape in c(0.4, 1, 2)) for (L in c(Lm, 0.07, 0.13)) {
    expect_identical(powerLawSmoothLink(x, pcS, 1.9, f, shape, L), powerLawLink(x, pcS, 1.9, f, shape))
    expect_equal(powerLawSmoothLink(c(-1e-14, 0), pcS, 1.9, f, shape, L), c(f, f), tolerance = 1e-3)
  }
  expect_identical(powerLawSmoothLink(c(-100, -1, -1e-9), pcS, 1.9, f, 1, lower = f), rep(f, 3))
  ## the x < 0 side, by the formula
  x <- c(-0.5, -3)
  expect_equal(powerLawSmoothLink(x, 0.25, 2, 0.13, 1, Lm), Lm + (0.13 - Lm) * exp(-((0.25 - 0.13) / (0.13 - Lm)) * abs(x) / 2))
  expect_equal(powerLawSmoothLink(x, 0.25, 2, 0.13, 2, 0.05),
               0.05 + 0.08 * exp(-((0.25 - 0.13) / 0.08) * (abs(x) / 2)^2))
  ## NA stays NA; a one-column matrix gives a plain vector
  expect_identical(powerLawSmoothLink(c(NA, -1), pcS, 1.9, f, 1, Lm)[1], NA_real_)
  expect_identical(powerLawSmoothLink(matrix(c(-1, 1), ncol = 1), pcS, 1.9, f, 1, Lm),
                   powerLawSmoothLink(c(-1, 1), pcS, 1.9, f, 1, Lm))
})

test_that("richards: nu = 0 is Gompertz, a tiny nu is within 1e-6 of it, nu = 1 is the logistic", {
  x <- seq(-30, 30, length.out = 121)
  B <- 1 / spreadPowerLawGamma
  for (L in c(Lm, 0.13)) {
    gomp <- L + (pcS - L) * exp(-exp(-B * x))
    expect_equal(richardsLink(x, pcS, L, 0), gomp)
    expect_lt(max(abs(richardsLink(x, pcS, L, 1e-8) - richardsLink(x, pcS, L, 0))), 1e-6)
    expect_equal(richardsLink(x, pcS, L, 1), L + (pcS - L) * plogis(B * x))
  }
  ## a vector of nu (one per cell) works, with zeros
  expect_equal(richardsLink(c(-1, 2, 3), pcS, 0.13, c(0, 0.5, 1)),
               c(richardsLink(-1, pcS, 0.13, 0), richardsLink(2, pcS, 0.13, 0.5), richardsLink(3, pcS, 0.13, 1)))
  expect_identical(richardsLink(c(NA, 1), pcS)[1], NA_real_)
})

test_that("richards does not use gamma, the calibration or shape: its cache key leaves them out", {
  expect_false(any(c("lowerAsymptote", "richardsNu") %in% powerLawOmitArgs("richards")) )
  expect_true(all(c("powerLawGamma", "powerLawShape") %in% powerLawOmitArgs("richards")))
  mat <- matOf("2004"); b <- c(clim = 2, fuel = 3, youngAge = -1)
  expect_identical(spreadProbFromCovs(mat, lp3, b, f, link = "richards", powerLawGamma = 0.3, powerLawShape = 3),
                   spreadProbFromCovs(mat, lp3, b, f, link = "richards"))
  expect_identical(powerLawOmitArgs("powerLawSmooth"), c("powerLawShape", "richardsNu"))
  expect_identical(powerLawOmitArgs("powerLawSmooth", 2), "richardsNu")
  expect_identical(powerLawOmitArgs("powerLaw", 2), c("lowerAsymptote", "richardsNu"))
})

test_that("both links go through logisticAll() and spreadProbFromCovs(), and the old links are unchanged", {
  mat <- matOf("2004"); b <- c(clim = 2, fuel = 3, youngAge = -1)
  x <- as.vector(mat %*% b)
  ps <- spreadProbFromCovs(mat, lp3, b, f, link = "powerLawSmooth", powerLawGamma = 1.5, powerLawShape = 2, lowerAsymptote = 0.05)
  expect_equal(as.numeric(ps), powerLawSmoothLink(x, lp3[1], 1.5, f, 2, 0.05))
  expect_equal(as.numeric(logisticAll(lp3, mat, b, f, link = "powerLawSmooth", gamma = 1.5, shape = 2, lowerAsymptote = 0.05)),
               as.numeric(ps))
  pr <- spreadProbFromCovs(mat, lp3, b, f, link = "richards", lowerAsymptote = 0.05, richardsNu = 2)
  expect_equal(as.numeric(pr), richardsLink(x, lp3[1], 0.05, 2))
  ## the default lower asymptote is the floor
  expect_identical(spreadProbFromCovs(mat, lp3, b, f, link = "richards"),
                   spreadProbFromCovs(mat, lp3, b, f, link = "richards", lowerAsymptote = spreadProbFloor,
                                      richardsNu = spreadRichardsNu))
  ## the old links ignore the new arguments, bit for bit
  expect_identical(spreadProbFromCovs(mat, lp3, b, f, link = "powerLaw", lowerAsymptote = Lm, richardsNu = 3),
                   spreadProbFromCovs(mat, lp3, b, f, link = "powerLaw"))
  expect_identical(spreadProbFromCovs(mat, lp3, b, f, lowerAsymptote = Lm, richardsNu = 3), logisticAll(lp3, mat, b, f))
  ## climate outside the ceiling does not go with a linear-predictor link
  expect_error(spreadProbFromCovs(mat, lp3, b, f, link = "richards", climateRef = c(clim = 0.3)), "richards")
  expect_true(all(c("powerLawSmooth", "richards") %in% spreadLinkOptions))
  expect_identical(spreadFitLowerAsymptoteTxt, "lowerAsymptote_spread")
  expect_identical(spreadFitRichardsNuTxt, "richardsNu_spread")
})

test_that("the R links are the C++ links to within 1e-12 (spreadCpp's linkPars from linearLinkPars())", {
  skipUnlessSmoothLinks()
  x <- c(-1e3, -50, -20, -5, -1, -1e-8, 0, 1e-8, 0.3, 1, 5, 20, 1e3)
  for (L in c(Lm, 0.13)) for (link in spreadLinearLinks) {
    pars <- linearLinkPars(link, 0.26, 1.9, 0.13, 0.4, L, 0.5)
    expect_identical(names(pars), spreadLinkParNames[[link]])
    ## to within 1e-12, not bit for bit: compilers and platforms differ in the last bit
    expect_equal(linearLinkProb(x, link, pars), SpaDES.tools:::linkProbCpp(match(link, spreadLinearLinks), x, pars), tolerance = 1e-12)
  }
  expect_identical(linearLinkProb(c(-1, 2), "richards", rbind(linearLinkPars("richards", 0.26), linearLinkPars("richards", 0.26))),
                   richardsLink(c(-1, 2), 0.26))
  expect_error(linearLinkProb(1, "probit", c(pc = 1)), "link")
})

## An escaped fire (it burns its first escapeMinPx cells whatever p, but cannot spread where no neighbour has p > 0)
## on a landscape whose linear predictor is -20 everywhere.
test_that("an escaped fire at x = -20 reaches its 9-cell minimum with the new links, and is stuck at 1 cell with powerLaw", {
  skipUnlessSmoothLinks()
  size <- function(link, L, x = -20, nu = spreadRichardsNu, jumpTries = 20, jumpMeanDist = 3, ...) {
    s <- wlObjective(whole = TRUE, Nreps = 3L, seed = 2, returnSims = TRUE, link = link, lowerSpreadProb = f,
                     lowerAsymptote = L, richardsNu = nu, escapeMinPx = 9, jumpTries = jumpTries, jumpMeanDist = jumpMeanDist,
                     cov = function(pix) rep(x, length(pix)), gateRefusal = "simulate", ...)
    s$sims$sim
  }
  expect_true(all(size("powerLaw", f) == 1))
  for (L in c(0.13, Lm)) for (lk in c("powerLawSmooth", "richards")) expect_true(all(size(lk, L) >= 9))
  ## Gompertz with the smallest lower asymptote, at a linear predictor far below any cell's
  for (x in c(-20, -200)) expect_true(all(size("richards", Lm, x = x, nu = 0) >= 9))
})

test_that("on the whole landscape the new links burn the cells of the full spreadProb vector, to the bit", {
  skipUnlessSmoothLinks()
  closed <- function() { fx <- wlFixture(11L); data.table::set(fx$buf, NULL, c("edge", "open"), list(FALSE, FALSE)); fx }
  for (link in c("powerLawSmooth", "richards")) for (L in c(Lm, 0.13)) {
    args <- list(yearSpreadSD = 0.5, escapeMinPx = 9, jumpTries = 20, jumpMeanDist = 3, link = link, lowerSpreadProb = f,
                 lowerAsymptote = L, returnSims = TRUE, returnBurned = TRUE, seed = 3, cov = function(pix) 0.5 + 1.5 * sin(pix))
    full <- do.call(wlObjective, c(list(fx = closed(), annPix = seq_len(3600L)), args))
    lazy <- do.call(wlObjective, c(list(whole = TRUE), args))
    expect_identical(lazy$burned, full$burned)
    expect_identical(lazy$sims, full$sims)
  }
})

test_that("the gates use the lower asymptote as their lower bound with the new links", {
  expect_identical(fireSenseUtils:::spreadLinkLower("richards", 0.13, Lm), Lm)
  expect_identical(fireSenseUtils:::spreadLinkLower("powerLawSmooth", 0.13, 0.05), 0.05)
  expect_identical(fireSenseUtils:::spreadLinkLower("powerLaw", 0.13, 0.05), 0.13)
  expect_identical(fireSenseUtils:::spreadLinkLower(NULL, 0.13, 0.05), 0.13)
  skipUnlessSmoothLinks()
  ## most pixels have x = 0, so p = 0.27 * exp(-1) = 0.099 (richards, nu = 0, L = the minimum): below the floor, but not below L
  par <- c(maxAsymptote = pcS, clim = 0, fuel = 0, youngAge = -1)
  g <- function(...) do.call(spreadProbGates, c(list(
    par = par, formulaToFit = formC, years = "2004", mutuallyExclusive = NULL, maxFireSpread = 0.28,
    lowerSpreadProb = f, lanscape1stQuantileThresh = 0.265, ...),
    climFix[c("annualDTx1000", "nonAnnualDTx1000", "historicalFires")]))
  expect_true(g(link = "richards", lowerAsymptote = Lm, richardsNu = 0)$medianOK)
  expect_lt(richardsLink(0, pcS, Lm, 0), f)        # so the floor as the gates' bound would refuse this year
})

test_that("runDEoptim() passes lowerAsymptote and richardsNu to the fit and the re-score, only with the links that use them", {
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
  testthat::local_mocked_bindings(Cache = function(FUN, ..., omitArgs = NULL, cachePath = NULL, .functionName = NULL,
                                                   .cacheExtra = NULL, useCache = NULL) {
    seen$omit <- c(seen$omit, list(omitArgs))
    if (is.function(FUN)) FUN(...) else FUN
  }, .package = "fireSenseUtils")
  run <- function(...) {
    seen$omit <- seen$rescore <- NULL
    suppressMessages(runDEoptim(
      landscape = NULL, annualDTx1000 = NULL, nonAnnualDTx1000 = NULL, fireBufferedListDT = NULL,
      historicalFires = NULL, itermax = 5, trace = FALSE, strategy = 2L, cores = c("hostA", "hostB"),
      paths = list(cachePath = withr::local_tempdir()), lower = lower, upper = lower + 1,
      mutuallyExclusive = NULL, formulaToFit = "~ 0 + clim + fuel", objFunCoresInternal = 1L, covMinMax = NULL,
      maxFireSpread = 0.3, Nreps = 1L, .verbose = FALSE, ...))
  }
  run(link = "richards", lowerAsymptote = 0.02, richardsNu = 2)
  expect_identical(seen$fit$lowerAsymptote, 0.02)
  expect_identical(seen$fit$richardsNu, 2)
  expect_identical(seen$rescore$lowerAsymptote, 0.02)
  expect_identical(seen$rescore$richardsNu, 2)
  expect_false(any(c("lowerAsymptote", "richardsNu") %in% seen$omit[[1]]))
  expect_false(any(c("powerLawGamma", "powerLawShape") %in% names(seen$rescore)))
  expect_true(all(c("powerLawGamma", "powerLawShape") %in% seen$omit[[1]]))
  run(link = "powerLawSmooth", powerLawGamma = 2.2, lowerAsymptote = Lm)
  expect_identical(seen$rescore$lowerAsymptote, Lm)
  expect_identical(seen$rescore$powerLawGamma, 2.2)
  expect_false("richardsNu" %in% names(seen$rescore))
  expect_true("richardsNu" %in% seen$omit[[1]])
  expect_false("lowerAsymptote" %in% seen$omit[[1]])
  ## the old links and the default: left out of the key and the re-score
  for (lk in list(NULL, "powerLaw")) {
    run(link = lk)
    expect_true(all(c("lowerAsymptote", "richardsNu") %in% seen$omit[[1]]))
    expect_false(any(c("lowerAsymptote", "richardsNu") %in% names(seen$rescore)))
  }
})

test_that("spreadFitValidationData() and plotSpreadFitResponse() work with the new links", {
  skipUnlessSmoothLinks()
  par <- c(maxAsymptote = pcS, clim = 1, fuel = 1.5, youngAge = -1)
  r <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20, vals = 1)
  fb <- lapply(climFix$historicalFires, function(h) dtC(pixelID = 1:nPixC, buffer = as.integer(1:nPixC <= 30), ids = 1L))
  for (L in c(Lm, 0.13)) {
    d <- spreadFitValidationData(par, r, climFix$annualDTx1000, climFix$nonAnnualDTx1000, formC,
                                 climFix$historicalFires, fb, mutuallyExclusive = NULL, Nreps = 2L,
                                 link = "richards", lowerAsymptote = L, richardsNu = 0.5, seed = 1)
    d4 <- d[d$year == "2004"]
    expect_equal(d4$p, richardsLink(as.vector(matOf("2004") %*% par[-1]), pcS, L, 0.5)[d4$pixelID])
    expect_identical(attr(d, "spreadFitValidation")$lowerAsymptote, L)
    expect_s3_class(plotSpreadFitResponse(d), "ggplot")
  }
})

test_that("the parameters of the links can be one per x", {
  x <- c(-3, -0.5, 0, 2, -8, 1)
  pc <- c(0.25, 0.3, 0.25, 0.26, 0.3, 0.27); L <- c(Lm, 0.13, 0.05, 0.1, 0.02, Lm)
  one <- function(i) powerLawSmoothLink(x[i], pc[i], 1.9, 0.13, 1.5, L[i])
  expect_identical(powerLawSmoothLink(x, pc, 1.9, 0.13, 1.5, L), vapply(seq_along(x), one, numeric(1)))
  oneR <- function(i) richardsLink(x[i], pc[i], L[i], c(0, 0.5, 1, 2, 0, 0.5)[i])
  expect_identical(richardsLink(x, pc, L, c(0, 0.5, 1, 2, 0, 0.5)), vapply(seq_along(x), oneR, numeric(1)))
})

test_that("checkLinkParams() is the one check of lowerAsymptote and richardsNu", {
  expect_true(checkLinkParams("richards", Lm, 0))
  expect_identical(Lm, 1e-3)
  ## 0 would let a fire under minSize retry for about 1 / p generations
  for (lk in c("powerLawSmooth", "richards")) {
    expect_error(checkLinkParams(lk, 0, 0.5), "at least 0.001")
    expect_error(checkLinkParams(lk, Lm / 2, 0.5), "minSize")
  }
  expect_true(checkLinkParams("powerLawSmooth", spreadProbFloor))
  expect_true(checkLinkParams("powerLaw", 5, -1))        # nothing to check
  expect_true(checkLinkParams(NULL))
  for (lk in c("powerLawSmooth", "richards")) for (bad in list(-0.01, 0.2, NA_real_, c(Lm, 0.1), NULL, "a"))
    expect_error(checkLinkParams(lk, bad, 0.5), "lowerAsymptote must be")
  for (bad in list(-1, NA_real_, NULL, c(1, 2)))
    expect_error(checkLinkParams("richards", 0.1, bad), "richardsNu must be")
  expect_identical(spreadRichardsB, 1 / spreadPowerLawGamma)
  expect_identical(formals(richardsLink)$B, quote(spreadRichardsB))
})

test_that("spreadLinkParNames is spreadCpp()'s list of link parameters, names and order", {
  skipUnlessSmoothLinks()
  expect_identical(spreadLinkParNames, SpaDES.tools:::.spreadCppLinks)
})
