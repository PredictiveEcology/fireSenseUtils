## An optional free intercept with centred covariates in the spread model (formula "~ 1 + ..."). Without
## one, the covariates are all >= 0 once rescaled, so the level of the linear predictor is set by the
## coefficients alone and they trade off against each other (ELF 13.1, 2026-10-05). With one, the
## coefficients describe variation about the centre and the intercept the level. The intercept is a
## column of 1s, named spreadInterceptTxt, first among the covariate coefficients; it is not rescaled or
## centred. Off (a "~ 0 + ..." formula), nothing changes.

dtI <- function(...) data.table::data.table(...)
nPixI <- 400L
set.seed(11)
icptFix <- list(
  annualDTx1000 = list(`2004` = dtI(pixelID = 1:nPixI, cov1 = sample(0:1000, nPixI, TRUE)),
                       `2005` = dtI(pixelID = 1:nPixI, cov1 = sample(0:1000, nPixI, TRUE)),
                       `2006` = dtI(pixelID = 1:nPixI, cov1 = sample(0:1000, nPixI, TRUE))),
  nonAnnualDTx1000 = list(`2000` = dtI(pixelID = 1:nPixI, cov2 = sample(0:1000, nPixI, TRUE))),
  historicalFires = list(`2004` = data.frame(size = 50, cells = 1L, ids = 1L),
                         `2005` = data.frame(size = 80, cells = 2L, ids = 1L),
                         `2006` = data.frame(size = 5, cells = 3L, ids = 1L)))
icptCentre <- list(cov1 = 0.4, cov2 = 0.55)

## objFunInner() for one year, with the real covariate pipeline and the real link; spread() is mocked and
## hands back the spreadProb it was given (pixels of the year, ignition cell excluded by the caller)
spreadProbSeen <- function(par, formula, covCentre = NULL, yr = "2004") {
  seen <- new.env()
  local_mocked_bindings(
    spreadCpp = function(...) {
      a <- list(...)
      seen$sp <- a$spreadProb
      dtI(initialLocus = a$loci, indices = a$loci, id = seq_along(a$loci), active = FALSE)
    }, .package = "SpaDES.tools", .env = parent.frame())
  colsToUse <- spreadDesignCols(formula)
  p <- fireSenseUtils:::splitSpreadPar(par)$par
  fireSenseUtils:::objFunInner(
    yr = yr, annDTx1000 = data.table::copy(icptFix$annualDTx1000[[yr]]), par = p,
    parsModel = length(colsToUse),
    annualFires = data.table::as.data.table(icptFix$historicalFires[[yr]]),
    nonAnnualDTx1000 = data.table::copy(icptFix$nonAnnualDTx1000),
    annualFireBufferedDT = dtI(ids = 1L, pixelID = 1:nPixI),
    indexNonAnnual = dtI(ind = 1L, date = "2000"), colsToUse = colsToUse, covMinMax = NULL,
    mutuallyExclusive = NULL, doAssertions = TRUE, maxFireSpread = 0.28, lowerSpreadProb = 0.13,
    cells = numeric(nPixI), lanscape1stQuantileThresh = Inf, weighted = TRUE,
    r = terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20), Nreps = 1L,
    doSNLL_FSTest = TRUE, doMADTest = FALSE, doADTest = FALSE, plot.it = FALSE, verbose = 0,
    returnSims = TRUE, covCentre = covCentre)
  seen$sp
}

## the spread probability from first principles: the plain logistic of intercept + centred covariates
handBuilt <- function(maxAsymptote, intercept, b1, b2, yr = "2004", centre = icptCentre) {
  x1 <- icptFix$annualDTx1000[[yr]]$cov1 / 1000
  x2 <- icptFix$nonAnnualDTx1000[[1]]$cov2 / 1000
  0.13 + (maxAsymptote - 0.13) * plogis(intercept + b1 * (x1 - centre$cov1) + b2 * (x2 - centre$cov2))
}

test_that("spreadDesignCols(): the intercept comes first, only when the formula has one", {
  expect_identical(spreadDesignCols("~ 0 + a + b"), c("a", "b"))
  expect_identical(spreadDesignCols("~ 1 + a + b"), c(spreadInterceptTxt, "a", "b"))
  expect_identical(spreadDesignCols(as.formula("~ a")), c(spreadInterceptTxt, "a"))
  expect_identical(spreadDesignCols(NULL), character(0))
  expect_identical(spreadInterceptTxt, "(Intercept)")
  expect_identical(spreadCovCols(c(spreadInterceptTxt, "a", "b")), c("a", "b"))
  expect_identical(spreadCovCols(c("a", "b")), c("a", "b"))
})

test_that("the objective's spreadProb with an intercept is the hand-built linear predictor", {
  skip_if_not_installed("SpaDES.tools")
  par <- c(maxAsymptote = 0.27, `(Intercept)` = -0.7, cov1 = 2.5, cov2 = -1.5)
  got <- spreadProbSeen(par, "~ 1 + cov1 + cov2", covCentre = icptCentre)
  expected <- handBuilt(0.27, -0.7, 2.5, -1.5)
  expected[1L] <- 1                                  # the ignition cell
  expect_equal(got[1:nPixI], expected)
  ## the intercept is the level: it moves every pixel by the same amount on the logit scale
  par2 <- par; par2[["(Intercept)"]] <- 0.3
  got2 <- spreadProbSeen(par2, "~ 1 + cov1 + cov2", covCentre = icptCentre)
  expect_equal(qlogis((got2[-1] - 0.13) / 0.14) - qlogis((got[-1] - 0.13) / 0.14), rep(1, nPixI - 1L))
  ## and an intercept with a centre is the intercept i - sum(b * centre) with none
  noCentre <- spreadProbSeen(c(maxAsymptote = 0.27, `(Intercept)` = -0.7 - 2.5 * 0.4 + 1.5 * 0.55,
                               cov1 = 2.5, cov2 = -1.5), "~ 1 + cov1 + cov2")
  expect_equal(noCentre, got)
})

test_that("off, the objective's spreadProb is what it was: no intercept, no centre", {
  skip_if_not_installed("SpaDES.tools")
  par <- c(maxAsymptote = 0.27, cov1 = 2.5, cov2 = -1.5)
  got <- spreadProbSeen(par, "~ 0 + cov1 + cov2")
  expected <- handBuilt(0.27, 0, 2.5, -1.5, centre = list(cov1 = 0, cov2 = 0))
  expected[1L] <- 1
  expect_equal(got[1:nPixI], expected)
})

test_that("spreadProbFromIntegerCovs() adds the intercept last, as 1s, never rescaled, centred or asserted", {
  run <- function(cols, covCentre = NULL) {
    spreadProbFromIntegerCovs(
      shortAnnDTx1000 = dtI(pixelID = 1:3, a = c(0L, 500L, 1000L), b = c(0L, 200L, 400L)), yr = 2000,
      covMinMax = NULL, mutuallyExclusive = NULL, colsToUse = cols, doAssertions = TRUE,
      logisticPars = c(0.2, 1, 1), covPars = rep(1, length(cols)), maxFireSpread = 0.28,
      lowerSpreadProb = 0.1, covCentre = covCentre)
  }
  x <- run(c(spreadInterceptTxt, "a", "b"), covCentre = list(a = 0.5, b = 0.2))
  expect_identical(x[[spreadInterceptTxt]], rep(1, 3))
  expect_equal(x$a, c(-0.5, 0, 0.5))
  expect_equal(x$b, c(-0.2, 0, 0.2))
  ## a centre named for the intercept is ignored: it is not a covariate
  expect_identical(run(c(spreadInterceptTxt, "a"), covCentre = list("(Intercept)" = 9))[[spreadInterceptTxt]], rep(1, 3))
  ## without it, no such column
  expect_false(spreadInterceptTxt %in% names(run(c("a", "b"))))
})

test_that("spreadProbGates() with an intercept agrees with the gates applied to the hand-built spreadProb", {
  pars <- list(c(maxAsymptote = 0.27, `(Intercept)` = -0.7, cov1 = 2.5, cov2 = -1.5),
               c(maxAsymptote = 0.27, `(Intercept)` = 9, cov1 = 0.1, cov2 = 0.1),    # saturated: too burny
               c(maxAsymptote = 0.27, `(Intercept)` = -9, cov1 = 0.1, cov2 = 0.1))   # at the floor: not spread
  g <- do.call(spreadProbGates, c(list(par = pars, formulaToFit = "~ 1 + cov1 + cov2",
                                       covCentre = icptCentre), icptFix))
  byHand <- vapply(pars, function(p) all(vapply(c("2004", "2005"), function(y) {
    sp <- handBuilt(0.27, p[["(Intercept)"]], p[["cov1"]], p[["cov2"]], yr = y)
    spreadProbGateTest(sp, 0.27, 0.13, 0.28, 0.265)$pass
  }, logical(1))), logical(1))
  expect_identical(g$pass, byHand)
  expect_true(any(g$pass) && any(!g$pass))
})

test_that("spreadCovCentre() is the mean of the rescaled covariates over every pixel-year", {
  cc <- spreadCovCentre(icptFix$annualDTx1000, icptFix$nonAnnualDTx1000, "~ 1 + cov1 + cov2")
  expect_named(cc, c("cov1", "cov2"))
  expect_equal(cc$cov1, mean(unlist(lapply(icptFix$annualDTx1000, `[[`, "cov1"))) / 1000)
  expect_equal(cc$cov2, mean(icptFix$nonAnnualDTx1000[[1]]$cov2) / 1000)
  ## the data are not changed by finding it
  expect_identical(icptFix$annualDTx1000[[1]]$cov1[1:3], icptFix$annualDTx1000[["2004"]]$cov1[1:3])
  ## rescaled with covMinMax, then made mutually exclusive, as the objective does it
  mm <- dtI(cov1 = c(0, 2), cov2 = c(0, 1))
  me <- spreadCovCentre(icptFix$annualDTx1000, icptFix$nonAnnualDTx1000, "~ 1 + cov1 + cov2",
                        covMinMax = mm, mutuallyExclusive = list(cov1 = "cov2"))
  expect_equal(me$cov1, cc$cov1 / 2)
  expect_lt(me$cov2, cc$cov2)
  ## without an intercept nothing is centred
  expect_null(spreadCovCentre(icptFix$annualDTx1000, icptFix$nonAnnualDTx1000, "~ 0 + cov1 + cov2"))
})

## ---- runDEoptim -----------------------------------------------------------------------------

callRunDEoptimI <- function(lower, ...) {
  withr::local_options(reproducible.useCache = FALSE)
  suppressMessages(runDEoptim(
    landscape = NULL, annualDTx1000 = NULL, nonAnnualDTx1000 = NULL,
    fireBufferedListDT = NULL, historicalFires = NULL,
    itermax = 5, trace = FALSE, strategy = 2L, cores = c("hostA", "hostB"),
    paths = list(cachePath = withr::local_tempdir()),
    lower = lower, upper = lower + 1, mutuallyExclusive = NULL, formulaToFit = "~ 1 + x",
    objFunCoresInternal = 1L, covMinMax = NULL, maxFireSpread = 0.3, Nreps = 1L,
    .verbose = FALSE, ...))
}
mockDEoptim <- function(seen, finalPop) {
  testthat::local_mocked_bindings(
    clusterSetup = function(...) list(itermax = 5, trace = FALSE, strategy = 2L, NP = 40L, cluster = NULL),
    DEoptimIterative = function(fn, lower, upper, control, ...) {
      seen$dots <- list(...)
      list(list(member = list(pop = finalPop)))
    },
    .package = "clusters", .env = parent.frame())
  testthat::local_mocked_bindings(
    rescorePopulation = function(pop, fn, reps, cl, seed = 1L, fnArgs = list()) {
      seen$fnArgs <- fnArgs
      data.table::data.table(member = 1L, rep = 1L, value = 1)
    },
    Cache = function(FUN, ..., omitArgs = NULL, cachePath = NULL, .functionName = NULL,
                     .cacheExtra = NULL, useCache = NULL) {
      seen$omitArgs <- c(seen$omitArgs, list(omitArgs))
      if (is.function(FUN)) FUN(...) else FUN
    }, .env = parent.frame())
}

test_that("runDEoptim passes covCentre to the objective in the fit and the re-score", {
  seen <- new.env()
  lowerI <- c(maxAsymptote = 0.25, `(Intercept)` = -1, x = -1)
  mockDEoptim(seen, matrix(c(0.26, 0, 1), nrow = 1))
  centre <- list(x = 0.37)
  callRunDEoptimI(lowerI, covCentre = centre)
  expect_identical(seen$dots$covCentre, centre)
  expect_identical(seen$fnArgs$covCentre, centre)
  ## the centre is part of the fit's key, so it is not omitted from it
  expect_false(any(vapply(seen$omitArgs, function(o) "covCentre" %in% o, logical(1))))
})

test_that("off, covCentre is omitted from the cache key and absent from the re-score's arguments", {
  seen <- new.env()
  lowerN <- c(maxAsymptote = 0.25, x = -1)
  mockDEoptim(seen, matrix(c(0.26, 1), nrow = 1))
  callRunDEoptimI(lowerN)
  expect_null(seen$dots$covCentre)
  expect_false("covCentre" %in% names(seen$fnArgs))
  ## so the call keeps the key it had before covCentre existed
  expect_true("covCentre" %in% seen$omitArgs[[1]])
  expect_identical(omitNullArgs(a = NULL, b = 1, c = NULL), c("a", "c"))
  expect_length(omitNullArgs(b = 1), 0L)
})

test_that("runDEoptim names the intercept a covariate, not a logistic parameter", {
  lowerI <- c(maxAsymptote = 0.25, `(Intercept)` = -1, x = -1)
  seen <- new.env()
  mockDEoptim(seen, matrix(c(0.26, 0, 1), nrow = 1))
  withr::local_options(reproducible.useCache = FALSE)
  msg <- testthat::capture_messages(runDEoptim(
    landscape = NULL, annualDTx1000 = NULL, nonAnnualDTx1000 = NULL,
    fireBufferedListDT = NULL, historicalFires = NULL, itermax = 5, trace = FALSE, strategy = 2L,
    cores = c("hostA", "hostB"), paths = list(cachePath = withr::local_tempdir()),
    lower = lowerI, upper = lowerI + 1, mutuallyExclusive = NULL, formulaToFit = "~ 1 + x",
    objFunCoresInternal = 1L, covMinMax = NULL, maxFireSpread = 0.3, Nreps = 1L, .verbose = FALSE))
  expect_match(paste(msg, collapse = ""), "logistic: maxAsymptote; covariates: (Intercept), x", fixed = TRUE)
})
