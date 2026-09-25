## The annual-area and area-weighted size-distribution terms of .objfunSpreadFit() (off by default),
## the jump pass-through and the minSize <= maxSize guard.

test_that("areaWeightedCvM() is 0 for identical area shares, matches a hand calculation, and caps an empty sim", {
  awc <- fireSenseUtils:::areaWeightedCvM
  expect_equal(awc(c(1, 3, 10), c(1, 3, 10)), 0)
  expect_equal(awc(rep(c(1, 3, 10), 5), c(1, 3, 10)), 0)
  ## obs 1 and 3: area shares at x = 1, 3 are 0.25, 1; sim all 1s: 1, 1 -> 2 * (0.25 * 0.75^2) = 0.28125
  expect_equal(awc(c(1, 1, 1, 1), c(1, 3)), 0.28125)
  obs <- c(2, 5, 20, 150, 900)
  expect_gt(awc(c(obs, 20000), obs), awc(c(obs, 1000), obs))
  expect_equal(awc(numeric(0), obs), length(obs))   # no simulated fires: the maximum, not an error
})

test_that("yearAreaNLL() is the t likelihood on annual totals, responds to a shift, and charges refused years", {
  set.seed(1)
  sims <- list(year2001 = round(rlnorm(50, log(400), 0.4)), year2002 = round(rlnorm(50, log(60), 0.4)))
  obs <- c(year2001 = 400, year2002 = 60)
  nll <- fireSenseUtils:::yearAreaNLL(obs, sims, sizeLik = "t", sizeLikDf = 5)
  byHand <- -log(fireSenseUtils:::sizeLikT(400, sims$year2001, 5)) - log(fireSenseUtils:::sizeLikT(60, sims$year2002, 5))
  expect_equal(nll, byHand)
  expect_gt(fireSenseUtils:::yearAreaNLL(c(year2001 = 1600, year2002 = 60), sims, "t", 5), nll)
  ## a year the objective refused to simulate costs -log(minLik), so refusing is never rewarded
  refused <- fireSenseUtils:::yearAreaNLL(obs, list(year2001 = sims$year2001), "t", 5)
  expect_equal(refused, -log(fireSenseUtils:::sizeLikT(400, sims$year2001, 5)) - log(1e-29))
  expect_true(is.finite(fireSenseUtils:::yearAreaNLL(obs, sims, sizeLik = "kde")))
})

test_that("the new arguments are off by default, in the objective and in runDEoptim()", {
  for (f in list(formals(fireSenseUtils:::.objfunSpreadFit), formals(fireSenseUtils::runDEoptim))) {
    expect_identical(f$jumpTries, 0); expect_identical(f$jumpMeanDist, 0)
    expect_identical(f$yearAreaWeight, 0); expect_identical(f$areaDistWeight, 0)
  }
  expect_false(formals(fireSenseUtils:::.objfunSpreadFit)$returnTerms)
  b <- paste(deparse(body(fireSenseUtils::runDEoptim)), collapse = "\n")
  expect_equal(lengths(regmatches(b, gregexpr("yearAreaWeight = yearAreaWeight", b))), 2L)  # fit and re-score
})

## objFunInner() with spread faked (the idiom of test-objFunInner-escapeMinSize.R): the arguments it gives spreadCpp()
r20 <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20)
spreadCppSeen <- function(capTo = NULL, ...) {
  seen <- list()
  mocks <- list(
    paramsSeparate = function(...) list(logisticPars = c(0.27, 1, 1), covPars = 1),
    spreadProbFromIntegerCovs = function(...) data.table::data.table(pixelID = 1:400, cov = 0),
    logisticAll = function(...) seq(0.15, 0.26, length.out = 400))
  if (!is.null(capTo)) mocks$multiplier <- function(size, ...) rep(capTo, length(size))
  do.call(local_mocked_bindings, c(mocks, .package = "fireSenseUtils"))
  local_mocked_bindings(
    spreadCpp = function(...) {
      seen[[length(seen) + 1L]] <<- list(...)
      data.table::data.table(initialLocus = list(...)$loci, indices = 1:20, id = 1L, active = FALSE)
    },
    .package = "SpaDES.tools")
  fireSenseUtils:::objFunInner(
    yr = "2004", annDTx1000 = NULL, par = 1, parsModel = NULL,
    annualFires = data.table::data.table(cells = 210L, size = 25L, ids = 1L), nonAnnualDTx1000 = NULL,
    annualFireBufferedDT = data.table::data.table(ids = 1L, pixelID = 1:400, buffer = 1L),
    indexNonAnnual = NULL, colsToUse = "cov", covMinMax = NULL, mutuallyExclusive = NULL, doAssertions = FALSE,
    maxFireSpread = 0.28, lowerSpreadProb = 0.13, cells = numeric(400), lanscape1stQuantileThresh = 0.25,
    weighted = FALSE, r = r20, Nreps = 3L, doSNLL_FSTest = TRUE, doMADTest = FALSE, doADTest = FALSE,
    plot.it = FALSE, verbose = 0, ...)
  seen
}

test_that("the per-fire size cap is raised to the escape size (spreadCpp needs minSize <= maxSize)", {
  skip_if_not_installed("SpaDES.tools")
  a <- spreadCppSeen(capTo = 3, escapeMinPx = 9L)
  expect_true(all(vapply(a, function(x) min(x$maxSize) >= 9 && identical(x$minSize, 9L), logical(1))))
  b <- spreadCppSeen(capTo = 3)                       # without an escape size the cap is left alone
  expect_true(all(vapply(b, function(x) identical(max(x$maxSize), 3), logical(1))))
})

test_that("jumpTries and jumpMeanDist reach spreadCpp() with an escape size", {
  skip_if_not_installed("SpaDES.tools")
  a <- spreadCppSeen(escapeMinPx = 9L, jumpTries = 20, jumpMeanDist = 3)
  expect_true(all(vapply(a, function(x) identical(x$jumpTries, 20) && identical(x$jumpMeanDist, 3), logical(1))))
})
