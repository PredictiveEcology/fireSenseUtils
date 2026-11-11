## Per-fire random effect (fitFireSpreadSD): in each replicate every fire draws its own eps ~ N(0, sd^2), added to
## the linear predictor of the cells it spreads to (spreadCpp()'s linkOffset, one per fire), on top of the year's
## draw. Only on the whole landscape (the powerLaw link computes p per cell there). See helper-wholeLandscape.R.

## two fires, 20 and 30 pixels observed, on the whole 60 x 60 landscape
twoFires <- function() {
  fx <- wlFixture(60L)
  fx$fires <- data.table::data.table(cells = c(wlCell(20L, 20L), wlCell(40L, 40L)), size = c(20L, 30L), ids = 1:2)
  fx
}
wholeTwo <- function(...) wlObjective(fx = twoFires(), whole = TRUE, Nreps = 30L, ...)

test_that("fireSpreadSD = 0 gives the objective and the fires of before it existed, to the bit", {
  skip_if_not_installed("SpaDES.tools")
  ## fixtures/objFunInnerWholeBase.rds: the same calls on fireSenseUtils 0.2.3.9156 (15ac7f3), Linux x86_64
  tol <- if (identical(Sys.info()[["sysname"]], "Linux") && identical(R.version$arch, "x86_64")) 0 else 1e-12
  before <- readRDS(test_path("fixtures", "objFunInnerWholeBase.rds"))
  expect_equal(wholeTwo(yearSpreadSD = 0.4, fireSpreadSD = 0, seed = 5), before$year, tolerance = tol)
  expect_equal(wholeTwo(fireSpreadSD = 0, seed = 6), before$noEffect, tolerance = tol)
  expect_equal(wholeTwo(yearSpreadSD = 0.4, fireSpreadSD = 0, seed = 5, returnSims = TRUE, returnBurned = TRUE),
               before$sims, tolerance = tol)
})

test_that("fireSpreadSD > 0 gives each fire its own offset, so per-fire sizes change", {
  skip_if_not_installed("SpaDES.tools")
  seen <- new.env()
  seen$offsets <- list()
  real <- SpaDES.tools::spreadCpp
  local_mocked_bindings(spreadCpp = function(...) {
    a <- list(...)
    seen$offsets[[length(seen$offsets) + 1L]] <- a$linkOffset
    real(...)
  }, .package = "SpaDES.tools")
  off <- wholeTwo(yearSpreadSD = 0.4, seed = 5, returnSims = TRUE)
  on <- wholeTwo(yearSpreadSD = 0.4, fireSpreadSD = 1, seed = 5, returnSims = TRUE)
  ## without it one offset per replicate (the year's); with it one per fire, and they differ
  expect_true(all(lengths(seen$offsets[1:30]) == 1L))
  withFire <- seen$offsets[31:60]
  expect_true(all(lengths(withFire) == 2L))
  expect_true(all(vapply(withFire, function(o) abs(diff(o)) > 0, logical(1))))
  expect_false(identical(on$sims$sim, off$sims$sim))
  ## per fire: the two fires' sizes vary independently, not only together
  bySim <- data.table::dcast(on$sims, rep ~ ids, value.var = "sim")
  expect_lt(abs(stats::cor(log(bySim[["1"]]), log(bySim[["2"]]))), 0.9)
})

test_that("fireSpreadSD needs the whole landscape", {
  expect_error(wlObjective(fireSpreadSD = 0.5), "needs `wholeLandscape = TRUE`")
  expect_silent(checkFireSpreadSD(FALSE, FALSE))
})

test_that("splitSpreadPar takes fireSpreadSD then yearSpreadSD from the end of par", {
  s <- splitSpreadPar(c(a = 0.26, cov1 = 2, fireSpreadSD = 0.3, yearSpreadSD = 0.7))
  expect_equal(c(s$fireSpreadSD, s$yearSpreadSD), c(0.3, 0.7))
  expect_identical(names(s$par), c("a", "hillSlope1", "inflectionPoint1", "cov1"))
  s <- splitSpreadPar(c(a = 0.26, cov1 = 2, fireSpreadSD = 0.3))
  expect_equal(c(s$fireSpreadSD, s$yearSpreadSD), c(0.3, 0))
  s <- splitSpreadPar(c(0.26, 2, 0.3, 0.7), fitYearSpreadSD = TRUE, fitFireSpreadSD = TRUE) # DEoptim: unnamed
  expect_equal(c(s$fireSpreadSD, s$yearSpreadSD), c(0.3, 0.7))
  expect_length(s$par, 4L)
  expect_error(splitSpreadPar(c(a = 0.26, yearSpreadSD = 0.7, fireSpreadSD = 0.3, cov1 = 2)),
               "`fireSpreadSD`, `yearSpreadSD` must be the last parameters, in that order")
  expect_identical(spreadRandomEffectTxt, c("fireSpreadSD", "yearSpreadSD"))
})

test_that("runDEoptim fits fireSpreadSD when the bounds name it, in the fit and the re-score, and only then", {
  seen <- new.env()
  lower <- stats::setNames(c(0.25, 0.1, 0, 0), c("maxAsymptote", "x", "fireSpreadSD", "yearSpreadSD"))
  testthat::local_mocked_bindings(
    clusterSetup = function(...) list(itermax = 5, trace = FALSE, strategy = 2L, NP = 40L, cluster = NULL),
    DEoptimIterative = function(fn, lower, ...) {
      seen$fit <- list(...)$fitFireSpreadSD
      list(list(member = list(pop = matrix(lower + 0.1, 1))))
    },
    .package = "clusters")
  testthat::local_mocked_bindings(
    rescorePopulation = function(pop, fn, reps, cl, seed = 1L, fnArgs = list()) {
      seen$rescore <- fnArgs$fitFireSpreadSD; data.table::data.table(member = 1L, rep = 1L, value = 1) })
  run <- function(lower, whole = TRUE) suppressMessages(runDEoptim(
    landscape = NULL, annualDTx1000 = NULL, nonAnnualDTx1000 = NULL, fireBufferedListDT = NULL,
    historicalFires = NULL, itermax = 5, trace = FALSE, strategy = 2L, cores = c("hostA", "hostB"),
    paths = list(cachePath = withr::local_tempdir()), lower = lower, upper = lower + 1,
    mutuallyExclusive = NULL, formulaToFit = NULL, objFunCoresInternal = 1L, covMinMax = NULL,
    maxFireSpread = 0.3, Nreps = 1L, .verbose = FALSE, link = "powerLaw", penaliseRunaways = FALSE,
    wholeLandscape = whole))
  withr::local_options(reproducible.useCache = FALSE)
  run(lower)
  expect_true(seen$fit); expect_true(seen$rescore)
  run(lower[-3])
  expect_null(seen$fit); expect_null(seen$rescore) # no argument: the fit's cache key is what it was
  expect_error(run(lower[c(1, 2, 4, 3)]), "`fireSpreadSD` then `yearSpreadSD` must be the last elements")
  expect_error(run(lower, whole = FALSE), "needs `wholeLandscape = TRUE`")
})
