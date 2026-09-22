## Fit diagnostics (?fitDiagnostics): the checks made by hand on the phase-2 fits (September 2026),
## now run by fireSense_SpreadFit after every fit. The held-out validation used trace() to pull the
## simulated fires out of the objective; .objfunSpreadFit(returnSims = TRUE) replaces it.

## objFunInner with only the spread mocked, as in test-logistic-upperTail.R: 400 pixels, one fire of
## observed size 10 at cell 210, every replicate burns 25 cells.
callInner <- function(seen, ...) {
  local_mocked_bindings(
    paramsSeparate = function(...) list(logisticPars = c(0.27, 1, 1), covPars = 1),
    spreadProbFromIntegerCovs = function(...) data.table::data.table(pixelID = 1:400, cov = 0),
    logisticAll = function(...) c(rep(0.2, 300), rep(0.269, 100)),
    .package = "fireSenseUtils", .env = parent.frame())
  local_mocked_bindings(
    spreadCpp = function(...) {
      seen$maxSize <- list(...)$maxSize
      data.table::data.table(initialLocus = list(...)$loci, indices = 1:25, id = 1L, active = FALSE)
    },
    .package = "SpaDES.tools", .env = parent.frame())
  r20 <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20)
  fireSenseUtils:::objFunInner(
    yr = "2004", annDTx1000 = NULL, par = 1, parsModel = NULL,
    annualFires = data.table::data.table(cells = 210L, size = 10, ids = 1L),
    nonAnnualDTx1000 = NULL,
    annualFireBufferedDT = data.table::data.table(ids = 1L, pixelID = 1:400, buffer = 1L),
    indexNonAnnual = NULL, colsToUse = "cov", covMinMax = NULL,
    mutuallyExclusive = NULL, doAssertions = FALSE, maxFireSpread = 0.28,
    lowerSpreadProb = 0.13, cells = numeric(400), lanscape1stQuantileThresh = Inf,
    weighted = FALSE, r = r20, Nreps = 3L,
    doSNLL_FSTest = TRUE, doMADTest = FALSE, doADTest = FALSE,
    plot.it = FALSE, verbose = 0, ...)
}

test_that("objFunInner(returnSims = TRUE) returns each replicate's simulated size beside the observed", {
  seen <- new.env()
  out <- callInner(seen, returnSims = TRUE)
  expect_equal(nrow(out$sims), 3L)
  expect_equal(out$sims$rep, 1:3)
  expect_equal(out$sims$sim, rep(25L, 3))
  expect_equal(out$sims$size, rep(10, 3))
  expect_identical(out$sims$yr, rep("2004", 3))
  ## the spread probabilities of the year's pixels come back summarised: 100 of 400 at the ceiling
  expect_equal(out$pSummary$n, 400)
  expect_equal(out$pSummary$atCeiling, 100)
})

test_that("capSizes = FALSE lifts the size cap the fit puts on each fire", {
  seen <- new.env()
  callInner(seen, returnSims = TRUE)
  expect_true(all(is.finite(seen$maxSize)))
  callInner(seen, returnSims = TRUE, capSizes = FALSE)
  expect_equal(seen$maxSize, Inf)
})

test_that(".objfunSpreadFit(returnSims = TRUE) stacks every year and combines the spread probabilities", {
  skip_if_not_installed("purrr")
  dt <- function(...) data.table::data.table(...)
  yrs <- c("2001", "2002", "2003")
  local_mocked_bindings(objFunInner = function(yr, ...) list(
    sims = dt(yr = yr, rep = 1L, initialLocus = 1L, ids = 1L, size = 5, sim = 7L),
    pSummary = spreadProbSummary(c(0.1, 0.26), ceiling = 0.26)))
  out <- fireSenseUtils:::.objfunSpreadFit(
    par = c(1, 1),
    landscape = terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10, vals = 1),
    annualDTx1000 = stats::setNames(lapply(yrs, function(y) dt(pixelID = 1:2, cov1 = c(100L, 200L))), yrs),
    nonAnnualDTx1000 = list(`2001_2003` = dt(pixelID = 1:2, cov1 = c(100L, 200L))),
    historicalFires = list(`2001` = data.frame(size = c(100, 200), cells = 1:2),
                           `2002` = data.frame(size = c(150, 250), cells = 3:4),
                           `2003` = data.frame(size = c(5, 6), cells = 5:6)),
    fireBufferedListDT = stats::setNames(lapply(yrs, function(y) dt(pixelID = 1:2, buffer = c(1L, 0L), ids = 1L)), yrs),
    formulaToFit = "~ cov1", tests = "snll_fs", Nreps = 1L, doAssertions = FALSE, verbose = 0,
    thresh = 1, returnSims = TRUE)
  ## thresh = 1 would stop an objective after the first batch; the simulations still cover every year
  expect_setequal(out$yr, yrs)
  expect_equal(attr(out, "spreadProb")$n, 6)
  expect_equal(attr(out, "spreadProb")$atCeiling, 3)
})

test_that("linkSaturation() reads the share at the ceiling and the quantiles from the summaries", {
  p <- c(rep(0.15, 50), rep(0.2, 20), rep(0.26, 30))
  s <- combineSpreadProbSummaries(list(spreadProbSummary(p[1:60], 0.26),
                                       spreadProbSummary(p[61:100], 0.26)))
  sims <- data.table::data.table(member = 1L)
  data.table::setattr(sims, "spreadProb", list(s, spreadProbSummary(rep(0.26, 10), 0.26)))
  ls <- linkSaturation(sims, probsP = c(0.25, 0.60, 0.90))
  expect_equal(ls$propAtCeiling, c(0.3, 1))
  expect_equal(ls$pixelYears, c(100, 10))
  expect_equal(unlist(ls[1, c("p_q25", "p_q60", "p_q90")]), c(p_q25 = 0.15, p_q60 = 0.2, p_q90 = 0.26),
               tolerance = 0.001)
})

test_that("scoreFireSizes(): a perfect simulation scores zero error; a 10x one scores 1", {
  d <- data.table::CJ(member = 1:2, yr = c("2001", "2002"), ids = 1:3, rep = 1:2)
  d[, size := ids * 10 + as.integer(yr) - 2000]
  d[, sim := size]
  s <- scoreFireSizes(d)
  expect_equal(s$fireBias, 0); expect_equal(s$fireRMSE, 0)
  expect_equal(s$yearBias, 0); expect_equal(s$totalAreaRatio, 1)
  expect_equal(s$qSim_q50, s$qObs_q50)
  expect_equal(s$years, 2L); expect_equal(s$fires, 6L)
  d[, sim := 10 * size]
  s <- scoreFireSizes(d)
  expect_equal(s$fireBias, 1); expect_equal(s$yearRMSE, 1)
  expect_equal(s$simsOver10xObs, 0) # exactly 10x is not over
  ## a year with no simulation is counted, and scored on the years that have one
  d[yr == "2002", sim := NA]
  s <- scoreFireSizes(d)
  expect_equal(s$memberYearsNotSimulated, 2L)
  expect_equal(s$years, 1L)
})

test_that("coefIdentifiability() scores only covariates, and zPop is the median in population sd", {
  set.seed(1)
  pop <- cbind(maxAsymptote = 0.26, hillSlope1 = runif(200, 0.5, 1), inflectionPoint1 = 1,
               agb_pine = seq(10, 30, length.out = 200),   # median 20, 5-95% range 18
               agb_fir = seq(-20, 20, length.out = 200))   # median 0
  lower <- c(maxAsymptote = 0.25, hillSlope1 = 0.1, inflectionPoint1 = 0.1, agb_pine = -60, agb_fir = -60)
  upper <- c(maxAsymptote = 0.27, hillSlope1 = 2, inflectionPoint1 = 4, agb_pine = 60, agb_fir = 60)
  id <- coefIdentifiability(pop, lower, upper)
  expect_equal(id$coef, c("agb_pine", "agb_fir"))
  expect_equal(id$zPop[1], 20 / (18 / 3.29), tolerance = 1e-6)
  expect_equal(id$spread90[1], 18 / 120, tolerance = 1e-6)
  expect_equal(id$signPinned, c(TRUE, FALSE))
  ## a coefficient every member holds at 0 has no sign to pin (zPop is 0/0)
  pop[, "agb_fir"] <- 0
  expect_false(coefIdentifiability(pop, lower, upper)$signPinned[2])
  ## and one every member holds at 5 is pinned
  pop[, "agb_fir"] <- 5
  expect_true(coefIdentifiability(pop, lower, upper)$signPinned[2])
})

test_that("profileCoefficients() pairs every point with the best member on the same seeds", {
  ## The objective's noise depends only on the seed, so paired differences have no noise at all:
  ## any deltaSE above zero would mean the points and the reference were drawn with different seeds.
  fn <- function(par, target) sum((par - target)^2) + stats::rnorm(1, sd = 5)
  target <- c(maxAsymptote = 0.26, a = 3, b = 0)
  pop <- rbind(target, target + 1, target - 1)
  pr <- profileCoefficients(best = target, pop = pop, fn = fn, reps = 4L,
                            fnArgs = list(target = target))
  expect_equal(pr$coef, c("", rep(c("a", "b"), each = 6)))
  expect_equal(pr$deltaSE, rep(0, 13))
  ## dropping `a` (optimum 3) costs 9; dropping `b` (optimum 0) costs nothing
  expect_equal(pr[isZero == TRUE]$delta, c(9, 0))
  id <- data.table::data.table(coef = c("a", "b"), signPinned = c(TRUE, TRUE))
  ii <- identifiedInIsolation(id, pr)
  expect_equal(ii$identified, c(TRUE, FALSE))
})

test_that("runDEoptim() profiles and simulates only when asked, on the re-score's arguments", {
  seen <- new.env()
  finalPop <- matrix(c(0.26, 1, 1, 0.5, 0.26, 1, 1, 0.7), nrow = 2, byrow = TRUE)
  lower <- stats::setNames(c(0.25, 0.2, 0.1, 0), c("maxAsymptote", "hillSlope1", "inflectionPoint1", "x"))
  testthat::local_mocked_bindings(
    clusterSetup = function(...) list(itermax = 5, trace = FALSE, strategy = 2L, NP = 40L, cluster = NULL),
    DEoptimIterative2 = function(...) list(list(member = list(pop = finalPop))),
    .package = "clusters")
  testthat::local_mocked_bindings(
    termsInDEoptim = function(...) invisible(NULL),
    rescorePopulation = function(pop, fn, reps, cl, seed = 1L, fnArgs = list()) {
      seen$rescoreArgs <- fnArgs
      data.table::data.table(member = 1:2, rep = 1L, value = c(2, 1))
    },
    profileCoefficients = function(best, pop, fn, reps, cl, fnArgs, ...) {
      seen$best <- best; seen$profileArgs <- fnArgs; "profile"
    },
    simulateFireSizes = function(pop, fn, fnArgs, cl, ...) {
      seen$simPop <- pop; "sims"
    })
  run <- function(...) suppressMessages(runDEoptim(
    landscape = NULL, annualDTx1000 = NULL, nonAnnualDTx1000 = NULL,
    fireBufferedListDT = NULL, historicalFires = NULL,
    itermax = 5, trace = FALSE, strategy = 2L, cores = c("hostA", "hostB"),
    paths = list(cachePath = withr::local_tempdir()),
    lower = lower, upper = lower + 1, mutuallyExclusive = NULL, formulaToFit = NULL,
    objFunCoresInternal = 1L, covMinMax = NULL, maxFireSpread = 0.3, Nreps = 1L,
    .verbose = FALSE, ...))
  withr::local_options(reproducible.useCache = FALSE)
  DE <- run()
  expect_null(attr(DE, "profile")); expect_null(attr(DE, "fitSims"))
  DE <- run(profileReps = 3L, simulateMembers = 2L, sizeLik = "t")
  expect_identical(attr(DE, "profile"), "profile")
  expect_identical(attr(DE, "fitSims"), "sims")
  ## the best member is the re-score's (member 2, mean 1), not DEoptim's
  expect_equal(unname(seen$best[["x"]]), 0.7)
  expect_equal(nrow(seen$simPop), 2L)
  expect_identical(seen$profileArgs, seen$rescoreArgs)
  expect_identical(seen$profileArgs$sizeLik, "t")
})

test_that("fitConvergence() gives one row per chunk, ignoring fail values", {
  DE <- list(list(member = list(bestvalit = c(50, 40), popval = c(60, 55, 1e6, 70, 65))),
             list(member = list(bestvalit = c(38, 35, 35), popval = c(40, 41, 39, 42, 1e7))))
  fc <- fitConvergence(DE)
  expect_equal(fc$generation, c(2L, 5L))
  expect_equal(fc$bestval, c(40, 35))
  expect_equal(fc$popMedian, c(62.5, 40.5))
})
