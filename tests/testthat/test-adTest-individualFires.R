## The Anderson-Darling term compares the DISTRIBUTION of simulated fire sizes with the observed
## one. The observed sample is single fires, so the simulated sample must be single fires too. It
## was each fire's MEAN over the Nreps replicates, which has a much shorter tail: on ELF 4.3 a model
## compared with its own replicate scored 50 * AD = 3260 with means and 42 with single fires.

test_that("objFunInner returns every replicate's fire size for the adTest, not the per-fire mean", {
  skip_if_not_installed("SpaDES.tools")
  r20 <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20)
  simSizes <- c(2, 10, 300)
  i <- 0L
  local_mocked_bindings(
    paramsSeparate = function(...) list(logisticPars = c(0.27, 1, 1), covPars = 1),
    spreadProbFromIntegerCovs = function(...) data.table::data.table(pixelID = 1:400, cov = 0),
    logisticAll = function(...) seq(0.15, 0.26, length.out = 400),
    .package = "fireSenseUtils"
  )
  local_mocked_bindings(
    spreadCpp = function(...) {
      i <<- i + 1L
      data.table::data.table(initialLocus = list(...)$loci, indices = seq_len(simSizes[i]), id = 1L, active = FALSE)
    },
    .package = "SpaDES.tools"
  )
  out <- fireSenseUtils:::objFunInner(
    yr = "2004", annDTx1000 = NULL, par = 1, parsModel = NULL,
    annualFires = data.table::data.table(cells = 210L, size = 100, ids = 1L),
    nonAnnualDTx1000 = NULL,
    annualFireBufferedDT = data.table::data.table(ids = 1L, pixelID = 1:400, buffer = 1L),
    indexNonAnnual = NULL, colsToUse = "cov", covMinMax = NULL,
    mutuallyExclusive = NULL, doAssertions = FALSE, maxFireSpread = 0.28,
    lowerSpreadProb = 0.13, cells = numeric(400), lanscape1stQuantileThresh = 0.25,
    weighted = FALSE, r = r20, Nreps = 3L,
    doSNLL_FSTest = FALSE, doMADTest = FALSE, doADTest = TRUE,
    plot.it = FALSE, verbose = 0
  )
  expect_identical(sort(as.numeric(out$allFireSizes)), c(2, 10, 300))
})

test_that(".objfunSpreadFit hands the adTest the individual simulated fires", {
  skip_if_not_installed("terra")
  dt <- function(...) data.table::data.table(...)
  yrs <- c("2001", "2002", "2003")
  seen <- NULL
  local_mocked_bindings(
    ## per-fire means of 104, but the individual fires are 2, 10 and 300
    objFunInner = function(...) list(SNLL_FS = 1, fireSizes = c(104, 104), allFireSizes = c(2, 10, 300, 2, 10, 300)),
    adStatistic = function(x, y) { seen <<- x; 1 }
  )
  fireSenseUtils:::.objfunSpreadFit(
    par = c(1, 1),
    landscape = terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10, vals = 1),
    annualDTx1000 = stats::setNames(lapply(yrs, function(y) dt(pixelID = 1:2, cov1 = c(100L, 200L))), yrs),
    nonAnnualDTx1000 = list(`2001_2003` = dt(pixelID = 1:2, cov1 = c(100L, 200L))),
    formulaToFit = "~ cov1",
    historicalFires = list(`2001` = data.frame(size = c(100, 200), cells = 1:2),
                           `2002` = data.frame(size = c(150, 250), cells = 3:4),
                           `2003` = data.frame(size = c(5, 6), cells = 5:6)),
    fireBufferedListDT = stats::setNames(lapply(yrs, function(y) dt(pixelID = 1:2, buffer = c(1L, 0L), ids = 1L)), yrs),
    tests = c("adTest", "snll_fs"), Nreps = 3L, doAssertions = FALSE, plot.it = FALSE, verbose = 0, thresh = Inf
  )
  expect_identical(sort(unique(seen)), c(2, 10, 300))
})

test_that("adStatistic() is the statistic kSamples::ad.test() reports, ties included", {
  skip_if_not_installed("kSamples")
  set.seed(2)
  for (shift in c(0, 1)) {
    obs <- round(rlnorm(300, 3, 2))
    sim <- round(rlnorm(300 * 25, 3 + shift, 2))
    expect_equal(adStatistic(sim, obs), kSamples::ad.test(sim, obs)$ad[1, 1], tolerance = 1e-4)
  }
})

test_that("adStatistic() errors on an empty sample, so the objective returns its fail value, not NaN", {
  ## every year can bail before simulating; DEoptim stops on a NaN objective
  expect_error(adStatistic(NULL, c(5, 10)))
  expect_error(adStatistic(c(5, 10), numeric(0)))
})

test_that("adStatistic() does not overflow on the sample sizes of a real fit", {
  ## 1800 fires x 50 replicates: N * M passes .Machine$integer.max
  set.seed(3)
  sim <- round(rlnorm(90000, 3, 2)); obs <- round(rlnorm(1800, 3, 2))
  expect_no_warning(ad <- adStatistic(sim, obs))
  expect_true(is.finite(ad))
})
