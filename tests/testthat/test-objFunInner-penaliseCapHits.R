## A simulated fire that reaches its cap (multiplier() of the observed size) is censored: "at least this
## big", a runaway. Before, it was scored as a fire of the capped size, about 5x too large, so a
## parameter set that burns the whole landscape was only mildly worse than one that burns a bit too much.
## `objFunInner()` runs the real spreadCpp() here; only the spreadProb surface is faked.

r60 <- terra::rast(nrows = 60, ncols = 60, xmin = 0, xmax = 60, ymin = 0, ymax = 60)

capHitObjective <- function(spLow, spHigh, seed = 1, ...) {
  local_mocked_bindings(
    paramsSeparate = function(...) list(logisticPars = c(0.27, 1, 1), covPars = 1),
    spreadProbFromIntegerCovs = function(...) data.table::data.table(pixelID = 1:3600, cov = 0),
    logisticAll = function(...) seq(spLow, spHigh, length.out = 3600),
    .package = "fireSenseUtils"
  )
  set.seed(seed)
  fireSenseUtils:::objFunInner(
    yr = "2004", annDTx1000 = NULL, par = 1, parsModel = NULL,
    annualFires = data.table::data.table(cells = c(400L, 1000L, 1800L, 2600L, 3300L),
                                         size = c(25L, 30L, 40L, 50L, 60L), ids = 1:5),
    nonAnnualDTx1000 = NULL,
    annualFireBufferedDT = data.table::data.table(ids = 1L, pixelID = 1:3600, buffer = 1L),
    indexNonAnnual = NULL, colsToUse = "cov", covMinMax = NULL,
    mutuallyExclusive = NULL, doAssertions = FALSE, maxFireSpread = 1,
    lowerSpreadProb = 0.05, cells = numeric(3600), lanscape1stQuantileThresh = 1,
    weighted = FALSE, r = r60, Nreps = 6L,
    doSNLL_FSTest = TRUE, doMADTest = TRUE, doADTest = TRUE, doYearArea = TRUE,
    plot.it = FALSE, verbose = 0, ...
  )
}

test_that("a parameter set that saturates spread is much worse with penaliseCapHits", {
  skip_if_not_installed("SpaDES.tools")
  off <- capHitObjective(0.9, 0.99)
  on <- capHitObjective(0.9, 0.99, runawaySize = 3600)
  expect_equal(on$capHits, on$nSims)            # every replicate hit its cap
  expect_gt(on$SNLL_FS, off$SNLL_FS + 50)
  expect_true(all(on$allFireSizes == 3600))     # censored sizes reach the AD term
  expect_true(all(off$allFireSizes < 3600))     # uncensored: the cap
  expect_gt(min(on$annualAreaByRep), max(off$annualAreaByRep))
})

test_that("with no fire at its cap, the objective is identical with and without the penalty", {
  skip_if_not_installed("SpaDES.tools")
  off <- capHitObjective(0.13, 0.16)
  on <- capHitObjective(0.13, 0.16, runawaySize = 3600)
  expect_equal(on$capHits, 0)
  expect_identical(on[setdiff(names(on), c("capHits", "nSims"))],
                   off[setdiff(names(off), c("capHits", "nSims"))])
})

test_that("returnSims keeps the simulated sizes, penalty or not", {
  skip_if_not_installed("SpaDES.tools")
  off <- capHitObjective(0.9, 0.99, returnSims = TRUE)
  on <- capHitObjective(0.9, 0.99, returnSims = TRUE, runawaySize = 3600)
  expect_identical(on$sims, off$sims)
  expect_true(all(on$sims$sim < 3600))
})

test_that(".objfunSpreadFit() defaults runawaySize to the landscape's non-NA pixels, and switches off", {
  seen <- list()
  local_mocked_bindings(
    objFunInner = function(runawaySize, ...) {
      seen[[length(seen) + 1L]] <<- runawaySize
      list(SNLL_FS = 0)
    },
    .package = "fireSenseUtils"
  )
  land <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20, vals = 1)
  land[1:40] <- NA
  run <- function(...) {
    seen <<- list()
    fireSenseUtils::.objfunSpreadFit(
      par = c(0.27, 1, 1, 1), landscape = land,
      annualDTx1000 = list(year2004 = data.table::data.table(pixelID = 1:400, cov = 0L)),
      nonAnnualDTx1000 = list(`year2004` = data.table::data.table(pixelID = 1:400)),
      formulaToFit = "~ 0 + cov",
      historicalFires = list(year2004 = data.frame(size = 40L, cells = 1L, ids = 1L)),
      fireBufferedListDT = list(year2004 = data.table::data.table(ids = 1L, pixelID = 1:400)),
      indexNonAnnual = data.table::data.table(ind = 1L, date = "2004"), doAssertions = FALSE,
      lrgSmallFireYears = list(1L, integer(0)), ...
    )
    seen
  }
  expect_equal(unlist(run()), 360)
  expect_equal(unlist(run(runawaySize = 99)), 99)
  expect_null(unlist(run(penaliseCapHits = FALSE)))
  expect_null(unlist(run(capSizes = FALSE)))
})
