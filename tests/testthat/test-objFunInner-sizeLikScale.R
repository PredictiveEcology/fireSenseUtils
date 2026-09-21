## The fire-size likelihood is the empirical density of the SIMULATED sizes, evaluated at the
## OBSERVED size. Both are square-rooted first. From 2021-02-10 (48f75ba) until this test, only the
## simulated sizes were: `demp(x = size, obs = sqrt(N))`. spread() caps a simulated fire at
## `multiplier(size)`, so sqrt(N) could never reach `size` for a fire above ~12 pixels; every such
## fire scored the `minLik` floor whatever the parameters, and simulations that matched a fire
## scored the same as simulations that burned its whole buffer.
##
## Same idiom as test-objFunInner-spreadProbVector.R: a fake spread() decides what burned.

## a 20 x 20 landscape, every cell is this year's; one fire of 100 pixels ignited at cell 210
r20 <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20)
floorSNLL <- -log(1e-29)

## SNLL_FS of that one fire when replicate i of spread() burns simSizes[i] cells
snllFor <- function(simSizes) {
  i <- 0L
  local_mocked_bindings(
    paramsSeparate = function(...) list(logisticPars = c(0.27, 1, 1), covPars = 1),
    spreadProbFromIntegerCovs = function(...) data.table::data.table(pixelID = 1:400, cov = 0),
    logisticAll = function(...) seq(0.15, 0.26, length.out = 400),
    .package = "fireSenseUtils"
  )
  local_mocked_bindings(
    spread = function(...) {
      i <<- i + 1L
      a <- list(...)
      data.table::data.table(initialLocus = a$loci, indices = seq_len(simSizes[i]), id = 1L, active = FALSE)
    },
    .package = "SpaDES.tools"
  )
  fireSenseUtils:::objFunInner(
    yr = "2004", annDTx1000 = NULL, par = 1, parsModel = NULL,
    annualFires = data.table::data.table(cells = 210L, size = 100, ids = 1L),
    nonAnnualDTx1000 = NULL,
    annualFireBufferedDT = data.table::data.table(ids = 1L, pixelID = 1:400, buffer = 1L),
    indexNonAnnual = NULL, colsToUse = "cov", covMinMax = NULL,
    mutuallyExclusive = NULL, doAssertions = FALSE, maxFireSpread = 0.28,
    lowerSpreadProb = 0.13, cells = numeric(400), lanscape1stQuantileThresh = 0.25,
    weighted = FALSE, r = r20, Nreps = length(simSizes),
    doSNLL_FSTest = TRUE, doMADTest = FALSE, doADTest = FALSE,
    plot.it = FALSE, verbose = 0
  )$SNLL_FS
}

test_that("simulations that match a 100-pixel fire score above the minLik floor", {
  skip_if_not_installed("SpaDES.tools")
  expect_lt(snllFor(c(90, 95, 98, 100, 102, 105, 110, 97, 103, 99)), floorSNLL)
})

test_that("matching a fire scores better than burning everything around it", {
  skip_if_not_installed("SpaDES.tools")
  matched <- snllFor(c(90, 95, 98, 100, 102, 105, 110, 97, 103, 99))
  burnAll <- snllFor(c(400, 398, 395, 400, 399, 397, 400, 396, 400, 398))
  expect_lt(matched, burnAll)
})
