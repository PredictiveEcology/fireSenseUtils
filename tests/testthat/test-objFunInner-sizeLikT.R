## `sizeLik = "t"`: the per-fire size likelihood is a Student-t on the square-root scale, centred on
## the mean of the simulated sizes with their standard deviation as its scale, instead of a kernel
## density of them. A kernel density is zero away from the simulated sizes, so an observed fire the
## simulations did not reach scores the `minLik` floor however far away it is. The t has no floor
## and a miss costs more the further it is.
##
## Same fake-spread() idiom as test-objFunInner-spreadProbVector.R.

r20 <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20)
floorSNLL <- -log(1e-29)

## SNLL_FS of one fire of `obs` pixels when replicate i of spread() burns simSizes[i] cells
snllFor <- function(simSizes, obs, ...) {
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
      data.table::data.table(initialLocus = list(...)$loci, indices = seq_len(simSizes[i]), id = 1L, active = FALSE)
    },
    .package = "SpaDES.tools"
  )
  fireSenseUtils:::objFunInner(
    yr = "2004", annDTx1000 = NULL, par = 1, parsModel = NULL,
    annualFires = data.table::data.table(cells = 210L, size = obs, ids = 1L),
    nonAnnualDTx1000 = NULL,
    annualFireBufferedDT = data.table::data.table(ids = 1L, pixelID = 1:400, buffer = 1L),
    indexNonAnnual = NULL, colsToUse = "cov", covMinMax = NULL,
    mutuallyExclusive = NULL, doAssertions = FALSE, maxFireSpread = 0.28,
    lowerSpreadProb = 0.13, cells = numeric(400), lanscape1stQuantileThresh = 0.25,
    weighted = FALSE, r = r20, Nreps = length(simSizes),
    doSNLL_FSTest = TRUE, doMADTest = FALSE, doADTest = FALSE,
    plot.it = FALSE, verbose = 0, ...
  )$SNLL_FS
}
sims <- c(20, 25, 30, 22, 28, 24, 26, 27, 23, 29) # simulated fires of about 25 pixels

test_that("the default is the kernel density: a fire the simulations never reach scores the floor", {
  skip_if_not_installed("SpaDES.tools")
  expect_equal(snllFor(sims, obs = 100), floorSNLL)
  expect_equal(snllFor(sims, obs = 100), snllFor(sims, obs = 100, sizeLik = "kde"))
})

test_that("with sizeLik = 't' a missed fire is above the floor, and a bigger miss costs more", {
  skip_if_not_installed("SpaDES.tools")
  near <- snllFor(sims, obs = 100, sizeLik = "t")
  far <- snllFor(sims, obs = 380, sizeLik = "t")
  hit <- snllFor(sims, obs = 25, sizeLik = "t")
  expect_lt(far, floorSNLL)
  expect_lt(near, far)
  expect_lt(hit, near)
})

test_that("sizeLikDf sets the tail: fewer degrees of freedom forgive a miss more", {
  skip_if_not_installed("SpaDES.tools")
  expect_lt(snllFor(sims, obs = 380, sizeLik = "t", sizeLikDf = 3), snllFor(sims, obs = 380, sizeLik = "t", sizeLikDf = 30))
})

test_that("sizeLikT() is the t density of sqrt(size), and survives identical simulated sizes", {
  N <- c(16, 25, 36, 49)
  m <- mean(sqrt(N)); s <- sd(sqrt(N))
  expect_equal(sizeLikT(30, N, df = 5), dt((sqrt(30) - m) / s, 5) / s)
  expect_true(is.finite(sizeLikT(30, c(25, 25, 25), df = 5))) # sd is 0: the scale has a lower bound
})
