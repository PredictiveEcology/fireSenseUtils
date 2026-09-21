## `weighted` is meant to give large fires more influence on the SNLL. It was
## log(pmax(minLik, lik * log(size))) = log(lik) + log(log(size)): an offset, not a weight -- every
## fire's likelihood moved the objective by the same amount whatever its size. A weight multiplies:
## w(size) * log(lik), with w = log(size) (`TRUE` or "log") or sqrt(size) ("sqrt"), divided by
## `sizeWeightMean` so that the weights of a fit average 1 and the SNLL keeps its scale.
##
## Same fake-spread() idiom as test-objFunInner-spreadProbVector.R.

r20 <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20)
fires <- data.table::data.table(cells = c(110L, 310L), size = c(16, 100), ids = 1:2)
simsOf <- list(`110` = c(12, 14, 16, 18, 20), `310` = c(80, 90, 100, 110, 120))

snllFor <- function(annualFires, ...) {
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
      data.table::rbindlist(lapply(as.character(annualFires$cells), function(cl)
        data.table::data.table(initialLocus = as.integer(cl), indices = seq_len(simsOf[[cl]][i]), id = 1L, active = FALSE)))
    },
    .package = "SpaDES.tools"
  )
  fireSenseUtils:::objFunInner(
    yr = "2004", annDTx1000 = NULL, par = 1, parsModel = NULL, annualFires = annualFires,
    nonAnnualDTx1000 = NULL,
    annualFireBufferedDT = data.table::data.table(ids = 1L, pixelID = 1:400, buffer = 1L),
    indexNonAnnual = NULL, colsToUse = "cov", covMinMax = NULL,
    mutuallyExclusive = NULL, doAssertions = FALSE, maxFireSpread = 0.28,
    lowerSpreadProb = 0.13, cells = numeric(400), lanscape1stQuantileThresh = 0.25,
    r = r20, Nreps = 5L, doSNLL_FSTest = TRUE, doMADTest = FALSE, doADTest = FALSE,
    plot.it = FALSE, verbose = 0, ...
  )$SNLL_FS
}

test_that("a size weight multiplies each fire's log-likelihood", {
  skip_if_not_installed("SpaDES.tools")
  small <- snllFor(fires[1], weighted = FALSE)
  large <- snllFor(fires[2], weighted = FALSE)
  expect_equal(snllFor(fires, weighted = FALSE), small + large)
  m <- mean(sqrt(fires$size))
  expect_equal(snllFor(fires, weighted = "sqrt", sizeWeightMean = m), (sqrt(16) * small + sqrt(100) * large) / m)
  m <- mean(log(fires$size))
  expect_equal(snllFor(fires, weighted = "log", sizeWeightMean = m), (log(16) * small + log(100) * large) / m)
  expect_equal(snllFor(fires, weighted = TRUE, sizeWeightMean = m), snllFor(fires, weighted = "log", sizeWeightMean = m))
})

test_that("sizeWeight() gives the weights, and rejects an unknown one", {
  expect_equal(sizeWeight(c(4, 100), FALSE), c(1, 1))
  expect_equal(sizeWeight(c(4, 100), "sqrt"), c(2, 10))
  expect_equal(sizeWeight(c(4, 100), TRUE), log(c(4, 100)))
  expect_error(sizeWeight(4, "size"))
})
