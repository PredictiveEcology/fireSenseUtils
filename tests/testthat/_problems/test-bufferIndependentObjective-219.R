# Extracted from test-bufferIndependentObjective.R:219

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "fireSenseUtils", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
fixture <- readRDS(test_path("fixtures", "bufferIndependentObjective.rds"))
h <- 0.25
shift <- log(dnorm(0) / h)
r3 <- terra::rast(nrows = 5, ncols = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 5, vals = 1)
r3[seq(1, 41, by = 10)] <- NA
cell3 <- function(rows, cols) as.vector(outer(rows, cols, function(r, c) (r - 1) * 10 + c))
toy3 <- data.table::data.table(ids = rep(1:3, c(15, 15, 10)), pixelID = c(cell3(1:5, 2:4), cell3(1:5, 5:7), cell3(1:5, 9:10)))
toyObjective <- function(burn, Nreps = 2L, ...) {
  seen <- new.env()
  fires <- data.table::data.table(cells = cell3(3, c(3, 6, 10)), size = c(10L, 10L, 5L), ids = 1:3)
  local_mocked_bindings(
    cropToCells = function(r, cells, ...) list(r = r, ncell = terra::ncell(r), toCrop = identity, toFull = identity),
    paramsSeparate = function(...) list(logisticPars = c(0.27, 1, 1), covPars = 1),
    spreadProbFromIntegerCovs = function(...) data.table::data.table(pixelID = toy3$pixelID, cov = 0),
    logisticAll = function(...) seq(0.06, 0.08, length.out = nrow(toy3)),
    .package = "fireSenseUtils")
  local_mocked_bindings(
    spreadCpp = function(landscape, loci, ...) {
      seen$args <- list(...)
      data.table::rbindlist(lapply(seq_along(loci), function(k)
        data.table::data.table(initialLocus = loci[k], indices = unique(c(loci[k], burn[[k]])))))
    },
    .package = "SpaDES.tools")
  args <- utils::modifyList(
    list(yr = "2001", annDTx1000 = NULL, par = 1, parsModel = NULL, annualFires = fires, nonAnnualDTx1000 = NULL,
         annualFireBufferedDT = addBufferEdge(list(toy3), r3)[[1]], indexNonAnnual = NULL, colsToUse = "cov",
         covMinMax = NULL, mutuallyExclusive = NULL, doAssertions = FALSE, maxFireSpread = 1, lowerSpreadProb = 0.05,
         cells = numeric(50), lanscape1stQuantileThresh = 1, weighted = FALSE, r = r3, Nreps = Nreps,
         doSNLL_FSTest = TRUE, doMADTest = FALSE, plot.it = FALSE, verbose = 0, sizeLik = "kernel"),
    list(...))
  out <- do.call(fireSenseUtils:::objFunInner, args)
  c(out, list(spreadArgs = seen$args))
}

# test -------------------------------------------------------------------------
skip_if_not_installed("SpaDES.tools")
burn <- list(cell3(2:4, 2:4), cell3(2:4, 5:7), cell3(3, 10))
out <- toyObjective(burn, doADTest = TRUE, doYearArea = TRUE, doMADTest = TRUE)
expect_equal(out$runaways, 2)
expect_equal(out$allFireSizes, rep(9, 4))
expect_equal(out$annualAreaByRep, c(18, 18))
expect_equal(out$fireSizes, c(9, 9, 1))
expect_equal(toyObjective(burn, sizeLik = "kde")$runaways, 2)
