## objFunInner() hands spread() the year's spreadProb: zero everywhere, this year's spreadProb at
## this year's pixels, 1 at the ignition cells. It once built that vector landscape-length before
## deciding whether the year bails, then scanned all of it for the values the bail tests summarise.
## Now the bail tests read the spreadProb column, and spread() gets the year's bounding box only
## (see test-cropToCells.R), with its result mapped back to landscape cell indices.
## These tests pin what reaches spread(), what comes back, and what the bail tests see.

## a 10 x 10 landscape; this year's pixels are 23, 24, 25 (row 3) and 34, 35 (row 4)
r10 <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10)
pix <- c(23L, 24L, 25L, 34L, 35L)

callInner <- function(sp, cells = numeric(100), loci = 24L, lanscape1stQuantileThresh = 0.25,
                      spreadReturns = NULL) {
  captured <- NULL
  local_mocked_bindings(
    paramsSeparate = function(...) list(logisticPars = c(0.27, 1, 1), covPars = 1),
    spreadProbFromIntegerCovs = function(...)
      data.table::data.table(pixelID = pix, cov = 0),
    logisticAll = function(...) sp,
    .package = "fireSenseUtils"
  )
  local_mocked_bindings(
    spread = function(...) {
      captured <<- list(...)
      if (!is.null(spreadReturns)) return(spreadReturns(list(...)))
      stop(structure(class = c("reachedSpread", "error", "condition"),
                     list(message = "reached spread()", call = NULL)))
    },
    .package = "SpaDES.tools"
  )
  out <- tryCatch(
    fireSenseUtils:::objFunInner(
      yr = "2004", annDTx1000 = NULL, par = 1, parsModel = NULL,
      annualFires = data.table::data.table(cells = loci, size = 50, ids = 1L),
      nonAnnualDTx1000 = NULL,
      annualFireBufferedDT = data.table::data.table(ids = 1L, pixelID = 24L, buffer = 1L),
      indexNonAnnual = NULL, colsToUse = "cov", covMinMax = NULL,
      mutuallyExclusive = NULL, doAssertions = FALSE, maxFireSpread = 0.28,
      lowerSpreadProb = 0.13, cells = cells,
      lanscape1stQuantileThresh = lanscape1stQuantileThresh,
      weighted = TRUE, r = r10, Nreps = 2L,
      doSNLL_FSTest = TRUE, doMADTest = FALSE, doADTest = FALSE,
      plot.it = FALSE, verbose = 0
    ),
    reachedSpread = function(e) "reachedSpread"
  )
  list(out = out, spreadArgs = captured)
}

sp <- c(0.15, 0.20, 0.22, 0.24, 0.26)

## what spread() was handed, put back on the 10 x 10 landscape using terra's coordinates only
onLandscape <- function(args) {
  full <- numeric(100)
  idx <- terra::cellFromXY(r10, terra::xyFromCell(args$landscape, seq_along(args$spreadProb)))
  full[idx] <- args$spreadProb
  list(spreadProb = full, loci = terra::cellFromXY(r10, terra::xyFromCell(args$landscape, args$loci)))
}

test_that("spread() receives this year's spreadProb at this year's pixels, 1 at the ignition", {
  skip_if_not_installed("SpaDES.tools")
  res <- callInner(sp)
  expect_identical(res$out, "reachedSpread")
  expected <- numeric(100); expected[pix] <- sp; expected[24L] <- 1
  got <- onLandscape(res$spreadArgs)
  expect_identical(got$spreadProb, expected)
  expect_identical(got$loci, 24)
})

test_that("spread() gets the year's bounding box plus a one-cell margin, not the landscape", {
  skip_if_not_installed("SpaDES.tools")
  a <- callInner(sp)$spreadArgs
  expect_identical(c(terra::nrow(a$landscape), terra::ncol(a$landscape)), c(2 + 2, 3 + 2))
  expect_length(a$spreadProb, 20L)
})

test_that("what spread() burned comes back as landscape cell indices", {
  skip_if_not_installed("SpaDES.tools")
  ## a fake spread() that burns the ignition and the crop cell to its right, every replicate
  burnTwo <- function(args) data.table::data.table(
    initialLocus = args$loci, indices = c(args$loci, args$loci + 1L), id = 1L, active = FALSE)
  out <- callInner(sp, spreadReturns = burnTwo)$out
  ## observed size 50 against two replicates that each burned 2 cells: what matters here is that
  ## the join on the ignition's LANDSCAPE cell (24) found the simulated fires at all
  expect_type(out, "list")
  expect_true(is.finite(out$SNLL_FS))
  missed <- callInner(sp, spreadReturns = function(args) data.table::data.table(
    initialLocus = 99L, indices = c(99L, 100L), id = 1L, active = FALSE))$out
  expect_false(isTRUE(all.equal(out$SNLL_FS, missed$SNLL_FS)))
})

test_that("a year that bails never reaches spread()", {
  skip_if_not_installed("SpaDES.tools")
  ## 1st quartile of sp is 0.20; a threshold below it fails `lowSPLowEnough`
  res <- callInner(sp, lanscape1stQuantileThresh = 0.10)
  expect_null(res$spreadArgs)
  expect_type(res$out, "list")
})

test_that("the bail tests summarise the spreadProb values, not the landscape's zeros", {
  skip_if_not_installed("SpaDES.tools")
  ## 5 of 100 cells carry a value. If the zeros leaked into the summary its 1st quartile would be
  ## 0 and ANY threshold would pass; 0.10 must still bail.
  expect_null(callInner(sp, lanscape1stQuantileThresh = 0.10)$spreadArgs)
  ## and with the threshold above the true 1st quartile (0.20) the year runs
  expect_identical(callInner(sp, lanscape1stQuantileThresh = 0.21)$out, "reachedSpread")
})
