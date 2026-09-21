## objFunInner() hands spread() a landscape-length spreadProb vector: zero everywhere, this
## year's spreadProb at this year's pixels, and 1 at the ignition cells. It used to build that
## vector before deciding whether the year bails, and then scanned all of it
## (`cells[cells > a | cells > b]`) to get the values the bail tests summarise -- four passes over
## 6.0M cells on ELF 5.3.1 to recover ~41k values that were already in the spreadProb column.
## Now the bail tests read the column, and the vector is filled only when spread() will run.
##
## Measured on ELF 5.3.1 with identical seeds: identical objective values, 1.1-1.4x faster.
## These tests pin the two things that change could break: what reaches spread(), and that the
## bail tests still see the same values.

callInner <- function(sp, cells, loci = 3L, lanscape1stQuantileThresh = 0.25) {
  captured <- NULL
  local_mocked_bindings(
    paramsSeparate = function(...) list(logisticPars = c(0.27, 1, 1), covPars = 1),
    spreadProbFromIntegerCovs = function(...)
      data.table::data.table(pixelID = c(2L, 3L, 5L, 7L, 8L), cov = 0),
    logisticAll = function(...) sp,
    .package = "fireSenseUtils"
  )
  local_mocked_bindings(
    spread = function(...) {
      captured <<- list(...)
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
      annualFireBufferedDT = data.table::data.table(ids = 1L, pixelID = 3L, buffer = 1L),
      indexNonAnnual = NULL, colsToUse = "cov", covMinMax = NULL,
      mutuallyExclusive = NULL, doAssertions = FALSE, maxFireSpread = 0.28,
      lowerSpreadProb = 0.13, cells = cells,
      lanscape1stQuantileThresh = lanscape1stQuantileThresh,
      weighted = TRUE, r = NULL, Nreps = 1L,
      doSNLL_FSTest = TRUE, doMADTest = FALSE, doADTest = FALSE,
      plot.it = FALSE, verbose = 0
    ),
    reachedSpread = function(e) "reachedSpread"
  )
  list(out = out, spreadArgs = captured)
}

sp <- c(0.15, 0.20, 0.22, 0.24, 0.26)

test_that("spread() receives this year's spreadProb at this year's pixels, 1 at the ignition", {
  skip_if_not_installed("SpaDES.tools")
  res <- callInner(sp, cells = numeric(10))
  expect_identical(res$out, "reachedSpread")
  ##          pixel:  1  2     3  4  5     6  7     8     9  10
  expect_identical(res$spreadArgs$spreadProb,
                   c(0, 0.15, 1, 0, 0.22, 0, 0.24, 0.26, 0, 0))
})

test_that("the vector spread() receives is the same whether `cells` arrives integer or numeric", {
  skip_if_not_installed("SpaDES.tools")
  expect_identical(callInner(sp, cells = integer(10))$spreadArgs$spreadProb,
                   callInner(sp, cells = numeric(10))$spreadArgs$spreadProb)
})

test_that(".objfunSpreadFit allocates `cells` as numeric, so filling it does not coerce it", {
  src <- paste(deparse(fireSenseUtils:::.objfunSpreadFit), collapse = "\n")
  expect_match(src, "cells <- numeric(ncells)", fixed = TRUE)
})

test_that("a year that bails never reaches spread()", {
  skip_if_not_installed("SpaDES.tools")
  ## 1st quartile of sp is 0.20; a threshold below it fails `lowSPLowEnough`
  res <- callInner(sp, cells = numeric(10), lanscape1stQuantileThresh = 0.10)
  expect_null(res$spreadArgs)
  expect_type(res$out, "list")
})

test_that("the bail tests summarise the spreadProb values, not the landscape's zeros", {
  skip_if_not_installed("SpaDES.tools")
  ## Here 5 of 1000 cells carry a value. If the zeros leaked into the summary
  ## its 1st quartile would be 0 and ANY threshold would pass; 0.10 must still bail.
  res <- callInner(sp, cells = numeric(1000), lanscape1stQuantileThresh = 0.10)
  expect_null(res$spreadArgs)
  ## and with the threshold above the true 1st quartile (0.20) the year runs
  expect_identical(callInner(sp, cells = numeric(1000), lanscape1stQuantileThresh = 0.21)$out,
                   "reachedSpread")
})
