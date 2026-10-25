## `escapeSizeHa`: the spread model is fitted as "given the fire escaped". Only observed fires of at least
## escapeSizePixels() pixels are fitted, and every simulated fire starts as an escaped fire, burning its
## first escapeSizePixels() cells whatever spreadProb (SpaDES.tools::spreadCpp(minSize =)). Before this,
## an escaped fire was any fire over 1 pixel, while 11-16% of simulated fires never left their first pixel.
##
## Same fake-spread() idiom as test-objFunInner-sizeLikT.R.

test_that("escapeSizePixels() is the fewest whole pixels covering escapeSizeHa", {
  r240 <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 2400, ymin = 0, ymax = 2400) # 5.76 ha
  r100 <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 1000, ymin = 0, ymax = 1000) # 1 ha
  expect_identical(escapeSizePixels(50, r240), 9L)   # 8 pixels are 46.1 ha, 9 are 51.8 ha
  expect_identical(escapeSizePixels(50, r100), 50L)  # exactly 50 pixels, not 51
  expect_identical(escapeSizePixels(0.5, r240), 1L)  # never below one pixel
})

r20 <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20)

## the arguments objFunInner() gives spreadCpp()
spreadCppArgs <- function(...) {
  seen <- list()
  local_mocked_bindings(
    paramsSeparate = function(...) list(logisticPars = c(0.27, 1, 1), covPars = 1),
    spreadProbFromIntegerCovs = function(...) data.table::data.table(pixelID = 1:400, cov = 0),
    logisticAll = function(...) seq(0.15, 0.26, length.out = 400),
    .package = "fireSenseUtils"
  )
  local_mocked_bindings(
    spreadCpp = function(...) {
      seen[[length(seen) + 1L]] <<- list(...)
      data.table::data.table(initialLocus = list(...)$loci, indices = 1:20, id = 1L, active = FALSE)
    },
    .package = "SpaDES.tools"
  )
  fireSenseUtils:::objFunInner(
    yr = "2004", annDTx1000 = NULL, par = 1, parsModel = NULL,
    annualFires = data.table::data.table(cells = 210L, size = 25L, ids = 1L),
    nonAnnualDTx1000 = NULL,
    annualFireBufferedDT = data.table::data.table(ids = 1L, pixelID = 1:400, buffer = 1L),
    indexNonAnnual = NULL, colsToUse = "cov", covMinMax = NULL,
    mutuallyExclusive = NULL, doAssertions = FALSE, maxFireSpread = 0.28,
    lowerSpreadProb = 0.13, cells = numeric(400), lanscape1stQuantileThresh = 0.25,
    weighted = FALSE, r = r20, Nreps = 3L,
    doSNLL_FSTest = TRUE, doMADTest = FALSE, doADTest = FALSE,
    plot.it = FALSE, verbose = 0, ...
  )
  seen
}

test_that("with escapeMinPx every replicate starts as an escaped fire of that many cells", {
  skip_if_not_installed("SpaDES.tools")
  a <- spreadCppArgs(escapeMinPx = 9L)
  expect_length(a, 3L)
  expect_true(all(vapply(a, function(x) identical(x$minSize, 9L), logical(1))))
})

test_that("without escapeMinPx spreadCpp() gets no minSize, as before", {
  skip_if_not_installed("SpaDES.tools")
  a <- spreadCppArgs()
  expect_false(any(vapply(a, function(x) "minSize" %in% names(x), logical(1))))
})

test_that("observed fires below escapeSizeHa are not fitted", {
  ## .objfunSpreadFit() hands each year's fires to objFunInner(); capture what it hands over
  handed <- NULL
  local_mocked_bindings(
    objFunInner = function(annualFires, ...) {
      handed <<- rbind(handed, annualFires)
      list(SNLL_FS = 0)
    },
    .package = "fireSenseUtils"
  )
  land <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 4800, ymin = 0, ymax = 4800) # 5.76-ha pixels
  fires <- list(year2004 = data.frame(size = c(3L, 8L, 9L, 40L), cells = c(1L, 2L, 3L, 4L), ids = 1:4))
  run <- function(...) {
    handed <<- NULL
    try(fireSenseUtils::.objfunSpreadFit(
      par = c(0.27, 1, 1, 1), landscape = land,
      annualDTx1000 = list(year2004 = data.table::data.table(pixelID = 1:400, cov = 0L)),
      nonAnnualDTx1000 = list(`year2004` = data.table::data.table(pixelID = 1:400)),
      formulaToFit = "~ 0 + cov", historicalFires = fires,
      fireBufferedListDT = list(year2004 = data.table::data.table(ids = 1:4, pixelID = 1:4, buffer = 1L)),
      tests = "snll_fs", Nreps = 2L, doAssertions = FALSE, thresh = Inf, verbose = 0, ...), silent = TRUE)
    sort(handed$size)
  }
  expect_identical(run(), c(3L, 8L, 9L, 40L))                 # minFireSize = 2
  expect_identical(run(escapeSizeHa = 50), c(9L, 40L))       # 9 pixels = 51.8 ha
})
