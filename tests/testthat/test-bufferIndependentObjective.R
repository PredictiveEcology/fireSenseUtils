## The spread objective scored a fire cut at its buffer at the size it had reached, which made the fit
## depend on the buffers (a biased year effect and area). Now each fire's observed size is scored by a
## kernel log score of its simulated sizes, a replicate that reached an open buffer edge (one next to
## flammable land outside every buffer) contributing nothing: kernelLogScore(), bufferOpen().
## The fixture is from the known-truth prototype (~/claudeSessions/2026-10-06-multiELF-predict/phase4/
## bufferCensor/tune, 2026-10-23): ELF 6.6.1, one dataset, its truth, production buffers, 50 replicates.

fixture <- readRDS(test_path("fixtures", "bufferIndependentObjective.rds"))
h <- 0.25
shift <- log(dnorm(0) / h) # the largest log density: the score of a perfect prediction is 0

test_that("kernelLogScore() is 0 for a perfect prediction, grows with the miss, and is finite and >= 0 when all are cut", {
  expect_equal(fireSenseUtils:::kernelLogScore(100, rep(100, 50), rep(FALSE, 50)), 0)
  near <- fireSenseUtils:::kernelLogScore(100, rep(110, 50), rep(FALSE, 50))
  far <- fireSenseUtils:::kernelLogScore(100, rep(200, 50), rep(FALSE, 50))
  expect_gt(near, 0)
  expect_gt(far, near)
  allCut <- fireSenseUtils:::kernelLogScore(100, rep(100, 50), rep(TRUE, 50))
  expect_true(is.finite(allCut))
  expect_equal(allCut, shift - log(1e-29))
  ## a cut replicate adds no density: cutting the ones at the observed size costs, cutting far ones does not
  sims <- c(rep(100, 25), rep(5000, 25))
  none <- fireSenseUtils:::kernelLogScore(100, sims, rep(FALSE, 50))
  expect_equal(fireSenseUtils:::kernelLogScore(100, sims, sims == 5000), none)
  expect_gt(fireSenseUtils:::kernelLogScore(100, sims, sims == 100), none)
  ## never negative, also where the Monte Carlo correction would take log(m) past its bound
  set.seed(1)
  s <- replicate(200, fireSenseUtils:::kernelLogScore(100, 100 * exp(rnorm(3, 0, 0.05)), rep(FALSE, 3)))
  expect_true(all(s >= 0))
})

test_that("kernelLogScore() is finite when the density is positive but its square underflows", {
  ## ELF 4.1, 2026-10-24: a 61577-pixel fire whose replicates burned 9-54 pixels. Every kernel value is
  ## below 1e-162, so k^2 and m^2 underflow to 0 and the Monte Carlo correction was 0 / 0: the score was
  ## NaN, and the objective's early-stop `if (over)` stopped the first DEoptim generation.
  sims <- c(11, 14, 14, 32, 10, 20, 9, 23, 9, 9, 12, 9, 12, 11, 54, 24, 10, 12, 18, 10, 34, 27, 31, 11, 10,
            14, 17, 9, 11, 25, 12, 10, 50, 10, 24, 32, 23, 32, 11, 10, 9, 11, 9, 20, 9, 10, 14, 13, 34, 22)
  s <- fireSenseUtils:::kernelLogScore(61577, sims, rep(FALSE, 50))
  expect_true(is.finite(s))
  ## m is far below minLik: the floor's score, less the correction (at most 1)
  expect_lte(s, shift - log(1e-29))
  expect_gte(s, shift - log(1e-29) - 1)
})

test_that("the objective's per-fire term equals the prototype's on its saved simulations", {
  ## objFunInner() with spreadCpp() replaced by the prototype's simulated fires: each fire gets its own run
  ## of cells, the first of which is its buffer's only open pixel; a replicate the prototype cut burns it.
  sc <- fixture$scores
  total <- 0
  for (y in unique(sc$obs$yr)) {
    o <- sc$obs[sc$obs$yr == y, ]
    s <- sc$sims[sc$sims$yr == y, ]
    len <- vapply(o$ids, function(i) max(s$sim[s$ids == i]) + 1L, numeric(1))
    first <- cumsum(c(1, head(len, -1)))
    cellsOf <- function(k) first[k] + seq_len(len[k]) - 1
    buf <- data.table::rbindlist(lapply(seq_along(o$ids), function(k)
      data.table::data.table(ids = o$ids[k], pixelID = cellsOf(k), buffer = 0L))) # pixel ids as numeric
    buf[, `:=`(edge = pixelID %in% first, open = pixelID %in% first)]
    fires <- data.table::data.table(cells = first + 1, size = o$size, ids = o$ids)
    r <- terra::rast(nrows = ceiling(sum(len) / 1000), ncols = 1000)
    nCall <- 0L
    local_mocked_bindings(
      cropToCells = function(r, cells, ...) list(r = r, ncell = terra::ncell(r), toCrop = identity, toFull = identity),
      paramsSeparate = function(...) list(logisticPars = c(0.27, 1, 1), covPars = 1),
      spreadProbFromIntegerCovs = function(...) data.table::data.table(pixelID = buf$pixelID, cov = 0),
      logisticAll = function(...) seq(0.06, 0.08, length.out = nrow(buf)),
      .package = "fireSenseUtils")
    local_mocked_bindings(
      spreadCpp = function(landscape, loci, ...) {
        nCall <<- nCall + 1L
        data.table::rbindlist(lapply(seq_along(o$ids), function(k) {
          z <- s[s$ids == o$ids[k] & s$rep == nCall, ]
          cells <- if (z$cut) cellsOf(k)[seq_len(z$sim)] else cellsOf(k)[1L + seq_len(z$sim)]
          data.table::data.table(initialLocus = loci[k], indices = cells)
        }))
      },
      .package = "SpaDES.tools")
    out <- fireSenseUtils:::objFunInner(
      yr = sub("year", "", y), annDTx1000 = NULL, par = 1, parsModel = NULL,
      annualFires = fires, nonAnnualDTx1000 = NULL, annualFireBufferedDT = buf,
      indexNonAnnual = NULL, colsToUse = "cov", covMinMax = NULL, mutuallyExclusive = NULL, doAssertions = FALSE,
      maxFireSpread = 1, lowerSpreadProb = 0.05, cells = numeric(terra::ncell(r)), lanscape1stQuantileThresh = 1,
      weighted = FALSE, r = r, Nreps = sc$Nreps, doSNLL_FSTest = TRUE, doMADTest = FALSE, plot.it = FALSE,
      verbose = 0)
    expect_equal(out$runaways, sum(s$cut))
    total <- total + out$SNLL_FS
  }
  ## the prototype's value is unshifted: -sum(log(m) + correction)
  expect_equal(total, sc$prototypeF + nrow(sc$obs) * shift, tolerance = 1e-10)
})

test_that("bufferOpen() reproduces the prototype's open edge pixels on a real year of production buffers", {
  op <- fixture$open
  r <- terra::rast(nrows = op$nrow, ncols = op$ncol, xmin = 0, xmax = op$ncol, ymin = 0, ymax = op$nrow)
  r[op$nonNA] <- 1
  made <- addBufferEdge(list(op$bufs[, c("pixelID", "ids")]), r)[[1]]
  expect_identical(made$edge, op$bufs$edge)
  expect_identical(made$open, op$bufs$open)
  expect_gt(sum(op$bufs$edge & !op$bufs$open), 1000) # the edge pixels that are not open: next to NA or another buffer
})

## Three buffers on a 5 x 10 landscape whose first column is non-flammable (NA):
## buffer 1, columns 2-4: its left side faces NA, its right side buffer 2, top and bottom the boundary
## buffer 2, columns 5-7: its left side faces buffer 1, its right side column 8, flammable and in no buffer
## buffer 3, columns 9-10: its left side faces column 8, its right side the boundary
r3 <- terra::rast(nrows = 5, ncols = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 5, vals = 1)
r3[seq(1, 41, by = 10)] <- NA
cell3 <- function(rows, cols) as.vector(outer(rows, cols, function(r, c) (r - 1) * 10 + c))
toy3 <- data.table::data.table(ids = rep(1:3, c(15, 15, 10)), pixelID = c(cell3(1:5, 2:4), cell3(1:5, 5:7), cell3(1:5, 9:10)))

test_that("open edge pixels: not next to NA, nor next to another buffer; next to free flammable land", {
  open <- bufferOpen(toy3, r3)
  expect_false(any(open[toy3$ids == 1L]))                                  # NA, a shared border, the boundary
  expect_setequal(toy3$pixelID[open & toy3$ids == 2L], cell3(1:5, 7))     # free border only, not the shared one
  expect_setequal(toy3$pixelID[open & toy3$ids == 3L], cell3(1:5, 9))     # not the raster boundary
  expect_true(all(bufferEdge(toy3, r3)[open]))                             # open pixels are edge pixels
  ## a raster without values: every cell is flammable, so buffer 1's side facing column 1 is open
  noValues <- terra::rast(r3)
  expect_setequal(toy3$pixelID[bufferOpen(toy3, noValues) & toy3$ids == 1L], cell3(1:5, 2))
  made <- addBufferEdge(list(year2001 = toy3), r3)$year2001
  expect_identical(made$open, open)
  expect_identical(addBufferEdge(list(made), r3)[[1]], made)               # a table with both columns is left as is
})

## objFunInner() on the toy landscape with spreadCpp() replaced: `burn` lists, per fire, the cells it burns
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
  out <- fireSenseUtils:::objFunInner(
    yr = "2001", annDTx1000 = NULL, par = 1, parsModel = NULL, annualFires = fires, nonAnnualDTx1000 = NULL,
    annualFireBufferedDT = addBufferEdge(list(toy3), r3)[[1]], indexNonAnnual = NULL, colsToUse = "cov",
    covMinMax = NULL, mutuallyExclusive = NULL, doAssertions = FALSE, maxFireSpread = 1, lowerSpreadProb = 0.05,
    cells = numeric(50), lanscape1stQuantileThresh = 1, weighted = FALSE, r = r3, Nreps = Nreps,
    doSNLL_FSTest = TRUE, doMADTest = FALSE, plot.it = FALSE, verbose = 0, ...)
  c(out, list(spreadArgs = seen$args))
}

test_that("a replicate is cut when it burns an open pixel of any buffer of its year, and only then", {
  skip_if_not_installed("SpaDES.tools")
  quiet <- list(cell3(2:4, 2:4), cell3(2:4, 5:6), cell3(2:4, 10))     # NA side, shared border, boundary: not cut
  expect_equal(toyObjective(quiet)$runaways, 0)
  expect_equal(toyObjective(list(cell3(2:4, 2:4), cell3(2:4, 5:7), cell3(3, 10)))$runaways, 2)  # fire 2 x 2 replicates
  ## fire 1 runs through buffer 2 to its open column: cut, though its own buffer has no open pixel
  expect_equal(toyObjective(list(cell3(3, 2:7), cell3(3, 6), cell3(3, 10)))$runaways, 2)
  ## a cut fire scores as one with no density at its observed size
  cutAll <- toyObjective(list(cell3(1:5, 2:7), cell3(1:5, 5:7), cell3(1:5, 9:10)))
  expect_equal(cutAll$SNLL_FS, 3 * fireSenseUtils:::kernelLogScore(1, 1, cut = TRUE))
})

test_that("stopCutFires stops a fire at the first open pixel of its own buffer", {
  skip_if_not_installed("SpaDES.tools")
  quiet <- list(cell3(2:4, 2:4), cell3(2:4, 5:6), cell3(2:4, 10))
  on <- toyObjective(quiet, stopCutFires = TRUE)$spreadArgs
  expect_equal(on$stopAt, 1L)
  expect_setequal(on$stopCells[on$stopEvent == 2L], cell3(1:5, 7))
  expect_setequal(on$stopCells[on$stopEvent == 3L], cell3(1:5, 9))
  expect_false(1L %in% on$stopEvent)                                    # buffer 1 has no open pixel
  off <- toyObjective(quiet, stopCutFires = FALSE)$spreadArgs
  expect_null(off$stopCells)
})

test_that("the first block's early stop is per fire, not per year", {
  ## two years in the first block with 3 and 1 fires: 400 in all is 100 per fire (it was 200 per year)
  mk <- function(thresh, pruneAbove = Inf) {
    local_mocked_bindings(
      objFunInner = function(yr, ...) list(SNLL_FS = c("2001" = 300, "2002" = 100, "2003" = 7)[[yr]],
                                           runaways = 0L, nSims = 0L, gated = FALSE),
      .package = "fireSenseUtils")
    land <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20, vals = 1)
    yrs <- c("2001", "2002", "2003")
    dt <- function(...) data.table::data.table(...)
    fireSenseUtils::.objfunSpreadFit(
      par = c(0.27, 1, 1, 1), landscape = land,
      annualDTx1000 = stats::setNames(lapply(yrs, function(y) dt(pixelID = 1:2, cov = 100L)), yrs),
      nonAnnualDTx1000 = list(`2001_2003` = dt(pixelID = 1:2)), formulaToFit = "~ 0 + cov",
      historicalFires = list(`2001` = data.frame(size = c(100, 200, 300), cells = 1:3, ids = 1:3),
                             `2002` = data.frame(size = 900, cells = 4, ids = 1),
                             `2003` = data.frame(size = c(5, 6), cells = 5:6, ids = 1:2)),
      fireBufferedListDT = stats::setNames(lapply(yrs, function(y) dt(pixelID = 1:6, ids = 1L)), yrs),
      tests = "snll_fs", Nreps = 1L, doAssertions = FALSE, verbose = 0, thresh = thresh, pruneAbove = pruneAbove,
      returnTerms = TRUE)
  }
  ok <- mk(Inf)
  expect_equal(ok[["firstBlockPerFire"]], 100)
  expect_equal(ok[["objective"]], 407)
  expect_equal(ok[["bailed"]], 0)
  expect_equal(mk(150)[["bailed"]], 0)       # 100 per fire: under 150 (200 per year was not)
  over <- mk(90)
  expect_equal(over[["bailed"]], 1)
  expect_equal(over[["objective"]], 1e6)
  expect_equal(mk(Inf, pruneAbove = 399)[["objective"]], 1e6) # pruneAbove keeps the total's unit
})

test_that("the removed terms' arguments warn once and change nothing", {
  local_mocked_bindings(objFunInner = function(...) list(SNLL_FS = 5, runaways = 0L, nSims = 0L, gated = FALSE),
                        .package = "fireSenseUtils")
  rm(list = ls(fireSenseUtils:::.deprecatedArgsWarned), envir = fireSenseUtils:::.deprecatedArgsWarned)
  land <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20, vals = 1)
  call1 <- function(...) fireSenseUtils::.objfunSpreadFit(
    par = c(0.27, 1, 1, 1), landscape = land,
    annualDTx1000 = list(year2004 = data.table::data.table(pixelID = 1:400, cov = 0L)),
    nonAnnualDTx1000 = list(`year2004` = data.table::data.table(pixelID = 1:400)), formulaToFit = "~ 0 + cov",
    historicalFires = list(year2004 = data.frame(size = 40L, cells = 1L, ids = 1L)),
    fireBufferedListDT = list(year2004 = data.table::data.table(ids = 1L, pixelID = 1:400)),
    doAssertions = FALSE, verbose = 0, ...)
  plain <- call1()
  expect_warning(withAD <- call1(tests = c("adTest", "SNLL_FS"), yearAreaWeight = "auto", areaDistWeight = "auto",
                                 runawayBufferMultiple = 2), "deprecated and ignored")
  expect_equal(withAD, plain)
  expect_no_warning(call1(tests = c("adTest", "SNLL_FS"), yearAreaWeight = "auto"))  # once per session
  expect_no_warning(call1(yearAreaWeight = 0, areaDistWeight = 0L))
})
