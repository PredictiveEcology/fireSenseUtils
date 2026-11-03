## A simulated fire that reaches an open edge pixel of the year's buffers (one next to flammable land outside
## every buffer) is cut: bufferOpen(), addBufferEdge(). The per-fire CRPS scores it as at least the observed
## size; the size-based terms (adTest, annual area, area-weighted size distribution) score the size it reached.
## The fixture is from the known-truth prototype (~/claudeSessions/2026-10-06-multiELF-predict/phase4/
## bufferCensor/tune, 2026-10-23): ELF 6.6.1, one dataset, its truth, production buffers, 50 replicates.

fixture <- readRDS(test_path("fixtures", "bufferIndependentObjective.rds"))

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
toyObjective <- function(burn, Nreps = 2L, buffers = toy3, ...) {
  seen <- new.env()
  fires <- data.table::data.table(cells = cell3(3, c(3, 6, 10)), size = c(10L, 10L, 5L), ids = 1:3)
  local_mocked_bindings(
    cropToCells = function(r, cells, ...) list(r = r, ncell = terra::ncell(r), toCrop = identity, toFull = identity),
    paramsSeparate = function(...) list(logisticPars = c(0.27, 1, 1), covPars = 1),
    spreadProbFromIntegerCovs = function(...) data.table::data.table(pixelID = buffers$pixelID, cov = 0),
    logisticAll = function(...) seq(0.06, 0.08, length.out = nrow(buffers)),
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
         annualFireBufferedDT = addBufferEdge(list(buffers), r3)[[1]], indexNonAnnual = NULL, colsToUse = "cov",
         covMinMax = NULL, mutuallyExclusive = NULL, doAssertions = FALSE, maxFireSpread = 1, lowerSpreadProb = 0.05,
         cells = numeric(50), lanscape1stQuantileThresh = 1, r = r3, Nreps = Nreps,
         doMADTest = FALSE, doADTest = FALSE, doCrps = TRUE, plot.it = FALSE, verbose = 0),
    list(...))
  out <- do.call(fireSenseUtils:::objFunInner, args)
  c(out, list(spreadArgs = seen$args))
}

test_that("a replicate is cut when it burns an open pixel of any buffer of its year, and only then", {
  skip_if_not_installed("SpaDES.tools")
  ## observed sizes 10, 10 and 5; a cut replicate is scored as at least the observed size, so one that
  ## stopped short of it adds 0 to the per-fire CRPS instead of log(observed / burned)
  crpsOf <- function(burn) toyObjective(burn)$crpsFire
  quiet <- list(cell3(2:4, 2:4), cell3(2:4, 5:6), cell3(2:4, 10))     # NA side, shared border, boundary: not cut
  expect_equal(crpsOf(quiet), log(10 / 9) + log(10 / 6) + log(5 / 3))
  expect_equal(crpsOf(list(cell3(2:4, 2:4), cell3(2:4, 5:7), cell3(3, 10))), log(10 / 9) + 0 + log(5))
  ## fire 1 runs through buffer 2 to its open column: cut, though its own buffer has no open pixel
  expect_equal(crpsOf(list(cell3(3, 2:7), cell3(3, 6), cell3(3, 10))), 0 + log(10) + log(5))
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

test_that("a cut replicate is scored at the size it reached by the size-based terms", {
  skip_if_not_installed("SpaDES.tools")
  ## fire 1 burns 9 pixels inside its buffer; fire 2 burns 9 pixels including its open column 7: cut
  burn <- list(cell3(2:4, 2:4), cell3(2:4, 5:7), cell3(3, 10))
  out <- toyObjective(burn, doADTest = TRUE, doYearArea = TRUE, doMADTest = TRUE)
  expect_equal(out$allFireSizes, rep(9, 4))           # fire 3 burns 1 pixel and is left out; the cut ones count 9
  expect_equal(out$annualAreaByRep, c(18, 18))
  expect_equal(out$fireSizes, c(9, 9, 1))
})

test_that("the Anderson-Darling, annual-area and area-weighted terms add to the objective", {
  local_mocked_bindings(
    objFunInner = function(yr, ...) list(allFireSizes = c(20, 30, 40, 55), annualAreaByRep = c(145, 150)),
    .package = "fireSenseUtils")
  land <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20, vals = 1)
  fit <- function(...) fireSenseUtils::.objfunSpreadFit(
    par = c(0.27, 1, 1, 1), landscape = land,
    annualDTx1000 = list(`2004` = data.table::data.table(pixelID = 1:400, cov = 0L)),
    nonAnnualDTx1000 = list(`2004` = data.table::data.table(pixelID = 1:400)), formulaToFit = "~ 0 + cov",
    historicalFires = list(`2004` = data.frame(size = c(25, 60), cells = 1:2, ids = 1:2)),
    fireBufferedListDT = list(`2004` = data.table::data.table(ids = 1L, pixelID = 1:400)),
    tests = "adTest", doAssertions = FALSE, verbose = 0, returnTerms = TRUE, ...)
  base <- fit()
  expect_true(is.finite(base[["adTest"]]) && base[["adTest"]] > 0)
  expect_equal(base[["objective"]], base[["adTest"]])
  expect_true(is.na(base[["yearArea"]]) && is.na(base[["areaDist"]]))      # off by default (weight 0)
  all <- fit(yearAreaWeight = "auto", areaDistWeight = "auto")
  expect_true(is.finite(all[["yearArea"]]) && is.finite(all[["areaDist"]]))
  expect_equal(all[["objective"]], all[["adTest"]] + all[["yearArea"]] + all[["areaDist"]])
})

test_that("runawayBufferMultiple and runawaySize warn once and change nothing", {
  local_mocked_bindings(objFunInner = function(...) list(allFireSizes = c(2, 10, 300)), .package = "fireSenseUtils")
  rm(list = ls(fireSenseUtils:::.deprecatedArgsWarned), envir = fireSenseUtils:::.deprecatedArgsWarned)
  land <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20, vals = 1)
  call1 <- function(...) fireSenseUtils::.objfunSpreadFit(
    par = c(0.27, 1, 1, 1), landscape = land,
    annualDTx1000 = list(year2004 = data.table::data.table(pixelID = 1:400, cov = 0L)),
    nonAnnualDTx1000 = list(`year2004` = data.table::data.table(pixelID = 1:400)), formulaToFit = "~ 0 + cov",
    historicalFires = list(year2004 = data.frame(size = 40L, cells = 1L, ids = 1L)),
    fireBufferedListDT = list(year2004 = data.table::data.table(ids = 1L, pixelID = 1:400)),
    tests = "adTest", doAssertions = FALSE, verbose = 0, ...)
  plain <- call1()
  expect_warning(withMult <- call1(runawayBufferMultiple = 2), "deprecated and ignored")
  expect_equal(withMult, plain)
  expect_no_warning(call1(runawayBufferMultiple = 2))  # once per session
})

## The rank-based Anderson-Darling term (tests = "pitAD"): each fitted fire's observed size is ranked among its
## replicates (a mid-rank probability integral transform, u), and the u of all fires are scored by the
## Anderson-Darling statistic against Uniform(0, 1). A rank does not change when a larger buffer would have let
## a replicate that is already larger than the observed fire burn on, so the term does not depend on the buffer.

test_that("pitMidRank(): all replicates smaller gives 1, all larger gives 0, ties count half", {
  pit <- fireSenseUtils:::pitMidRank
  no <- rep(FALSE, 4)
  expect_equal(pit(10, c(1, 2, 3, 4), no, 4L), 1)
  expect_equal(pit(10, c(11, 12, 13, 14), no, 4L), 0)
  expect_equal(pit(10, c(10, 10, 10, 10), no, 4L), 0.5)
  expect_equal(pit(10, c(5, 10, 20, 30), no, 4L), (1 + 0.5) / 4)
  ## a replicate with no row (it burned nothing) is smaller than any observed fire
  expect_equal(pit(10, c(20, 30), c(FALSE, FALSE), 4L), 2 / 4)
  expect_equal(pit(10, c(5, 20), c(FALSE, FALSE), 4L), 3 / 4)
})

test_that("pitMidRank(): a cut replicate is larger if it reached the observed size, half a count if it did not", {
  pit <- fireSenseUtils:::pitMidRank
  ## reached >= observed: the true size is at least that, so it is larger (a tie included)
  expect_equal(pit(10, c(10, 50, 5, 5), c(TRUE, TRUE, FALSE, FALSE), 4L), 2 / 4)
  ## reached < observed: its order is unknown, counted 0.5 (not 0 or 1)
  expect_equal(pit(10, c(3, 5, 5, 5), c(TRUE, FALSE, FALSE, FALSE), 4L), (3 + 0.5) / 4)
  expect_equal(pit(10, c(3, 20), c(TRUE, FALSE), 2L), 0.5 / 2)
})

test_that("adUniform(): small for evenly spread u, large for u piled in the middle or at an end", {
  n <- 100
  even <- (seq_len(n) - 0.5) / n
  expect_lt(fireSenseUtils:::adUniform(even, 50L), 0.1)
  expect_gt(fireSenseUtils:::adUniform(rep(0.5, n), 50L), 30)  # exactly 2 * log(2) * n - n = 38.6
  expect_gt(fireSenseUtils:::adUniform(rep(0, n), 50L), 50)       # finite: u is kept off 0 and 1
  expect_true(is.finite(fireSenseUtils:::adUniform(rep(1, n), 50L)))
  ## the mean of the statistic for u ~ Uniform is 1
  set.seed(3)
  expect_equal(mean(replicate(2000, fireSenseUtils:::adUniform(runif(50), 50L))), 1, tolerance = 0.1)
})

test_that("objFunInner() returns one PIT per fire, scoring a cut replicate as larger once it reached the observed size", {
  skip_if_not_installed("SpaDES.tools")
  ## fire 1 (observed 10): 9 pixels both replicates, inside its buffer: both smaller, u = 1
  ## fire 2 (observed 10): 9 pixels, the open column 7 included: cut, reached < 10: 0.5 each, u = 0.5
  ## fire 3 (observed 5): 1 pixel: smaller, u = 1
  out <- toyObjective(list(cell3(2:4, 2:4), cell3(2:4, 5:7), cell3(3, 10)), doPitADTest = TRUE)
  expect_equal(out$pit, c(1, 0.5, 1))
  ## a fire of 20 pixels that burns the open column is cut at a size above the observed one: larger, u = 0
  big <- toyObjective(list(cell3(2:4, 2:4), cell3(1:5, 5:7), cell3(3, 10)), doPitADTest = TRUE)
  expect_equal(big$pit[2], 0) # 15 pixels >= 10
})

test_that("the PITs do not change when the buffers grow and no replicate was cut below the observed size", {
  skip_if_not_installed("SpaDES.tools")
  ## fire 2 burns columns 5-7 (15 pixels, observed 10): cut by the small buffers (column 7 is open), exact with
  ## buffers that extend to column 8. Same simulations, same u.
  burn <- list(cell3(2:4, 2:4), cell3(1:5, 5:7), cell3(3, 10))
  small <- toyObjective(burn, doPitADTest = TRUE)
  grown <- rbind(toy3, data.table::data.table(ids = 2L, pixelID = cell3(1:5, 8)))
  withBig <- toyObjective(burn, doPitADTest = TRUE, buffers = grown)
  expect_equal(withBig$pit, small$pit)
})

test_that("tests = \"pitAD\" adds a finite, weighted term, beside or instead of the pooled adTest", {
  local_mocked_bindings(
    objFunInner = function(yr, ...) list(allFireSizes = c(20, 30, 40, 55), pit = c(0.2, 0.9)),
    .package = "fireSenseUtils")
  land <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20, vals = 1)
  fit <- function(tests, ...) fireSenseUtils::.objfunSpreadFit(
    par = c(0.27, 1, 1, 1), landscape = land,
    annualDTx1000 = list(`2004` = data.table::data.table(pixelID = 1:400, cov = 0L)),
    nonAnnualDTx1000 = list(`2004` = data.table::data.table(pixelID = 1:400)), formulaToFit = "~ 0 + cov",
    historicalFires = list(`2004` = data.frame(size = c(25, 60), cells = 1:2, ids = 1:2)),
    fireBufferedListDT = list(`2004` = data.table::data.table(ids = 1L, pixelID = 1:400)),
    tests = tests, Nreps = 10L, doAssertions = FALSE, verbose = 0, returnTerms = TRUE, ...)
  p <- fit("pitAD")
  expect_true(is.finite(p[["pitAD"]]) && p[["pitAD"]] > 0)
  expect_true(is.na(p[["adTest"]]))
  expect_equal(p[["objective"]], p[["pitAD"]])
  expect_equal(p[["pitAD"]], fireSenseUtils:::adUniform(c(0.2, 0.9), 10L) * fireSenseUtils::adWeightAuto(2, "kde", FALSE))
  expect_equal(fit("pitAD", adWeight = 7)[["pitAD"]], fireSenseUtils:::adUniform(c(0.2, 0.9), 10L) * 7)
  both <- fit(c("adTest", "pitAD"))
  expect_equal(both[["objective"]], both[["adTest"]] + both[["pitAD"]])
  expect_true(is.na(fit("adTest")[["pitAD"]]))
})

test_that("tests = \"crps\": objFunInner() returns the per-fire CRPS and every simulated fire", {
  skip_if_not_installed("SpaDES.tools")
  ## fires burn 9, 6 and 1 pixels in both replicates; observed 10, 10 and 5
  burn <- list(cell3(2:4, 2:4), cell3(2:4, 5:6), cell3(3, 10))
  out <- toyObjective(burn)
  expect_equal(out$crpsFire, log(10 / 9) + log(10 / 6) + log(5))
  expect_equal(sort(out$crpsSizes$N), sort(rep(c(9, 6, 1), 2)))
  expect_equal(sort(unique(out$crpsSizes$rep)), 1:2)
  ## fire 2 burns its open column 7 (cut) and stops at 4 pixels of a 10-pixel observed fire: scored as 10
  cutShort <- toyObjective(list(cell3(2:4, 2:4), cell3(3, 5:7), cell3(3, 10)))
  expect_equal(cutShort$crpsFire, log(10 / 9) + 0 + log(5))
})

test_that("gateRefusal = \"simulate\" scores a year the gates would refuse", {
  skip_if_not_installed("SpaDES.tools")
  burn <- list(cell3(2:4, 2:4), cell3(2:4, 5:6), cell3(3, 10))
  ## the first quartile of the spread probabilities (0.06-0.08) is above 0.01: "too burny", refused
  refused <- toyObjective(burn, lanscape1stQuantileThresh = 0.01)
  expect_true(all(refused$crpsSizes$N == 1L))
  expect_equal(refused$crpsFire, sum(log(c(10, 10, 5))))                 # every fire left at one pixel
  simulated <- toyObjective(burn, lanscape1stQuantileThresh = 0.01, gateRefusal = "simulate")
  expect_equal(simulated$crpsFire, log(10 / 9) + log(10 / 6) + log(5))
  expect_lt(simulated$crpsFire, refused$crpsFire)
})
