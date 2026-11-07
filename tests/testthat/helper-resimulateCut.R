## Shared by test-resimulateCutFires.R: objFunInner() on a 60 x 60 landscape with no values (every cell flammable),
## one fire ignited at the centre of an 11 x 11 buffer. Only the spreadProb surface is faked (a flat `p`); spreadCpp()
## is the real one. With `whole` the whole-landscape tables are given, and cut fires are re-simulated on them.

resimLandscape <- terra::rast(nrows = 60, ncols = 60, xmin = 0, xmax = 60, ymin = 0, ymax = 60)
resimCell <- function(row, col) (row - 1L) * 60L + col

## the fixture: a buffer of `side` x `side` cells around cell (30, 30), or the whole landscape if side = 60
resimFixture <- function(side = 11L) {
  rows <- max(30L - side %/% 2L, 1L) + seq_len(side) - 1L
  list(buf = data.table::data.table(ids = 1L, pixelID = as.vector(outer(rows, rows, resimCell)), buffer = 0L),
       fires = data.table::data.table(cells = resimCell(30L, 30L), size = 20L, ids = 1L))
}

resimWholeTables <- function() {
  list(ann = data.table::data.table(pixelID = seq_len(3600L)),
       nonAnn = stats::setNames(list(data.table::data.table(pixelID = seq_len(3600L))), "2003"))
}

resimObjective <- function(p = 0.2, side = 11L, whole = FALSE, seed = 1, Nreps = 20L, fx = resimFixture(side),
                           env = parent.frame(), ...) {
  pix <- fx$buf$pixelID
  local_mocked_bindings(
    paramsSeparate = function(...) list(logisticPars = c(0.27, 1, 1), covPars = 1),
    spreadProbFromIntegerCovs = function(shortAnnDTx1000 = NULL, annDTx1000, ...)
      data.table::data.table(pixelID = annDTx1000$pixelID, cov = 0),
    logisticAll = function(logisticPars, mat, ...) rep(p, nrow(mat)),
    .package = "fireSenseUtils", .env = env)
  wt <- if (isTRUE(whole)) resimWholeTables()
  args <- list(
    yr = "2004", annDTx1000 = data.table::data.table(pixelID = pix), par = 1, parsModel = NULL,
    annualFires = fx$fires, nonAnnualDTx1000 = list(`2003` = data.table::data.table(pixelID = pix)),
    annualFireBufferedDT = data.table::copy(fx$buf), indexNonAnnual = nonAnnualIndex(list(`2003` = 1)),
    colsToUse = "cov", covMinMax = NULL, mutuallyExclusive = NULL, doAssertions = FALSE, maxFireSpread = 1,
    lowerSpreadProb = 0.05, cells = numeric(3600), lanscape1stQuantileThresh = 1,
    r = resimLandscape, Nreps = Nreps, doMADTest = TRUE, doADTest = TRUE, doYearArea = TRUE, doCrps = TRUE,
    plot.it = FALSE, verbose = 0, gateRefusal = "simulate")
  if (isTRUE(whole)) args <- c(args, list(annDTx1000Whole = wt$ann, nonAnnualDTx1000Whole = wt$nonAnn))
  set.seed(seed)
  do.call(fireSenseUtils:::objFunInner, utils::modifyList(args, list(...)))
}
