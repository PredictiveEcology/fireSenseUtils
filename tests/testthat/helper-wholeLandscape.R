## Shared by test-objFunInner-wholeLandscape.R: objFunInner() on a 60 x 60 landscape with no values (every cell
## flammable), one fire ignited at the centre of a `side` x `side` buffer. Only the covariates are faked; the
## spreadProb link and spreadCpp() are the real ones. With `whole` the tables cover the whole landscape and
## `wholeLandscape = TRUE`: the fire spreads on all of it.

wlLandscape <- terra::rast(nrows = 60, ncols = 60, xmin = 0, xmax = 60, ymin = 0, ymax = 60)
wlCell <- function(row, col) (row - 1L) * 60L + col

## the fixture: a buffer of `side` x `side` cells around cell (30, 30), or the whole landscape if side = 60
wlFixture <- function(side = 11L) {
  rows <- max(30L - side %/% 2L, 1L) + seq_len(side) - 1L
  list(buf = data.table::data.table(ids = 1L, pixelID = as.vector(outer(rows, rows, wlCell)), buffer = 0L),
       fires = data.table::data.table(cells = wlCell(30L, 30L), size = 20L, ids = 1L))
}

wlTables <- function(pix) {
  list(ann = data.table::data.table(pixelID = pix),
       nonAnn = stats::setNames(list(data.table::data.table(pixelID = pix)), "2003"))
}

## `p = NULL`: one covariate that varies from pixel to pixel, through the real "powerLaw" link (pc 0.4, gamma 1.9).
## A number `p`: a flat spreadProb `p` and no link, as before the whole-landscape tables existed (the fixture).
## `annPix`: the pixels of the fit's tables (default: the buffer's; with `whole`, every pixel).
wlObjective <- function(p = NULL, side = 11L, whole = FALSE, seed = 1, Nreps = 20L, fx = wlFixture(side),
                        annPix = if (whole) seq_len(3600L) else fx$buf$pixelID, env = parent.frame(), ...) {
  if (is.null(p)) {
    local_mocked_bindings(
      paramsSeparate = function(...) list(logisticPars = c(0.4, 1, 1), covPars = 1),
      spreadProbFromIntegerCovs = function(shortAnnDTx1000 = NULL, annDTx1000, ...)
        data.table::data.table(pixelID = annDTx1000$pixelID, cov = 1.6 + 1.2 * sin(annDTx1000$pixelID)),
      .package = "fireSenseUtils", .env = env)
  } else {
    local_mocked_bindings(
      paramsSeparate = function(...) list(logisticPars = c(0.27, 1, 1), covPars = 1),
      spreadProbFromIntegerCovs = function(shortAnnDTx1000 = NULL, annDTx1000, ...)
        data.table::data.table(pixelID = annDTx1000$pixelID, cov = 0),
      logisticAll = function(logisticPars, mat, ...) rep(p, nrow(mat)),
      .package = "fireSenseUtils", .env = env)
  }
  fit <- wlTables(annPix)
  args <- list(
    yr = "2004", annDTx1000 = fit$ann, par = 1, parsModel = NULL,
    annualFires = fx$fires, nonAnnualDTx1000 = fit$nonAnn,
    annualFireBufferedDT = data.table::copy(fx$buf), indexNonAnnual = nonAnnualIndex(list(`2003` = 1)),
    colsToUse = "cov", covMinMax = NULL, mutuallyExclusive = NULL, doAssertions = FALSE, maxFireSpread = 1,
    lowerSpreadProb = 0.05, cells = numeric(3600), lanscape1stQuantileThresh = 1,
    r = wlLandscape, Nreps = Nreps, doMADTest = TRUE, doADTest = TRUE, doYearArea = TRUE, doCrps = TRUE,
    plot.it = FALSE, verbose = 0, gateRefusal = "simulate")
  if (is.null(p)) args <- c(args, list(doPitADTest = TRUE, link = "powerLaw", powerLawGamma = 1.9))
  if (isTRUE(whole)) args <- c(args, list(wholeLandscape = TRUE))
  set.seed(seed)
  do.call(fireSenseUtils:::objFunInner, utils::modifyList(args, list(...)))
}
