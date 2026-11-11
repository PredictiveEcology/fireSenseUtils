## One definition of an escaped fire, and one set of fire starts, for every part of fireSense: the spread fit, the
## escape and ignition fits, the burn and the hindcast (Eliot, 2026-11-10). The buffer of an ELF only removes edge
## effects: a fire belongs to the ELF if and only if its start is in the core, and it then counts at its whole extent,
## core and buffer. Its start is the polygon centroid moved to the nearest burnable pixel of the fire, as the spread
## fit has always had it ([harmonizeBufferAndPoints()], [cleanUpSpreadFirePoints()]); never the NFDB point.

utils::globalVariables(".poly")

#' The definition of an escaped fire
#'
#' An escape is a burned-area polygon of one of `fireCauses` whose whole size is at least `escapeSizeHa` (and
#' larger than one pixel). The same definition must be used to fit the spread, ignition and escape models and to
#' burn and score with them; a fit stores it in its ledger row, and [checkEscapeDefinition()] compares a
#' ledger row with the settings of a later run.
#'
#' @param escapeSizeHa Size (ha) of the whole polygon at and above which a fire escaped. Default
#'   [fireSenseEscapeSizeHa]. `NA` is the definition of a fit made before one was recorded: it is stored as such and
#'   [checkEscapeDefinition()] does not check it.
#' @param fireCauses Cause groups of [fireCauseCodes()]. Default [fireSenseFireCauses].
#'
#' @return A one-row `data.table` with `escapeSizeHa` (numeric) and `fireCauses` (the groups, sorted and
#'   joined by `","`), the columns a ledger row stores.
#' @export
#' @seealso [escapeFireSet()], [checkEscapeDefinition()], [escapeDefinitionOf()]
escapeDefinition <- function(escapeSizeHa = fireSenseEscapeSizeHa, fireCauses = fireSenseFireCauses) {
  if (length(escapeSizeHa) != 1L || !(is.na(escapeSizeHa) || (is.numeric(escapeSizeHa) && escapeSizeHa >= 0)))
    stop("escapeSizeHa must be one number, 0 or more, or NA")
  fireCauseCodes(fireCauses) # stops on an unknown group
  data.table::data.table(escapeSizeHa = as.numeric(escapeSizeHa),
                         fireCauses = paste(sort(unique(fireCauses)), collapse = ","))
}

#' The escape definition of a SpaDES run
#'
#' [escapeDefinition()] of the current module's `escapeSizeHa` and `fireCauses` parameters, after
#' `SpaDES.core::paramCheckOtherMods()` has checked that every other module of `sim` that has them has the same
#' values (an error otherwise). Call it in an event of a module that defines both parameters.
#'
#' @param sim A `simList`.
#'
#' @return As [escapeDefinition()].
#' @export
#' @importFrom SpaDES.core paramCheckOtherMods
escapeDefinitionOf <- function(sim) {
  do.call(escapeDefinition, lapply(stats::setNames(nm = names(formals(escapeDefinition))), function(p)
    SpaDES.core::paramCheckOtherMods(sim, p, ifSetButDifferent = "error")))
}

## the cause groups of an escape definition
escapeCauses <- function(definition) strsplit(definition$fireCauses[1], ",", fixed = TRUE)[[1]]

#' Stop if a stored escape definition differs from the current one
#'
#' A fit's ledger row records the [escapeDefinition()] it was fitted with. Burning or scoring with another
#' definition would compare fires of one kind with a model of another, so this is an error.
#'
#' @param recorded Ledger rows (`data.frame` or list), e.g. one per ELF, with the columns of [escapeDefinition()].
#'   Rows without them, or with `escapeSizeHa` `NA` (fitted before the definition was recorded), are not checked.
#'   With a [polygonIDTxt] column, the message names the ELF of the row that differs.
#' @param current The definition in use, from [escapeDefinition()] or [escapeDefinitionOf()]. One with
#'   `escapeSizeHa` `NA` (a run with an old fit's definition) checks nothing.
#' @param what Where `recorded` comes from, for the message.
#'
#' @return `TRUE` invisibly if every row checked agrees, `NA` invisibly if none was checked; otherwise an error.
#' @export
checkEscapeDefinition <- function(recorded, current = escapeDefinition(), what = "the fit") {
  cols <- names(escapeDefinition())
  if (!all(cols %in% names(recorded)) || is.na(current$escapeSizeHa)) return(invisible(NA))
  sizes <- as.numeric(unlist(recorded[["escapeSizeHa"]]))
  causes <- as.character(unlist(recorded[["fireCauses"]]))
  ids <- if (polygonIDTxt %in% names(recorded)) as.character(unlist(recorded[[polygonIDTxt]]))
  rows <- which(!is.na(sizes))
  for (i in rows) {
    rec <- escapeDefinition(sizes[i], strsplit(causes[i], ",", fixed = TRUE)[[1]])
    if (!isTRUE(all.equal(rec$escapeSizeHa, current$escapeSizeHa)) || !identical(rec$fireCauses, current$fireCauses))
      stop("the escape definition of ", what, if (!is.null(ids)) paste0(" of ELF ", ids[i]), " (escapeSizeHa ",
           rec$escapeSizeHa, ", fireCauses ", rec$fireCauses, ") differs from this run's (escapeSizeHa ",
           current$escapeSizeHa, ", fireCauses ", current$fireCauses, ")", call. = FALSE)
  }
  invisible(if (length(rows)) TRUE else NA)
}

#' One start per fire polygon, and whether it is in the core
#'
#' `fireStarts()` gives every polygon of `polys` its start: the start in `starts` with the same `YEAR` and `idCol`
#' (the starts the spread fit was given, e.g. `spreadFirePoints` from [harmonizeFireData()]), or, with
#' `starts = NULL`, the polygon's centroid moved to the nearest burnable pixel of the fire by the functions
#' [harmonizeFireData()] uses ([harmonizeBufferAndPoints()], then [cleanUpSpreadFirePoints()]), year by year.
#' A polygon with no burnable pixel has no start and is left out. `escapeFireSet()` does the same for the
#' escape polygons of `definition` only.
#'
#' A fire's pixels are the spread fit's ([burnableFireHa()]): the cells of `flammableRTM` whose centre is in the
#' polygon, one fire per cell where polygons of a year overlap, less the `NA` cells; a polygon covering no cell centre
#' has none, so no start, as [harmonizeFireData()] drops it.
#'
#' @param polys Fire polygons (`SpatVector` or `sf`, or a list of them, one per year) with column `YEAR`.
#' @param flammableRTM `SpatRaster` of the burnable pixels (1; 0 or `NA` cannot burn) of the whole landscape,
#'   core and buffer. Its grid gives `pixelID`, the fires' pixels and the burnable area.
#' @param burnablePixels `NULL`, or the `pixelID`s the simulation can burn: those of the fire year's `landcoverDT`, as
#'   the spread fit counts them ([burnableFireHa()]). `NULL`: the cells of `flammableRTM` above 0.
#' @param coreRaster `SpatRaster` of the core, read by [inELFCore()] (with `coreIs`), e.g. `rasterToMatchLargeELF`
#'   (2 core, 1 buffer). `NULL` makes every start a core start.
#' @param starts `NULL`, or start points (`SpatVector` or `sf`, or a list of them) sharing `YEAR` and `idCol` with
#'   `polys`.
#' @param idCol The fire ID column of `polys` and `starts`; needed with `starts`, and copied to the result if given.
#' @param sizeCol Column of `polys` with the whole fire's size (ha); default [polySizeColTxt]. A `polys` without it
#'   is measured (`terra::expanse()`).
#' @param polyCauseCol Column of `polys` holding the cause; default [polyCauseColTxt].
#' @param coreIs What `coreRaster` is, as in [inELFCore()]: `"ELF"` (default; core [fireSenseELFCoreValue], buffer 1)
#'   or `"mask"`.
#'
#' @return A `SpatVector` of points in the CRS of `flammableRTM`, one per polygon with a start: `YEAR`, `.poly`
#'   (the polygon's row in `polys`), `idCol` (if given), `polyHa` (the whole polygon), `burnableHa` (its burnable
#'   pixels, core and buffer: the area the simulation can burn; [burnableFireHa()]), `cause` (if `polys` has `polyCauseCol`),
#'   `pixelID` (the start's cell of `flammableRTM`) and `core` (logical: the start is in the core).
#' @export
#' @importFrom data.table data.table
#' @importFrom terra cellFromXY cells centroids crds crs extract project res same.crs vect
fireStarts <- function(polys, flammableRTM, coreRaster = NULL, starts = NULL, idCol = NULL,
                       sizeCol = polySizeColTxt, polyCauseCol = polyCauseColTxt, coreIs = "ELF",
                       burnablePixels = NULL) {
  crsR <- terra::crs(flammableRTM)
  polys <- asOneVect(polys)
  if (!NROW(polys)) return(fireStartsOf(NULL, matrix(numeric(), ncol = 2), crsR, idCol, polyCauseCol))
  if (!terra::same.crs(polys, flammableRTM)) polys <- terra::project(polys, crsR)
  if (!".poly" %in% names(polys)) polys$.poly <- seq_len(NROW(polys))
  if (!".ha" %in% names(polys)) polys$.ha <- polySizeHa(polys, sizeCol)

  ## the pixels of each polygon, as harmonizeFireData() has them, year by year
  cl <- data.table::rbindlist(lapply(split(seq_len(NROW(polys)), as.integer(polys$YEAR)), function(i)
    removeBufferedFiresOutsideRTM(firePixels(polys[i, ".poly"], flammableRTM, ".poly"), flammableRTM)))
  data.table::setnames(cl, "ids", ".poly")
  if (is.null(burnablePixels)) {
    flam <- terra::extract(flammableRTM, cl$pixelID)[[1]]
    burnablePixels <- cl$pixelID[!is.na(flam) & flam > 0]
  }

  if (is.null(starts)) {
    xy <- snappedCentroids(polys, cl, flammableRTM)
    at <- match(polys$.poly, xy$.poly)
    xy <- as.matrix(xy[, c("x", "y")])
  } else {
    if (is.null(idCol)) stop("fireStarts(): `idCol` is needed to match `starts` to `polys`")
    starts <- terra::project(asOneVect(starts), crsR)
    at <- startOf(polys, starts, idCol)
    xy <- terra::crds(starts)
  }
  has <- !is.na(at)
  polys <- polys[has, ]
  xy <- xy[at[has], , drop = FALSE]
  pts <- fireStartsOf(polys, xy, crsR, idCol, polyCauseCol)
  burnableHa <- burnableFireHa(data.table::data.table(pixelID = cl$pixelID, ids = cl$.poly), burnablePixels,
                               prod(terra::res(flammableRTM)) / 1e4)
  burnableHa <- as.vector(burnableHa[as.character(polys$.poly)])
  burnableHa[is.na(burnableHa)] <- 0
  pts$burnableHa <- burnableHa
  pts$pixelID <- terra::cellFromXY(flammableRTM, xy)
  pts$core <- inELFCore(pts, coreRaster, coreIs)
  pts
}

#' @rdname fireStarts
#' @param definition The escape definition, from [escapeDefinition()].
#' @export
escapeFireSet <- function(polys, flammableRTM, coreRaster = NULL, definition = escapeDefinition(), starts = NULL,
                          idCol = NULL, sizeCol = polySizeColTxt, polyCauseCol = polyCauseColTxt, coreIs = "ELF",
                          burnablePixels = NULL) {
  polys <- escapesOf(polys, definition, flammableRTM, sizeCol, polyCauseCol)
  fireStarts(polys, flammableRTM, coreRaster, starts, idCol, sizeCol, polyCauseCol, coreIs, burnablePixels)
}

## the polygons of `polys` that are escapes of `definition`, with `.poly` their row in `polys` (as [fireStarts()]'s)
escapesOf <- function(polys, definition, flammableRTM, sizeCol, polyCauseCol) {
  escapePolys(polys, definition$escapeSizeHa, prod(terra::res(flammableRTM)) / 1e4, sizeCol,
              fireCauseCodes(escapeCauses(definition), "nbac"), polyCauseCol)
}

#' The ignitions the ignition and escape models are fitted to, and the hindcast counts
#'
#' The NFDB points of `definition`'s causes that start in the core, with the escapes of `fireSet`: the point
#' [escapedIgnitions()] matches to an escape polygon (in it or within `maxDistM`, one per polygon) is replaced by
#' that fire's start, so a fire is one ignition and one escape, at the start the spread fit uses. A core-start
#' escape polygon that no point matches gets an ignition at its start too, as [polygonIgnitions()] adds. Matched
#' points of escape polygons that start in the buffer are dropped with their fire.
#'
#' @param points NFDB ignition points (`SpatVector` or `sf`) with `YEAR`, `CAUSE` and `pointSizeCol`.
#' @param polys The fire polygons `fireSet` was made from.
#' @param fireSet The escape polygons' starts, from [escapeFireSet()] with the same `definition`, `polys` and
#'   `sizeCol`.
#' @param coreRaster,coreIs As in [escapeFireSet()]: the core the points must start in.
#' @param definition The escape definition, from [escapeDefinition()].
#' @param pixSizeHa Pixel area (ha), as in [escapedIgnitions()].
#' @param maxDistM As in [ignitionPolygonPairs()].
#' @param sizeCol,polyCauseCol As in [escapeFireSet()].
#' @param pointSizeCol The column of `points` that picks one point per polygon ([escapedIgnitions()]).
#'
#' @return A `SpatVector` of points in the CRS of `points`: the core points that are not an escape polygon's
#'   (`escaped` `FALSE`), then one per core-start escape of `fireSet` (`escaped` `TRUE`, `SIZE_HA` the polygon's
#'   whole size, `CAUSE` the NFDB code of its cause group).
#' @export
escapeIgnitionPoints <- function(points, polys, fireSet, coreRaster, definition = escapeDefinition(), pixSizeHa,
                                 coreIs = "ELF", maxDistM = 1000, sizeCol = polySizeColTxt,
                                 pointSizeCol = "SIZE_HA", polyCauseCol = polyCauseColTxt) {
  causes <- escapeCauses(definition)
  points <- asOneVect(points)
  points <- points[points$CAUSE %in% fireCauseCodes(causes, "nfdb"), ]
  m <- ignitionPolygonPairs(points, polys, definition$escapeSizeHa, pixSizeHa, maxDistM, sizeCol = sizeCol,
                            causes = fireCauseCodes(causes, "nbac"), polyCauseCol = polyCauseCol)
  points <- m$points[!pickEscapePoints(m, pointSizeCol), ]
  points <- points[inELFCore(points, coreRaster, coreIs), ]
  points$escaped <- rep(FALSE, NROW(points))
  esc <- fireSet[fireSet$core, ]
  esc <- terra::project(esc, terra::crs(points))
  added <- ignitionsAt(terra::crds(esc), terra::crs(points), esc$YEAR, nfdbCauseOf(esc$cause), esc$polyHa)
  added$escaped <- rep(TRUE, NROW(added))
  if (!NROW(points)) return(added)
  if (!NROW(added)) return(points)
  rbind(points, added)
}

#' Functions whose code affects what [escapeFireSet()] and [escapeIgnitionPoints()] return
#'
#' Pass as `reproducible::Cache()`'s `.cacheExtra` (see [harmonizeFireDataDeps()]).
#'
#' @return A list of functions.
#' @export
escapeFireSetDeps <- function() {
  list(escapeDefinition, escapeCauses, escapeFireSet, fireStarts, escapeIgnitionPoints, escapePolys, polySizeHa,
       isEscapeSize, asOneVect, startOf, snappedCentroids, fireStartsOf, inELFCore, escapesOf, ignitionsAt, nfdbCauseOf,
       ignitionPolygonPairs, pickEscapePoints, fireCauseCodes, harmonizeBufferAndPoints, cleanUpSpreadFirePoints,
       firePixels, removeBufferedFiresOutsideRTM, burnableFireHa)
}

## the whole fire's size (ha): `sizeCol` if `polys` has it, else the area given
polySizeHa <- function(polys, sizeCol) {
  if (!is.null(sizeCol) && sizeCol %in% names(polys)) as.numeric(polys[[sizeCol]][, 1]) else
    terra::expanse(polys, unit = "ha")
}

## each polygon's centroid moved to the nearest burnable pixel of the fire, as harmonizeFireData() does, year by
## year: `cl` holds the polygons' pixels (`.poly`, `pixelID`). A data.frame of `.poly`, `x`, `y`.
snappedCentroids <- function(polys, cl, flammableRTM) {
  out <- lapply(split(seq_len(NROW(polys)), as.integer(polys$YEAR)), function(i) {
    p <- polys[i, ]
    cent <- sf::st_as_sf(terra::centroids(p)[, ".poly"])
    buff <- cl[.poly %in% p$.poly, list(pixelID, ids = .poly, buffer = 1L)]
    h <- harmonizeBufferAndPoints(list(cent), list(buff), flammableRTM, idCol = ".poly")[[1]]
    if (is.null(h)) return(NULL)
    h <- cleanUpSpreadFirePoints(h, buff, flammableRTM, idCol = ".poly")$SpatialPoints
    if (!NROW(h)) return(NULL)
    data.frame(.poly = h$.poly, sf::st_coordinates(h)[, c("X", "Y"), drop = FALSE])
  })
  out <- do.call(rbind, c(list(data.frame(.poly = integer(), X = numeric(), Y = numeric())), out))
  data.frame(.poly = out$.poly, x = out$X, y = out$Y)
}

## the start points of `polys` at `xy`, with the polygons' YEAR, .poly, idCol, polyHa and cause
fireStartsOf <- function(polys, xy, crs, idCol, polyCauseCol) {
  pts <- terra::vect(xy, crs = crs)
  n <- NROW(xy)
  pts$YEAR <- if (n) as.integer(polys$YEAR) else integer()
  pts$.poly <- if (n) polys$.poly else integer()
  if (!is.null(idCol)) pts[[idCol]] <- if (n) polys[[idCol]][, 1] else integer()
  pts$polyHa <- if (n) polys$.ha else numeric()
  if (n && polyCauseCol %in% names(polys)) pts$cause <- polys[[polyCauseCol]][, 1] else
    pts$cause <- rep(NA_character_, n)
  pts$burnableHa <- numeric(n)
  pts$pixelID <- integer(n)
  pts$core <- logical(n)
  pts
}

#' Which points or cells are in an ELF's core
#'
#' The one rule for the core of an ELF, for every part of fireSense: the fit's fires, the escape and ignition data,
#' the burn and the hindcast. Three kinds of raster are in use. The values of the first two cannot tell them apart
#' (a 1 is buffer in one and core in the other), so `coreIs` says which a numeric `coreRaster` is:
#' - `"ELF"`: a one-ELF raster, `rasterToMatchLargeELF`: the core is [fireSenseELFCoreValue] (2), the buffer 1.
#' - `"mask"`: the core alone, as `rasterToMatchELF` and the `coreRaster` of [hindcastObserved()] and
#'   [hindcastSimulated()]: the core is every value above 0; `NA` or 0 elsewhere.
#' - A categorical raster (`terra::is.factor()`), as `fireSense_ELFs` makes `rasterToMatchLargeELF` of a study
#'   area over several ELFs: each cell is labelled with its ELF (the codes start at 0) and the study area has no
#'   buffer, so every labelled cell is core. `coreIs` is not used.
#'
#' @param x Points (`SpatVector` or `sf`), or cell numbers of `grid`.
#' @param coreRaster `SpatRaster` of the core, as above. `NULL`: everything is core.
#' @param coreIs `"ELF"` (default) or `"mask"`: what a numeric `coreRaster` is.
#' @param grid `SpatRaster` whose cells `x` numbers; default `coreRaster`. On another grid each cell is placed by its
#'   centre.
#'
#' @return A logical vector, one per point or cell: `TRUE` in the core.
#' @export
#' @importFrom terra compareGeom crds extract project same.crs xyFromCell
inELFCore <- function(x, coreRaster, coreIs = c("ELF", "mask"), grid = coreRaster) {
  coreIs <- match.arg(coreIs)
  cells <- is.numeric(x)
  n <- if (cells) length(x) else NROW(x)
  if (is.null(coreRaster)) return(rep(TRUE, n))
  if (!n) return(logical())
  v <- if (cells) {
    if (terra::compareGeom(grid, coreRaster, stopOnError = FALSE, messages = FALSE)) terra::extract(coreRaster, x)[[1]]
    else terra::extract(coreRaster, terra::xyFromCell(grid, x))[[1]]
  } else {
    if (!inherits(x, "SpatVector")) x <- terra::vect(x)
    if (!terra::same.crs(x, coreRaster)) x <- terra::project(x, terra::crs(coreRaster))
    terra::extract(coreRaster, terra::crds(x))[[1]]
  }
  if (terra::is.factor(coreRaster)[1]) return(!is.na(v))
  !is.na(v) & if (coreIs == "ELF") v %in% fireSenseELFCoreValue else v > 0
}

## the NFDB CAUSE of NBAC cause values: the cause group's main NFDB code (N, U, H, RE)
nfdbCauseOf <- function(nbac) {
  out <- rep(NA_character_, length(nbac))
  for (g in c("natural", "unknown", "human", "other")) {
    code <- setdiff(fireCauseCodes(g, "nfdb"), c("L", "", NA, "H-PB"))[1]
    out[nbac %in% fireCauseCodes(g, "nbac") & is.na(out)] <- code
  }
  out
}
