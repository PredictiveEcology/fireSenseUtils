## An escaped fire is one that reached `escapeSizeHa` (the same threshold in the escape model, the spread fit
## and fireSense's burning). Fires below it, though recorded, are not spread; their sizes are kept so a
## forecast can give its non-escaped ignitions a size.

#' Fires that escaped
#'
#' @param fires Fire points or polygons (`sf`, `SpatVector` or `data.frame`) with a size column in ha.
#' @param escapeSizeHa Size (ha) a fire must reach to count as escaped.
#' @param pixSizeHa Pixel area (ha). A fire must also be larger than one pixel.
#' @param sizeCol The size column.
#' @return The rows of `fires` that escaped.
#' @export
escapedFires <- function(fires, escapeSizeHa, pixSizeHa, sizeCol = "SIZE_HA") {
  s <- as.data.frame(fires)[[sizeCol]]
  fires[which(s > pixSizeHa & s >= escapeSizeHa), ]
}

#' Ignition points and the fire polygons they match
#'
#' The matching [escapedIgnitions()] and [polygonIgnitions()] share. A polygon counts if its area is at least
#' `escapeSizeHa` and larger than one pixel (the polygons [escapedFires()] keeps, the fires the spread model is
#' fitted to; area as given, by `terra::expanse()`, so a polygon clipped to a study area is judged on the part
#' inside it). A point matches a polygon of its year that it lies in, or is within `maxDistM` of: NFDB locations
#' are often rounded to 0.01 degree (about 1 km), and on ELF 6.6.1 20 of the 31 natural points of 50 ha or
#' more outside every polygon of their year were within 1 km of one.
#'
#' @param points Ignition points (`SpatVector` or `sf`) with column `YEAR`.
#' @param polys Fire polygons (`SpatVector` or `sf`) with column `YEAR`, or a list of them (one per year, as
#'   `firePolys` is in fireSense; `NULL` elements are dropped). Any CRS: they are projected to that of `points`.
#' @param escapeSizeHa Size (ha) a polygon must reach.
#' @param pixSizeHa Pixel area (ha). A polygon must also be larger than one pixel.
#' @param maxDistM Distance (m) from a polygon's edge within which a point outside it still matches it.
#' @return A list: `points` and `polys` (the polygons that count, with their area in `.ha`), both
#'   `SpatVector`, and `pairs`, a matrix with columns `point` and `poly` (rows of those) and logical `inside`.
#' @export
#' @importFrom terra buffer crs expanse project relate same.crs vect
ignitionPolygonPairs <- function(points, polys, escapeSizeHa, pixSizeHa, maxDistM = 1000) {
  asVect <- function(x) if (inherits(x, "SpatVector")) x else terra::vect(x)
  points <- asVect(points)
  if (is.list(polys) && !inherits(polys, "sf")) {
    polys <- Filter(Negate(is.null), polys)
    polys <- if (length(polys)) do.call(rbind, unname(lapply(polys, asVect)))
  }
  none <- list(points = points, polys = NULL, pairs = cbind(point = integer(), poly = integer(), inside = integer()))
  if (!NROW(polys)) return(none)
  polys <- asVect(polys)
  if (!terra::same.crs(polys, points)) polys <- terra::project(polys, terra::crs(points))
  polys$.ha <- terra::expanse(polys, unit = "ha")
  polys <- escapedFires(polys, escapeSizeHa, pixSizeHa, sizeCol = ".ha")
  none$polys <- polys
  if (!NROW(polys) || !NROW(points)) return(none)

  pairs <- function(y) {
    h <- terra::relate(points, y, "intersects", pairs = TRUE)
    h[as.integer(points$YEAR[h[, 1]]) == as.integer(polys$YEAR[h[, 2]]), , drop = FALSE]
  }
  inside <- pairs(polys)
  hits <- if (maxDistM > 0) pairs(terra::buffer(polys, maxDistM)) else inside
  isIn <- paste(hits[, 1], hits[, 2]) %in% paste(inside[, 1], inside[, 2])
  list(points = points, polys = polys, pairs = cbind(point = hits[, 1], poly = hits[, 2], inside = isIn))
}

#' Ignitions that escaped: those that became a burned-area polygon of at least `escapeSizeHa`
#'
#' The escape model's response. An ignition point escaped if it matches a fire polygon of its year (see
#' [ignitionPolygonPairs()]: at least `escapeSizeHa`, larger than one pixel, the point inside it or within
#' `maxDistM`). A polygon matching several points is one escape: points inside it come before points only near
#' it, then the largest `sizeCol` (the first of ties); a point is marked at most once. A natural-cause polygon
#' that matches no point is added as an ignition by [polygonIgnitions()], and that point then escapes here.
#'
#' @inheritParams ignitionPolygonPairs
#' @param points Ignition points (`SpatVector` or `sf`) with columns `YEAR` and `sizeCol`.
#' @param sizeCol The column of `points` that picks one point per polygon.
#' @return Logical, one per row of `points`: `TRUE` for an escape.
#' @export
escapedIgnitions <- function(points, polys, escapeSizeHa, pixSizeHa, sizeCol = "SIZE_HA", maxDistM = 1000) {
  m <- ignitionPolygonPairs(points, polys, escapeSizeHa, pixSizeHa, maxDistM)
  out <- logical(NROW(m$points))
  hits <- m$pairs
  ## per polygon: points inside first, then the largest; each polygon takes its first point not already taken
  size <- as.data.frame(m$points)[[sizeCol]][hits[, "point"]]
  hits <- hits[order(hits[, "poly"], !hits[, "inside"], -size, hits[, "point"]), , drop = FALSE]
  for (p in unique(hits[, "poly"])) {
    free <- hits[hits[, "poly"] == p & !out[hits[, "point"]], "point"]
    if (length(free)) out[free[1]] <- TRUE
  }
  out
}

#' Ignitions for natural-cause fire polygons that no ignition point matches
#'
#' Eliot, 2026-10-25: a natural-cause burned-area polygon of at least `escapeSizeHa` that no natural ignition
#' point matches (see [ignitionPolygonPairs()]) was still an ignition, and it escaped. It is added at the
#' fire's start point, the one the spread fit uses: `starts`, the polygon centroids that
#' [harmonizeFireData()] moved to the nearest burnable pixel of the fire ([harmonizeBufferAndPoints()],
#' [cleanUpSpreadFirePoints()]). Bind the result to the ignition points: each adds one ignition, and
#' [escapedIgnitions()] then counts it as the polygon's escape.
#'
#' @inheritParams ignitionPolygonPairs
#' @param points Natural-cause ignition points (`SpatVector` or `sf`) with column `YEAR`.
#' @param starts One start point per fire polygon (`SpatVector` or `sf`, or a list of them), sharing `idCol`
#'   with `polys`, e.g. `spreadFirePoints` after `fireSense_dataPrepFit`'s `prepSpreadFitData`.
#' @param idCol The fire ID column of `polys` and `starts`.
#' @param polyCauseCol,naturalPolyCauses Column of `polys` holding the cause, and its natural values.
#' @return `SpatVector` of points in the CRS of `points`, one per unmatched natural-cause polygon that has a
#'   start point, with `YEAR`, `CAUSE` (`"N"`) and `SIZE_HA` (the polygon's area), the columns of `points`
#'   that the ignition and escape counts read.
#' @export
polygonIgnitions <- function(points, polys, starts, idCol, escapeSizeHa, pixSizeHa, maxDistM = 1000,
                             polyCauseCol = "FIRECAUS", naturalPolyCauses = "Natural") {
  m <- ignitionPolygonPairs(points, polys, escapeSizeHa, pixSizeHa, maxDistM)
  asPoints <- function(xy, p) {
    out <- terra::vect(xy, crs = terra::crs(m$points))
    out$YEAR <- as.integer(p$YEAR)
    out$CAUSE <- rep("N", NROW(p))
    out$SIZE_HA <- p$.ha
    out
  }
  p <- m$polys
  if (NROW(p) && !polyCauseCol %in% names(p)) {
    warning("polys has no column '", polyCauseCol, "': no polygon ignitions added")
    p <- NULL
  }
  if (NROW(p)) p <- p[p[[polyCauseCol]][, 1] %in% naturalPolyCauses & !seq_len(NROW(p)) %in% m$pairs[, "poly"], ]
  if (!NROW(p)) return(asPoints(matrix(numeric(), ncol = 2), data.frame(YEAR = integer(), .ha = numeric())))
  asVect <- function(x) if (inherits(x, "SpatVector")) x else terra::vect(x)
  if (is.list(starts) && !inherits(starts, "sf"))
    starts <- do.call(rbind, unname(lapply(Filter(Negate(is.null), starts), asVect)))
  starts <- terra::project(asVect(starts), terra::crs(m$points))
  key <- function(x) paste(as.integer(x$YEAR), x[[idCol]][, 1])
  at <- match(key(p), key(starts))
  p <- p[!is.na(at), ]
  asPoints(terra::crds(starts)[at[!is.na(at)], , drop = FALSE], p)
}
