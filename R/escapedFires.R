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

#' Ignitions that escaped: those that became a burned-area polygon of at least `escapeSizeHa`
#'
#' The escape model's response. An ignition point escaped if it lies in a fire polygon of the same year whose
#' area is at least `escapeSizeHa` and larger than one pixel (the polygons [escapedFires()] keeps, the same
#' fires the spread model is fitted to). A polygon holding several points is one escape: only its point with
#' the largest `sizeCol` (the first of ties) is marked, and a point is marked at most once. A polygon with no
#' point is not an escape, as it has no ignition to count it against. The polygon's area is its area as given
#' (`terra::expanse()`), so polygons clipped to a study area are judged on the part inside it.
#'
#' @param points Ignition points (`SpatVector` or `sf`) with columns `YEAR` and `sizeCol`.
#' @param polys Fire polygons (`SpatVector` or `sf`) with column `YEAR`, or a list of them (one per year, as
#'   `firePolys` is in fireSense; `NULL` elements are dropped). Any CRS: they are projected to that of `points`.
#' @param escapeSizeHa Size (ha) a polygon must reach.
#' @param pixSizeHa Pixel area (ha). A polygon must also be larger than one pixel.
#' @param sizeCol The column of `points` that picks one point per polygon.
#' @return Logical, one per row of `points`: `TRUE` for an escape.
#' @export
#' @importFrom terra crs expanse project relate same.crs vect
escapedIgnitions <- function(points, polys, escapeSizeHa, pixSizeHa, sizeCol = "SIZE_HA") {
  asVect <- function(x) if (inherits(x, "SpatVector")) x else terra::vect(x)
  points <- asVect(points)
  if (is.list(polys) && !inherits(polys, "sf")) {
    polys <- Filter(Negate(is.null), polys)
    polys <- if (length(polys)) do.call(rbind, unname(lapply(polys, asVect)))
  }
  out <- logical(NROW(points))
  if (!length(out) || !NROW(polys)) return(out)
  polys <- asVect(polys)
  if (!terra::same.crs(polys, points)) polys <- terra::project(polys, terra::crs(points))
  polys$.ha <- terra::expanse(polys, unit = "ha")
  polys <- escapedFires(polys, escapeSizeHa, pixSizeHa, sizeCol = ".ha")
  if (!NROW(polys)) return(out)

  hits <- terra::relate(points, polys, "intersects", pairs = TRUE)
  hits <- hits[as.integer(points$YEAR[hits[, 1]]) == as.integer(polys$YEAR[hits[, 2]]), , drop = FALSE]
  ## per polygon, its largest point first; then each polygon takes its first point not already taken
  hits <- hits[order(hits[, 2], -as.data.frame(points)[[sizeCol]][hits[, 1]], hits[, 1]), , drop = FALSE]
  for (p in unique(hits[, 2])) {
    free <- hits[hits[, 2] == p & !out[hits[, 1]], 1]
    if (length(free)) out[free[1]] <- TRUE
  }
  out
}
