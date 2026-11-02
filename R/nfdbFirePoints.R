## NFDB fire points of a study area, for the ignition fits (fireSense_dataPrepFit) and for the sizes of the
## fires that did not escape, which fireSense_burn gives its non-escaped ignitions. Moved from
## fireSense_dataPrepFit (fireRecords.R, escape.R), so fireSense_burn can make them in a run without it.

#' NFDB fire points in a study area, of every size and cause
#'
#' @param shp path of the NFDB point shapefile.
#' @param years integer years to keep.
#' @param studyArea polygon to crop to; the result is in its CRS.
#' @export
#' @return spatial points, with `YEAR` and `SIZE_HA` harmonized by fireregimetools.
nfdbFirePoints <- function(shp, years, studyArea) {
  fireregimetools::load_nfdb_points(shp, study_area = studyArea, fire_years = years, min_size_ha = 0)
}

#' NFDB fire points of the given cause groups
#'
#' @param points NFDB fire points (`SpatVector`) with column `CAUSE`.
#' @param groups Cause groups, see `fireCauseCodes()`.
#' @export
#' @return The rows of `points` whose `CAUSE` is in `groups`.
pointsOfCauses <- function(points, groups) {
  points[points$CAUSE %in% fireCauseCodes(groups, "nfdb"), ]
}

#' NFDB fire points with one row per fire
#'
#' Some NFDB rows repeat another row's date (`YEAR`, `MONTH`, `DAY`) and location exactly, nearly always under
#' a different `FIRE_ID`: 554 of the 134,007 natural fires of 1986-2025 in release 20260811, 319 of them in BC.
#' Each such group is taken as one fire recorded more than once, and the row with the largest `SIZE_HA` is kept.
#'
#' @param pts `SpatVector` of NFDB points with `YEAR`, `MONTH`, `DAY` and `SIZE_HA`.
#' @export
#' @return `pts` without the repeated rows, in their original order.
dropRepeatedFirePoints <- function(pts) {
  xy <- terra::crds(pts)
  key <- paste(pts$YEAR, pts$MONTH, pts$DAY, xy[, 1], xy[, 2])
  ord <- order(pts$SIZE_HA, decreasing = TRUE)
  pts[sort(ord[!duplicated(key[ord])]), ]
}

#' Sizes of the fires that did not escape
#'
#' @param fires Fire points or polygons (`sf`, `SpatVector` or `data.frame`) with a size column in ha.
#' @param escapeSizeHa Size (ha) a fire must reach to count as escaped.
#' @param sizeCol The size column.
#' @export
#' @return Numeric vector (ha) of the positive fire sizes below `escapeSizeHa`.
nonEscapedFireSizes <- function(fires, escapeSizeHa, sizeCol = "SIZE_HA") {
  s <- as.data.frame(fires)[[sizeCol]]
  s[is.finite(s) & s > 0 & s < escapeSizeHa]
}

#' The study area's NFDB fire points of the given causes, one row per fire
#'
#' The fires of `years` in `studyArea` ([nfdbFirePoints()]), of `causes` ([pointsOfCauses()]), one row per
#' fire ([dropRepeatedFirePoints()]): `fireSense_dataPrepFit`'s `ignitionFirePoints`, and the fires whose
#' sizes `fireSense_burn` draws for the ignitions that do not escape.
#'
#' @param shp path of the NFDB point shapefile, e.g.
#'   `fireRecordShapefile(nfdbPointUrl(), destinationPath)`; its name carries the release, so a caller
#'   can key a `Cache()` on it.
#' @param studyArea polygon (`SpatVector` or `sf`) to crop to; the result is in its CRS.
#' @param years integer years to keep.
#' @param causes cause groups, see [fireCauseCodes()].
#' @return `SpatVector` of points with `YEAR`, `SIZE_HA` and `CAUSE`.
#' @export
nfdbFirePointsOfCauses <- function(shp, studyArea, years, causes) {
  dropRepeatedFirePoints(pointsOfCauses(nfdbFirePoints(shp, years = years, studyArea = studyArea), causes))
}
