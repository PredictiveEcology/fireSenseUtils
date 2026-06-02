#' Discover the URL of the most recent NBAC shapefile archive
#'
#' Parses the Apache directory listing at `base` and returns the URL of the
#' newest file matching `NBAC_<startYear>to<endYear>_<YYYYMMDD>_shp.zip`,
#' selected by the YYYYMMDD stamp embedded in the filename. If the listing
#' cannot be read or no matching file is found, `fallback` is returned, so
#' callers can rely on getting *some* URL back.
#'
#' @param base Character. Directory listing URL (must be the parent folder
#'   that serves the autoindex page).
#' @param fallback Character. URL to return when discovery fails.
#'
#' @return A single character URL.
#'
#' @export
latestNBACUrl <- function(
  base = "https://cwfis.cfs.nrcan.gc.ca/downloads/nbac/",
  fallback = paste0(base, "NBAC_1972to2025_20260513_shp.zip")
) {
  if (!endsWith(base, "/")) base <- paste0(base, "/")
  tryCatch({
    html <- paste(readLines(base, warn = FALSE), collapse = "\n")
    hits <- regmatches(html, gregexpr("NBAC_\\d+to\\d+_\\d{8}_shp\\.zip", html))[[1L]]
    hits <- unique(hits)
    if (!length(hits)) return(fallback)
    stamp <- sub(".*_(\\d{8})_shp\\.zip$", "\\1", hits)
    paste0(base, hits[which.max(stamp)])
  }, error = function(e) fallback)
}

#' Download and prepare fire data from National Fire Database
#'
#' @param url Character (or `missing`). URL to the fire polygons archive
#'   (typically an NBAC `*_shp.zip` on the CFS server). When `missing`,
#'   [latestNBACUrl()] is queried so the newest NBAC release available on
#'   the CFS download index is used; this means re-running the same code
#'   later may pick up a newer dataset. Pass an explicit URL to pin to a
#'   specific release when reproducibility matters.
#' @param years years to filter fire polygons by
#' @param ... additional arguments passed to [reproducible::prepInputs()]
#' @param useInnerCache logical indicating whether to cache the `prepInputs` call
#'
#' @return list of fire polygons by year
#'
#' @seealso [latestNBACUrl()] for the discovery helper used when `url` is missing.
#'
#' @export
#' @importFrom reproducible prepInputs
#' @importFrom sf st_area
getFirePolygons <- function(url, years, useInnerCache = FALSE, ...) {
  if (missing(url))
    # alt: "https://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_poly/current_version/NFDB_poly.zip"
    url <- latestNBACUrl()
  firePolys <- prepInputs(
    url = url,
    useCache = useInnerCache,
    ...
  ) ## this object will cache several gigabytes of cached for a small object

  firePolys$YEAR <- as.numeric(firePolys$YEAR) ## it has been character before...
  if (is(firePolys, "SpatVector"))
    firePolys$POLY_HA <- round(terra::expanse(firePolys, unit = "ha"), 2)
  else
    firePolys$POLY_HA <- round(sf::st_area(firePolys, unit = "ha"), 2)

  firePolygonsList <- lapply(years, FUN = function(x, polys = firePolys) {
    firePoly <- polys[polys$YEAR == x, ]
    if (nrow(firePoly) > 0) {
      return(firePoly)
    } else {
      return(NULL)
    }
  })

  names(firePolygonsList) <- paste0("year", years)
  return(firePolygonsList)
}
