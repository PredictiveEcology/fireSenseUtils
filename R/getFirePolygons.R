#' Download and prepare fire data from National Fire Database
#'
#' @param years years to filter fire polygons by
#' @param ... additional arguments passed to [reproducible::prepInputs()]
#' @param useInnerCache logical indicating whether to cache the `prepInputs` call
#'
#' @return list of fire polygons by year
#'
#' @export
#' @importFrom reproducible prepInputs
#' @importFrom sf st_area
getFirePolygons <- function(url, years, useInnerCache = FALSE, ...) {
  # This is stuck at 2021 max year
  if (missing(url))
    # url <- "https://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_poly/current_version/NFDB_poly.zip"
    url <- "https://cwfis.cfs.nrcan.gc.ca/downloads/nbac/NBAC_1972to2024_20250506_shp.zip"
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
