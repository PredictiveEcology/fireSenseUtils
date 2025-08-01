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
getFirePolygons <- function(years, useInnerCache = FALSE, ...) {
  currentURL <- "https://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_poly/current_version/NFDB_poly.zip"
  firePolys <- prepInputs(
    url = currentURL,
    useCache = useInnerCache,
    ...
  ) ## this object will cache several gigabytes of cached for a small object

  firePolys$YEAR <- as.numeric(firePolys$YEAR) ## it has been character before...

  firePolygonsList <- lapply(years, FUN = function(x, polys = firePolys) {
    firePoly <- polys[polys$YEAR == x, ]
    if (nrow(firePoly) > 0) {
      firePoly$POLY_HA <- round(st_area(firePoly, unit = "ha") / 1e4, digits = 2)
      # firePoly <- firePoly[!duplicated(firePoly$FIRE_ID), ] ## TODO: what was this for?
      return(firePoly)
    } else {
      return(NULL)
    }
  })

  names(firePolygonsList) <- paste0("year", years)
  return(firePolygonsList)
}
