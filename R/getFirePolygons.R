#' Download and prepare fire data from National Fire Database
#'
#' @param years years to filter fire polygons by
#' @param ... additional arguments passed to prepInputs
#' @param useInnerCache cache the prepInputs call
#'
#' @return DESCRIPTION NEEDED
#'
#' @export
#' @importFrom reproducible prepInputs
#' @importFrom terra expanse
getFirePolygons <- function(years, useInnerCache = FALSE, ...) {

  currentURL <- "https://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_poly/current_version/NFDB_poly.zip"
  firePolys <- prepInputs(
    url = currentURL,
    useCache = useInnerCache,
    ...
  ) # this object will cache several gigabytes of cached for a small object

  firePolys$YEAR <- as.numeric(firePolys$YEAR) #it has been character before.. 

  firePolygonsList <- lapply(years, FUN = function(x, polys = firePolys) {
    firePoly <- polys[polys$YEAR == x, ]
    if (nrow(firePoly) > 0) {
      firePoly$POLY_HA <- round(expanse(firePoly, unit = "ha"), digits = 2)
      # firePoly <- firePoly[!duplicated(firePoly$FIRE_ID), ] what was this for? 
      return(firePoly)
    } else {
      return(NULL)
    }
  })

  names(firePolygonsList) <- paste0("year", years)
  return(firePolygonsList)
}
