#' Download and prepare fire data from National Fire Database
#'
#' @param years years to filter fire polygons by
#' @param studyArea study area to subset fire polygons to
#' @param destinationPath the directory where the fire polys will be stored
#' @param useInnerCache cache the prepInputs call
#'
#' @return DESCRIPTION NEEDED
#'
#' @export
#' @importFrom reproducible prepInputs
getFirePolygons <- function(years, studyArea, destinationPath,
                            useInnerCache = FALSE) {
  currentURL <- "https://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_poly/current_version/NFDB_poly.zip"
  firePolys <- prepInputs(url = currentURL,
                          studyArea = studyArea,
                          destinationPath = destinationPath,
                          useCache = useInnerCache) #this object will cache several gigabytes of cached for a small object

  firePolys$YEAR <- as.numeric(firePolys$YEAR) #why is it character?
  firePolygonsList <- lapply(years, FUN = function(x, polys = firePolys){
    firePoly <- polys[polys$YEAR == x,]
    if (nrow(firePoly) > 0) {
    firePoly$POLY_HA <- rgeos::gArea(firePoly, byid = TRUE)
    firePoly <- firePoly[!duplicated(firePoly$FIRE_ID),]
    return(firePoly)
    } else {
      return(NULL)
    }
  })

  names(firePolygonsList) <- paste0("year", years)
  return(firePolygonsList)
}
