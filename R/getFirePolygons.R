#' Download and prepare fire data from National Fire Database
#'
#' @param years years to filter fire polygons by
#' @template studyArea
#' @param destinationPath the directory where the fire polys will be stored
#' @param useInnerCache cache the prepInputs call
#'
#' @return DESCRIPTION NEEDED
#'
#' @export
#' @importFrom reproducible prepInputs
#' @importFrom sf st_area st_zm st_transform st_is_longlat
getFirePolygons <- function(years, studyArea, destinationPath, useInnerCache = FALSE) {

  currentURL <- "https://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_poly/current_version/NFDB_poly.zip"
  firePolys <- prepInputs(
    url = currentURL,
    studyArea = studyArea,
    useSAcrs = TRUE,
    fun = "terra::vect",
    destinationPath = destinationPath,
    useCache = useInnerCache
  ) # this object will cache several gigabytes of cached for a small object

  #consider terra compareGeom for two vectors... unsure

  firePolys$YEAR <- as.numeric(firePolys$YEAR) # why is it character?

  if (!is(firePolys, "sf")) {
    firePolys <- sf::st_as_sf(firePolys)
  }
  firePolys <- st_zm(firePolys)

  firePolygonsList <- lapply(years, FUN = function(x, polys = firePolys) {
    firePoly <- polys[polys$YEAR == x, ]
    if (nrow(firePoly) > 0) {
      if (st_is_longlat(firePoly)) {
        stop("please use a study area that is projected in metres")
      }
      firePoly$POLY_HA <- round(sf::st_area(firePoly, by_element = TRUE) / 1e4, digits = 2)
      firePoly <- firePoly[!duplicated(firePoly$FIRE_ID), ]
      return(firePoly)
    } else {
      return(NULL)
    }
  })

  names(firePolygonsList) <- paste0("year", years)
  return(firePolygonsList)
}
