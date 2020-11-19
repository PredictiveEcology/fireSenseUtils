#' preparing terrain covariates for fireSense modules
#'
#' @param rasterToMatch template raster
#' @param studyArea \code{SpatialPolygonsDataFrame} file of \code{studyArea}
#' @param destinationPath directory where data is downloaded
#'
#' @return a raster stack of terrain indices
#'
#' @export
#' @importFrom raster raster stack terrain
#' @importFrom reproducible prepInputs
#' @importFrom spatialEco hli
#' @rdname prepTerrainCovariates
prepTerrainCovariates <- function(rasterToMatch, studyArea, destinationPath) {
  DEM <- prepInputs(url = 'https://drive.google.com/file/d/121x_CfWy2XP_-1av0cYE7sxUfb4pmsup/view?usp=sharing',
                    destinationPath = destinationPath,
                    studyArea = studyArea,
                    rasterToMatch = rasterToMatch,
                    userTags = c("DEM", "fireSenseUtils"))
  otherIndices <- raster::terrain(x = DEM,
                                  opt = c("TPI"),
                                  unit = "degrees")
  hli <- spatialEco::hli(DEM)

  #if we add saga wetness index or dynatop's TWI or catchment area, here is the place

  terrainCovariates <- stack(otherIndices, hli)
  names(terrainCovariates) <- c("TPI", "HLI")
  return(terrainCovariates)
}
