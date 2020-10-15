#' preparing terrain covariates for fireSense modules
#'
#' @param rasterToMatch template raster
#' @param studyArea SpatialPolygonsDataFrame file of studyArea
#' @param destinationPath directory where data is downloaded
#' @return a raster stack of terrain indices
#' @importFrom raster raster stack terrain
#' @importFrom reproducible prepInputs
#' @importFrom spatialEco hli
#' @export
#' @rdname prepTerrainCovariates
prepTerrainCovariates <- function(rasterToMatch, studyArea, dPath) {

  DEM <- prepInputs(url = 'https://drive.google.com/file/d/121x_CfWy2XP_-1av0cYE7sxUfb4pmsup/view?usp=sharing',
                    destinationPath = dPath,
                    studyArea = studyArea,
                    rasterToMatch = rasterToMatch,
                    userTags = c("DEM", "fireSenseUtils"))
  otherIndices <- raster::terrain(x = DEM,
                                  opt = c('slope', 'TPI'),
                                  unit = 'degrees')
  hli <- spatialEco::hli(DEM)
  terrainCovariates <- stack(otherIndices, hli)
  names(terrainCovariates) <- c('TPI', 'slope', 'HLI')
  return(terrainCovariates)
}
