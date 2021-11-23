globalVariables(c(
 "flammable"
))

#' preparing terrain covariates for fireSense modules
#'
#' @template rasterToMatch
#' @template studyArea
#' @param destinationPath directory where data is downloaded
#'
#' @return a raster stack of terrain indices
#'
#' @export
#' @importFrom raster raster stack terrain
#' @importFrom reproducible prepInputs
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
  if (requireNamespace("spatialEco", quietly = TRUE)) {
    hli <- spatialEco::hli(DEM)
  } else {
    stop("Please install.packages('spatialEco')")
  }
  #if we add saga wetness index or dynatop's TWI or catchment area, here is the place

  terrainCovariates <- stack(otherIndices, hli)
  names(terrainCovariates) <- c("TPI", "HLI")
  return(terrainCovariates)
}

#' Convert terrain raster into \code{data.table}
#'
#' @param terrainCovariates raster stack with terrain values for PCA
#' @template flammableRTM
#'
#' @return a data.table with terrain values for each flammable pixel
#'
#' @export
#' @importFrom data.table set setDT
#' @importFrom raster getValues nlayers
#' @rdname buildTerrainDT
buildTerrainDT <- function(terrainCovariates, flammableRTM){
  layers <- seq(nlayers(terrainCovariates))
  names(layers) <- names(terrainCovariates)
  # Seems to be the fastest way to get data out of a RasterStack
  terrainDT <- setDT(lapply(layers, FUN = function(x)
    getValues(terrainCovariates[[x]])
  ))

  set(terrainDT, j = "pixelID", value = 1:ncell(flammableRTM))
  set(terrainDT, j = "flammable", value = getValues(flammableRTM))
  terrainDT <- terrainDT[flammable == 1, ] %>%
    set(., NULL, "flammable", NULL) %>%
    na.omit(.)

  return(terrainDT)
}
