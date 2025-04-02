#' Retrieve NTEMS LCC with relevant nonflammable cover corrected to flammable classes during
#' reprojection from native resolution to that fo rasterToMatch.
#'
#' @param neededYear year of NTEMS LCC
#' @param writeTo optional filename for output raster
#' @template studyArea
#' @template rasterToMatch
#' @param nonflammableLCC nonflammable cover values
#' @param destinationPath passed to `prepInputs`. Directory to store downloaded file
#' @param flammabilityThreshold the flammability threshold, below which pixels become NA
#'
#' @importFrom terra mask
#' @returns a landcover raster where pixels below mean flammabilityThreshold are converted to 0,
#' otherwise the nearest flammable landcover is returned at the resolution of rasterToMatch
#' @export
#'
makeFireSenseLCC <- function(neededYear, rasterToMatch, studyArea,
                             nonflammableLCC = c(20, 31, 32, 33),
                             flammabilityThreshold = 0.1, writeTo = NULL,
                             destinationPath) {

  #1. NA the non-flammable
  #2. reproject the landcover for fireSense
  #3. make a flammableMap from the original resolution
  #4. burn the non-flammable pixels (below thresh) into the LCC from 2


  rstLCC <- prepInputs_NTEMS_LCC_FAO(year = neededYear,
                                     disturbedCode = 240,
                                     overwrite = TRUE,
                                     destinationPath = destinationPath,
                                     cropTo = rasterToMatch,
                                     maskTo = studyArea)
  #we want the mode of non-NA pixels but that is not an option
  allFlam <- terra::mask(rstLCC, rstLCC, maskvalue = nonflammableLCC, updatevalue = NA) |>
    project(y = rasterToMatch, method = "mode")
  #these are the values that will be assigned

  #you can probably make the inverse of allFlam, focal that to get the amount of non-flammable cover, and then burn in the
  #non flam into allFlam
  nonFlam <- terra::mask(rstLCC, rstLCC, updatevalue = 0, maskvalues = nonflammableLCC)
  nonFlam[nonFlam > 0] <- 1
  nonFlam <- terra::project(nonFlam, rasterToMatch, method = "average")

  allFlam[nonFlam < flammabilityThreshold] <- 0 #for nonflammable
  if (!is.null(writeTo)) {
    allFlam <- writeRaster(allFlam, filename = file.path(destinationPath, writeTo), overwrite = TRUE)
  }

  return(allFlam)

}
