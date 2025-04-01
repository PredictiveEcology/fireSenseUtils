#' Retrieve NTEMS LCC with relevant nonflammable cover corrected to flammable classes during
#' reprojection from native resolution to that fo rasterToMatch.
#'
#' @param neededYear
#' @param writeTo
#' @template studyArea
#' @template rasterToMatch
#' @param nonflammableLCC nonflammable cover values
#' @param destinationPath passed to `prepInputs`. Directory to store downloaded file
#' @param flammabilityThreshold
#'
#' @returns
#' @export
#'
makeFireSenseLCC <- function(neededYear, writeTo, studyArea, rasterToMatch,
                             nonflammableLCC, destinationPath, flammabilityThreshold = 0.1) {
  #1. NA the non-flammable
  #2. reproject the landcover for fireSense
  #3. make a flammableMap from the original resolution
  #4. burn the non-flammable pixels (below thresh) into the LCC from 2
  #5. write to disk
  browser()
  rstLCC <- prepInputs_NTEMS_LCC_FAO(year = neededYear, disturbedCode = 240,
                                     cropTo = rasterToMatch, maskTo = studyArea,
                                     destinationPath = destinationPath)


  return(LCC)

}
