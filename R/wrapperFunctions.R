
#' preparing cohortData and nonFireAgeMap
#'
#' @template cohortData
#' @template rasterToMatch
#' @param terrainDT optional data table of terrain covariates with pixelIndex
#' @param lcc data.table of dummified landcover
#' @param pixelGroupMap pixelGroupMap to join terrain with pixelGroup
#' @param firePolysForAge a list of sf objects comprising annual fire polygons dated to 1986
#' @param missingLCC Pixels on forested LCC but absent from cohortData will be assigned this LCC.
#' must be a character matching a nonForestedLCC group, e.g. 'nonForest_highFlam'
#' @param fireYears vector of fire years representing range used to determine time since disturbance
#'
#' @return a trimmed cohortData with wide-layout and rows for every pixel in lcc
#'
#' @export
#' @importFrom fasterize fasterize
#' @importFrom magrittr %>%
#' @importFrom raster setValues getValues
#' @importFrom sf st_collection_extract
#' @rdname cohortsToCovariates

cohortsToCovariates <- function(cohortData, terrainDT, lcc, pixelGroupMap, firePolysForAge, missingLCC, rasterToMatch, fireYears) {

  nonForestAgeMap <- firePolysForAge[firePolysForAge$YEAR > min(fireYears) & firePolysForAge$YEAR < max(fireYears),] %>%
    sf::st_collection_extract(., "POLYGON") %>%
    fasterize(sf = ., raster = rasterToMatch, background = min(fireYears) - 1, field = "YEAR", fun = "max") %>%
    setValues(., values = max(fireYears) + 1 - getValues(.))
  #a pixel burning in 2000 should have time since disturbance 1.

  cohorts <- castCohortData(cohortData = cohortData,
                            pixelGroupMap = pixelGroupMap,
                            ageMap = nonForestAgeMap,
                            lcc = lcc,
                            terrainDT = terrainDT,
                            missingLCC = missingLCC)

  return(cohorts)
}
