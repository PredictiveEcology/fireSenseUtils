globalVariables(c(
  ".", "age", "B", "pixelGroup", "standAge", "nonforest", "youngAge"
))

#' preparing covariates for fitting modules
#'
#' @template cohortData
#' @param terrainDT optional data table of terrain covariates with pixelIndex
#' @param lcc data.table of dummified landcover
#' @param pixelGroupMap pixelGroupMap to join terrain with pixelGroup
#' @param ageMap a stand age map to assign ages to non-forest LCC used during predict
#' @param missingLCC LCC class to assign forested pixels absent from cohortData
#'   must be a character matching a nonForestedLCC group, e.g. 'nonForest_highFlam'
#' @param year numeric representing the year represented by cohortData
#' @param cutoffForYoungAge Numeric. Default is 15. This is the age below which the pixel is considered
#'   "young" --> youngAge column will be 1 if age <= 15
#'
#' @return a trimmed cohortData with wide-layout and rows for every pixel in lcc
#'
#' @export
#' @importFrom data.table copy dcast set setnafill
#' @rdname castCohortData
castCohortData <- function(cohortData, terrainDT = NULL, pixelGroupMap, lcc, ageMap = NULL,
                           missingLCC, year = NULL,
                           cutoffForYoungAge = 15) {

  #need stand age for predictions but it won't be included in PCAer
  # cohortData[, standAge := sum(B * age) / sum(B), .(pixelGroup)] # don't chain with magrittr
  #     Eliot removed the above line because it caused NAs to occur wherever there was sum(B) == 0
  #           We also don't want biomass-weighted stand age here, I believe. This should be "time since disturbance"
  cohortData[, standAge := max(age, na.rm = TRUE), .(pixelGroup)] # don't chain with magrittr
  cohortData <- dcast(cohortData, pixelGroup + standAge ~ speciesCode,
                      value.var = c("B"), fun.aggregate = sum, fill = 0)

  cohortDataLong <- data.table('pixelID' = 1:ncell(pixelGroupMap), 'pixelGroup' = getValues(pixelGroupMap))
  cohortData <- cohortData[cohortDataLong, on = c("pixelGroup")]
  rm(cohortDataLong)
  cohortData <- cohortData[lcc, on = c("pixelID")]

  #reclassify forested lcc that is missing from cohortData as missingLCC group
  lccNames <- colnames(lcc)[!colnames(lcc) %in% c("pixelID")]

   #this should be faster
  cohortData[, nonforest := rowSums(cohortData[, .SD, .SDcols = lccNames])]
  cohortData[is.na(pixelGroup) & nonforest == 0, eval(missingLCC) := 1] #these are forest by LCC2005 and not LandR
  set(cohortData, ,'nonforest', NULL)

  if (!is.null(ageMap)) { #should not be NULL in predict
    cohortData[is.na(standAge), standAge := ageMap[cohortData[is.na(standAge)]$pixelID]]
  }
  if (!is.null(terrainDT)) {
    cohortData <- cohortData[terrainDT, on = c("pixelID")]
    setnafill(cohortData, fill = 0, cols = colnames(cohortData)[!colnames(cohortData) %in%
                                                                  c('pixelGroup', colnames(terrainDT))])
  }

  set(cohortData, NULL, "youngAge", as.integer(cohortData$standAge <= cutoffForYoungAge))
  # cohortData[, youngAge := ifelse(cohortData$standAge <= cutoffForYoungAge, 1, 0)]
  set(cohortData, NULL, 'standAge', NULL)
  #we only care if stand age is young
  #we don't care at all about this during fit as we will reclassify it anyway
  if (!is.null(year)) {
    set(cohortData, NULL, 'year', year)
  }
  return(cohortData)
}
