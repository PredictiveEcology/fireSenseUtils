globalVariables(c(
  ".", "age", "B", "pixelGroup", "standAge"
))

#' preparing covariates for fitting modules
#'
#' @param cohortData DESCRIPTION NEEDED
#' @param terrainDT optional data table of terrain covariates with pixelIndex
#' @param lcc data.table of dummified landcover
#' @param pixelGroupMap pixelGroupMap to join terrain with pixelGroup
#'
#' @return a trimmed cohortData with wide-layout
#'
#' @export
#' @importFrom data.table copy dcast set setnafill
#' @rdname castCohortData
castCohortData <- function(cohortData, terrainDT = NULL, pixelGroupMap, lcc) {
  cohortData <- copy(cohortData)
  #this can be made more efficient if the terrain data.table is constructed outside the function and passed in...
  #need stand age for predictions
  cohortData[, standAge := sum(B * age) / sum(B), .(pixelGroup)] # don't chain with magrittr
  cohortData <- dcast(cohortData, pixelGroup + standAge ~ speciesCode, value.var = c("B"), fill = 0)

  cohortDataLong <- data.table('pixelID' = 1:ncell(pixelGroupMap), 'pixelGroup' = getValues(pixelGroupMap))
  cohortData <- cohortData[cohortDataLong, on = c("pixelGroup")]
  rm(cohortDataLong)
  cohortData <- cohortData[lcc, on = c("pixelID")]

  if (!is.null(terrainDT)) {
    cohortData <- cohortData[terrainDT, on = c("pixelID")]
    setnafill(cohortData, fill = 0, cols = colnames(cohortData)[!colnames(cohortData) %in%
                                                                  c('pixelGroup', colnames(terrainDT))])
  }

  return(cohortData)
}
