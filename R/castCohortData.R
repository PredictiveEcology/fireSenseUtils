#' preparing covariates for fitting modules
#'
#' @param cohortData DESCRIPTION NEEDED
#' @param terrainDT optional data table of terrain covariates with pixelIndex
#' @param pixelGroupMap pixelGroupMap to join terrain with pixelGroup
#' @param optional index list of relevant pixel IDs for subsetting, if applicable
#' @return a trimmed cohortData with wide-layout
#' @importFrom data.table dcast set setnafill
#' @importFrom dplyr bind_cols
#' @export
#' @rdname castCohortData
castCohortData <- function(cohortData, terrainDT = NULL, pixelGroupMap, index = NULL) {

  cohortData <- copy(cohortData)
  #this can be made more efficient if the terrain data.table is constructed outside the function and passed in...
  #need stand age for predictions
  cohortData[, standAge := sum(B * age)/sum(B), .(pixelGroup)] #don't chain with magrittr
  cohortData <- dcast(cohortData, pixelGroup + standAge ~ speciesCode, value.var = c("B"), fill = 0)

  cohortDataLong <- data.table('pixelID' = 1:ncell(pixelGroupMap), 'pixelGroup' = getValues(pixelGroupMap)) %>%
    na.omit(.)
  cohortData <- cohortData[cohortDataLong, on = c("pixelGroup")]

  rm(cohortDataLong)
  if (!is.null(terrainDT)){
    cohortData <- cohortData[terrainDT, on = c("pixelID")]
    setnafill(cohortData, fill = 0, cols = colnames(cohortData)[!colnames(cohortData) %in% c('pixelGroup', colnames(terrainDT))])
  }
  #we need pixelIndex due to terrain - this way we don't have to store long version of sim$terrainDT

  return(cohortData)
}
