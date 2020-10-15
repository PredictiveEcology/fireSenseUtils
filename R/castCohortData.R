#' preparing covariates for fitting modules
#'
#' @param cohortData DESCRIPTION NEEDED
#' @param terrainDT rasterStack of terrain covariates converted to data.table
#' @param pixelGroupMap pixelGroupMap to join terrain with pixelGroup
#' @param climateRasters if predicting, a single raster layer (for now), else raster stack of historical climate
#' @param index list of relevant pixel IDs for subsetting, if applicable
#' @return a trimmed cohortData with wide-layout
#' @importFrom data.table dcast set setnafill
#' @importFrom dplyr bind_cols
#' @export
#' @rdname castCohortData
castCohortData <- function(cohortData, terrainDT, pixelGroupMap, index = NULL, climateRasters) {
  #this can be made more efficient if the terrain data.table is constructed outside the function and passed in...
  terrainDT <- copy(terrainDT)
  terrainNames <- colnames(terrainDT)
  #need stand age for predictions
  cohortData[, standAge := sum(B * age)/sum(B), .(pixelGroup)]

  cohortDataWide <- dcast(cohortData, pixelGroup + standAge ~ speciesCode, value.var = c("B"),
                          fill = 0)

  set(terrainDT, j = 'pixelGroup', value = getValues(pixelGroupMap))

  #if index is present, we need pixelIndices too
  if (!is.null(index)) {
    terrainDT<- terrainDT[!is.na(get(terrainNames))]
    terrainDT <- cohortDataWide[terrainDT, on = c("pixelGroup")]
    terrainDT <- data.table::setnafill(terrainDT, fill = 0, cols = names(cohortDataWide))
    set(terrainDT, , 'pixelGroup', NULL)
    #NAs from non-treed data

    #need to grab climate layers for each year
    fullCovariates <- lapply(names(index), FUN = function(x, ind = index, clim = climateRasters, ter = terrainDT) {
      ind <- ind[[x]]
      clim <- clim[[x]]
      climDat <- data.table(pixelID = 1:ncell(clim), MDC = getValues(clim))
      ind <- climDat[ind, on = c("pixelID")]
      fullYear <- ter[ind, on = c("pixelID")]
      return(fullYear)
    })
    fullCovariates <- rbindlist(fullCovariates) #I don't think these need to be separated by year anymore..

  } else {
    set(terrainDT, j = 'MDC', value = getValues(climateRasters))
    #this will be for predicting- when we get there
    #do stuff and make fulLCovariates
  }

  return(fullCovariates)
}
