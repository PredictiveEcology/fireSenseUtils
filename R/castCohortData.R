#' preparing covariates for fitting modules
#'
#' @param cohortData DESCRIPTION NEEDED
#' @param terrainRasters rasterStack of terrain covariates
#' @param pixelGroupMap pixelGroupMap to join terrain with pixelGroup
#' @param climateRasters if predicting, a single raster layer (for now), else raster stack of historical climate
#' @param index relevant pixel IDs for subsetting, if applicable. list okay
#' @return a trimmed cohortData with wide-layout
#' @importFrom data.table dcast set setnafill
#' @importFrom dplyr bind_cols
#' @export
#' @rdname castCohortData
castCohortData <- function(cohortData, terrainRasters, pixelGroupMap, index = NULL, climateRasters) {
  browser()
  #need stand age for predictions
  cohortData[, standAge := sum(B * age)/sum(B), .(pixelGroup)]

  cohortDataWide <- dcast(cohortData, pixelGroup + standAge ~ speciesCode, value.var = c("B"),
                          fill = 0)
  #get terrain of pixelGroups
  terrain <- lapply(names(terrainRasters), FUN = function(x){
    y <- data.table(getValues(terrainRasters[[x]]))
  }) %>%
    dplyr::bind_cols(.)
  setnames(terrain, names(terrainRasters))

  setDT(terrain) #needed to pre-allocate space for new columns
  set(terrain, j = 'pixelGroup', value = getValues(pixelGroupMap))

  #if index is present, we need pixelIndices too
  if (!is.null(index)) {
    set(terrain, j = 'pixelIndex', value = 1:ncell(pixelGroupMap))
    terrain <- terrain[!is.na(get(names(terrainRasters)))]
    terrain <- cohortDataWide[terrain, on = c("pixelGroup")]
    terrain <- data.table::setnafill(terrain, fill = 0, cols = names(cohortDataWide))
    set(terrain, , 'pixelGroup', NULL)
    #NAs from non-treed data

    #need to grab climate layers for each year
    climate <- lapply()


 } else {
   set(terrain, j = 'MDC', value = getValues(climateRasters))

 }


  #if index is non null, this function is being used for fitting. need to get years
}
