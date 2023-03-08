utils::globalVariables(c(
  "sumRows"
))
#' preparing a time since disturbance map from stand age and fire data
#'
#' @param standAgeMap initial stand age map
#' @param firePolys list of \code{spatialPolygon} objects comprising annual fires.
#' fireRaster will supersede firePolys if provided
#' @param fireRaster a \code{RasterLayer} with values representing fire years
#' @param year the year represented by \code{standAge}
#' @param lcc \code{data.table} with landcover values - \code{landcoverDT}
#' @inheritParams castCohortData
#'
#' @return a raster layer with values representing time since disturbance
#'
#' @export
#' @importFrom data.table data.table
#' @importFrom terra values rast setValues rasterize vect set.names
#' @importFrom sf %>% st_as_sf st_collection_extract
makeTSD <- function(year, firePolys = NULL, fireRaster = NULL,
                    standAgeMap, lcc, cutoffForYoungAge = 15) {

  if (!is.null(fireRaster)){
    baseYear <- raster(fireRaster)
    baseYear[] <- year
    initialTSD <- baseYear - fireRaster
    initialTSD[initialTSD < 0] <- cutoffForYoungAge + 1
    #these pixels burn in the future - can't infer prior disturbance
  } else if (!is.null(firePolys)) {
    ## get particular fire polys in format that can be fasterized
    polysNeeded <- firePolys[names(firePolys) %in% paste0("year", c(year - cutoffForYoungAge - 1):year - 1)]
    polysNeeded <- polysNeeded[sapply(polysNeeded, length) > 0]
    polysNeeded <- do.call(rbind, polysNeeded)

    # create background raster with TSD
    initialTSD <- rasterize(polysNeeded, y = standAgeMap,
                            background = year - cutoffForYoungAge - 1,
                            field = "YEAR", fun = "max")
    initialTSD <- year - initialTSD
  } else {
    stop("Please provide either firePolys or fireRaster")
  }

  nfLCC <- names(lcc)[!names(lcc) %in% "pixelID"]
  lcc[, sumRows := rowSums(.SD), .SDcol = nfLCC]
  pixToUpdate <- lcc[sumRows > 0]$pixelID
  lcc[, sumRows := NULL]

  standAgeVals <- values(standAgeMap, mat = FALSE)
  #these have no disturbance history but are apparently young
  falseYoungs <- standAgeVals[pixToUpdate] <= cutoffForYoungAge | is.na(standAgeVals[pixToUpdate])
  #disturbance history suggests young
  trueYoungs <- initialTSD[pixToUpdate] <= cutoffForYoungAge

  standAgeVals[pixToUpdate[falseYoungs]] <- cutoffForYoungAge + 1
  #note that by doing this second, pixels in both groups are correctly set to trueYoung
  standAgeVals[pixToUpdate[trueYoungs]] <- values(initialTSD, mat = FALSE)[pixToUpdate[trueYoungs]]
  standAgeMap <- setValues(standAgeMap, standAgeVals)
  set.names(standAgeMap, paste0("timeSinceDisturbance", year))

  return(standAgeMap)
}

#' Iteratively calculate \code{youngAge} column  in FS covariates
#'
#' @param standAgeMap template raster
#' @param years the years over which to iterate
#' @param fireBufferedListDT data.table containing non-annual burn and buffer pixelIDs
#' @param annualCovariates list of data.table objects with pixelID
#' @inheritParams castCohortData
#'
#' @return a raster layer with unified standAge and time-since-disturbance values
#'
#' @export
#' @importFrom terra rast setValues values
#' @importFrom data.table data.table
calcYoungAge <- function(years, annualCovariates, standAgeMap, fireBufferedListDT,
                         cutoffForYoungAge = 15) {
  # this is safest way to subset given the NULL year
  for (year in years) {
    yearChar <- paste0("year", year)
    ann <- annualCovariates[[yearChar]] # no copy made
    fires <- fireBufferedListDT[[yearChar]]
    if (!is.null(fires)) {
      ageVals <- values(standAgeMap, mat = FALSE)
      set(ann, NULL, "youngAge", as.integer(ageVals[ann$pixelID] <= cutoffForYoungAge) |
        is.na(ageVals[ann$pixelID])) ## cannot have NAs
      burnedPix <- fires$pixelID[fires$buffer == 1]
      ageVals[burnedPix] <- 0
      standAgeMap <- setValues(standAgeMap, values = ageVals)
    }
    standAgeMap <- setValues(standAgeMap, values = values(standAgeMap, mat = FALSE) + 1)
  }
  return(annualCovariates)
}
