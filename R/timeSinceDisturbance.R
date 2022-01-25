utils::globalVariables(c(
  "sumRows"
))
#' preparing a time since disturbance map from stand age and fire data
#'
#' @param standAgeMap initial stand age map
#' @param firePolys list of \code{spatialPolygon} objects comprising annual fires
#' @param year the year represented by \code{standAge}
#' @param lcc \code{data.table} with landcover values - \code{landcoverDT}
#' @inheritParams castCohortData
#'
#' @return a raster layer with values representing time since disturbance
#'
#' @export
#' @importFrom data.table data.table
#' @importFrom fasterize fasterize
#' @importFrom magrittr %>%
#' @importFrom raster getValues raster setValues
#' @importFrom sf st_as_sf st_collection_extract
makeTSD <- function(year, firePolys, standAgeMap, lcc, cutoffForYoungAge = 15) {
  ## get particular fire polys in format that can be fasterized
  polysNeeded <- names(firePolys) %in% paste0("year", c(year - cutoffForYoungAge - 1):year - 1) %>%
    firePolys[.] %>%
    .[lengths(.) > 0] %>% # gets rid of missing years that break function
    lapply(., FUN = sf::st_as_sf) %>%
    lapply(., FUN = function(x) {
      x <- x[, "YEAR"]
    }) %>%
    do.call(rbind, .)

  # create background raster with TSD
  initialTSD <- sf::st_collection_extract(polysNeeded, "POLYGON") %>%
    fasterize(
      sf = ., raster = standAgeMap, background = year - cutoffForYoungAge - 1,
      field = "YEAR", fun = "max"
    ) %>%
    setValues(., values = year - getValues(.))

  nfLCC <- names(lcc)[!names(lcc) %in% "pixelID"]
  lcc[, sumRows := rowSums(.SD), .SDcol = nfLCC]
  pixToUpdate <- lcc[sumRows > 0]$pixelID
  lcc[, sumRows := NULL]


  falseYoungs <- standAgeMap[pixToUpdate] <= cutoffForYoungAge | is.na(standAgeMap[pixToUpdate])
  trueYoungs <- initialTSD[pixToUpdate] <= cutoffForYoungAge
  standAgeMap[pixToUpdate[falseYoungs]] <- cutoffForYoungAge + 1
  standAgeMap[pixToUpdate[trueYoungs]] <- initialTSD[pixToUpdate[trueYoungs]]

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
#' @importFrom raster raster setValues
#' @importFrom data.table data.table
calcYoungAge <- function(years, annualCovariates, standAgeMap, fireBufferedListDT,
                         cutoffForYoungAge = 15) {
  # this is safest way to subset given the NULL year
  for (year in years) {
    yearChar <- paste0("year", year)
    ann <- annualCovariates[[yearChar]] # no copy made
    fires <- fireBufferedListDT[[yearChar]]
    if (!is.null(fires)) {
      set(ann, NULL, "youngAge", as.integer(standAgeMap[ann$pixelID] <= cutoffForYoungAge) |
        is.na(standAgeMap[ann$pixelID])) ## cannot have NAs

      # pix <- annualCovariates[[yearChar]]$pixelID
      # ages <- standAgeMap[pix]
      # young <- ifelse(ages <= cutoffForYoungAge, 1, 0)
      # annualCovariates[[yearChar]][, youngAge := as.integer(standAgeMap[][pixelID] <= cutoffForYoungAge)]
      burnedPix <- fires$pixelID[fires$buffer == 1]
      standAgeMap[burnedPix] <- 0
    }
    standAgeMap <- setValues(standAgeMap, getValues(standAgeMap) + 1)
  }
  return(annualCovariates)
}
