globalVariables(c(
  "nonForest_highFlam", "nonForest_lowFlam"
))

#' preparing a time since disturbance map from stand age and fire data
#'
#' @param standAgeMap initial stand age map
#' @param firePolys list of \code{spatialPolygon} objects comprising annual fires
#' @param year the year represented by \code{standAge}
#' @param lcc \code{data.table} with landcover values - \code{landcoverDT}
#'
#' @return a raster layer with values representing time since disturbance
#'
#' @export
#' @importFrom data.table data.table
#' @importFrom fasterize fasterize
#' @importFrom magrittr %>%
#' @importFrom raster getValues raster setValues
#' @importFrom sf st_collection_extract
makeTSD <- function(year, firePolys, standAgeMap, lcc) {
  ## TODO: Make this work with lcc values that aren't hardcoded ie line 33
  ## get particular fire polys in format that can be fasterized
  polysNeeded <- names(firePolys) %in% paste0("year", c(year - 16):year - 1) %>%
    firePolys[.] %>%
    .[lengths(.) > 0] %>% #gets rid of missing years that break function
    lapply(., FUN = sf::st_as_sf) %>%
    lapply(., FUN = function(x){
      x <- x[, "YEAR"]
    }) %>%
    do.call(rbind, .)

  #create background raster with TSD
  initialTSD <-  sf::st_collection_extract(polysNeeded, "POLYGON") %>%
    fasterize(sf = ., raster = standAgeMap, background = year - 16, field = "YEAR", fun = "max") %>%
    setValues(., values = year - getValues(.))

  pixToUpdate <- lcc[nonForest_highFlam > 0 | nonForest_lowFlam > 0,]$pixelID

  standAgeMap[pixToUpdate] <- initialTSD[pixToUpdate]

  return(standAgeMap)
}

#' Iteratively calculate \code{youngAge} column  in FS covariates
#'
#' @param standAgeMap template raster
#' @param years the years over which to iterate
#' @param fireBufferedListDT data.table containing non-annual burn and buffer pixelIDs
#' @param annualCovariates list of data.table objects with pixelID
#'
#' @return a raster layer with unified standAge and time-since-disturbance values
#'
#' @export
#' @importFrom raster raster setValues
#' @importFrom data.table data.table
calcYoungAge <- function(years, annualCovariates, standAgeMap, fireBufferedListDT) {
  #this is safest way to subset given the NULL year
  for (year in years) {
    fires <- fireBufferedListDT[[paste0("year", year)]]
    if (!is.null(fires)) {
      pix <- annualCovariates[[paste0("year", year)]]$pixelID
      ages <- standAgeMap[pix]
      young <- ifelse(ages < 16, 1, 0)
      annualCovariates[[paste0("year", year)]][, youngAge := young]
      burnedPix <- fires[buffer == 1]$pixelID
      standAgeMap[burnedPix] <- 0
    }
    standAgeMap <- setValues(standAgeMap, getValues(standAgeMap) + 1)
  }
  return(annualCovariates)
}
