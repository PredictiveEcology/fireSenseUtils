utils::globalVariables(c(
  "nFires", "cells", "simSize"
))



#' prepare covariate table with ignition year, fuel class, climate value, and land cover
#'
#' @param years character vector of fire years with FS notation e.g. year2002
#' @param fuel raster brick of aggregated fuel classes
#' @param LCC raster brick of aggregated LCC classes
#' @param fires list of spatial points representing annual ignitions
#' @param climate raster stack of climate layers with names matching \code{years}
#' @param climVar the name of the variables in \code{climate}
#' @param flamIndex index of flammable cells - necessary when pixels are not aggregated
#'
#' @return a data.frame with cell numbers, ignitions, and covariates for each year
#'
#' @export
#' @importFrom magrittr %>%
#' @importFrom data.table rbindlist data.table
#' @importFrom raster extract brick
#' @importFrom stats na.omit
stackAndExtract <- function(years, climate, fuel, LCC, fires, climVar, flamIndex) {

  ignitionYears <- lapply(years, FUN = function(year) {

    fires <- fires[paste0("year", fires$YEAR) %in% year,]
    climate <- climate[[year]]
    names(climate) <- climVar
    dtNames <- c(names(climate), names(fuel), names(LCC))
    yearCovariates <- raster::brick(climate, fuel, LCC)
    fireDT <- raster::extract(x = yearCovariates, y = fires, cellnumbers = TRUE) %>%
      as.data.table(.)
    fireDT <- fireDT[, ignition := 1]
    noFireDT <- raster::getValues(yearCovariates) %>%
      as.data.table(.) %>%
      .[, cells := 1:ncell(yearCovariates)] %>%
      na.omit(.)
    noFireDT <- noFireDT[cells %in% flamIndex]
    ignitionYear <- fireDT[noFireDT, on = c('cells', dtNames)]
    ignitionYear[is.na(ignition), ignition := 0]
    return(ignitionYear)
  })
  ignitionYears <- rbindlist(ignitionYears)
  return(ignitionYears)
}
