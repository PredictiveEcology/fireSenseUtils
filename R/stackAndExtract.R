utils::globalVariables(c(
  "nFires", "cells"))



#' prepare covariate table with ignition year, fuel class, climate value, and land cover
#'
#' @param years character vector of fire years with FS notation e.g. year2002
#' @param fuel raster brick of aggregated fuel classes
#' @param LCC raster brick of aggregated LCC classes
#' @param climate raster stack of climate layers with names matching \code{years}
#' @param climVar the name of the climate variable
#' @param fires list of spatial points representing annual ignitions
#' @return a data.frame with cell numbers, ignitions, and covariates for each year
#'
#' @export
#' @importFrom data.table rbindlist data.table setnames as.data.table
#' @importFrom raster stack extract
#' @importFrom stats na.omit
#' @importFrom magrittr %>%
stackAndExtract <- function(years, fuel, LCC, climate, climVar, fires) {

  ignitionYears <- lapply(years, FUN = function(year, climateRas = climate, climvar = climVar,
                                                LCCras = LCC, fuelRas = fuel, ignitions = fires) {

    climate <- climate[[year]] #get single climate layer
    ignitions <- ignitions[paste0("year", ignitions$YEAR) %in% year,] #get annual ignitions
    names(climate) <- climVar

    yearCovariates <- raster::stack(climateRas, LCCras, fuelRas)

    #get cells with ignitions and aggregate by repeated ignitions
    ignitionDT <- extract(x = yearCovariates, y = ignitions, cellnumbers = TRUE) %>%
     as.data.table(.) %>%
      .[, .(ignitions = .N), .(cells)]

    #get covariate values of all cells
    noIgnitionsDT <- raster::getValues(yearCovariates) %>%
      as.data.table(.) %>%
      .[, cells := 1:ncell(yearCovariates)] %>%
      na.omit(.)

    #join and assign 0 to non-ignited pixels
    ignitionYear <- ignitionDT[noIgnitionsDT, on = c('cells')]
    ignitionYear[is.na(ignitions), ignitions := 0]

    return(ignitionYear)
  })
  ignitionYears <- rbindlist(ignitionYears)
  return(ignitionYears)
}
