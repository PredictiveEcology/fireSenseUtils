utils::globalVariables(c(
  "nFires", "cells", "simSize"
))



#' prepare covariate table with ignition year, fuel class, climate value, and land cover
#'
#' @param years character vector of fire years with FS notation e.g. year2002
#' @param fuel raster brick of aggregated fuel classes
#' @param LCC raster brick of aggregated LCC classes
#' @param indices a list of data.tables output by \code{bufferIgnitionPoints}
#' @param climate raster stack of climate layers with names matching \code{years}
#' @param climVar the name of the climate variable
#' @return a data.frame with cell numbers, ignitions, and covariates for each year
#'
#' @export
#' @importFrom data.table rbindlist data.table setnames
#' @importFrom raster stack
#' @importFrom stats na.omit
stackAndExtract <- function(years, climate, fuel, LCC, indices, climVar) {
  ignitionYears <- lapply(years, FUN = function(year, climateRas = climate, climvar = climVar,
                                                LCCras = LCC, fuelRas = fuel, igIndex = indices) {
    climate <- climate[[year]]
    igIndex <- igIndex[[year]]
    dtNames <- c(climvar, names(LCCras), names(fuelRas))
    yearCovariates <- raster::stack(climateRas, LCCras, fuelRas)
    annDT <- raster::getValues(yearCovariates)
    annDT <- as.data.table(annDT)
    setnames(annDT, dtNames)
    annDT[, pixelID := 1:ncell(yearCovariates)]
    annDT <- annDT[igIndex, on = c("pixelID")]
    annDT <- na.omit(annDT)
    return(annDT)
  })
  ignitionYears <- rbindlist(ignitionYears)
  return(ignitionYears)
}
