utils::globalVariables(c(
  "nFires", "cell"
))



#' prepare covariate table with ignition year, fuel class, climate value, and land cover
#'
#' @param years character vector of fire years with FS notation e.g. year2002
#' @param fuel raster brick of aggregated fuel classes
#' @param LCC raster brick of aggregated LCC classes
#' @param climate raster stack of climate layers with names matching `years`
#' @param climVar the name of the climate variable
#' @param fires list of spatial points representing annual ignitions
#' @return a data.frame with cell numbers, ignitions, and covariates for each year
#'
#' @export
#' @importFrom data.table rbindlist data.table setnames as.data.table
#' @importFrom terra extract rast values ncell
#' @importFrom sf %>%
#' @importFrom stats na.omit
stackAndExtract <- function(years, fuel, LCC, climate, climVar, fires) {
  ignitionYears <- lapply(years, FUN = function(year, climateRas = climate, climvar = climVar,
                                                LCCras = LCC, fuelRas = fuel, ignitions = fires) {

    climate <- climate[[year]] # get single climate layer
    ignitions <- ignitions[paste0("year", ignitions$YEAR) %in% year, ] # get annual ignitions
    names(climate) <- climVar

    yearCovariates <- c(climateRas, LCCras, fuelRas)

    # get cells with ignitions and aggregate by repeated ignitions
    ignitionDT <- extract(x = yearCovariates, y = ignitions, cells = TRUE) %>%
      as.data.table(.) %>% #some ignitions are not in cells
      .[, .(ignitions = .N), .(cell)]

    # get covariate values of all cells
    noIgnitionsDT <- values(yearCovariates) %>%
      as.data.table(.) %>%
      .[, cell := 1:ncell(yearCovariates)] %>%
      na.omit(.)
    # join and assign 0 to non-ignited pixels
    ignitionYear <- ignitionDT[noIgnitionsDT, on = c("cell")]
    ignitionYear[is.na(ignitions), ignitions := 0]

    ignitionYear[, year := gsub("year", "", year)]

    return(ignitionYear)
  })
  ignitionYears <- rbindlist(ignitionYears)
  return(ignitionYears)
}
