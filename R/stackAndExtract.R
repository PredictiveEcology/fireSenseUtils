utils::globalVariables(c(
  "nFires", "cell"
))


#' prepare covariate table with ignition year, fuel class, climate value, and land cover
#'
#' @param years character vector of fire years with FS notation e.g. year2002
#' @param fuel raster brick of aggregated fuel classes
#' @param LCC raster brick of aggregated LCC classes
#' @param climate list of raster layers named by climate variable
#'  with raster layer names matching `years`
#' @param fires list of spatial points representing annual ignitions
#' @return a data.frame with cell numbers, ignitions, and covariates for each year
#'
#' @export
#' @importFrom data.table rbindlist data.table setnames as.data.table
#' @importFrom terra extract rast values ncell
#' @importFrom sf %>%
#' @importFrom stats na.omit
stackAndExtract <- function(years, fuel, LCC, climate, fires) {

  ignitionYears <- lapply(years, FUN = function(year, climateRas = climate, fuelRas = fuel,
                                                LCCras = LCC, ignitions = fires) {
    climateVariables <- names(climateRas)
    thisYearsClimate <- lapply(climateRas,
                               FUN = function(x){x[[which(names(x) == year, useNames = TRUE)]]}) |>
      rast()

    ignitions <- ignitions[paste0("year", ignitions$YEAR) %in% year, ] # get annual ignitions
    #TODO: check this works with length = 1 and length > 1. str of lapply changes..

    #rename from year to climate variable
    names(thisYearsClimate) <- climateVariables

    yearCovariates <- c(thisYearsClimate, LCCras, fuelRas)

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
