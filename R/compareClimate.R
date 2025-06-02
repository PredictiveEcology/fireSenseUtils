utils::globalVariables(c(
  "stat"
))

#' Plot the historical and projected climate values for flammable land surface
#'
#' @param historicalClimate SpatRaster of historical climate variable
#' @param projectedClimate SpatRaster of projected climate variable
#' @param flammableRTM an optional raster of flammable pixels to subset data
#' @param Ylimits the upper and lower MDC range for the plot
#' @param firstHistoricalYear the earliest year of historical data
#' @param firstProjectedYear the earliest year of projected data
#' @param climVar the name of the climate variable - for axis label
#' @return a ggplot object
#'
#' @export
#' @importFrom data.table data.table rbindlist
#' @importFrom stats median
#' @importFrom terra ncell nlyr rast values
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#' compareClimate(
#'   historicalClimate = simOutPreamble$historicalClimateRasters$Climate,
#'   projectedClimate = simOutPreamble$projectedClimateRasters$Climate,
#'   flammableRTM = fSsimDataPrep$flammableRTM
#' )
#' }
compareClimate <- function(historicalClimate, projectedClimate, flammableRTM = NULL,
                           Ylimits = NULL, firstHistoricalYear = 2001, firstProjectedYear = 2011,
                           climVar = "climate variable") {

  if (is(historicalClimate, "Raster")) {
    historicalClimate <- rast(historicalClimate)
  }

  if (is(projectedClimate, "Raster")) {
    projectedClimate <- rast(projectedClimate)
  }

  if (is(flammableRTM, "Raster")) {
    flammableRTM <- rast(flammableRTM)
  }

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    valfun <- function(x, flamMap = NULL) {
      years <- 1:nlyr(x)
      if (!is.null(flamMap)) {
        isFlam <- as.vector(values(flamMap))
        flamMap <- c(1:ncell(flamMap))[!is.na(isFlam) & isFlam > 0]
      }
      medVal <- lapply(years, FUN = function(year, index = flamMap) {
        dt <- data.table("year" = year, "Climate" = as.vector(values(x[[year]])))
        if (!is.null(index)) {
          dt <- dt[index]
        }
        dt <- dt[, .(Climate = median(Climate, na.rm = TRUE)), .(year)]
        return(dt)
      })
      medVal <- rbindlist(medVal)
      return(medVal)
    }

    proj <- valfun(projectedClimate, flamMap = flammableRTM)
    proj$year <- proj$year + firstProjectedYear - 1
    proj$year
    proj$stat <- "projected"
    hist <- valfun(historicalClimate, flamMap = flammableRTM)
    hist$year <- hist$year + firstHistoricalYear - 1
    hist$stat <- "historical"
    Climate <- rbind(hist, proj)

    if (!is.null(Ylimits)) {
      if (min(Climate$Climate) == 0) {
        Ymin <- 0
      } else {
        Ymin <- round(min(Climate$Climate) * 0.95)
      }
      Ymax <- round(max(Climate$Climate) * 1.05)
      Ylimits <- c(Ymin, Ymax)
    }

    ggplot(data = Climate, aes(y = Climate, x = year, color = stat)) +
      geom_line() +
      geom_smooth() +
      coord_cartesian(ylim = Ylimits) +
      ggplot2::labs(y = climVar)
  }
}
