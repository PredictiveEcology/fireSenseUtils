utils::globalVariables(c(
  "stat"
))

#' Download and prepare fire data from National Fire Database
#'
#' @param historicalMDC raster stack of historical MDC
#' @param projectedMDC raster stack of projected MDC
#' @param flammableRTM an optional raster of flammable pixels to subset data
#' @param Ylimits the upper and lower MDC range for the plot
#' @param firstHistoricalYear the earliest year of historical data
#' @param firstProjectedYear the earliest year of projected data
#' @return a ggplot object
#'
#' @export
#' @importFrom data.table data.table rbindlist
#' @importFrom stats median
#' @importFrom terra ncell nlyr rast values
#'
#' @examples
#' \dontrun{
#' compareMDC(
#'   historicalMDC = simOutPreamble$historicalClimateRasters$MDC,
#'   projectedMDC = simOutPreamble$projectedClimateRasters$MDC,
#'   flammableRTM = fSsimDataPrep$flammableRTM
#' )
#' }
compareMDC <- function(historicalMDC, projectedMDC, flammableRTM = NULL, Ylimits = c(80, 220),
                       firstHistoricalYear = 2001, firstProjectedYear = 2011) {
  if (is(historicalMDC, "Raster")) {
    historicalMDC <- rast(historicalMDC)
  }

  if (is(projectedMDC, "Raster")) {
    projectedMDC <- rast(projectedMDC)
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
        dt <- data.table("year" = year, "MDC" = as.vector(values(x[[year]])))
        if (!is.null(index)) {
          dt <- dt[index]
        }
        dt <- dt[, .(MDC = median(MDC, na.rm = TRUE)), .(year)]
        return(dt)
      })
      medVal <- rbindlist(medVal)
      return(medVal)
    }

    proj <- valfun(projectedMDC, flamMap = flammableRTM)
    proj$year <- proj$year + firstProjectedYear - 1
    proj$year
    proj$stat <- "projected"
    hist <- valfun(historicalMDC, flamMap = flammableRTM)
    hist$year <- hist$year + firstHistoricalYear - 1
    hist$stat <- "historical"
    MDC <- rbind(hist, proj)
    ggplot2::ggplot(data = MDC, ggplot2::aes(y = MDC, x = year, color = stat)) +
      ggplot2::geom_line() +
      ggplot2::geom_smooth() +
      ggplot2::coord_cartesian(ylim = Ylimits)
  }
}
