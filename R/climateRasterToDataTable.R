#' Converts stacks of climate rasters to data.table and optionally subsets to index
#'
#' @param historicalClimateRasters named list of `SpatRaster` objects
#' @param Index optional list of `data.table` objects named by fireYear and
#' containing fire buffer indices
#'
#' @return a long-layout `data.table` of climate values in each pixel and year
#'
#' @export
#' @importFrom data.table as.data.table setnames melt.data.table
#' @importFrom LandR asInteger
#' @importFrom terra values ncell
#' @rdname climateRasterToDataTable
climateRasterToDataTable <- function(historicalClimateRasters, Index = NULL) {

  varName <- names(historicalClimateRasters)
  temp <- as.data.table(values(historicalClimateRasters[[1]]))
  temp[, pixelID := 1:ncell(historicalClimateRasters[[1]][[1]])]
  temp <- temp[pixelID %in% Index]
  temp <- melt.data.table(data = temp, id.vars = "pixelID", value.name = varName,
                          variable.name = "year")
  temp[, (varName) := lapply(.SD, asInteger), .SDcols = varName]
  gc()
  return(temp)
}
