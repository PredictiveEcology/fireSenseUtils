#' Converts stacks of climate rasters to data.table and optionally subsets to index
#'
#' @param historicalClimateRasters named list of `SpatRaster` objects
#' @param Index optional list of `data.table` objects named by `fireYear` and
#' containing fire buffer indices
#'
#' @return a long-layout `data.table` of climate values in each pixel and year, at full precision: the
#' x1000 integer storage ([toX1000()]) is where they are rounded, as `fireSense_spreadPredict` does, so a
#' value is stored as the same integer in the fit and in the prediction
#'
#' @export
#' @importFrom data.table as.data.table melt.data.table merge.data.table setnames
#' @importFrom terra as.data.frame
#' @rdname climateRasterToDataTable
climateRasterToDataTable <- function(historicalClimateRasters, Index = NULL) {
  climVar <- names(historicalClimateRasters)
  wideClimate <- lapply(historicalClimateRasters, terra::as.data.frame, cells = TRUE)
  wideClimate <- lapply(wideClimate, as.data.table)
  out <- lapply(climVar, FUN = function(var) {
    longClimate <- wideClimate[[var]]
    longClimate <- melt.data.table(longClimate,
      id.vars = "cell",
      value.name = var, variable.name = "year"
    )
  })

  if (length(out) > 1) {
    out$by <- c("cell", "year")
    out <- do.call(data.table::merge.data.table, out)
  } else {
    out <- out[[1]]
  }
  setnames(out, old = "cell", new = "pixelID")

  if (!is.null(Index)) {
    out <- out[pixelID %in% Index]
  }
  gc()
  return(out)
}
