#' Converts stacks of climate rasters to data.table and optionally subsets to index
#'
#' @param historicalClimateRasters named list of `SpatRaster` objects
#' @param Index optional. Either a list of `data.table` objects named by `fireYear` (the raster layer
#' names, e.g. `year1990`), each with a column `pixelID`: the fire buffer indices, so only those
#' year-and-pixel pairs are returned; or a vector of `pixelID`s kept in every year. Only the cells
#' of `Index` are read from the rasters, not every cell, so a large raster is never made into a table.
#' `NULL` (the default) returns every cell.
#'
#' @return a long-layout `data.table` of climate values in each pixel and year, at full precision: the
#' x1000 integer storage ([toX1000()]) is where they are rounded, as `fireSense_spreadPredict` does, so a
#' value is stored as the same integer in the fit and in the prediction
#'
#' @export
#' @importFrom data.table as.data.table melt.data.table merge.data.table rbindlist setnames
#' @importFrom terra as.data.frame
#' @rdname climateRasterToDataTable
climateRasterToDataTable <- function(historicalClimateRasters, Index = NULL) {
  climVar <- names(historicalClimateRasters)
  pairs <- NULL
  cells <- NULL
  if (is.list(Index)) {
    pairs <- unique(rbindlist(lapply(Index, function(d) as.data.table(d)[, "pixelID"]), idcol = "year"))
    cells <- unique(pairs$pixelID)
  } else if (!is.null(Index)) {
    cells <- unique(Index)
  }
  wideClimate <- lapply(historicalClimateRasters, function(r) {
    if (is.null(cells)) {
      return(terra::as.data.frame(r, cells = TRUE))
    }
    inR <- sort(cells[!is.na(cells) & cells >= 1 & cells <= terra::ncell(r)])
    d <- r[inR] # only these cells are read
    keep <- rowSums(!is.na(d)) > 0 # as as.data.frame(): no row with every layer NA
    data.frame(cell = inR[keep], d[keep, , drop = FALSE], check.names = FALSE)
  })
  wideClimate <- lapply(wideClimate, as.data.table)
  out <- lapply(climVar, FUN = function(var) {
    longClimate <- wideClimate[[var]]
    longClimate <- melt.data.table(longClimate,
      id.vars = "cell",
      value.name = var, variable.name = "year"
    )
  })

  out <- Reduce(function(a, b) data.table::merge.data.table(a, b, by = c("cell", "year")), out)
  setnames(out, old = "cell", new = "pixelID")

  if (!is.null(pairs)) {
    pairs$year <- factor(pairs$year, levels = levels(out$year))
    out <- out[pairs, on = c("year", "pixelID"), nomatch = NULL]
  }
  gc()
  return(out)
}
