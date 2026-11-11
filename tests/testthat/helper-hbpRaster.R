## 20 x 20 cells of 100 m, and the cells of given columns and rows; used by the harmonizeBufferAndPoints and snappedCentroids tests
hbpRaster <- function(vals = 1L) {
  terra::rast(nrows = 20, ncols = 20, extent = c(0, 2000, 0, 2000), crs = "EPSG:3978", vals = vals)
}
fireCells <- function(ras, cols, rows) {
  terra::cellFromRowColCombine(ras, rows, cols)
}
