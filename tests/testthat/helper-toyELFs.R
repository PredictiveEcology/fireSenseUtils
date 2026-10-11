## Toy ELF maps for the tests: nothing is downloaded.
##
## One common grid of 4 rows x 12 columns of 5 km cells (one cell = 2500 ha) in EPSG:3978.
## Each ELF's core (value 2) is three columns wide, with a one-column buffer (value 1) on each
## side that exists, and 0 elsewhere:
##
##   ELF     core cols   buffer cols
##   3.1.1   1-3         4
##   3.1.2   4-6         3, 7
##   5.1     7-9         6, 10
##   1.1     10-12       9          (ecozone 1: permanently left out of fireSense)

toyGrid <- function() {
  terra::rast(nrows = 4, ncols = 12, xmin = 0, xmax = 60000, ymin = 0, ymax = 20000,
              crs = "EPSG:3978", vals = 0)
}

toyELFraster <- function(name, coreCols) {
  r <- toyGrid()
  col <- terra::colFromCell(r, seq_len(terra::ncell(r)))
  r[col %in% c(min(coreCols) - 1, max(coreCols) + 1)] <- 1
  r[col %in% coreCols] <- 2
  names(r) <- name
  r
}

toyELFs <- function() {
  cores <- list("3.1.1" = 1:3, "3.1.2" = 4:6, "5.1" = 7:9, "1.1" = 10:12)
  rasWhole <- Map(toyELFraster, names(cores), cores)
  ## rasCentered is the ELF alone, 0 -> NA, trimmed (the real one is also reprojected)
  rasCentered <- lapply(rasWhole, function(r) terra::trim(terra::classify(r, cbind(0, NA))))
  list(rasWhole = rasWhole, rasCentered = rasCentered)
}

## Centre of the cell at (row, col) of the toy grid
toyXY <- function(row, col) terra::xyFromCell(toyGrid(), terra::cellFromRowCol(toyGrid(), row, col))

## SIZE_HA 5000 is above the 2500 ha toy pixel and the 50 ha escape size, so a point escapes by default.
toyPoints <- function(row, col, YEAR, CAUSE, SIZE_HA = 5000) {
  terra::vect(toyXY(row, col), type = "points", crs = "EPSG:3978",
              atts = data.frame(YEAR = YEAR, CAUSE = CAUSE, SIZE_HA = SIZE_HA))
}

## A rectangle covering whole cells: rows r1..r2, cols c1..c2
toyPoly <- function(r1, r2, c1, c2, YEAR) {
  g <- toyGrid()
  e <- terra::ext(terra::xFromCol(g, c1) - 2500, terra::xFromCol(g, c2) + 2500,
                  terra::yFromRow(g, r2) - 2500, terra::yFromRow(g, r1) + 2500)
  p <- terra::as.polygons(e, crs = "EPSG:3978")
  p$YEAR <- YEAR
  p
}

toyFirePoints <- function() {
  rbind(
    toyPoints(2, 2, 2001L, "L"),
    toyPoints(3, 2, 2001L, "H"),          # human-caused
    toyPoints(1, 5, 2001L, "N"),
    toyPoints(2, 5, 2002L, "L"),
    toyPoints(3, 5, 2003L, "L", SIZE_HA = 100),   # counted, but too small to escape
    toyPoints(4, 5, 1999L, "L"),          # outside fireYears
    toyPoints(2, 8, 2002L, "L"),
    toyPoints(2, 11, 2002L, "L")          # in arctic ELF 1.1
  )
}

toyFirePolys <- function() {
  half <- toyPoly(1, 1, 8, 8, 2002L)
  half <- terra::crop(half, terra::ext(terra::xmin(half), terra::xmin(half) + 2500,
                                       terra::ymin(half), terra::ymax(half))) # 2.5 x 5 km = 1250 ha
  rbind(toyPoly(1, 1, 1, 2, 2001L),      # 2 cells = 5000 ha, in 3.1.1
        toyPoly(3, 4, 5, 5, 2002L),      # 2 cells = 5000 ha, in 3.1.2 (col 5 is in no buffer)
        half)
}
