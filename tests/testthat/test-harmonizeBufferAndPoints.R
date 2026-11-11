## harmonizeBufferAndPoints() moves an ignition point that is not in its fire to a nearby fire cell
## with enough flammable neighbours, drops points whose fire has no buffer, and returns NULL for a
## year with no buffered cells. 20 x 20 cells of 100 m.

hbpRaster <- function(vals = 1L) {
  terra::rast(nrows = 20, ncols = 20, extent = c(0, 2000, 0, 2000), crs = "EPSG:3978", vals = vals)
}
hbpPoints <- function(x, y, id, crs = 3978) {
  sf::st_as_sf(data.frame(FIRE_ID = id, x = x, y = y), coords = c("x", "y"), crs = 3978) |>
    sf::st_transform(crs)
}
fireCells <- function(ras, cols, rows) {
  terra::cellFromRowColCombine(ras, rows, cols)
}

test_that("a point outside its fire moves into it; a point without a fire is dropped", {
  ras <- hbpRaster()
  inFire <- fireCells(ras, 5:10, 5:10)
  buff <- data.table::data.table(pixelID = c(inFire, fireCells(ras, 3:12, 3)), ids = 1L,
                                 buffer = c(rep(1L, length(inFire)), rep(0L, 10)))
  ## fire 1's point is 1 km east of the fire; fire 2 has no buffer; points given in lat/long
  cent <- hbpPoints(c(1850, 150), c(1250, 150), c(1L, 2L), crs = 4326)

  out <- harmonizeBufferAndPoints(list(year2001 = cent), list(year2001 = buff), ras)[[1]]

  expect_identical(out$FIRE_ID, 1L)
  expect_true(sf::st_crs(out) == sf::st_crs(ras))
  expect_true(terra::cellFromXY(ras, sf::st_coordinates(out)) %in% inFire)
})

test_that("in a fire with no well-connected cell, the point goes to the best one found", {
  ## a one-cell-wide burn in a non-flammable landscape: no cell has 7 flammable neighbours
  line <- fireCells(hbpRaster(), 2:19, 10)
  ras <- hbpRaster(0L)
  ras[line] <- 1L
  buff <- data.table::data.table(pixelID = line, ids = 1L, buffer = 1L)
  cent <- hbpPoints(1950, 150, 1L)

  out <- suppressWarnings(capture.output(
    res <- harmonizeBufferAndPoints(list(cent), list(buff), ras)[[1]]))

  expect_true(terra::cellFromXY(ras, sf::st_coordinates(res)) %in% line)
})

test_that("a year with no buffered cells gives NULL; points already in their fire are kept", {
  ras <- hbpRaster()
  inFire <- fireCells(ras, 5:10, 5:10)
  buff <- data.table::data.table(pixelID = inFire, ids = 1L, buffer = 1L)
  cent <- hbpPoints(750, 1250, 1L)

  out <- harmonizeBufferAndPoints(list(a = cent, b = cent), list(a = buff, b = buff[0]), ras)

  expect_null(out$b)
  expect_equal(sf::st_coordinates(out$a), sf::st_coordinates(cent))
})

test_that("an id column ending in y or starting with x (e.g. `.poly`) is kept when a point moves", {
  ras <- hbpRaster()
  inFire <- fireCells(ras, 5:10, 5:10)
  buff <- data.table::data.table(pixelID = inFire, ids = 1L, buffer = 1L)
  for (idCol in c(".poly", "xid", "FIRE_ID", "city")) {
    cent <- hbpPoints(1850, 1250, 1L)
    names(cent)[names(cent) == "FIRE_ID"] <- idCol
    out <- harmonizeBufferAndPoints(list(cent), list(buff), ras, idCol = idCol)[[1]]
    expect_identical(out[[idCol]], 1L)
    expect_true(terra::cellFromXY(ras, sf::st_coordinates(out)) %in% inFire)
  }
})

test_that("snappedCentroids() moves a centroid that is off its fire and keeps the `.poly` id", {
  ras <- hbpRaster()
  ## a U-shaped fire: its centroid is in the gap, which is not burnable
  uCells <- c(fireCells(ras, 5:10, 5:10))
  gap <- fireCells(ras, 7:8, 5:8)
  ras[gap] <- NA
  uCells <- setdiff(uCells, gap)
  r <- terra::rast(ras); r[] <- NA; r[uCells] <- 1L
  fire <- terra::as.polygons(r, dissolve = TRUE) |> terra::aggregate()
  fire$.poly <- 1L
  fire$YEAR <- 2001L
  cl <- data.table::data.table(.poly = 1L, pixelID = uCells)

  out <- snappedCentroids(fire, cl, ras)

  expect_identical(out$.poly, 1L)
  expect_true(terra::cellFromXY(ras, as.matrix(out[, c("x", "y")])) %in% uCells)
})
