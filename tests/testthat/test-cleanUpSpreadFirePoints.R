## cleanUpSpreadFirePoints() moves an ignition point on a non-flammable pixel to the nearest
## flammable pixel of its fire, and drops a fire with no flammable pixel at all. It used
## terra::extract()'s ID -- the row of the point -- as the fire ID, so with real fire IDs neither
## happened.

test_that("cleanUpSpreadFirePoints fixes points by fire ID, not by row", {
  ## 10 x 10 cells of 100 m; columns 6-10 are a lake
  rtm <- terra::rast(nrows = 10, ncols = 10, extent = c(0, 1000, 0, 1000), crs = "EPSG:3978",
                     vals = rep(c(1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L), 10))
  cell <- function(x, y) terra::cellFromXY(rtm, cbind(x, y))
  bufferDT <- data.table::data.table(
    pixelID = c(cell(150, 950), cell(450, 550), cell(650, 550), cell(850, 150)),
    ids = c(11L, 12L, 12L, 13L), buffer = 1L)
  ## fire 11 fine; fire 12's point is on its lake cell; fire 13 is all lake
  pts <- sf::st_as_sf(data.frame(FIRE_ID = c(11L, 12L, 13L), x = c(150, 650, 850),
                                 y = c(950, 550, 150)), coords = c("x", "y"), crs = 3978)

  out <- cleanUpSpreadFirePoints(pts, bufferDT, rtm, idCol = "FIRE_ID")

  p <- out$SpatialPoints
  expect_setequal(p$FIRE_ID, c(11L, 12L))
  expect_equal(unname(sf::st_coordinates(p[p$FIRE_ID == 12L, ])[1, ]), c(450, 550))
  expect_equal(unname(sf::st_coordinates(p[p$FIRE_ID == 11L, ])[1, ]), c(150, 950))
  expect_false(13L %in% out$FireBuffered$ids)
})
