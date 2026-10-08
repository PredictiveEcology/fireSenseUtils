## A fire whose target buffer does not fit in the landscape is kept, with the largest buffer the
## landscape allows. It used to vanish from bufferToArea() because it never outgrew its target.

circleFire <- function(x, y, r, id) {
  sf::st_buffer(sf::st_sf(FIRE_ID = id, geometry = sf::st_sfc(sf::st_point(c(x, y)), crs = 3978)), r)
}
firePoint <- function(x, y, id) {
  sf::st_sf(FIRE_ID = id, geometry = sf::st_sfc(sf::st_point(c(x, y)), crs = 3978))
}
## 60 x 60 cells of 1 km; east of x = 54 km is outside the study area (NA)
toyRTM <- function() {
  rtm <- terra::rast(nrows = 60, ncols = 60, extent = c(0, 60000, 0, 60000), crs = "EPSG:3978", vals = 1L)
  rtm[terra::xFromCell(rtm, seq_len(terra::ncell(rtm))) > 54000] <- NA
  rtm
}
## small central fire, and a large fire near the west edge whose 10 x buffer cannot fit
toyPolys <- function() {
  rbind(circleFire(30000, 30000, 1500, 1), circleFire(6000, 30000, 4000, 2))
}
toyPoints <- function() rbind(firePoint(30000, 30000, 1), firePoint(6000, 30000, 2))
am <- function(x) 100 * x ## the large fire's target (about 5000 cells) exceeds the 3600-cell landscape

test_that("bufferToArea keeps a fire whose buffer cannot reach its target size", {
  rtm <- toyRTM()
  set.seed(1)
  b <- bufferToArea(toyPolys(), rtm, areaMultiplier = am, field = "FIRE_ID", minSize = 10)
  expect_setequal(b$ids, 1:2)
  big <- b[ids == 2]
  expect_gt(sum(big$buffer == 1L), 0)
  expect_true(all(big$pixelID >= 1 & big$pixelID <= terra::ncell(rtm)))
  expect_equal(anyDuplicated(big$pixelID), 0)
  expect_lt(nrow(big), am(sum(big$buffer == 1L))) ## target not reached: whole landscape used
})

test_that("harmonizeFireData keeps the large fire, buffer inside the study area, small fire unchanged", {
  rtm <- toyRTM()
  pts <- list(year2001 = toyPoints())
  set.seed(1)
  res <- suppressMessages(capture.output(
    out <- harmonizeFireData(list(year2001 = toyPolys()), rtm, pts, areaMultiplier = am, minSize = 10)))
  d <- out$fireBufferedListDT$year2001
  expect_setequal(d$ids, 1:2)
  expect_false(anyNA(rtm[d$pixelID][[1]]))
  expect_gt(sum(d$buffer[d$ids == 2] == 1L), 0)
  ## the small fire's buffer is what it was before the change (recorded from the base code)
  small <- sort(d[ids == 1]$pixelID)
  expect_equal(c(length(small), sum(as.numeric(small)^2)), c(400, 1344613400))
})

test_that("harmonizeFireData counts a fire lost to the buffer step as removed", {
  rtm <- toyRTM()
  pts <- list(year2001 = toyPoints())
  set.seed(1)
  msgs <- testthat::capture_messages(capture.output(
    harmonizeFireData(list(year2001 = toyPolys()), rtm, pts, areaMultiplier = am, minSize = 10)))
  expect_match(paste(msgs, collapse = ""), "there were 2 escaped fires and 0 were removed")
})
