## Tests for makeTSD
##
## NOTE: makeTSD currently has a fragile path when fireRaster has NA values
## in pixels that lcc flags as nonForested (`pixToUpdate`). The NA propagates
## into the `pixToUpdate[initialTSD[pixToUpdate] <= cutoffForYoungAge]` subset
## and produces a length mismatch on `:= trueAges`. These tests therefore
## supply fireRasters without NAs over the active `pixToUpdate` set.

library(data.table)

test_that("makeTSD: fireRaster path returns 2020 - fireRaster for recent burns", {
  withr::local_package("terra")

  ## 3 pixels burned in 2005, 2010, 2018 (no NAs)
  fireRaster   <- terra::rast(nrows = 1, ncols = 3, vals = c(2005, 2010, 2018))
  standAgeMap  <- terra::rast(nrows = 1, ncols = 3, vals = c(50, 30, 10))
  lcc <- data.table(pixelID = 1:3,
                    nonForest_lowFlam  = 1L,
                    nonForest_highFlam = 0L)
  out <- makeTSD(year = 2020, fireRaster = fireRaster, standAgeMap = standAgeMap,
                 lcc = lcc, cutoffForYoungAge = 15)
  expect_s4_class(out, "SpatRaster")

  vals <- terra::values(out, mat = FALSE)
  ## pixel 3 (burned 2018) -> TSD = 2; "young" so age comes from initialTSD
  expect_equal(vals[3], 2)
  ## pixel 2 (burned 2010) -> TSD = 10; "young" so age comes from initialTSD
  expect_equal(vals[2], 10)
})

test_that("makeTSD: errors when neither firePolys nor fireRaster is provided", {
  withr::local_package("terra")

  standAgeMap <- terra::rast(nrows = 1, ncols = 1, vals = 10)
  lcc <- data.table(pixelID = 1L,
                    nonForest_lowFlam  = 0L,
                    nonForest_highFlam = 0L)
  expect_error(
    makeTSD(year = 2020, standAgeMap = standAgeMap, lcc = lcc),
    regexp = "firePolys or fireRaster"
  )
})

test_that("makeTSD: names output layer 'timeSinceDisturbance<year>'", {
  withr::local_package("terra")

  fireRaster   <- terra::rast(nrows = 1, ncols = 1, vals = 2018)
  standAgeMap  <- terra::rast(nrows = 1, ncols = 1, vals = 10)
  lcc <- data.table(pixelID = 1L,
                    nonForest_lowFlam  = 1L,
                    nonForest_highFlam = 0L)
  out <- makeTSD(year = 2020, fireRaster = fireRaster, standAgeMap = standAgeMap,
                 lcc = lcc, cutoffForYoungAge = 15)
  expect_equal(names(out), "timeSinceDisturbance2020")
})

test_that("makeTSD: future burns (year < fireRaster) are clamped to cutoff + 1", {
  withr::local_package("terra")

  ## Pixel disturbed in 2025 (the future relative to neededYear 2020) -> initialTSD
  ## is negative, then clamped to cutoff + 1. Standage is preserved (not "young").
  fireRaster   <- terra::rast(nrows = 1, ncols = 1, vals = 2025)
  standAgeMap  <- terra::rast(nrows = 1, ncols = 1, vals = 50)
  lcc <- data.table(pixelID = 1L,
                    nonForest_lowFlam  = 1L,
                    nonForest_highFlam = 0L)
  out <- makeTSD(year = 2020, fireRaster = fireRaster, standAgeMap = standAgeMap,
                 lcc = lcc, cutoffForYoungAge = 15)
  expect_false(terra::values(out, mat = FALSE)[1] <= 15)
})

test_that("makeTSD: pixels not in lcc are set to NA in output", {
  withr::local_package("terra")

  ## Pixel 2 is missing from lcc -> should be NA in output (non-flammable)
  fireRaster   <- terra::rast(nrows = 1, ncols = 2, vals = c(2018, 2018))
  standAgeMap  <- terra::rast(nrows = 1, ncols = 2, vals = c(10, 10))
  lcc <- data.table(pixelID = 1L,
                    nonForest_lowFlam  = 1L,
                    nonForest_highFlam = 0L)
  out <- makeTSD(year = 2020, fireRaster = fireRaster, standAgeMap = standAgeMap,
                 lcc = lcc, cutoffForYoungAge = 15)
  expect_true(is.na(terra::values(out, mat = FALSE)[2]))
})
