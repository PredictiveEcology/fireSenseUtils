## Tests for castCohortData

library(data.table)

test_that("castCohortData: produces one row per pixel of pixelGroupMap", {
  withr::local_package("terra")

  pixelGroupMap <- terra::rast(nrows = 3, ncols = 3, vals = c(1, 1, 2,
                                                              2, NA, 3,
                                                              3, 3, NA))
  cohortData <- data.table(
    pixelGroup   = c(1L, 1L, 2L, 3L),
    speciesCode  = c("Pice_mar", "Pinu_con", "Pice_mar", "Popu_tre"),
    age          = c(10L, 20L, 50L, 5L),
    B            = c(100, 200, 300, 50)
  )
  lcc <- data.table(
    pixelID            = 1:9,
    nonForest_lowFlam  = 0L,
    nonForest_highFlam = 0L
  )
  out <- castCohortData(
    cohortData    = cohortData,
    pixelGroupMap = pixelGroupMap,
    lcc           = lcc,
    missingLCC    = "nonForest_lowFlam"
  )
  expect_equal(nrow(out), terra::ncell(pixelGroupMap))
  expect_true("pixelID" %in% colnames(out))
  expect_true("youngAge" %in% colnames(out))
  expect_false("standAge" %in% colnames(out))
})

test_that("castCohortData: youngAge is 1 where standAge <= cutoffForYoungAge", {
  withr::local_package("terra")

  pixelGroupMap <- terra::rast(nrows = 1, ncols = 2, vals = c(1, 2))
  cohortData <- data.table(
    pixelGroup  = c(1L, 2L),
    speciesCode = c("sp1", "sp1"),
    age         = c(5L, 100L),
    B           = c(100, 100)
  )
  lcc <- data.table(pixelID = 1:2,
                    nonForest_lowFlam  = 0L,
                    nonForest_highFlam = 0L)
  out <- castCohortData(
    cohortData    = cohortData,
    pixelGroupMap = pixelGroupMap,
    lcc           = lcc,
    missingLCC    = "nonForest_lowFlam",
    cutoffForYoungAge = 15
  )
  setorder(out, pixelID)
  expect_equal(out$youngAge, c(1L, 0L))
})

test_that("castCohortData: assigns missingLCC to pixels with NA pixelGroup but forested in LCC", {
  withr::local_package("terra")

  ## Pixel 2 has NA pixelGroup (so it never appears in cohortData), and lcc
  ## flags it as forested (all nonForest cols are 0); it should be reassigned
  ## to the missingLCC group.
  pixelGroupMap <- terra::rast(nrows = 1, ncols = 2, vals = c(1, NA))
  cohortData <- data.table(
    pixelGroup  = 1L,
    speciesCode = "sp1",
    age         = 5L,
    B           = 100
  )
  lcc <- data.table(pixelID = 1:2,
                    nonForest_lowFlam  = 0L,
                    nonForest_highFlam = 0L)
  out <- castCohortData(
    cohortData    = cohortData,
    pixelGroupMap = pixelGroupMap,
    lcc           = lcc,
    missingLCC    = "nonForest_highFlam"
  )
  setorder(out, pixelID)
  expect_equal(out$nonForest_highFlam[2], 1)
})

test_that("castCohortData: adds 'year' column when year is supplied", {
  withr::local_package("terra")

  pixelGroupMap <- terra::rast(nrows = 1, ncols = 1, vals = 1)
  cohortData <- data.table(pixelGroup = 1L, speciesCode = "sp1",
                           age = 5L, B = 100)
  lcc <- data.table(pixelID = 1L,
                    nonForest_lowFlam  = 0L,
                    nonForest_highFlam = 0L)
  out <- castCohortData(
    cohortData    = cohortData,
    pixelGroupMap = pixelGroupMap,
    lcc           = lcc,
    missingLCC    = "nonForest_lowFlam",
    year          = 2020
  )
  expect_true("year" %in% colnames(out))
  expect_true(all(out$year == 2020))
})

test_that("castCohortData: ageMap fills standAge where cohortData lacks it", {
  withr::local_package("terra")

  pixelGroupMap <- terra::rast(nrows = 1, ncols = 2, vals = c(1, NA))
  ageMap        <- terra::rast(nrows = 1, ncols = 2, vals = c(NA, 80))
  cohortData <- data.table(pixelGroup = 1L, speciesCode = "sp1",
                           age = 5L, B = 100)
  ## Pixel 2 (NA pixelGroup, non-forest LCC) gets age from ageMap
  lcc <- data.table(pixelID = 1:2,
                    nonForest_lowFlam  = c(0L, 1L),
                    nonForest_highFlam = 0L)
  out <- castCohortData(
    cohortData    = cohortData,
    pixelGroupMap = pixelGroupMap,
    lcc           = lcc,
    ageMap        = ageMap,
    missingLCC    = "nonForest_lowFlam",
    cutoffForYoungAge = 15
  )
  setorder(out, pixelID)
  expect_equal(out$youngAge, c(1L, 0L))
})
