## makeFireSenseLCC()'s base land-cover step (LandR's SCANFI/NTEMS crop, mask and forest-land update,
## ~10 min per study area and year) was only cached as part of the whole function, so any change to
## the later dominant-flammable step (its inputs or its code) redid it. It is now cached on its own,
## keyed on its own arguments and the code of the functions it reaches.

skip_if_not_installed("terra")

test_that("the base land-cover step is reused when only the later step changes", {
  withr::local_options(reproducible.cachePath = withr::local_tempdir())
  dest <- withr::local_tempdir()
  ## a 60 x 60 source grid (30 m) aggregating to 6 x 6 target cells (300 m): conifer (210) everywhere,
  ## with every other source row water (20) in the right half, so those target cells are half flammable
  v <- matrix(210L, 60, 60)
  v[seq(2, 60, 2), 31:60] <- 20L
  src <- terra::rast(nrows = 60, ncols = 60, xmin = 0, xmax = 1800, ymin = 0, ymax = 1800, crs = "EPSG:3978",
                     vals = as.vector(t(v)))
  to <- terra::rast(nrows = 6, ncols = 6, xmin = 0, xmax = 1800, ymin = 0, ymax = 1800, crs = "EPSG:3978",
                    vals = 1L)
  calls <- 0L
  scratch <- withr::local_tempdir()
  testthat::local_mocked_bindings(
    .scanfiLCC = function(year, ...) {
      calls <<- calls + 1L
      ## like LandR's builder, a raster on disk in a scratch file that does not outlive the run
      terra::writeRaster(src, file.path(scratch, paste0("base_", year, "_", calls, ".tif")),
                         overwrite = TRUE, datatype = "INT2S")
    },
    .cwimWetland = function(to, ...) terra::rast(to, vals = 0L)
  )
  run <- function(flammabilityThreshold, neededYear = 2020)
    makeFireSenseLCC(neededYear = neededYear, to = to, destinationPath = dest, lccSource = "SCANFI",
                     scanfiVersion = "V2", flammabilityThreshold = flammabilityThreshold)

  first <- run(0.1)
  expect_equal(calls, 1L)
  expect_equal(terra::values(first$lcc, mat = FALSE), rep(210, 36))
  ## the base step's own file can vanish (terra scratch, another job's cleanup): the cache keeps a copy
  unlink(list.files(scratch, full.names = TRUE))

  ## a later-step input changes: the base step is not redone, and the result is still right
  second <- run(0.9)
  expect_equal(calls, 1L)
  expect_equal(terra::values(second$flammableProp, mat = FALSE), terra::values(first$flammableProp, mat = FALSE))
  expect_equal(terra::values(second$lcc, mat = FALSE), rep(c(210, 210, 210, 0, 0, 0), 6))

  ## one of its own arguments changes: it is redone
  run(0.1, neededYear = 2015)
  expect_equal(calls, 2L)
})
