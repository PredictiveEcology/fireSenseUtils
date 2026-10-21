## The national CWIM wetland layer (5.4 GB) is downloaded by LandR::prepInputs_CWIM() into its
## destinationPath; makeFireSenseLCC() must pass its own destinationPath on, or the file lands in the
## working directory (the same defect as PredictiveEcology/Biomass_borealDataPrep#144).

skip_if_not_installed("terra")

test_that("makeFireSenseLCC() passes destinationPath on to the CWIM wetland layer", {
  target <- terra::rast(nrows = 6, ncols = 6, xmin = 0, xmax = 1800, ymin = 0, ymax = 1800,
                        crs = "EPSG:3978", vals = 1L)
  dPath <- withr::local_tempdir()
  seen <- new.env()
  testthat::local_mocked_bindings(
    .scanfiLCC = function(...) terra::rast(target, vals = 210L),
    .cwimWetland = function(to, ...) { seen$args <- list(...); terra::rast(to, vals = 0L) },
    .wetlandToLCC = function(lcc, wet, ...) lcc
  )
  makeFireSenseLCC(neededYear = 2020, to = target, destinationPath = dPath,
                   lccSource = "SCANFI", scanfiVersion = "V2")
  expect_identical(seen$args$destinationPath, dPath)
})
