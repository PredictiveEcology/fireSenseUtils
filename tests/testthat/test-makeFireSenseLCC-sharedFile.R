## makeFireSenseLCC() used to hand its own `writeTo` to the land-cover builder, so the builder's
## native-resolution raster lived at the path of the final result, and the function kept reading it
## there block by block. A second job sharing `destinationPath` that restored the cached final file
## (or wrote its own) replaced it mid-read: "[classify] too few values for writing: 0 < N".

skip_if_not_installed("terra")

test_that("makeFireSenseLCC() does not read its source from the shared final path", {
  dest <- withr::local_tempdir()
  finalName <- "rstLCC_shared.tif"
  finalPath <- file.path(dest, finalName)
  ## 240 x 240 source cells (30 m) over 7200 m, so the final 240 m grid is 30 x 30
  ext <- c(0, 7200, 0, 7200)
  mk <- function(n, v) terra::rast(nrows = n, ncols = n, xmin = ext[1], xmax = ext[2], ymin = ext[3],
                                   ymax = ext[4], crs = "EPSG:3978", vals = v)
  fine <- mk(240, 210L)
  to <- mk(30, 1L)

  testthat::local_mocked_bindings(
    .scanfiLCC = function(..., writeTo = NULL, destinationPath = NULL) {
      ## like LandR::prepInputs_SCANFI_LCC_FAO(): write to destinationPath/writeTo when given
      f <- if (is.null(writeTo)) tempfile(fileext = ".tif") else file.path(destinationPath, writeTo)
      handle <- terra::writeRaster(fine, f, overwrite = TRUE, datatype = "INT2S")
      ## another job restores the cached 240 m result over the shared final file
      coarse <- mk(30, 210L)
      unlink(finalPath)
      terra::writeRaster(coarse, finalPath, overwrite = TRUE, datatype = "INT2S")
      handle
    },
    .cwimWetland = function(to, ...) terra::rast(to, vals = 0L)
  )
  out <- makeFireSenseLCC(neededYear = 2020, to = to, destinationPath = dest, writeTo = finalName,
                          lccSource = "SCANFI", scanfiVersion = "V2")
  expect_true(all(terra::values(out$lcc) == 210))
  expect_true(file.exists(finalPath))
  expect_true(file.exists(reproducible::.suffix(finalPath, "_propFlam")))
  expect_equal(dim(terra::rast(finalPath))[1:2], c(30, 30))
  ## nothing but the two final files is left in destinationPath
  expect_setequal(list.files(dest), c(finalName, basename(reproducible::.suffix(finalPath, "_propFlam"))))
})
