## makeFireSenseLCC()'s `lcc` is a land-cover class map, but it came back as a floating-point raster
## (terra::ifel() and the wetland recoding make doubles). LandR::defineFlammable() refuses such a map
## ("LandCoverClassifiedMap must be an integer"), which stopped fireSense_dataPrepPredict's
## .inputObjects() in a run without fireSense_dataPrepFit (2026-11-02); dataPrepFit converts it first.

skip_if_not_installed("terra")

test_that("makeFireSenseLCC() returns an integer land-cover map that defineFlammable() accepts", {
  dest <- withr::local_tempdir()
  ext <- c(0, 7200, 0, 7200)
  mk <- function(n, v) terra::rast(nrows = n, ncols = n, xmin = ext[1], xmax = ext[2], ymin = ext[3],
                                   ymax = ext[4], crs = "EPSG:3978", vals = v)
  fine <- mk(240, rep(c(210L, 220L, 20L), length.out = 240 * 240))   # 20: water, never flammable
  to <- mk(30, 1L)
  testthat::local_mocked_bindings(
    .scanfiLCC = function(...) fine,
    .cwimWetland = function(to, ...) terra::rast(to, vals = 0L)
  )
  for (writeTo in list(NULL, "rstLCC_int.tif")) {
    out <- makeFireSenseLCC(neededYear = 2020, to = to, destinationPath = dest, writeTo = writeTo,
                            lccSource = "SCANFI", scanfiVersion = "V2")
    expect_true(terra::is.int(out$lcc))
    expect_no_error(LandR::defineFlammable(out$lcc, nonFlammClasses = fireSenseNonflammableLCC, to = to))
  }
})
