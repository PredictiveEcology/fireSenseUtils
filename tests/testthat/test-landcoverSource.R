## fireSense land cover from SCANFI + CWIM by default, with NTEMS still available.
##
## Biomass_borealDataPrep's default land cover is moving from NTEMS to SCANFI, with the wetland classes 80
## (wetland) and 81 (treed wetland) added from the Canadian Wetland Inventory Map (SCANFI has none). The fire
## models must be fitted on the same land cover the simulation predicts with, so makeFireSenseLCC() follows.
## NTEMS must stay reachable (`lccSource = "NTEMS"`, or options(fireSense.lccSource = "NTEMS")) and unchanged.

skip_if_not_installed("terra")

## a 60 x 60 source grid (30 m) aggregating to 6 x 6 target cells (300 m)
srcGrid <- function(vals) terra::rast(nrows = 60, ncols = 60, xmin = 0, xmax = 1800, ymin = 0, ymax = 1800,
                                      crs = "EPSG:3978", vals = vals)
target <- terra::rast(nrows = 6, ncols = 6, xmin = 0, xmax = 1800, ymin = 0, ymax = 1800, crs = "EPSG:3978",
                      vals = 1L)
## source classes by target cell: rows of target cells are conifer (210), shrub (50), water (20)
srcVals <- function() {
  rowOfTarget <- rep(rep(1:6, each = 10), each = 60)        # 10 source rows (of 60 cells) per target row
  c(210L, 210L, 50L, 50L, 20L, 20L)[rowOfTarget]
}

runLCC <- function(source, seen, wetRows = c(1, 3, 5)) {
  testthat::local_mocked_bindings(
    prepInputs_NTEMS_LCC_FAO = function(...) { seen$ntems <- TRUE; srcGrid(srcVals()) },
    .scanfiLCC = function(...) { seen$scanfi <- TRUE; srcGrid(srcVals()) },
    .cwimWetland = function(to, ...) {
      seen$cwim <- TRUE
      w <- terra::rast(to); terra::values(w) <- as.integer(terra::rowFromCell(w, seq_len(terra::ncell(w))) %in% wetRows)
      w
    }
  )
  makeFireSenseLCC(neededYear = 2020, to = target, destinationPath = withr::local_tempdir(),
                   lccSource = source)$lcc
}
rowsOf <- function(r) matrix(terra::values(r, mat = FALSE), nrow = 6, byrow = TRUE)[, 1]

test_that("NTEMS: the NTEMS layer is used as before, and no wetland layer is applied", {
  seen <- new.env()
  lcc <- runLCC("NTEMS", seen)
  expect_true(isTRUE(seen$ntems)); expect_null(seen$scanfi); expect_null(seen$cwim)
  expect_equal(rowsOf(lcc), c(210, 210, 50, 50, 0, 0))       # water -> 0 (non-flammable)
})

test_that("SCANFI: the SCANFI layer is used, and wet ground becomes 81 (treed) or 80 (not treed)", {
  ## the recoding is LandR's own wetlandToLCC() (PredictiveEcology/LandR#228), not a mock
  skip_if_not("wetlandToLCC" %in% getNamespaceExports("LandR"), "LandR without wetlandToLCC() (LandR#228)")
  seen <- new.env()
  lcc <- runLCC("SCANFI", seen)
  expect_true(isTRUE(seen$scanfi)); expect_null(seen$ntems); expect_true(isTRUE(seen$cwim))
  ## wet rows 1, 3, 5: conifer -> 81, shrub -> 80, and a non-flammable 0 stays 0
  expect_equal(rowsOf(lcc), c(81, 210, 80, 50, 0, 0))
})

test_that("the source defaults to SCANFI and can be switched with options(fireSense.lccSource)", {
  expect_identical(formals(makeFireSenseLCC)$lccSource, quote(getOption("fireSense.lccSource", "SCANFI")))
  seen <- new.env()
  withr::local_options(fireSense.lccSource = "NTEMS")
  testthat::local_mocked_bindings(
    prepInputs_NTEMS_LCC_FAO = function(...) { seen$ntems <- TRUE; srcGrid(srcVals()) },
    .scanfiLCC = function(...) { seen$scanfi <- TRUE; srcGrid(srcVals()) },
    .cwimWetland = function(to, ...) { seen$cwim <- TRUE; terra::rast(to, vals = 0L) }
  )
  makeFireSenseLCC(neededYear = 2020, to = target, destinationPath = withr::local_tempdir())
  expect_true(isTRUE(seen$ntems)); expect_null(seen$scanfi)
  expect_error(makeFireSenseLCC(neededYear = 2020, to = target, destinationPath = tempdir(), lccSource = "EOSD"),
               "SCANFI.*NTEMS")
})
