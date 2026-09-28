## plotHistoricFires(), plotCumulativeBurns() and plotBurnSummary() read per-replicate outputs
## (`simFiles`) and write figures under <outputDir>/<studyAreaName>/figures/. The replicates here
## are made up: three reps of ten years of fires, and a 5 x 5 burn map per rep.

makeRepOutputs <- function(dir, nReps = 3, years = 2011:2020, lastYear = max(years)) {
  set.seed(1)
  files <- unlist(lapply(seq_len(nReps), function(rep) {
    repDir <- file.path(dir, paste0("rep", rep))
    dir.create(repDir, recursive = TRUE)
    nFires <- sample(2:6, length(years), replace = TRUE)
    burn <- data.frame(year = rep(years, nFires), N = 1L,
                       areaBurnedHa = round(stats::rexp(sum(nFires), 1 / 500)) + 1)
    csv <- file.path(repDir, "fireSense_burnSummary.csv")
    data.table::fwrite(burn, csv)
    tif <- file.path(repDir, paste0("burnMap_year", lastYear, ".tif"))
    r <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 500, ymin = 0, ymax = 500,
                     crs = "EPSG:3978", vals = stats::rbinom(25, 1, 0.5))
    terra::writeRaster(r, tif)
    c(csv, tif)
  }))
  files
}

test_that("plotBurnSummary writes one figure named by study area and scenario", {
  skip_if_not_installed("cowplot")
  withr::local_options(mc.cores = 1L)
  dir <- withr::local_tempdir()
  simFiles <- makeRepOutputs(dir)

  f <- suppressMessages(plotBurnSummary("CanESM5_ssp370", "SA", dir, Nreps = 3,
                                        years = c(2011, 2020), pixelSize = 250,
                                        simFiles = simFiles))

  expect_identical(f, file.path(dir, "SA", "figures", "burnSummary_SA_CanESM5_ssp370.png"))
  expect_true(file.exists(f))
})

test_that("plotBurnSummary stops with the reading errors instead of failing later", {
  skip_if_not_installed("cowplot")
  withr::local_options(mc.cores = 1L)
  dir <- withr::local_tempdir()
  simFiles <- makeRepOutputs(dir)
  local_mocked_bindings(mclapply = function(X, FUN, ...) {
    list(structure("Error : cannot open file\n", class = "try-error"))
  }, .package = "parallel")

  expect_error(plotBurnSummary(NA, "SA", dir, 3, c(2011, 2020), 250, simFiles),
               "failed to read 1 of 1 burn summaries:\nError : cannot open file")
})

test_that("plotCumulativeBurns writes the cumulative burn map, labelling no scenario as NRV", {
  skip_if_not_installed("raster")
  skip_if_not_installed("rasterVis")
  withr::local_options(mc.cores = 1L)
  dir <- withr::local_tempdir()
  simFiles <- makeRepOutputs(dir)
  rtm <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 500, ymin = 0, ymax = 500,
                     crs = "EPSG:3978", vals = 1)

  expect_message(
    f <- plotCumulativeBurns(NA, "SA", dir, Nreps = 3, years = c(2011, 2020),
                             rasterToMatch = rtm, simFiles = simFiles),
    "declaring it to be NRV")

  expect_identical(f, file.path(dir, "SA", "figures", "cumulBurnMap_SA_NRV.png"))
  expect_true(file.exists(f))
})

test_that("plotHistoricFires writes ignition, escape and area-burned figures", {
  withr::local_options(mc.cores = 1L)
  dir <- withr::local_tempdir()
  simFiles <- makeRepOutputs(dir, years = 2001:2010)
  pts <- sf::st_as_sf(data.frame(x = 1:20, y = 1:20, YEAR = rep(2001:2010, 2)),
                      coords = c("x", "y"), crs = 3978)
  polys <- terra::vect(sf::st_buffer(pts, 1))
  polys$SIZE_HA <- polys$POLY_HA <- seq(100, 2000, by = 100)

  f <- suppressMessages(plotHistoricFires("CanESM5_ssp370", "SA", dir, pixelSize = 250,
                                          firePolys = polys, ignitionPoints = pts,
                                          simFiles = simFiles))

  expect_named(f, c("ignition", "escape", "spread"))
  expect_identical(basename(f), paste0(c("simulated_Ignitions", "simulated_Escapes",
                                         "simulated_burnArea"), "_SA_CanESM5_ssp370.png"))
  expect_true(all(file.exists(f)))
})

test_that(".attributes reads SpatVector, sf, Spatial and lists of them", {
  df <- data.frame(YEAR = 2001:2002, x = 1:2, y = 1:2)
  pts <- sf::st_as_sf(df, coords = c("x", "y"), crs = 3978)
  expect_identical(.attributes(pts), data.frame(YEAR = 2001:2002))
  expect_identical(.attributes(terra::vect(pts)), data.frame(YEAR = 2001:2002))
  ## a Spatial*DataFrame keeps its attributes in @data
  fakeSpatial <- methods::setClass("fakeSpatial", representation(data = "data.frame"),
                                   where = environment())
  sp <- fakeSpatial(data = data.frame(YEAR = 2001:2002))
  expect_identical(.attributes(sp), data.frame(YEAR = 2001:2002))
  expect_identical(.attributes(list(pts, terra::vect(pts)))$YEAR, c(2001:2002, 2001:2002))
  expect_identical(.attributes(data.frame(YEAR = 1L)), data.frame(YEAR = 1L))
})
