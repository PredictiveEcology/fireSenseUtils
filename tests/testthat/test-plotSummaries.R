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

## the observed rows of the area-burned figure's data, with ggsave() captured
burnedData <- function(polys, pts, simFiles, dir, ...) {
  plots <- list()
  local_mocked_bindings(ggsave = function(plot, filename, ...) plots[[basename(filename)]] <<- plot,
                        .package = "ggplot2")
  suppressMessages(plotHistoricFires("CanESM5_ssp370", "SA", dir, pixelSize = 250, firePolys = polys,
                                     ignitionPoints = pts, simFiles = simFiles, ...))
  d <- plots[["simulated_burnArea_SA_CanESM5_ssp370.png"]]$data
  d[d$stat == "observed", ]
}

test_that("plotHistoricFires needs only SIZE_HA, and gives the same burns as when POLY_HA is there too", {
  withr::local_options(mc.cores = 1L)
  dir <- withr::local_tempdir()
  simFiles <- makeRepOutputs(dir, years = 2001:2010)
  pts <- sf::st_as_sf(data.frame(x = 1:20, y = 1:20, YEAR = rep(2001:2010, 2)),
                      coords = c("x", "y"), crs = 3978)
  both <- terra::vect(sf::st_buffer(pts, 1))
  both$SIZE_HA <- both$POLY_HA <- seq(100, 2000, by = 100)
  only <- both
  only$POLY_HA <- NULL

  expect_no_error(obsOnly <- burnedData(only, pts, simFiles, dir))
  expect_equal(obsOnly, burnedData(both, pts, simFiles, dir))
  expect_gt(nrow(obsOnly), 0L)
})

test_that("plotHistoricFires filters and sums the one column named by sizeCol", {
  withr::local_options(mc.cores = 1L)
  dir <- withr::local_tempdir()
  simFiles <- makeRepOutputs(dir, years = 2001:2010)
  pts <- sf::st_as_sf(data.frame(x = 1:20, y = 1:20, YEAR = rep(2001:2010, 2)),
                      coords = c("x", "y"), crs = 3978)
  polys <- terra::vect(sf::st_buffer(pts, 1))
  polys$ADJ_HA <- seq(100, 2000, by = 100)

  obs <- burnedData(polys, pts, simFiles, dir, sizeCol = "ADJ_HA")

  expect_equal(sum(obs$sumBurn), sum(polys$ADJ_HA))
})

## every figure of plotHistoricFires(), with ggsave() captured
historicFigs <- function(simFiles, dir) {
  plots <- list()
  local_mocked_bindings(ggsave = function(plot, filename, ...) plots[[basename(filename)]] <<- plot,
                        .package = "ggplot2")
  pts <- sf::st_as_sf(data.frame(x = 1:20, y = 1:20, YEAR = rep(2001:2010, 2)),
                      coords = c("x", "y"), crs = 3978)
  polys <- terra::vect(sf::st_buffer(pts, 1))
  polys$SIZE_HA <- seq(100, 2000, by = 100)
  suppressMessages(plotHistoricFires("CanESM5_ssp370", "SA", dir, pixelSize = 250, firePolys = polys,
                                     ignitionPoints = pts, simFiles = simFiles))
  plots
}

test_that("plotHistoricFires counts the same fires on both sides of each figure", {
  withr::local_options(mc.cores = 1L)
  dir <- withr::local_tempdir()
  simFiles <- makeRepOutputs(dir, nReps = 2, years = 2001:2004)
  ## observed fires of mixed sizes, some below pixelSize (250), two per year
  pts <- sf::st_as_sf(data.frame(x = 1:8, y = 1:8, YEAR = rep(2001:2004, each = 2)),
                      coords = c("x", "y"), crs = 3978)
  polys <- terra::vect(sf::st_buffer(pts, 1))
  polys$SIZE_HA <- c(10, 1000, 200, 300, 5, 6, 400, 20)
  plots <- list()
  local_mocked_bindings(ggsave = function(plot, filename, ...) plots[[basename(filename)]] <<- plot,
                        .package = "ggplot2")
  suppressMessages(plotHistoricFires("CanESM5_ssp370", "SA", dir, pixelSize = 250, firePolys = polys,
                                     ignitionPoints = pts, simFiles = simFiles))
  dat <- function(nm) plots[[paste0("simulated_", nm, "_SA_CanESM5_ssp370.png")]]$data
  obs <- function(d) d[d$stat == "observed", ]
  proj <- function(d) d[d$stat == "projected", ]
  byYear <- function(d, col) tapply(d[[col]], d$year, sum)
  sim <- data.table::rbindlist(lapply(simFiles[grepl("csv$", simFiles)], data.table::fread))

  ## area burned: all fires, observed and simulated
  expect_equal(as.numeric(byYear(obs(dat("burnArea")), "sumBurn")),
               as.numeric(tapply(polys$SIZE_HA, polys$YEAR, sum)))
  expect_equal(sum(proj(dat("burnArea"))$sumBurn), sum(sim$areaBurnedHa))

  ## escapes: fires larger than pixelSize, both sides
  expect_equal(as.numeric(byYear(obs(dat("Escapes")), "nFires")),
               as.numeric(tapply(polys$SIZE_HA > 250, polys$YEAR, sum)))
  expect_equal(sum(proj(dat("Escapes"))$nFires), sum(sim$areaBurnedHa > 250))

  ## ignitions: every fire, both sides
  expect_equal(as.numeric(byYear(obs(dat("Ignitions")), "N")), rep(2, 4))
  expect_equal(sum(proj(dat("Ignitions"))$N), nrow(sim))
})

test_that("plotHistoricFires draws every replicate, and the median across replicates", {
  withr::local_options(mc.cores = 1L)
  dir <- withr::local_tempdir()
  simFiles <- makeRepOutputs(dir, nReps = 3, years = 2001:2010)

  plots <- historicFigs(simFiles, dir)
  d <- plots[["simulated_burnArea_SA_CanESM5_ssp370.png"]]$data
  proj <- d[d$stat == "projected", ]
  expect_setequal(proj$rep, 1:3)

  med <- d[d$stat == "projected median", ]
  expect_equal(med$sumBurn[order(med$year)],
               as.numeric(tapply(proj$sumBurn, proj$year, stats::median)))
  expect_false(any(grepl("rep", plots[["simulated_burnArea_SA_CanESM5_ssp370.png"]]$labels$title)))
  for (f in plots) expect_s3_class(f, "ggplot")
})

test_that("plotHistoricFires still works with a single replicate", {
  withr::local_options(mc.cores = 1L)
  dir <- withr::local_tempdir()
  simFiles <- makeRepOutputs(dir, nReps = 1, years = 2001:2010)
  d <- historicFigs(simFiles, dir)[["simulated_burnArea_SA_CanESM5_ssp370.png"]]$data
  expect_equal(d$sumBurn[d$stat == "projected median"][order(d$year[d$stat == "projected median"])],
               d$sumBurn[d$stat == "projected"][order(d$year[d$stat == "projected"])])
})

test_that("plotBurnSummary plots the median fire size, and titles the x axis Simulation year", {
  skip_if_not_installed("cowplot")
  withr::local_options(mc.cores = 1L)
  dir <- withr::local_tempdir()
  simFiles <- makeRepOutputs(dir)
  panels <- NULL
  local_mocked_bindings(
    plot_grid = function(...) {
      args <- list(...)
      if (is.null(panels) && length(args) >= 3L) panels <<- args[1:3]
      ggplot2::ggplot()
    }, .package = "cowplot")
  local_mocked_bindings(ggsave = function(...) NULL, .package = "ggplot2")
  suppressMessages(plotBurnSummary("CanESM5_ssp370", "SA", dir, Nreps = 3, years = c(2011, 2020),
                                   pixelSize = 250, simFiles = simFiles))

  ## the median per rep and year, of the fires larger than a pixel (250 m pixel = 6.25 ha)
  fires <- do.call(rbind, lapply(1:3, function(r)
    cbind(as.data.frame(data.table::fread(simFiles[grepl(paste0("rep", r, "/.*csv"), simFiles)])), rep = r)))
  fires <- fires[fires$areaBurnedHa > 250^2 / 1e4, ]
  expected <- tapply(fires$areaBurnedHa, list(fires$year, fires$rep), stats::median)
  got <- panels[[3]]$data
  expect_equal(sort(got$val), sort(as.numeric(expected[!is.na(expected)])))
  expect_match(panels[[3]]$labels$y, "median fire size")
  expect_identical(panels[[3]]$labels$x, "Simulation year")
  expect_s3_class(panels[[1]]$theme$axis.title.x, "element_blank")
  expect_s3_class(panels[[2]]$theme$axis.title.x, "element_blank")
  expect_null(panels[[3]]$theme$axis.title.x)
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
