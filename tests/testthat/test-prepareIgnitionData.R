## The ignition-fit data chain: climate, fuel and land-cover rasters stacked per year, with the
## ignitions counted per cell (stackAndExtract, mergePreparedCovs), then rescaled
## (rescaleCovariates). Small made-up rasters: 4 x 4 cells of 1 km, two fire years.

ignRas <- function(vals, names = NULL) {
  r <- terra::rast(nrows = 4, ncols = 4, extent = c(0, 4000, 0, 4000), crs = "EPSG:3978",
                   nlyrs = NCOL(vals), vals = vals)
  if (!is.null(names)) names(r) <- names
  r
}
ignClimate <- function() {
  list(MDC = ignRas(cbind(100 + 1:16, 200 + 1:16), c("year2001", "year2002")))
}
ignFires <- function() {
  ## two ignitions in the first cell in 2001, one in the last cell in 2002
  v <- terra::vect(cbind(c(500, 500, 3500), c(3500, 3500, 500)), crs = "EPSG:3978")
  v$YEAR <- c(2001, 2001, 2002)
  v
}

test_that("stackAndExtract gives each year's covariates per cell with its ignition count", {
  fuel <- ignRas(rep(c(1, 0), 8), "class2")
  lcc <- ignRas(rep(c(0, 1), 8), "nf1")

  out <- stackAndExtract(years = c("year2001", "year2002"), fuel = fuel, LCC = lcc,
                         climate = ignClimate(), fires = ignFires())

  expect_setequal(names(out), c("MDC", "nf1", "class2", "cell", "ignitions", "year"))
  expect_identical(nrow(out), 32L)
  expect_identical(out[year == "2001" & cell == 1, ignitions], 2L)
  expect_identical(out[year == "2002" & cell == 16, ignitions], 1L)
  expect_identical(sum(out$ignitions), 3L)
  expect_identical(out[year == "2002" & cell == 3, MDC], 203)

  noFires <- stackAndExtract("year2001", fuel, lcc, ignClimate())
  expect_true(all(noFires$ignitions == 0))
})

test_that("mergePreparedCovs drops cells with no fuel or land cover and adds lightning", {
  withr::local_options(reproducible.useCache = FALSE)
  fuelCov <- c(ignRas(c(1, 1, 0, rep(1, 13)), "class2"), ignRas(c(0, 0, 0, rep(0, 13)), "nf1"))
  lightning <- list(lightningDensity = ignRas(1:16 / 10))

  out <- mergePreparedCovs(years = list(c("year2001", "year2002")), fuelCovsCoarse = list(fuelCov),
                           ignitionFirePoints = ignFires(), nonForestedLCCGroups = list(nf1 = 1),
                           ignitionClimateCoarse = ignClimate(), lightningMap = lightning,
                           digest = NULL, useCache = FALSE)

  expect_false(3 %in% out$pixelID)                         # no cover at all in cell 3
  expect_identical(names(out)[1:3], c("pixelID", "ignitions", "MDC"))
  expect_true(is.numeric(out$year))
  expect_equal(out[pixelID == 16 & year == 2002, lightningDensity], 1.6)
})

test_that("prepare_ignitionClimate aggregates each climate variable to the coarse grid", {
  out <- prepare_ignitionClimate(ignClimate(), fact = 2, useCache = FALSE)
  expect_equal(dim(out$MDC), c(2, 2, 2))
  expect_equal(unname(terra::values(out$MDC)[1, "year2001"]), mean(c(101, 102, 105, 106)))
})

test_that("prepare_FuelCovsCoarse rasterizes the fuel covariates and aggregates them", {
  local_mocked_bindings(fireSenseCovariatesCreate = function(...)
    data.table::data.table(pixelID = 1:16, class2 = rep(c(1, 0), 8)))
  out <- prepare_FuelCovsCoarse(rasTemplate = ignRas(0), fact = 2)
  expect_identical(names(out), "class2")
  expect_equal(unname(terra::values(out)[, 1]), rep(0.5, 4))
})

test_that("rescaleCovariates brings variables above 10 to [0, 10] by their magnitude", {
  covs <- data.frame(ignitions = c(0, 2, 1), MDC = c(150, 320, 90), class2 = c(0.2, 0.5, 1),
                     year = 2001:2003)
  out <- suppressMessages(rescaleCovariates(ignitions ~ MDC + class2, covs, rescaleVars = TRUE,
                                            modelAlgorithm = "glm"))
  expect_equal(out$covariates$MDC, c(1.5, 3.2, 0.9))
  expect_equal(out$covariates$class2, covs$class2)
  expect_equal(unname(out$ignitionRescalers), 100)
  expect_identical(out$xvar, "year")

  kept <- rescaleCovariates(ignitions ~ MDC, covs, rescaleVars = FALSE, modelAlgorithm = "glm")
  expect_null(kept$ignitionRescalers)
  expect_equal(kept$covariates$MDC, covs$MDC)

  expect_error(rescaleCovariates(~ MDC, covs, rescaleVars = TRUE, modelAlgorithm = "glm"),
               "the LHS is missing")
})

test_that("prepareCovariatesOuter centres and scales for xgboost, keeping the counts", {
  covs <- data.frame(ignitions = c(0, 2, 1), MDC = c(150, 320, 90), year = 2001:2003)
  out <- prepareCovariatesOuter(covs, algorithm = "xgb", rescaleVars = TRUE, useCache = FALSE)
  expect_equal(mean(out$covariates$MDC), 0)
  expect_identical(out$covariates$ignitions, c(0, 2, 1))
  expect_false(is.null(attr(out$covariates, "scaleData")))
  expect_null(out$digestOfData)
})

test_that("igOrEscNames builds object names in the asked case", {
  expect_identical(igOrEscNames("ignition", post = "Covariates"), "fireSense_ignitionCovariates")
  expect_identical(igOrEscNames("escape", post = "Fit", case = "camel"), "fireSense_EscapeFit")
  expect_identical(igOrEscNames("escape", pre = "", post = "", case = "title"), "Escape")
})

test_that("climateRasterToDataTable makes one integer column per variable, by pixel and year", {
  clim <- list(MDC = ignRas(cbind(1:16 + 0.4, 1:16), c("year2001", "year2002")),
               CMD = ignRas(cbind(rep(5, 16), rep(6, 16)), c("year2001", "year2002")))
  out <- climateRasterToDataTable(clim)
  expect_identical(nrow(out), 32L)
  expect_true(is.integer(out$MDC))
  expect_identical(out[pixelID == 2 & year == "year2001", MDC], 2L)
  expect_identical(nrow(climateRasterToDataTable(clim["MDC"], Index = 1:3)), 6L)
})

test_that("predictIgnition scales the model's predicted rate by both factors", {
  d <- data.frame(ignitions = c(0, 1, 2, 4), MDC = 1:4)
  m <- stats::glm(ignitions ~ MDC, family = stats::poisson(), data = d)
  expect_equal(predictIgnition(m, d, rescaleFactor = 2, lambdaRescaleFactor = 3),
               6 * stats::fitted(m))
})

test_that("prepare_LightningData reads each lightning product and aggregates it", {
  withr::local_options(reproducible.useCache = FALSE)
  rtm <- ignRas(0)
  local_mocked_bindings(getRemoteMetadata = function(...) list(remoteHash = "h"),
                        googledriveIDtoHumanURL = function(x) x, .package = "reproducible")
  local_mocked_bindings(prepInputs = function(url, fun, destinationPath, useCache) ignRas(1:16))

  out <- prepare_LightningData(rtm, igAggFactor = 2, dPath = withr::local_tempdir())

  expect_named(out, c("lightningDays", "lightningDensity", "positiveCG", "positiveCGdensity"))
  expect_equal(dim(out$lightningDays), c(2, 2, 1))
})

test_that("readLightningData rasterizes the Lat/Long/density triples at 10 km", {
  ## the file repeats Lat, Long, LightningDensity across columns
  csv <- withr::local_tempfile(fileext = ".csv")
  ## a 0.5-degree grid of points, 2 x 2 degrees, split across two triples of columns
  grid <- expand.grid(lat = seq(49.05, 51.05, by = 0.5), long = seq(-101.05, -99.05, by = 0.5))
  half <- seq_len(nrow(grid) / 2 + 0.5)
  a <- cbind(grid[half, ], ld = 1.5)
  b <- cbind(grid[-half, ], ld = 4.5)
  data.table::fwrite(cbind(a[seq_len(nrow(b)), ], b), csv, col.names = FALSE)

  out <- readLightningData(csv)

  expect_s4_class(out, "SpatRaster")
  expect_identical(terra::res(out), c(1e4, 1e4))
  expect_equal(range(terra::values(out), na.rm = TRUE), c(1.5, 4.5))

  ## with a target grid: back on that grid, in its projection
  to <- terra::rast(terra::ext(terra::project(out, "EPSG:3978")), res = 5000, crs = "EPSG:3978",
                    vals = 1)   # it is also the mask
  ## sf warns about point attributes when reproducible crops and masks them
  onTo <- suppressWarnings(readLightningData(csv, to = to))
  expect_true(terra::compareGeom(onTo, to, stopOnError = FALSE))
  expect_true(any(!is.na(terra::values(onTo))))
})

test_that("prepare_FuelCovsCoarseByYear: a pixel burned after the data year is young, and cleared, from the next year", {
  ## covariates as fireSenseCovariatesCreate(youngAge = FALSE) gives them: nothing zeroed
  local_mocked_bindings(fireSenseCovariatesCreate = function(..., youngAge) {
    expect_false(youngAge)
    data.table::data.table(pixelID = 1:16, spruce = log(3000), nfLCC_40 = 1, treedWetland = 1)
  })
  tsd <- ignRas(100)                           # every pixel old at the data year, 2000
  ## coarse cell 1 is pixels 1, 2, 5, 6; they burn in 2002
  fires <- list("2002" = c(1L, 2L, 5L, 6L))
  out <- prepare_FuelCovsCoarseByYear(years = c("year2001", "year2003"), dataYear = 2000,
                                      nonForest_timeSinceDisturbance = tsd,
                                      firePixelsByYear = fires, cutoffForYoungAge = 15,
                                      rasTemplate = ignRas(0), fact = 2)
  expect_named(out, c("year2001", "year2003"))
  v1 <- terra::values(out$year2001)
  v3 <- terra::values(out$year2003)
  expect_equal(unname(v1[, "youngAge"]), rep(0, 4))
  expect_equal(unname(v1[1, "spruce"]), log(3000))
  expect_equal(unname(v3[, "youngAge"]), c(1, 0, 0, 0))
  expect_equal(unname(v3[1, "spruce"]), logMinB(0))
  expect_equal(unname(v3[2, "spruce"]), log(3000))
  expect_equal(unname(v3[1, c("nfLCC_40", "treedWetland")]), c(0, 0))

  ## the per-year rasters feed stackAndExtract: ignition covariates at 2003 see the burn
  clim <- list(MDC = ignRas(cbind(1:16, 11:26), c("year2001", "year2003")))
  ex <- stackAndExtract(c("year2001", "year2003"),
                        fuel = lapply(out, function(r) r[[c("youngAge", "spruce")]]),
                        LCC = lapply(out, function(r) r[[c("nfLCC_40")]]), climate = lapply(clim, function(r) terra::aggregate(r, 2, mean)))
  expect_equal(ex[year == "2003" & cell == 1, youngAge], 1)
  expect_equal(ex[year == "2001" & cell == 1, youngAge], 0)
})

test_that("per-year ignition covariates run through fireSenseCovariatesCreate and mergePreparedCovs", {
  withr::local_package("terra")
  withr::local_package("data.table")
  ## 4 x 4 pixels of 1 km; pixels 1, 2, 5, 6 (coarse cell 1) burn in 2002. No cohorts: the only
  ## cover is non-forest, nfLCC_40 everywhere, so a young pixel is cleared to nothing
  pgm <- ignRas(0L)
  noCohorts <- data.table(pixelGroup = integer(), speciesCode = character(),
                          age = integer(), B = integer())
  landcoverDT <- data.table(pixelID = 1:16, nfLCC_40 = 1)
  fuel <- suppressWarnings(prepare_FuelCovsCoarseByYear(
    cohortData = noCohorts, pixelGroupMap = pgm, flammableRTM = ignRas(1),
    sppEquiv = data.table(LandR = character(), FuelClass = character()),
    landcoverDT = landcoverDT, fuelClassCol = "FuelClass", sppEquivCol = "LandR",
    missingLCCgroup = "nfLCC_40", nonForestedLCCGroups = c(nfLCC_40 = 40),
    nonForestCanBeYoungAge = TRUE, studyAreaName = "test", useCache = FALSE,
    nonForest_timeSinceDisturbance = ignRas(100), cutoffForYoungAge = 15,
    years = c("year2001", "year2002"), dataYear = 2000,
    firePixelsByYear = list("2001" = c(1L, 2L, 5L, 6L)),
    rasTemplate = ignRas(0), fact = 2))
  expect_named(fuel, c("year2001", "year2002"))
  expect_equal(unname(terra::values(fuel$year2002)[, "youngAge"]), c(1, 0, 0, 0))
  expect_equal(unname(terra::values(fuel$year2002)[, "nfLCC_40"]), c(0, 1, 1, 1))
  expect_equal(unname(terra::values(fuel$year2001)[, "youngAge"]), rep(0, 4))

  out <- mergePreparedCovs(years = list(c("year2001", "year2002")), fuelCovsCoarse = list(fuel),
                           ignitionFirePoints = ignFires(), nonForestedLCCGroups = list(nfLCC_40 = 1),
                           ignitionClimateCoarse = prepare_ignitionClimate(
                             list(MDC = ignRas(cbind(1:16, 11:26), c("year2001", "year2002"))), 2,
                             useCache = FALSE),
                           lightningMap = list(lightningDays = ignRas(1:16)),
                           digest = NULL, useCache = FALSE)
  ## coarse cell 1: all four pixels burned, so in 2002 (not 2001) it is young, without land cover
  expect_equal(out[pixelID == 1 & year == 2002, youngAge], 1)
  expect_equal(out[pixelID == 1 & year == 2002, nfLCC_40], 0)
  expect_equal(out[pixelID == 1 & year == 2001, youngAge], 0)
  expect_equal(out[pixelID == 1 & year == 2001, nfLCC_40], 1)
})
