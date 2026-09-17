## A treed-wetland site covariate for the spread and ignition models.
##
## Land-cover class 81 (treed wetland) is a forested class, so its pixels enter the fire models only through
## their fuel biomass and look exactly like upland forest. With `rstLCC` given, fireSenseCovariatesCreate()
## adds `treedWetland` (1 on class 81). It is a SITE attribute, not a fuel state, so it is added after the
## youngAge exclusivity (a burned bog is still wet) and it is not "cover" for the ignition data's
## all-cover-is-zero filter. Without `rstLCC` the output is unchanged.

skip_if_not_installed("terra")

fuelRas <- function() {
  r <- terra::rast(nrows = 2, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 2, crs = "EPSG:3978", nlyrs = 2)
  names(r) <- c("Pice_mar", "youngAge")
  terra::values(r) <- cbind(c(500, 800, 0, 900, 0, 700), c(0, 0, 0, 1, 0, 0))
  r
}
lccRas <- function() terra::rast(fuelRas()[[1]], vals = c(81, 210, 50, 81, 81, 20))
landcover <- function() data.table::data.table(pixelID = 1:6, nfLCC_50 = c(0, 0, 1, 0, 0, 0))

covs <- function(rstLCC = NULL) {
  testthat::local_mocked_bindings(cohortsToFuelClasses = function(...) fuelRas())
  args <- list(cohortData = NULL, pixelGroupMap = NULL, flammableRTM = fuelRas()[[1]], sppEquiv = NULL,
               landcoverDT = landcover(), fuelClassCol = "FuelClass", sppEquivCol = "LandR",
               missingLCCgroup = "nfLCC_50", nonForestedLCCGroups = list(nfLCC_50 = 50),
               nonForest_timeSinceDisturbance = NULL, cutoffForYoungAge = 15, nonForestCanBeYoungAge = FALSE,
               studyAreaName = "test", useCache = FALSE)
  if (!is.null(rstLCC)) args$rstLCC <- rstLCC
  do.call(fireSenseCovariatesCreate, args)
}

test_that("without rstLCC the covariates are unchanged: no treedWetland column", {
  expect_false(treedWetlandTxt %in% names(covs()))
})

test_that("with rstLCC, treedWetland is 1 exactly on class 81, and youngAge does not clear it", {
  out <- covs(lccRas())[order(pixelID)]
  expect_identical(out[[treedWetlandTxt]], c(1, 0, 0, 1, 1, 0))
  expect_identical(out$youngAge[4], 1)             # pixel 4: young, and still a treed wetland
  expect_identical(out$Pice_mar[4], logMinB(0))    # its fuel was cleared by youngAge as before
})

test_that("treedWetland is a site layer, not cover, in the ignition data's zero-cover filter", {
  src <- paste(deparse(mergePreparedCovs), collapse = "\n")
  expect_match(src, "treedWetlandTxt", fixed = TRUE)
})
