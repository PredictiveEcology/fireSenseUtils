## fuelCovariates = "domSecOther": fireSenseCovariatesCreate() collapses per-species fuel-class
## columns to dom_agb_<class>, sec_agb_<class> (the two classes with the most total treed AGB),
## other_agb (the rest, pooled) and treedWetland_agb (all tree AGB on treed-wetland pixels,
## removed from the other three there). Prediction must be able to force the same dom/sec classes
## the fit chose, even when a different class dominates the prediction area.

skip_if_not_installed("terra")

## 2x3 landscape, 3 tree fuel classes.
##            px1  px2  px3  px4  px5  px6
## Pice_mar:  500  800    0  900    0  700   -> total 2900 (dominant)
## Pinu_ban:  200  100    0  300    0  200   -> total  800 (secondary)
## Popu_tre:   50   50    0  100    0   50   -> total  250 (pooled into other_agb)
## youngAge:    0    0    0    1    0    0   -> px4 young
## lcc:        81  210   50   81   81   20   -> px1, px4, px5 treed wetland; px3 non-forest (nfLCC_50)
fuelRas <- function() {
  r <- terra::rast(nrows = 2, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 2, crs = "EPSG:3978", nlyrs = 4)
  names(r) <- c("Pice_mar", "Pinu_ban", "Popu_tre", "youngAge")
  terra::values(r) <- cbind(
    c(500, 800, 0, 900, 0, 700),
    c(200, 100, 0, 300, 0, 200),
    c(50, 50, 0, 100, 0, 50),
    c(0, 0, 0, 1, 0, 0)
  )
  r
}
lccRas <- function() terra::rast(fuelRas()[[1]], vals = c(81, 210, 50, 81, 81, 20))
landcover <- function() data.table::data.table(pixelID = 1:6, nfLCC_50 = c(0, 0, 1, 0, 0, 0))

covs <- function(rstLCC = NULL, domClass = NULL, secClass = NULL, fuelCovariates = "domSecOther") {
  testthat::local_mocked_bindings(cohortsToFuelClasses = function(...) fuelRas())
  args <- list(cohortData = NULL, pixelGroupMap = NULL, flammableRTM = fuelRas()[[1]], sppEquiv = NULL,
               landcoverDT = landcover(), fuelClassCol = "FuelClass", sppEquivCol = "LandR",
               missingLCCgroup = "nfLCC_50", nonForestedLCCGroups = list(nfLCC_50 = 50),
               nonForest_timeSinceDisturbance = NULL, cutoffForYoungAge = 15, nonForestCanBeYoungAge = FALSE,
               studyAreaName = "test", useCache = FALSE, fuelCovariates = fuelCovariates,
               domClass = domClass, secClass = secClass)
  if (!is.null(rstLCC)) args$rstLCC <- rstLCC
  do.call(fireSenseCovariatesCreate, args)
}

test_that("dom/sec are chosen by total AGB and named after the class; other_agb is the rest", {
  out <- covs()[order(pixelID)]
  expect_true(all(c("dom_agb_Pice_mar", "sec_agb_Pinu_ban", "other_agb") %in% names(out)))
  expect_false(any(c("Pice_mar", "Pinu_ban", "Popu_tre") %in% names(out)))

  ## pixel 2: not wetland, not young -- straightforward collapse
  expect_equal(out[pixelID == 2, dom_agb_Pice_mar], logMinB(800))
  expect_equal(out[pixelID == 2, sec_agb_Pinu_ban], logMinB(100))
  expect_equal(out[pixelID == 2, other_agb], logMinB(50)) # Popu_tre alone

  expect_equal(out[pixelID == 6, dom_agb_Pice_mar], logMinB(700))
  expect_equal(out[pixelID == 6, sec_agb_Pinu_ban], logMinB(200))
  expect_equal(out[pixelID == 6, other_agb], logMinB(50))

  ## pixel 3: non-forest, no tree AGB at all
  expect_equal(out[pixelID == 3, dom_agb_Pice_mar], logMinB(0))
  expect_equal(out[pixelID == 3, sec_agb_Pinu_ban], logMinB(0))
  expect_equal(out[pixelID == 3, other_agb], logMinB(0))
})

test_that("a treed-wetland pixel's AGB is only in treedWetland_agb, removed from dom/sec/other", {
  out <- covs(rstLCC = lccRas())[order(pixelID)]
  expect_true("treedWetland_agb" %in% names(out))

  ## pixel 1: wetland, not young -- all AGB (500 + 200 + 50 = 750) moves to treedWetland_agb
  expect_equal(out[pixelID == 1, treedWetland_agb], logMinB(750))
  expect_equal(out[pixelID == 1, dom_agb_Pice_mar], logMinB(0))
  expect_equal(out[pixelID == 1, sec_agb_Pinu_ban], logMinB(0))
  expect_equal(out[pixelID == 1, other_agb], logMinB(0))

  ## pixel 2: not wetland -- treedWetland_agb is 0, dom/sec/other unaffected
  expect_equal(out[pixelID == 2, treedWetland_agb], logMinB(0))
  expect_equal(out[pixelID == 2, dom_agb_Pice_mar], logMinB(800))
})

test_that("a young pixel has all four AGB terms (and nfLCC) at 0, and youngAge = 1", {
  out <- covs(rstLCC = lccRas())[order(pixelID)]
  ## pixel 4: young AND wetland -- youngAge wins: nothing survives, including treedWetland_agb
  expect_equal(out[pixelID == 4, youngAge], 1)
  expect_equal(out[pixelID == 4, dom_agb_Pice_mar], logMinB(0))
  expect_equal(out[pixelID == 4, sec_agb_Pinu_ban], logMinB(0))
  expect_equal(out[pixelID == 4, other_agb], logMinB(0))
  expect_equal(out[pixelID == 4, treedWetland_agb], logMinB(0))
})

test_that("prediction can force the fit's dom/sec classes even when a different class dominates here", {
  ## Popu_tre is the smallest class in this fixture (total 250), but a prediction is told it was
  ## the fit's dominant class: the columns built must follow that, not this area's own ranking.
  out <- covs(domClass = "Popu_tre", secClass = "Pice_mar")[order(pixelID)]
  expect_true(all(c("dom_agb_Popu_tre", "sec_agb_Pice_mar", "other_agb") %in% names(out)))
  expect_equal(out[pixelID == 2, dom_agb_Popu_tre], logMinB(50))
  expect_equal(out[pixelID == 2, sec_agb_Pice_mar], logMinB(800))
  expect_equal(out[pixelID == 2, other_agb], logMinB(100)) # Pinu_ban alone
})

test_that("a forced class absent from this area's fuel classes gets a zero column, not an error", {
  out <- covs(domClass = "Betu_pap", secClass = "Pice_mar")[order(pixelID)]
  expect_true("dom_agb_Betu_pap" %in% names(out))
  expect_equal(out[pixelID == 2, dom_agb_Betu_pap], logMinB(0))
  expect_equal(out[pixelID == 2, sec_agb_Pice_mar], logMinB(800))
  ## Pinu_ban and Popu_tre both fall into other_agb
  expect_equal(out[pixelID == 2, other_agb], logMinB(150))
})

test_that("fuelCovariates = \"species\" (the default) is unchanged: per-species columns, no dom/sec/other", {
  out <- covs(fuelCovariates = "species")
  expect_true(all(c("Pice_mar", "Pinu_ban", "Popu_tre") %in% names(out)))
  expect_false(any(c("dom_agb_Pice_mar", "sec_agb_Pinu_ban", "other_agb") %in% names(out)))
})

test_that("chooseDomSecFuelClasses() ranks by total AGB, ties broken alphabetically", {
  testthat::local_mocked_bindings(cohortsToFuelClasses = function(...) fuelRas())
  out <- chooseDomSecFuelClasses(
    cohortData = NULL, pixelGroupMap = NULL, flammableRTM = fuelRas()[[1]], landcoverDT = landcover(),
    sppEquiv = NULL, fuelClassCol = "FuelClass", sppEquivCol = "LandR", cutoffForYoungAge = 15
  )
  expect_identical(out, list(domClass = "Pice_mar", secClass = "Pinu_ban"))
})

test_that("chooseDomSecFuelClasses() with no tree fuel classes returns NA for both", {
  noTrees <- function() {
    r <- terra::rast(nrows = 2, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 2, crs = "EPSG:3978")
    terra::values(r) <- rep(0, 6)
    names(r) <- "youngAge"
    r
  }
  testthat::local_mocked_bindings(cohortsToFuelClasses = function(...) noTrees())
  out <- chooseDomSecFuelClasses(
    cohortData = NULL, pixelGroupMap = NULL, flammableRTM = noTrees(), landcoverDT = landcover(),
    sppEquiv = NULL, fuelClassCol = "FuelClass", sppEquivCol = "LandR", cutoffForYoungAge = 15
  )
  expect_identical(out, list(domClass = NA_character_, secClass = NA_character_))
})
