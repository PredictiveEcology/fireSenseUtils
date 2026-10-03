## fuelCovariates = "domSecWetland": fireSenseCovariatesCreate() collapses per-species fuel-class
## columns to dom_agb_<class>, sec_agb_<class> (the two classes with the most total treed AGB),
## and treedWetland_agb (all tree AGB on treed-wetland pixels, removed from dom/sec there).
## The remaining classes are not covariates: there is no other_agb.
## Prediction must be able to force the same dom/sec classes
## the fit chose, even when a different class dominates the prediction area.

skip_if_not_installed("terra")

## 2x3 landscape, 3 tree fuel classes.
##            px1  px2  px3  px4  px5  px6
## Pice_mar:  500  800    0  900    0  700   -> total 2900 (dominant)
## Pinu_ban:  200  100    0  300    0  200   -> total  800 (secondary)
## Popu_tre:   50   50    0  100    0   50   -> total  250 (neither dom nor sec: dropped)
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

covs <- function(rstLCC = NULL, domClass = NULL, secClass = NULL, fuelCovariates = "domSecWetland", treedWetland = TRUE) {
  testthat::local_mocked_bindings(cohortsToFuelClasses = function(...) fuelRas())
  args <- list(cohortData = NULL, pixelGroupMap = NULL, flammableRTM = fuelRas()[[1]], sppEquiv = NULL,
               landcoverDT = landcover(), fuelClassCol = "FuelClass", sppEquivCol = "LandR",
               missingLCCgroup = "nfLCC_50", nonForestedLCCGroups = list(nfLCC_50 = 50),
               nonForest_timeSinceDisturbance = NULL, cutoffForYoungAge = 15, nonForestCanBeYoungAge = FALSE,
               studyAreaName = "test", useCache = FALSE, fuelCovariates = fuelCovariates,
               domClass = domClass, secClass = secClass, treedWetland = treedWetland)
  if (!is.null(rstLCC)) args$rstLCC <- rstLCC
  do.call(fireSenseCovariatesCreate, args)
}

test_that("dom/sec are chosen by total AGB and named after the class; nothing else is kept", {
  out <- covs()[order(pixelID)]
  expect_true(all(c("dom_agb_Pice_mar", "sec_agb_Pinu_ban") %in% names(out)))
  expect_false(any(c("Pice_mar", "Pinu_ban", "Popu_tre") %in% names(out)))

  ## pixel 2: not wetland, not young -- straightforward collapse
  expect_equal(out[pixelID == 2, dom_agb_Pice_mar], logMinB(800))
  expect_equal(out[pixelID == 2, sec_agb_Pinu_ban], logMinB(100))

  expect_equal(out[pixelID == 6, dom_agb_Pice_mar], logMinB(700))
  expect_equal(out[pixelID == 6, sec_agb_Pinu_ban], logMinB(200))

  ## pixel 3: non-forest, no tree AGB at all
  expect_equal(out[pixelID == 3, dom_agb_Pice_mar], logMinB(0))
  expect_equal(out[pixelID == 3, sec_agb_Pinu_ban], logMinB(0))
})

test_that("the fuel columns produced are exactly dom_agb_*, sec_agb_* and treedWetland_agb", {
  fuelCols <- function(out) setdiff(names(out), c("pixelID", "nfLCC_50", "youngAge"))
  expect_setequal(fuelCols(covs(rstLCC = lccRas())), c("dom_agb_Pice_mar", "sec_agb_Pinu_ban", "treedWetland_agb"))
  expect_setequal(fuelCols(covs()), c("dom_agb_Pice_mar", "sec_agb_Pinu_ban"))
  expect_false("other_agb" %in% names(covs(rstLCC = lccRas())))
})

test_that("a treed-wetland pixel's AGB is only in treedWetland_agb, removed from dom/sec", {
  out <- covs(rstLCC = lccRas())[order(pixelID)]
  expect_true("treedWetland_agb" %in% names(out))

  ## pixel 1: wetland, not young -- all AGB (500 + 200 + 50 = 750) moves to treedWetland_agb
  expect_equal(out[pixelID == 1, treedWetland_agb], logMinB(750))
  expect_equal(out[pixelID == 1, dom_agb_Pice_mar], logMinB(0))
  expect_equal(out[pixelID == 1, sec_agb_Pinu_ban], logMinB(0))

  ## pixel 2: not wetland -- treedWetland_agb is 0, dom/sec unaffected
  expect_equal(out[pixelID == 2, treedWetland_agb], logMinB(0))
  expect_equal(out[pixelID == 2, dom_agb_Pice_mar], logMinB(800))
})

test_that("a young pixel has all its AGB terms (and nfLCC) at 0, and youngAge = 1", {
  out <- covs(rstLCC = lccRas())[order(pixelID)]
  ## pixel 4: young AND wetland -- youngAge wins: nothing survives, including treedWetland_agb
  expect_equal(out[pixelID == 4, youngAge], 1)
  expect_equal(out[pixelID == 4, dom_agb_Pice_mar], logMinB(0))
  expect_equal(out[pixelID == 4, sec_agb_Pinu_ban], logMinB(0))
  expect_equal(out[pixelID == 4, treedWetland_agb], logMinB(0))
})

test_that("prediction can force the fit's dom/sec classes even when a different class dominates here", {
  ## Popu_tre is the smallest class in this fixture (total 250), but a prediction is told it was
  ## the fit's dominant class: the columns built must follow that, not this area's own ranking.
  out <- covs(domClass = "Popu_tre", secClass = "Pice_mar")[order(pixelID)]
  expect_true(all(c("dom_agb_Popu_tre", "sec_agb_Pice_mar") %in% names(out)))
  expect_equal(out[pixelID == 2, dom_agb_Popu_tre], logMinB(50))
  expect_equal(out[pixelID == 2, sec_agb_Pice_mar], logMinB(800))
})

test_that("a forced class absent from this area's fuel classes gets a zero column, not an error", {
  out <- covs(domClass = "Betu_pap", secClass = "Pice_mar")[order(pixelID)]
  expect_true("dom_agb_Betu_pap" %in% names(out))
  expect_equal(out[pixelID == 2, dom_agb_Betu_pap], logMinB(0))
  expect_equal(out[pixelID == 2, sec_agb_Pice_mar], logMinB(800))
})

test_that("fuelCovariates = \"species\" (the default) is unchanged: per-species columns, no dom/sec", {
  out <- covs(fuelCovariates = "species")
  expect_true(all(c("Pice_mar", "Pinu_ban", "Popu_tre") %in% names(out)))
  expect_false(any(c("dom_agb_Pice_mar", "sec_agb_Pinu_ban") %in% names(out)))
  expect_false("other_agb" %in% names(out))
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

## An ELF with no tree fuel class at all (only youngAge in the fuel raster): collapseFuelClassesToDomSec()
## takes its no-dom/no-sec branch, so there are no fuel columns to build.
noTreeRas <- function() {
  r <- terra::rast(nrows = 2, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 2, crs = "EPSG:3978")
  terra::values(r) <- rep(0, 6)
  names(r) <- "youngAge"
  r
}

covsNoTrees <- function(rstLCC = NULL) {
  testthat::local_mocked_bindings(cohortsToFuelClasses = function(...) noTreeRas())
  args <- list(cohortData = NULL, pixelGroupMap = NULL, flammableRTM = noTreeRas(), sppEquiv = NULL,
               landcoverDT = landcover(), fuelClassCol = "FuelClass", sppEquivCol = "LandR",
               missingLCCgroup = "nfLCC_50", nonForestedLCCGroups = list(nfLCC_50 = 50),
               nonForest_timeSinceDisturbance = NULL, cutoffForYoungAge = 15, nonForestCanBeYoungAge = FALSE,
               studyAreaName = "test", useCache = FALSE, fuelCovariates = "domSecWetland")
  if (!is.null(rstLCC)) args$rstLCC <- rstLCC
  do.call(fireSenseCovariatesCreate, args)
}

test_that("an ELF with no tree fuel class has no dom/sec columns, and no other_agb", {
  out <- covsNoTrees()
  expect_setequal(setdiff(names(out), c("pixelID", "nfLCC_50", "youngAge")), character())
  expect_false(any(grepl("^(dom|sec)_agb_|^other_agb$", names(out))))
  expect_identical(attr(out, "fuelClassRoles"), list(domClass = NA_character_, secClass = NA_character_))
})

test_that("with no tree fuel class, treedWetland_agb is the only fuel column and is all zero", {
  out <- covsNoTrees(rstLCC = lccRas())
  expect_setequal(setdiff(names(out), c("pixelID", "nfLCC_50", "youngAge")), "treedWetland_agb")
  expect_true(all(out$treedWetland_agb == logMinB(0)))
})

test_that("a pixel whose only AGB is in a dropped class is not flagged as a missing land-cover group", {
  ## Popu_tre is neither dom nor sec; px5's only AGB is Popu_tre (not wetland here: lcc 210)
  testthat::local_mocked_bindings(cohortsToFuelClasses = function(...) {
    r <- fuelRas()
    terra::values(r) <- cbind(
      c(500, 800, 0, 900, 0, 700),
      c(200, 100, 0, 300, 0, 200),
      c(50, 50, 0, 100, 40, 50),
      c(0, 0, 0, 1, 0, 0)
    )
    r
  })
  out <- fireSenseCovariatesCreate(
    cohortData = NULL, pixelGroupMap = NULL, flammableRTM = fuelRas()[[1]], sppEquiv = NULL,
    landcoverDT = landcover(), fuelClassCol = "FuelClass", sppEquivCol = "LandR",
    missingLCCgroup = "nfLCC_50", nonForestedLCCGroups = list(nfLCC_50 = 50),
    nonForest_timeSinceDisturbance = NULL, cutoffForYoungAge = 15, nonForestCanBeYoungAge = FALSE,
    studyAreaName = "test", useCache = FALSE, fuelCovariates = "domSecWetland"
  )
  out <- out[order(pixelID)]
  expect_equal(out[pixelID == 5, nfLCC_50], 0)
  expect_equal(out[pixelID == 5, dom_agb_Pice_mar], logMinB(0))
})

test_that("forcing domClass = NA drops every tree class, even when this area has them", {
  out <- covs(domClass = NA_character_, rstLCC = lccRas())[order(pixelID)]
  expect_setequal(setdiff(names(out), c("pixelID", "nfLCC_50", "youngAge")), "treedWetland_agb")
  ## wetland pixel 1 still carries all its tree AGB (500 + 200 + 50) in treedWetland_agb
  expect_equal(out[pixelID == 1, treedWetland_agb], logMinB(750))
  expect_equal(out[pixelID == 2, treedWetland_agb], logMinB(0))
  expect_identical(attr(out, "fuelClassRoles"), list(domClass = NA_character_, secClass = NA_character_))
})

test_that("treedWetland = FALSE: no treedWetland_agb, and wetland AGB stays in dom/sec", {
  ## an ELF with too little treed wetland to estimate it (fireSense_dataPrepFit's minCovariateProp)
  out <- covs(rstLCC = lccRas(), treedWetland = FALSE)[order(pixelID)]
  expect_false(treedWetlandAgbTxt %in% names(out))
  ## pixel 1 is class 81: its AGB is ordinary forest, exactly as on pixel 2's upland
  expect_equal(out[pixelID == 1, dom_agb_Pice_mar], logMinB(500))
  expect_equal(out[pixelID == 1, sec_agb_Pinu_ban], logMinB(200))
  ## and rstLCC then changes nothing at all
  expect_equal(as.data.frame(out), as.data.frame(covs(treedWetland = FALSE)[order(pixelID)]))
})
