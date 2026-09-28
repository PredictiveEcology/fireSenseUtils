## fireSenseCovariatesCreate(): youngAge must be mutually exclusive with every other non-climate
## covariate, including non-forest land-cover columns for pixels that only become youngAge through
## the non-forest time-since-disturbance path (nonForestCanBeYoungAge = TRUE). Before the fix,
## makeMutuallyExclusive() ran before that path set youngAge, so those pixels kept their nfLCC_*
## value.

test_that("a young non-forest pixel ends with nfLCC_* = 0 and youngAge = 1", {
  withr::local_package("terra")
  withr::local_package("data.table")

  ## 2x2 landscape, no tree cohorts at all: fuelClassesRas has only a youngAge layer (all 0),
  ## exactly like "cohortsToFuelClasses with no tree species and no required classes gives
  ## youngAge only" in test-fuelClasses.R.
  pixelGroupMap <- rast(nrows = 2, ncols = 2, vals = 0L)
  flammableRTM <- rast(pixelGroupMap, vals = 1)
  sppEquiv <- data.table(LandR = character(), FuelClass = character())
  noCohorts <- data.table(pixelGroup = integer(), speciesCode = character(),
                          age = integer(), B = integer())

  ## pixel 1: non-forest (nfLCC_40), disturbed 5 years ago -- young
  ## pixel 2: non-forest (nfLCC_50), disturbed 30 years ago -- not young
  ## pixel 3: forested LCC absent from cohortData (missingLCCgroup)
  ## pixel 4: non-forest (nfLCC_40), disturbed 200 years ago -- not young
  landcoverDT <- data.table(pixelID = 1:4,
                            nfLCC_40 = c(1, 0, 0, 1),
                            nfLCC_50 = c(0, 1, 0, 0))
  nonForestedLCCGroups <- c(nfLCC_40 = 40, nfLCC_50 = 50)
  nonForest_timeSinceDisturbance <- c(5, 30, NA, 200)

  covs <- suppressWarnings( ## max(age) over zero rows, as in the cohortsToFuelClasses tests
    fireSenseCovariatesCreate(
      cohortData = noCohorts,
      pixelGroupMap = pixelGroupMap,
      flammableRTM = flammableRTM,
      sppEquiv = sppEquiv,
      landcoverDT = landcoverDT,
      fuelClassCol = "FuelClass",
      sppEquivCol = "LandR",
      missingLCCgroup = "nfLCC_40",
      nonForestedLCCGroups = nonForestedLCCGroups,
      nonForest_timeSinceDisturbance = nonForest_timeSinceDisturbance,
      cutoffForYoungAge = 15,
      nonForestCanBeYoungAge = TRUE,
      studyAreaName = "test",
      useCache = FALSE
    )
  )
  setkey(covs, pixelID)

  expect_equal(covs[pixelID == 1, youngAge], 1)
  expect_equal(covs[pixelID == 1, nfLCC_40], 0)   ## zeroed: the pixel is young
  expect_equal(covs[pixelID == 1, nfLCC_50], 0)

  expect_equal(covs[pixelID == 2, youngAge], 0)
  expect_equal(covs[pixelID == 2, nfLCC_50], 1)   ## not young: kept

  expect_equal(covs[pixelID == 4, youngAge], 0)
  expect_equal(covs[pixelID == 4, nfLCC_40], 1)   ## not young: kept
})
