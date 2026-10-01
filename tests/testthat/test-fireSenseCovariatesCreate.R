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

test_that("youngAge = FALSE builds no youngAge column and leaves fuels and nfLCC unzeroed", {
  withr::local_package("terra")
  withr::local_package("data.table")

  pixelGroupMap <- rast(nrows = 2, ncols = 2, vals = c(1L, 2L, 0L, 0L))
  flammableRTM <- rast(pixelGroupMap, vals = 1)
  sppEquiv <- data.table(LandR = c("Pice_mar", "Pinu_ban"), FuelClass = c("spruce", "pine"))
  ## pixel 1: old spruce; pixel 2: a 5-year-old pine stand (young at the cutoff)
  cohortData <- data.table(pixelGroup = 1:2, speciesCode = c("Pice_mar", "Pinu_ban"),
                           age = c(100L, 5L), B = c(3000L, 3000L))
  landcoverDT <- data.table(pixelID = 1:4, nfLCC_40 = c(0, 0, 1, 1))
  args <- list(cohortData = cohortData, pixelGroupMap = pixelGroupMap, flammableRTM = flammableRTM,
               sppEquiv = sppEquiv, landcoverDT = landcoverDT, fuelClassCol = "FuelClass",
               sppEquivCol = "LandR", missingLCCgroup = "nfLCC_40",
               nonForestedLCCGroups = c(nfLCC_40 = 40),
               nonForest_timeSinceDisturbance = c(100, 5, 5, 100),
               cutoffForYoungAge = 15, nonForestCanBeYoungAge = TRUE,
               studyAreaName = "test", useCache = FALSE, fuelCovariates = "species")

  yes <- do.call(fireSenseCovariatesCreate, args)
  setkey(yes, pixelID)
  expect_true("youngAge" %in% names(yes))   ## default unchanged: young pixels are zeroed
  expect_equal(yes[pixelID == 3, nfLCC_40], 0)
  expect_equal(yes[pixelID == 2, youngAge], 1)

  no <- do.call(fireSenseCovariatesCreate, c(args, list(youngAge = FALSE)))
  setkey(no, pixelID)
  expect_false("youngAge" %in% names(no))
  expect_equal(no[pixelID == 3, nfLCC_40], 1)              ## young non-forest kept
  expect_true(no[pixelID == 2, pine] > no[pixelID == 1, pine]) ## young pine stand keeps its biomass
  expect_equal(no[pixelID == 2, pine], log(3000))
})
