## A spread-ledger row can have FuelClass = NA for a species with no cohorts at fit time, and
## fireSense_dataPrepPredict passes that NA in `requiredFuelClasses`. An NA fuel class is "not a
## fuel": no column named NA, and the other classes are exactly what they are without the species.

skip_if_not_installed("terra")

pgm <- function() terra::rast(nrows = 2, ncols = 2, vals = c(1, 1, 2, 2))
flam <- function() terra::rast(pgm(), vals = 1)
seNA <- function() data.table::data.table(LandR = c("Pice_mar", "Pinu_res", "Popu_tre"),
                                          FuelClass = c("Pice_mar", NA, "Popu_tre"))
seNoNA <- function() seNA()[!is.na(FuelClass)]
cohorts <- function(withNAspp) {
  sp <- c("Pice_mar", "Popu_tre", if (withNAspp) "Pinu_res")
  data.table::data.table(pixelGroup = c(1L, 2L, if (withNAspp) 2L), speciesCode = factor(sp),
                         age = 50L, B = c(100L, 200L, if (withNAspp) 300L))
}
fuel <- function(cd, se, required) {
  cohortsToFuelClasses(cd, pgm(), flam(), sppEquiv = se, sppEquivCol = "LandR", cutoffForYoungAge = 15,
                       fuelClassCol = "FuelClass", requiredFuelClasses = required, asTable = TRUE)
}

test_that("an NA-class species adds no NA column and leaves the other classes unchanged", {
  for (withNAspp in c(FALSE, TRUE)) {
    ref <- fuel(cohorts(FALSE), seNoNA(), seNoNA()$FuelClass)
    out <- fuel(cohorts(withNAspp), seNA(), seNA()$FuelClass)
    expect_false(anyNA(names(out)))
    expect_false("NA" %in% names(out))
    expect_equal(out, ref)
  }
})

test_that("fireSenseCovariatesCreate reaches the row check with an NA required class", {
  se <- seNA()
  pg <- pgm()
  lc <- data.table::data.table(pixelID = 1:4, nfLCC_50 = 0)
  for (withNAspp in c(FALSE, TRUE)) {
    out <- fireSenseCovariatesCreate(
      cohortData = cohorts(withNAspp), pixelGroupMap = pg, flammableRTM = flam(), sppEquiv = se,
      landcoverDT = lc, fuelClassCol = "FuelClass", requiredFuelClasses = se$FuelClass,
      sppEquivCol = "LandR", missingLCCgroup = "nfLCC_50", nonForestedLCCGroups = list(nfLCC_50 = 50),
      nonForest_timeSinceDisturbance = NULL, cutoffForYoungAge = 15, nonForestCanBeYoungAge = FALSE,
      studyAreaName = "test", useCache = FALSE)
    expect_false(any(c("NA", NA) %in% names(out)))
    expect_true(all(c("Pice_mar", "Popu_tre") %in% names(out)))
  }
})
