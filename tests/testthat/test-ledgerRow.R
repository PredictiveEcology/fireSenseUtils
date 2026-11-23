## The ledger-row readers (ledger*()), which fireSense_spreadPredict and spreadFuelEffect() share. Rows are small synthetic ones;
## a column a row lacks is absent or NA, as a ledger bound from older and newer rows gives.

test_that("a row without link, climate or centre columns is read as it always was", {
  sa <- data.frame(polygonID = "A")
  expect_null(ledgerSpreadLink(sa, 1L)); expect_null(ledgerClimateRef(sa, 1L)); expect_null(ledgerCovCentre(sa, 1L))
  expect_null(ledgerPowerLawGamma(sa, 1L)); expect_null(ledgerPowerLawPc50(sa, 1L))
  expect_identical(ledgerSpreadClimate(sa, 1L), "predictor")
  expect_identical(ledgerCovMinMax(sa, 1L, "default"), "default")
  expect_identical(ledgerLinkArgs(sa, 1L), list(powerLawShape = spreadPowerLawShape))
  ## NA, as when rows with and without the column are bound
  sa$covCentre_spread <- NA; sa$spreadLink <- NA_character_
  expect_null(ledgerCovCentre(sa, 1L)); expect_null(ledgerLinkArgs(sa, 1L)$link)
})

test_that("a power-law row gives its link, gamma and shape; a missing gamma stops, naming the row", {
  sa <- data.frame(polygonID = "A", spreadLink = "powerLaw", powerLawGamma_spread = 2.2, powerLawShape_spread = 0.5)
  expect_identical(ledgerLinkArgs(sa, 1L), list(link = "powerLaw", powerLawGamma = 2.2, powerLawShape = 0.5))
  sa$powerLawGamma_spread <- NA_real_
  expect_error(ledgerPowerLawGamma(sa, 1L), "ledger row 1 has spreadLink 'powerLaw' but no valid")
})

test_that("the climate mode of a row without spreadClimate follows its climateRef and link", {
  sa <- data.frame(polygonID = "A", spreadLink = "powerLawSmooth", lowerAsymptote_spread = 0.1)
  sa$climateRef_spread <- list(list(MDC = 0.4))
  expect_identical(ledgerClimateRef(sa, 1L), list(MDC = 0.4))
  expect_identical(ledgerSpreadClimate(sa, 1L), "ceiling")
  expect_identical(ledgerLowerAsymptote(sa, 1L, floor = 0.1), 0.1)
  sa$spreadClimate <- "bad"
  expect_error(ledgerSpreadClimate(sa, 1L), "unknown")
})
