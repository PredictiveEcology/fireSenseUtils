test_that("defaultDataYears() is every SCANFI V2 year before the last default fire year", {
  testthat::local_mocked_bindings(latestHistoricalYear = function() 2024L, .package = "climateData")
  expect_identical(defaultDataYears(), seq(1985L, 2020L, 5L))
})

test_that("defaultDataYears() gains 2025 once historical climate passes it", {
  testthat::local_mocked_bindings(latestHistoricalYear = function() 2026L, .package = "climateData")
  expect_identical(defaultDataYears(), seq(1985L, 2025L, 5L))
})

test_that("a data year is never the last fire year", {
  testthat::local_mocked_bindings(latestHistoricalYear = function() 2020L, .package = "climateData")
  expect_identical(max(defaultDataYears()), 2015L)
})

test_that("scanfiV2Years are the SCANFI V2 years", {
  expect_identical(scanfiV2Years, seq(1985L, 2025L, 5L))
})
