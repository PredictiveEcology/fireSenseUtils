test_that("defaultFireYears() runs from 1985 to the latest year with historical climate", {
  skip_if_not_installed("climateData")
  fy <- defaultFireYears()
  expect_type(fy, "integer")
  expect_identical(fy, 1985L:climateData::latestHistoricalYear())
  expect_identical(min(fy), 1985L)
})
