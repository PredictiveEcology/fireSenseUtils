test_that("defaultFireYears() runs from 1986 to the latest year with historical climate", {
  skip_if_not_installed("climateData")
  fy <- defaultFireYears()
  expect_type(fy, "integer")
  expect_identical(fy, 1986L:climateData::latestHistoricalYear())
  expect_identical(min(fy), 1986L) # 1985, the first SCANFI year, has no earlier vegetation
})
