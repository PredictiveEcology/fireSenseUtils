## A fire year uses the latest vegetation table strictly before it (fireSense_dataPrepFit builds
## table "1990" from the pixels of fire years 1991-1995), so fire year 1990 reads table "1985".

library(data.table)

vegYearArgs <- function(yr, nonAnnual, pix) {
  list(annDTx1000 = data.table(pixelID = pix), nonAnnualDTx1000 = nonAnnual, indexNonAnnual = NULL,
       yr = yr, covMinMax = NULL, mutuallyExclusive = NULL, colsToUse = "a", doAssertions = FALSE,
       logisticPars = c(0.2, 1, 1), covPars = 1, maxFireSpread = 0.28, lowerSpreadProb = 0.1)
}

test_that("dataYearBefore is the latest data year strictly before each year", {
  dataYears <- c(1985L, 1990L, 2000L)
  expect_identical(dataYearBefore(dataYears, c(1986L, 1990L, 1991L, 2000L, 2001L)),
                   c(1985L, 1985L, 1990L, 1990L, 2000L))
  expect_identical(dataYearBefore(rev(dataYears), 1995L), 1990L) # any order
  expect_identical(dataYearFor(dataYears, 1990L), 1990L)         # dataYearFor keeps <=
  expect_error(dataYearBefore(dataYears, 1985L), "no vegetation data")
})

test_that("dataYearElement(before = TRUE) names the element of the data year strictly before", {
  x <- list(year2010 = "a", year2020 = "b")
  expect_identical(dataYearElement(x, 2020L, before = TRUE), "year2010")
  expect_identical(dataYearElement(x, 2021L, before = TRUE), "year2020")
  expect_identical(dataYearElement(x, 2020L), "year2020")
  expect_error(dataYearElement(x, 2010L, before = TRUE), "before the first dataYear")
})

test_that("spreadProbFromIntegerCovs reads the table strictly before the fire year", {
  nonAnnual <- list(`1985` = data.table(pixelID = 1:3, a = c(100L, 200L, 300L)),
                    `1990` = data.table(pixelID = 4L, a = 999L))
  got <- do.call(spreadProbFromIntegerCovs, vegYearArgs(1990L, nonAnnual, 1:3))
  expect_false(anyNA(got$a))
  expect_equal(got$a, c(0.1, 0.2, 0.3))
  ## 1991 reads table "1990"
  got <- do.call(spreadProbFromIntegerCovs, vegYearArgs(1991L, nonAnnual, 4L))
  expect_equal(got$a, 0.999)
})

test_that("spreadProbFromIntegerCovs stops, naming the fire year and table, when pixels are missing", {
  nonAnnual <- list(`1985` = data.table(pixelID = 1:2, a = c(100L, 200L)))
  expect_error(do.call(spreadProbFromIntegerCovs, vegYearArgs(1990L, nonAnnual, 1:3)),
               "fire year 1990.*table \"1985\"")
})
