## dataYearFor() is the one rule for matching a year to a vegetation data year: the data year
## closest to, and not after, it. fireSense_dataPrepFit (fire years, the fit's structure decisions)
## and fireSense_dataPrepPredict (the vegetation of time(sim)) both call it.

test_that("dataYearFor is the data year closest to and not after each year", {
  dataYears <- c(1985L, 1990L, 2000L, 2010L, 2020L)
  expect_identical(dataYearFor(dataYears, c(1985L, 1989L, 1990L, 1999L, 2000L, 2015L, 2020L, 2024L)),
                   c(1985L, 1985L, 1990L, 1990L, 2000L, 2010L, 2020L, 2020L))
  expect_identical(dataYearFor(rev(dataYears), 2015L), 2010L) # any order
  expect_identical(dataYearFor(2020, 2021.5), 2020L)        # a decimal time is in its calendar year
  expect_error(dataYearFor(dataYears, 1984L), "before the first dataYear")
})

test_that("dataYearElement names the element a year uses", {
  x <- list(year2010 = "a", year2020 = "b")
  expect_identical(dataYearElement(x, 2010L), "year2010")
  expect_identical(dataYearElement(x, 2015L), "year2010")
  expect_identical(dataYearElement(x, 2024L), "year2020")
  expect_identical(dataYearElement(c(`2020` = 1, `1985` = 2), 1990L), "1985")
  expect_error(dataYearElement(x, 2005L), "before the first dataYear")
})
