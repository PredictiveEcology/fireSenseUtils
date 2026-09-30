## youngAge is resolved per fire year: time since disturbance at year t is the smaller of the data
## year's TSD aged to t and the years since the last fire before t (from ALL fires, not just the
## fit buffers). NA time since disturbance is never young.

test_that("youngAgeAtYear: a fire between the data year and t makes an old pixel young at t", {
  tsd <- data.table::data.table(pixelID = 1:3, tsd = c(100, 100, 100))
  fires <- list("2008" = 1L, "2009" = 2L)
  ## data year 2005; at t = 2010 pixel 1 burned 2 years ago, pixel 2 one year ago, pixel 3 never
  ya <- youngAgeAtYear(tsd, dataYear = 2005, year = 2010, firePixelsByYear = fires,
                       cutoffForYoungAge = 15)
  expect_equal(ya, c(1L, 1L, 0L))
  ## before the fires it was old
  expect_equal(youngAgeAtYear(tsd, 2005, 2008, fires, 15), c(0L, 0L, 0L))
  ## a fire in year t itself is not yet in t's age
  expect_equal(youngAgeAtYear(tsd, 2005, 2009, fires, 15), c(1L, 0L, 0L))
})

test_that("youngAgeAtYear: a pixel young at the data year ages out after the cutoff", {
  tsd <- data.table::data.table(pixelID = 1:2, tsd = c(5, 20))
  ## pixel 1 has TSD 5 in 2000: 5 + (t - 2000) <= 15 until t = 2010
  expect_equal(youngAgeAtYear(tsd, 2000, 2000, list(), 15)[1], 1L)
  expect_equal(youngAgeAtYear(tsd, 2000, 2010, list(), 15)[1], 1L)
  expect_equal(youngAgeAtYear(tsd, 2000, 2011, list(), 15)[1], 0L)
  ## pixel 2 was already old
  expect_equal(youngAgeAtYear(tsd, 2000, 2000, list(), 15)[2], 0L)
})

test_that("youngAgeAtYear: NA time since disturbance stays not young, even when fires cover it", {
  tsd <- data.table::data.table(pixelID = 1:2, tsd = c(NA, 3))
  ya <- youngAgeAtYear(tsd, 2000, 2003, list("2002" = 1L), 15)
  expect_equal(ya, c(0L, 1L))
  expect_false(anyNA(ya))
})

test_that("youngAgeAtYear: fires outside the requested pixels still reset, and pixelID selects rows", {
  tsd <- data.table::data.table(pixelID = 1:5, tsd = 100)
  ## pixels 4 and 5 burned in 2001; only pixels 2 and 4 are asked for
  ya <- youngAgeAtYear(tsd, 2000, 2003, list("2001" = c(4L, 5L, 99L)), 15, pixelID = c(2L, 4L))
  expect_equal(ya, c(0L, 1L))
})

test_that("youngAgeAtYear: young means at or below the cutoff", {
  tsd <- data.table::data.table(pixelID = 1:2, tsd = 100)
  fires <- list("1990" = 1L, "1989" = 2L)
  ## 2005 - 1990 = 15 (young); 2005 - 1989 = 16 (not)
  expect_equal(youngAgeAtYear(tsd, 1985, 2005, fires, 15), c(1L, 0L))
})

test_that("firePixelsByYear: pixel IDs burned in each year, from polygons or a fire-year raster", {
  withr::local_package("terra")
  tmpl <- rast(nrows = 2, ncols = 2, vals = 1)
  ## polygon covering the left column (cells 1 and 3)
  poly <- as.polygons(ext(xmin(tmpl), 0, ymin(tmpl), ymax(tmpl)), crs = crs(tmpl))
  out <- firePixelsByYear(firePolys = list(year2001 = poly, year2002 = NULL), template = tmpl)
  expect_equal(out[["2001"]], c(1L, 3L))
  expect_false("2002" %in% names(out))

  fr <- rast(nrows = 2, ncols = 2, vals = c(2001, NA, 2003, 2001))
  out2 <- firePixelsByYear(fireRaster = fr)
  expect_equal(out2[["2001"]], c(1L, 4L))
  expect_equal(out2[["2003"]], 3L)
})
