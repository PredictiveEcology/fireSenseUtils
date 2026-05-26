## Tests for updateStackYearNames

test_that("updateStackYearNames: prefixes 'year' to layer names matching desiredYears", {
  withr::local_package("terra")
  r <- terra::rast(nrows = 2, ncols = 2, nlyrs = 3)
  names(r) <- c("layer_2001", "layer_2002", "layer_2003")
  out <- updateStackYearNames(r, desiredYears = c(2001L, 2002L, 2003L))
  expect_equal(names(out), c("year2001", "year2002", "year2003"))
})

test_that("updateStackYearNames: errors when stack names lack a 4-digit year suffix", {
  withr::local_package("terra")
  r <- terra::rast(nrows = 2, ncols = 2, nlyrs = 2)
  names(r) <- c("layer_one", "layer_two")
  expect_error(
    updateStackYearNames(r, desiredYears = c(2001L, 2002L)),
    regexp = "4 digit year"
  )
})

test_that("updateStackYearNames: errors when desiredYears lack a 4-digit year", {
  withr::local_package("terra")
  r <- terra::rast(nrows = 2, ncols = 2, nlyrs = 2)
  names(r) <- c("layer_2001", "layer_2002")
  expect_error(
    updateStackYearNames(r, desiredYears = c("yr1", "yr2")),
    regexp = "4 digit year"
  )
})

test_that("updateStackYearNames: errors when stack years do not match desiredYears", {
  withr::local_package("terra")
  r <- terra::rast(nrows = 2, ncols = 2, nlyrs = 2)
  names(r) <- c("layer_2001", "layer_2002")
  expect_error(
    updateStackYearNames(r, desiredYears = c(2003L, 2004L))
  )
})
