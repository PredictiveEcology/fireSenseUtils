## The ignition table has a row for every coarse pixel, buffer included, but the ignitions it counts start in
## the core. The Poisson fit therefore spread the core's ignitions over core + buffer, and the hindcast
## simulated about 0.8 of the observed ignitions. Each coarse pixel now carries the share of its flammable
## fine cells that are core (`coreFraction`) so fireSense_ignitionFit can weight the pixel by it.

coreFractionRasters <- function() {
  ## 8 x 8 fine cells, fact = 4: four coarse pixels. rasterToMatchLargeELF: core (2) in columns 1-6, buffer (1)
  ## in columns 7-8, so coarse pixel 1 (columns 1-4) is all core, pixel 2 (columns 5-8) is half, and the
  ## bottom row of coarse pixels the same.
  elf <- terra::rast(nrows = 8, ncols = 8, xmin = 0, xmax = 8, ymin = 0, ymax = 8, vals = 1)
  terra::values(elf) <- rep(c(rep(2, 6), 1, 1), 8)
  flam <- terra::rast(elf, vals = 1)
  list(elf = elf, flam = flam)
}

## coreFraction of a fire year, in pixelID order (data.table's [ is not used: a test is not data.table-aware)
fractionOf <- function(out, yr) {
  d <- as.data.frame(out)
  d <- d[d$year == yr, ]
  stats::setNames(d$coreFraction[order(d$pixelID)], sort(d$pixelID))
}

test_that("coreFraction is 1 in the core, 0 in pure buffer, fractional at the edge", {
  r <- coreFractionRasters()
  years <- list(`2010` = c("year2011", "year2012"))
  out <- coreFractionByYear(years, list(year2010 = r$flam), r$elf, fact = 4)
  expect_setequal(names(out), c("pixelID", "year", "coreFraction"))
  expect_identical(sort(unique(out$year)), c(2011, 2012))
  expect_equal(fractionOf(out, 2011), c(`1` = 1, `2` = 0.5, `3` = 1, `4` = 0.5))
  expect_equal(fractionOf(out, 2012), fractionOf(out, 2011))

  ## pure buffer: columns 5-8 become buffer
  terra::values(r$elf) <- rep(c(rep(2, 4), rep(1, 4)), 8)
  out <- coreFractionByYear(years, list(year2010 = r$flam), r$elf, fact = 4)
  expect_equal(unname(fractionOf(out, 2011)), c(1, 0, 1, 0))
})

test_that("coreFraction counts flammable cells only, and follows the data year's flammableRTM", {
  r <- coreFractionRasters()
  flam2 <- r$flam
  flam2[, 5:6] <- 0    # the core part of coarse pixel 2 is not flammable (0, as LandR::defineFlammable() gives): only buffer cells are left
  years <- list(`2010` = "year2011", `2020` = "year2021")
  out <- coreFractionByYear(years, list(year2010 = r$flam, year2020 = flam2), r$elf, fact = 4)
  expect_equal(unname(fractionOf(out, 2011)), c(1, 0.5, 1, 0.5))
  expect_equal(unname(fractionOf(out, 2021)), c(1, 0, 1, 0))
  flam2[, 5:6] <- NA   # NA is not flammable either
  out <- coreFractionByYear(years, list(year2010 = r$flam, year2020 = flam2), r$elf, fact = 4)
  expect_equal(unname(fractionOf(out, 2021)), c(1, 0, 1, 0))
})

test_that("a categorical (multi-ELF, no buffer) raster and no raster are all core", {
  r <- coreFractionRasters()
  cat <- terra::rast(r$elf, vals = rep(c(rep(0L, 4), rep(1L, 4)), 8))
  levels(cat) <- data.frame(id = 0:1, ELF = c("a", "b"))
  years <- list(`2010` = "year2011")
  expect_true(all(coreFractionByYear(years, list(year2010 = r$flam), cat, fact = 4)$coreFraction == 1))
  expect_true(all(coreFractionByYear(years, list(year2010 = r$flam), NULL, fact = 4)$coreFraction == 1))
})

test_that("addCoreFraction adds the column and ignitions are untouched", {
  r <- coreFractionRasters()
  dt <- data.table::data.table(pixelID = c(1L, 2L, 2L), year = c(2011, 2011, 2012), ignitions = c(1L, 0L, 2L))
  years <- list(`2010` = c("year2011", "year2012"))
  out <- addCoreFraction(dt, years, list(year2010 = r$flam), r$elf, fact = 4)
  expect_equal(out$coreFraction, c(1, 0.5, 0.5))
  expect_identical(out$ignitions, dt$ignitions)
  noFlammable <- data.table::copy(dt)
  data.table::set(noFlammable, NULL, "pixelID", 99L)
  expect_error(addCoreFraction(noFlammable, years, list(year2010 = r$flam), r$elf, 4), "coreFraction")
})
