test_that("climateRasterToDataTable merges any number of climate variables", {
  mkStack <- function(offset) {
    r <- terra::rast(nrows = 2, ncols = 2, nlyr = 2, vals = offset + 1:8)
    names(r) <- c("1985", "1986")
    r
  }
  allVars <- c("CMD", "CMDsm", "cumMDC", "MDC")
  for (n in 1:4) {
    vars <- allVars[seq_len(n)]
    rasters <- setNames(lapply(seq_len(n), function(i) mkStack(100 * i)), vars)
    out <- climateRasterToDataTable(rasters)
    expect_setequal(names(out), c("pixelID", "year", vars))
    expect_equal(nrow(out), 8L) # 4 cells x 2 years
    for (i in seq_len(n)) {
      expect_equal(sort(out[[vars[i]]]), 100 * i + 1:8)
    }
    # the variables of one pixel-year row come from the same pixel-year
    expect_true(all(out[[vars[n]]] - out[[vars[1]]] == 100 * (n - 1)))
  }
})

test_that("climateRasterToDataTable with Index returns the full result's rows for those years and pixels", {
  yrs <- c("year1985", "year1986", "year1987")
  mkStack <- function(offset) {
    r <- terra::rast(nrows = 5, ncols = 5, nlyr = length(yrs), vals = offset + seq_len(75) + 0.123)
    names(r) <- yrs
    r
  }
  rasters <- list(CMD = mkStack(0), MDC = mkStack(1000))
  full <- climateRasterToDataTable(rasters)

  # buffer pixels differ by year; a year with no entry (year1987) has no rows
  Index <- list(year1985 = data.table::data.table(pixelID = c(3L, 7L, 20L)),
                year1986 = data.table::data.table(pixelID = c(7L, 11L)))
  out <- climateRasterToDataTable(rasters, Index = Index)

  wanted <- data.table::rbindlist(Index, idcol = "year")
  want <- full[wanted, on = c("year", "pixelID"), nomatch = NULL]
  want$year <- factor(want$year, levels = levels(full$year))
  expect_identical(nrow(out), 5L)
  data.table::setkeyv(out, c("year", "pixelID"))
  data.table::setkeyv(want, c("year", "pixelID"))
  expect_equal(out, want, ignore_attr = TRUE)
  expect_identical(names(out), names(full))
  expect_identical(class(out$year), class(full$year))

  # a plain vector of pixelIDs still keeps those pixels in every year
  expect_equal(nrow(climateRasterToDataTable(rasters, Index = c(3L, 7L))), 6L)
})

test_that("climateRasterToDataTable with no raster (a '~ 1' spread model) has no climate column", {
  noClim <- list()
  out <- climateRasterToDataTable(noClim)
  expect_s3_class(out, "data.table")
  expect_identical(names(out), c("pixelID", "year"))
  expect_identical(nrow(out), 0L)

  Index <- list(year1985 = data.table::data.table(pixelID = c(3L, 7L)),
                year1986 = data.table::data.table(pixelID = 7L))
  out <- climateRasterToDataTable(noClim, Index = Index)
  expect_identical(names(out), c("pixelID", "year"))
  expect_identical(nrow(out), 3L)

  # the join that builds a year's spread covariates keeps every pixel of the year
  fbl <- data.table::rbindlist(Index, idcol = "year")
  expect_identical(nrow(out[fbl, on = c("year", "pixelID"), nomatch = NULL]), 3L)

  # `year` is the factor the non-empty path gives (the layer names are its levels), so joins behave the same
  r <- terra::rast(nrows = 2, ncols = 2, nlyr = 2, vals = 1:8)
  names(r) <- names(Index)
  withClim <- climateRasterToDataTable(list(CMD = r), Index = Index)
  expect_identical(class(out$year), class(withClim$year))
  expect_identical(levels(out$year), levels(withClim$year))
})
