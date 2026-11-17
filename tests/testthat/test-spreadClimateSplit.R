## spreadClimateSplit: the climate covariate as a between-year term (the year's ELF median) and a
## within-year term (the departure from it)

yearTable <- function() {
  data.table::data.table(year = rep(c("year1990", "year1991"), each = 4), pixelID = rep(1:4, 2),
                         CMD = c(10, 20, 30, 40, 100, 110, 120, 130))
}

test_that("the two terms are named <var>_yr and <var>_dev, per variable", {
  expect_identical(climateSplitCols(c("CMD", "Tmaxsm")), c("CMD_yr", "CMD_dev", "Tmaxsm_yr", "Tmaxsm_dev"))
  expect_identical(climateSplitBase(c("CMD_yr", "Tmaxsm_dev", "youngAge", "CMD")), c("CMD", "Tmaxsm", "youngAge", "CMD"))
  expect_identical(isClimateSplitCol(c("CMD_yr", "CMD_dev", "CMD", "class1")), c(TRUE, TRUE, FALSE, FALSE))
})

test_that("the between term is the year's median over the reference pixels, constant within a year; within = C - between", {
  dt <- yearTable()
  ## the reference is not the table's rows: the burnable pixels of the ELF (here two more per year)
  ref <- rbind(dt, data.table::data.table(year = c("year1990", "year1991"), pixelID = 5L, CMD = c(50, 140)))
  splitSpreadClimate(dt, "CMD", ref, by = "year")
  expect_false("CMD" %in% names(dt))
  expect_identical(dt$CMD_yr, rep(c(30, 120), each = 4))                 # median of 10,20,30,40,50 and 100..140
  expect_identical(dt$CMD_dev, c(10, 20, 30, 40, 100, 110, 120, 130) - dt$CMD_yr)
  expect_identical(tapply(dt$CMD_yr, dt$year, data.table::uniqueN)[[1]], 1L)
})

test_that("fit and prediction agree: a year's terms for the fitted (buffer) pixels equal those of the prediction's whole-ELF table", {
  burnable <- yearTable()                                  # the ELF's burnable pixels, two years
  fit <- burnable[burnable$pixelID %in% 1:2, ]                        # the fitted pixels: a buffer
  splitSpreadClimate(fit, "CMD", burnable, by = "year")
  pred <- burnable[burnable$year == "year1991", c("pixelID", "CMD")]    # prediction: the year's layer over the ELF
  splitSpreadClimate(pred, "CMD", data.table::copy(pred))
  fit91 <- fit[fit$year == "year1991", ]
  expect_identical(pred$CMD_yr[pred$pixelID %in% fit91$pixelID], fit91$CMD_yr)
  expect_identical(pred$CMD_dev[pred$pixelID %in% fit91$pixelID], fit91$CMD_dev)
  expect_identical(fit91$CMD_yr, rep(115, 2))              # the ELF's median, not the buffer's (100 and 110: 105)
})

test_that("NA climate is left out of the median", {
  dt <- data.table::data.table(pixelID = 1:4, CMD = c(1, 2, 3, NA))
  splitSpreadClimate(dt, "CMD", data.table::copy(dt))
  expect_identical(dt$CMD_yr, rep(2, 4))
  expect_identical(dt$CMD_dev, c(-1, 0, 1, NA))
})

test_that("the fixed ranges: _yr as the variable, _dev the same width centred on 0; given ranges win", {
  r <- climateSplitRanges(list(CMD = c(0, 100), Tmaxsm = c(0, 10)))
  expect_identical(r$CMD_yr, c(0, 100))
  expect_identical(r$CMD_dev, c(-50, 50))
  expect_identical(r$Tmaxsm_dev, c(-5, 5))
  expect_setequal(names(climateSplitRanges(climateCovRanges)),
                  c(names(climateCovRanges), climateSplitCols(names(climateCovRanges))))
  expect_identical(climateSplitRanges(list()), list())
  kept <- climateSplitRanges(list(CMD = c(0, 100), CMD_dev = c(-1, 1)))
  expect_identical(kept$CMD_dev, c(-1, 1))
})

test_that("a split fit's climate variable is found from its parameter names (spreadFitLedgerInputs)", {
  skip_if_not_installed("climateData")
  pars <- data.frame(maxAsymptote = 0.3, hillSlope1 = 1, inflectionPoint1 = 1, CMD_yr = 0.1, CMD_dev = 0.1,
                     youngAge = 0.1, class1 = 0.1)
  row <- data.frame(numIterations = 10L)
  row[[polygonIDTxt]] <- "6.6.2"
  row$params <- list(pars)
  row$sppEquiv <- list(data.table::data.table(LandR = "Pice_mar", FuelClass = "class1"))
  row$nonForestedLCCGroups <- list(list(nf_A = 40))
  row$missingLCCgroup <- "nf_A"
  out <- suppressWarnings(suppressMessages(
    spreadFitLedgerInputs(row, sppEquivCol = "LandR", fuelClassCol = "FuelClass", studyArea = NULL)))
  expect_identical(out[["6.6.2"]]$theseClimVars, "CMD")
  expect_identical(out[["6.6.2"]]$theseClimVarsNoUnderscore, "CMD")
})

test_that("the reference pixels are the flammable ones of the ELF: neither non-flammable nor blend-zone rows", {
  flam <- terra::rast(nrows = 2, ncols = 5, vals = c(1, 1, 1, 0, 1, 1, NA, 1, 1, 1))
  expect_identical(climateSplitRefPixels(1:10), 1:10)
  expect_identical(climateSplitRefPixels(1:10, flam), c(1L, 2L, 3L, 5L, 6L, 8L, 9L, 10L))
  expect_identical(climateSplitRefPixels(1:10, flam, elfPixels = c(2L, 4L, 7L, 8L, 9L)), c(2L, 8L, 9L))
})

test_that("a year's terms for the pixels of the fit equal those of a prediction table with extra rows, once filtered", {
  flam <- terra::rast(nrows = 1, ncols = 6, vals = c(1, 1, 1, 1, 0, 1))
  clim <- c(10, 20, 30, 40, 1000, 50)                      # cell 5 cannot burn and has an extreme value
  burnable <- c(1L, 2L, 3L, 4L, 6L)                        # what the fit has (its data year's landcoverDT)
  fit <- data.table::data.table(pixelID = 1:2, CMD = clim[1:2])
  splitSpreadClimate(fit, "CMD", data.table::data.table(CMD = clim[burnable]))
  tab <- data.table::data.table(pixelID = 1:6, CMD = clim)  # a prediction table that holds the extra row
  keep <- climateSplitRefPixels(tab$pixelID, flam)
  expect_identical(keep, burnable)
  pred <- data.table::copy(tab)
  splitSpreadClimate(pred, "CMD", tab[tab$pixelID %in% keep, "CMD", with = FALSE])
  expect_identical(pred$CMD_yr[1:2], fit$CMD_yr)
  expect_identical(pred$CMD_dev[1:2], fit$CMD_dev)
  expect_identical(fit$CMD_yr, rep(30, 2))
})

test_that("elfCellLabels() gives a categorical raster's labels and a plain raster's values", {
  r <- terra::rast(nrows = 1, ncols = 4, vals = c(2, 1, 2, NA))
  expect_identical(elfCellLabels(r), c("2", "1", "2", NA))
  levels(r) <- data.frame(id = 1:2, ELFind = c("A", "B"))
  expect_identical(elfCellLabels(r), c("B", "A", "B", NA))
})
