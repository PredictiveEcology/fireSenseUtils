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
