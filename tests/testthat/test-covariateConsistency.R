## Covariate consistency between the spread fit and prediction: fixed climate ranges (provisional),
## one climate precision, 0/1 indicator ranges, empty fuel columns, one spread-probability ceiling.

library(data.table)

test_that("climateCovRanges is one table of c(min, max) per climate variable, provisionally 0-100", {
  expect_type(climateCovRanges, "list")
  expect_setequal(names(climateCovRanges), c("CMD", "CMDsm", "CMDsp", "cumMDC"))
  for (r in climateCovRanges) {
    expect_length(r, 2L)
    expect_true(r[2] > r[1])
  }
  expect_true(all(vapply(climateCovRanges, identical, logical(1), c(0, 100))))
})

test_that("a climate value is stored as the same integer by the fit and the predict path", {
  skip_if_not_installed("terra")
  r <- terra::rast(nrows = 2, ncols = 2, vals = c(104.6, 0.4, 99.9996, 297.123))
  names(r) <- "year2001"
  fit <- climateRasterToDataTable(list(CMD = r))
  fitX1000 <- toX1000(list(data.frame(pixelID = fit$pixelID, CMD = fit$CMD)))[[1]]
  predX1000 <- toX1000(list(data.frame(pixelID = 1:4, CMD = c(104.6, 0.4, 99.9996, 297.123))))[[1]]
  expect_identical(fitX1000$CMD, predX1000$CMD)
  expect_identical(fitX1000$CMD[1], 104600L)
})

test_that("spreadIndicatorCols finds youngAge, nfLCC_* and treedWetland, and no fuel or climate", {
  nms <- c("pixelID", "CMD", "youngAge", "nfLCC_40_50", "nfLCC_100", "treedWetland", "treedWetland_agb",
           "BlkSprc")
  expect_setequal(spreadIndicatorCols(nms), c("youngAge", "nfLCC_40_50", "nfLCC_100", "treedWetland"))
  expect_identical(spreadIndicatorCols(c("CMD", "BlkSprc")), character(0))
  ## the three groups are one definition: youngAge excludes the other two and fuel
  expect_setequal(unlist(youngAgeExclusiveCols(nms, fuelCols = "BlkSprc")),
                  c("BlkSprc", "nfLCC_40_50", "nfLCC_100", "treedWetland"))
})

test_that("spreadIndicatorRanges gives every indicator the fixed range c(0, 1)", {
  out <- spreadIndicatorRanges(c("CMD", "youngAge", "nfLCC_100", "BlkSprc"))
  expect_identical(out, list(nfLCC_100 = c(0, 1), youngAge = c(0, 1)))
})

test_that("a rescale with the fixed indicator range does not give Inf or NaN for a constant column", {
  x <- rep(0, 5) # a constant youngAge column: no young pixels
  rg <- spreadIndicatorRanges("youngAge")$youngAge
  expect_true(all(is.finite(rescaleKnown2(x * 1000, 0, 1000, rg[1] * 1000, rg[2] * 1000))))
  ## with the data's range it was 0 / 0
  expect_true(all(is.nan(rescaleKnown2(x * 1000, 0, 1000, 0, 0))))
})

test_that("emptySpreadCovariates finds all-zero indicators and fuel columns all on the logMinB floor", {
  dt <- data.table(
    nfLCC_100 = c(0, 0, 0),
    nfLCC_50 = c(0, 1, 0),
    youngAge = c(0, 0, 0),
    BlkSprc = logMinB(c(0, 0, 0)),          # no biomass anywhere: 3.6 everywhere
    WhtSprc = logMinB(c(0, 500, 12000)),
    Pine = logMinB(c(NA, 0, 0))
  )
  expect_setequal(emptySpreadCovariates(dt, names(dt)), c("nfLCC_100", "youngAge", "BlkSprc", "Pine"))
  ## an indicator with some 1s is not "on the floor": 0/1 values are below it, whatever the column is called
  expect_identical(emptySpreadCovariates(data.table(nf = c(0, 1, 0), nfLCC_50 = c(0, 1, 0)), c("nf", "nfLCC_50")),
                   character(0))
  ## a fuel column with one pixel off the floor stays
  dt[2, BlkSprc := logMinB(80)]
  expect_false("BlkSprc" %in% emptySpreadCovariates(dt, names(dt)))
})

test_that("the spread-probability ceiling and floor are defined once and are the defaults", {
  expect_identical(spreadProbCeiling, 0.276)
  expect_identical(spreadProbFloor, 0.13)
  for (f in list(.objfunSpreadFit, spreadFitValidationData, spreadProbGates)) {
    expect_identical(eval(formals(f)$maxFireSpread), spreadProbCeiling)
    expect_identical(eval(formals(f)$lowerSpreadProb), spreadProbFloor)
  }
})
