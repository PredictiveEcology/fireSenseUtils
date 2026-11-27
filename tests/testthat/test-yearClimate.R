test_that("yearClimateColumn accepts the name with or without underscores, and stops on a missing one", {
  expect_identical(yearClimateColumn("CMD_sm", c("youngAge", "CMDsm")), "CMDsm")
  expect_identical(yearClimateColumn("CMDsm", c("youngAge", "CMD_sm")), "CMD_sm")
  expect_error(yearClimateColumn("CMDsm", c("youngAge", "x")), "not among the covariates")
})

test_that("yearClimateMeans is the unweighted mean per year over the rows with the variable and with core", {
  d <- data.frame(year = c(2, 2, 2, 1, 1, 1), CMDsm = c(1, 3, 100, 10, 20, 30), other = c(1, 1, 1, NA, 1, 1),
                  coreFraction = c(1, 0.2, 0, 1, 1, 1))
  m <- yearClimateMeans(d, "CMD_sm")
  expect_equal(m$year, c(1, 2))
  expect_equal(m$mean, c(20, 2))          # year 2 drops coreFraction 0; no weighting; NA in `other` does not drop a row
  expect_equal(yearClimateMeans(data.frame(CMDsm = c(1, 3)), "CMDsm")$mean, 2)   # no year column: one year
})

test_that("only the variable and coreFraction decide which rows count, not other columns", {
  d <- data.frame(year = 1, CMDsm = c(1, 3, NA, 5), coreFraction = c(1, 1, 1, NA), other = c(NA, 1, 1, 1))
  expect_equal(yearClimateMeans(d, "CMDsm")$mean, 2)
  expect_equal(yearClimateMeans(d[, c("year", "CMDsm", "other")], "CMDsm")$mean, 3)    # no coreFraction: every row with the variable
})

test_that("a year with no rows left is dropped, not NA", {
  d <- data.frame(year = c(1, 1, 2), CMDsm = c(1, 3, NA))
  expect_equal(yearClimateMeans(d, "CMDsm"), data.frame(year = 1, mean = 2))
  expect_equal(nrow(yearClimateMeans(data.frame(year = 1, CMDsm = NA_real_), "CMDsm")), 0L)
})

test_that("yearClimateMultiplier is exp(beta z) / norm, and 1 without a term", {
  yc <- list(variable = "CMDsm", center = 5, scale = 2, beta = 0.5, norm = 2)
  expect_equal(yearClimateMultiplier(yc, c(3, 5, 7)), exp(0.5 * c(-1, 0, 1)) / 2)
  expect_equal(yearClimateMultiplier(yc[names(yc) != "norm"], 5), 1)   # a fit stored without norm
  expect_equal(yearClimateMultiplier(NULL, c(3, 7)), c(1, 1))
  expect_equal(yearClimateMultiplier(modifyList(yc, list(beta = 0, norm = 1)), 9), 1)
})

test_that("the logistic shape is an S-curve through z = 0 that saturates; a term without shape is unchanged", {
  yc <- list(variable = "CMDsm", center = 5, scale = 2, beta = 0.5, norm = 2)
  z <- c(-1, 0, 1, 2)
  expect_identical(yearClimateMultiplier(yc, 5 + 2 * z), exp(0.5 * standardiseYearClimate(5 + 2 * z, 5, 2)) / 2)
  expect_equal(yearClimateMultiplier(yc, 5 + 2 * z), exp(0.5 * (5 + 2 * z - 5) / 2) / 2, tolerance = 1e-14)  # old arithmetic, to rounding
  lg <- list(variable = "CMDsm", center = 5, scale = 2, shape = "logistic", k = 3, ratio = 4, norm = 1)
  expect_equal(yearClimateMultiplier(lg, 5 + 2 * z), 1 + 3 * plogis(3 * z))
  expect_identical(yearClimateMultiplier(lg, 5 + 2 * z), logisticYearCurve(z, 3, 4))   # the one formula
  expect_equal(yearClimateMultiplier(lg, 5), 2.5)                       # midpoint: (floor + upper) / 2
  expect_equal(yearClimateMultiplier(lg, c(-1e3, 1e3)), c(1, 4))        # floor 1, upper `ratio`
  expect_true(all(diff(yearClimateMultiplier(lg, seq(-5, 15, 0.5))) >= 0))
  expect_equal(yearClimateMultiplier(modifyList(lg, list(k = 0)), c(1, 9)), c(2.5, 2.5))   # k = 0: flat
  expect_equal(yearClimateMultiplier(modifyList(lg, list(norm = 2)), 5), 1.25)
})

test_that("logisticYearCurve is an S-curve from 1 to ratio through its midpoint at z = 0", {
  expect_equal(logisticYearCurve(c(-1e3, 0, 1e3), k = 2, ratio = 4), c(1, 2.5, 4))
  expect_equal(logisticYearCurve(0.7, k = 0, ratio = 4), 2.5)
  expect_equal(logisticYearCurve(1, k = 2, ratio = 1), 1)
  expect_true(all(diff(logisticYearCurve(seq(-4, 4, 0.5), 1.5, 3)) > 0))
})

test_that("the capped shape is log-linear up to zCap and flat above; with zCap above the data it is the log-linear term", {
  cp <- list(variable = "CMDsm", center = 5, scale = 2, shape = "capped", beta = 0.5, zCap = 1, norm = 2)
  z <- c(-2, 0, 1, 2, 50)
  expect_identical(yearClimateMultiplier(cp, 5 + 2 * z), cappedYearCurve(z, 0.5, 1) / 2)   # the one formula
  expect_equal(yearClimateMultiplier(cp, 5 + 2 * z), exp(0.5 * c(-2, 0, 1, 1, 1)) / 2)
  ll <- list(variable = "CMDsm", center = 5, scale = 2, beta = 0.5, norm = 2)
  zz <- c(-2, 0, 0.9)
  expect_equal(yearClimateMultiplier(cp, 5 + 2 * zz), yearClimateMultiplier(ll, 5 + 2 * zz))
  expect_equal(cappedYearCurve(c(-1, 3), beta = 0, zCap = 0), c(1, 1))
})
