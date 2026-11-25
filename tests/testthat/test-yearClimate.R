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

test_that("beta2 adds a hinge on the dry half; a term without beta2 is unchanged", {
  yc <- list(variable = "CMDsm", center = 5, scale = 2, beta = 0.5, norm = 2)
  z <- c(-1, 0, 1, 2)
  expect_identical(yearClimateMultiplier(yc, 5 + 2 * z), exp(0.5 * (5 + 2 * z - 5) / 2) / 2)  # the arithmetic of old
  expect_identical(yearClimateMultiplier(modifyList(yc, list(beta2 = 0)), 5 + 2 * z), yearClimateMultiplier(yc, 5 + 2 * z))
  expect_equal(yearClimateMultiplier(modifyList(yc, list(beta2 = 0.3)), 5 + 2 * z), exp(0.5 * z + 0.3 * pmax(0, z)) / 2)
})
