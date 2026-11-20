test_that("yearClimateColumn accepts the name with or without underscores, and stops on a missing one", {
  expect_identical(yearClimateColumn("CMD_sm", c("youngAge", "CMDsm")), "CMDsm")
  expect_identical(yearClimateColumn("CMDsm", c("youngAge", "CMD_sm")), "CMD_sm")
  expect_error(yearClimateColumn("CMDsm", c("youngAge", "x")), "not among the covariates")
})

test_that("yearClimateMeans is the unweighted mean per year over complete rows with core", {
  d <- data.frame(year = c(2, 2, 2, 1, 1, 1), CMDsm = c(1, 3, 100, 10, 20, 30), other = c(1, 1, 1, NA, 1, 1),
                  coreFraction = c(1, 0.2, 0, 1, 1, 1))
  m <- yearClimateMeans(d, "CMD_sm")
  expect_equal(m$year, c(1, 2))
  expect_equal(m$mean, c(25, 2))          # year 1 drops the NA row; year 2 drops coreFraction 0; no weighting
  expect_equal(yearClimateMeans(data.frame(CMDsm = c(1, 3)), "CMDsm")$mean, 2)   # no year column: one year
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
