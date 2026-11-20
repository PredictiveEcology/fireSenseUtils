test_that("yearClimateColumn accepts the name with or without underscores, and stops on a missing one", {
  expect_identical(yearClimateColumn("CMD_sm", c("youngAge", "CMDsm")), "CMDsm")
  expect_identical(yearClimateColumn("CMDsm", c("youngAge", "CMD_sm")), "CMD_sm")
  expect_error(yearClimateColumn("CMDsm", c("youngAge", "x")), "not among the covariates")
})

test_that("yearClimateMeans is the weighted mean per year", {
  m <- yearClimateMeans(c(1, 3, 10, 20), year = c(2, 2, 1, 1), weight = c(1, 3, 1, 0))
  expect_equal(m$year, c(1, 2))
  expect_equal(m$mean, c(10, 2.5))
  expect_equal(yearClimateMeans(c(1, 3), year = c(1, 1))$mean, 2)
})

test_that("yearClimateMultiplier is exp(beta z), and 1 without a term", {
  yc <- list(variable = "CMDsm", center = 5, scale = 2, beta = 0.5)
  expect_equal(yearClimateMultiplier(yc, c(3, 5, 7)), exp(0.5 * c(-1, 0, 1)))
  expect_equal(yearClimateMultiplier(NULL, c(3, 7)), c(1, 1))
  expect_equal(yearClimateMultiplier(modifyList(yc, list(beta = 0)), 9), 1)
})
