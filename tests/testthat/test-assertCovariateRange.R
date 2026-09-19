## Regression tests for assertCovariateRange().
##
## The previous assertion stopped when any covariate exceeded 1. That blocked fuel
## covariates expressed as biomass/1e4, which legitimately reach ~2.3. Nothing in the
## objective requires an upper bound of 1 (see the function's documentation), so only
## non-negativity and finiteness are enforced now.

library(data.table)

test_that("assertCovariateRange: passes covariates within [0, 1]", {
  dt <- data.table(a = c(0, 0.5, 1), b = c(0.2, 0.2, 0.2))
  expect_true(assertCovariateRange(dt, c("a", "b")))
})

test_that("assertCovariateRange: allows values above 1 (biomass/1e4 fuel covariates)", {
  dt <- data.table(fuel = c(0, 1.5, 2.26), fuelSq = c(0, 2.25, 5.11))
  expect_true(assertCovariateRange(dt, c("fuel", "fuelSq")))
})

test_that("assertCovariateRange: still rejects negative covariates", {
  dt <- data.table(a = c(0, 0.5), b = c(-0.566, 0.2))
  expect_error(assertCovariateRange(dt, c("a", "b")), "non-negative and finite")
})

test_that("assertCovariateRange: tolerates negatives smaller than the rounding tolerance", {
  # the stored log floor rescales to -2.7e-05 rather than exactly 0
  dt <- data.table(a = c(-2.671e-05, 0.5))
  expect_true(assertCovariateRange(dt, "a"))
})

test_that("assertCovariateRange: rejects non-finite covariates", {
  expect_error(assertCovariateRange(data.table(a = c(0, NA_real_)), "a"), "non-negative and finite")
  expect_error(assertCovariateRange(data.table(a = c(0, Inf)), "a"), "non-negative and finite")
})

test_that("assertCovariateRange: names only the offending columns", {
  dt <- data.table(ok = c(0, 2.3), bad = c(0, -1))
  expect_error(assertCovariateRange(dt, c("ok", "bad")), "bad")
})
