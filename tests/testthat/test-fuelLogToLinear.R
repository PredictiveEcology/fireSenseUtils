## Fuel biomass reaches the spread model logged by logMinB(). fuelLogToLinear() undoes that for
## fireSense_SpreadFit and fireSense_SpreadPredict alike, and fuelLinearRange is the covMinMax that
## divides it by 1e4. These tests pin the values both modules depend on.

test_that("fuelLogToLinear() inverts logMinB() above the floor", {
  b <- c(40, 300, 7803, 22601, 76308)
  expect_equal(fuelLogToLinear(logMinB(b)), b)
})

test_that("biomass on or below the logMinB() floor becomes exactly 0", {
  floorB <- exp(log(100) - 1) # 36.79
  expect_identical(fuelLogToLinear(logMinB(c(0, 1, 36, floorB))), c(0, 0, 0, 0))
  ## the floor as it is stored in the x1000 integer tables: 3605, not 3605.170
  expect_identical(fuelLogToLinear(3605 / 1000), 0)
  ## just above the tolerance it is biomass again
  expect_equal(fuelLogToLinear(log(38)), 38)
})

test_that("it reproduces the transform the model-selection experiment used", {
  ## runArm.R: m <- col / 1000; b <- exp(m); b[m <= logMinB + 1e-3] <- 0
  x1000 <- c(3605L, 3605L, 5704L, 8962L, 10026L)
  m <- x1000 / 1000; b <- exp(m); b[m <= (log(100) - 1) + 1e-3] <- 0
  expect_identical(fuelLogToLinear(x1000 / 1000), b)
})

test_that("fuelLinearRange makes rescaleKnown2() a division by 1e4, without clamping", {
  b <- c(0, 5000, 10000, 22601)
  r <- rescaleKnown2(b, 0, 1, fuelLinearRange[1], fuelLinearRange[2])
  expect_equal(r, b / 1e4)
  expect_gt(max(r), 1)
})

test_that("isLinearFuelRange() tells a linear fit from a log fit", {
  expect_true(isLinearFuelRange(c(0, 1e4)))
  expect_true(isLinearFuelRange(c(0L, 10000L)))
  expect_false(isLinearFuelRange(c(3.605, 10.06)))   # a log-scale fit
  expect_false(isLinearFuelRange(c(0, 1)))           # an indicator
  expect_false(isLinearFuelRange(c(0, 22601)))       # a data maximum is not the fixed divisor
  expect_false(isLinearFuelRange(NULL))
})
