## Tests for multiplier() in makeBufferedFires.R

test_that("multiplier: respects minSize floor for small fires", {
  expect_equal(multiplier(1, minSize = 1000, baseMultiplier = 5), 1000)
  expect_equal(multiplier(10, minSize = 1000, baseMultiplier = 5), 1000)
})

test_that("multiplier: uses baseMultiplier when logIntercept - logSlope * log(size) drops below it", {
  ## For large size, 21 - 2 * log(size) becomes small; baseMultiplier wins via pmax
  size <- 1e6  ## log(1e6) ~ 13.8; 21 - 27.6 < 5
  expect_equal(multiplier(size, minSize = 1, baseMultiplier = 5),
               round(5 * size, 0))
})

test_that("multiplier: uses (logIntercept - logSlope * log(size)) when above baseMultiplier", {
  ## At size = 100, log(100) ~ 4.6; 21 - 2 * 4.6 = 11.8 > 7.5 (baseMultiplier); 14 - 4.6 = 9.4 > 5
  size <- 100
  expect_equal(multiplier(size, minSize = 1), round(pmax(7.5, 21 - 2 * log(size)) * size, 0))
  expect_equal(multiplier(size, minSize = 1, baseMultiplier = 5, logIntercept = 14, logSlope = 1),
               round(pmax(5, 14 - log(size)) * size, 0))
})

test_that("multiplier: is vectorised over size", {
  out <- multiplier(c(1, 100, 1e6), minSize = 1, baseMultiplier = 5)
  expect_length(out, 3)
  expect_true(all(out > 0))
})

test_that("multiplier: returns rounded integers", {
  out <- multiplier(c(1, 100, 1e6), minSize = 1, baseMultiplier = 5)
  expect_equal(out, round(out, 0))
})

test_that("multiplier: minSize floor still applies element-wise", {
  out <- multiplier(c(1, 1e6), minSize = 1000, baseMultiplier = 5)
  expect_gte(out[1], 1000)
  expect_gte(out[2], 1000)
})

test_that("multiplier: the defaults keep the spread objective buffer-independent (see test-bufferIndependentObjective.R)", {
  expect_equal(multiplier(c(9, 100, 3000, 18000)), c(500, 1179, 22500, 135000))
  ## with the module's factor 10: 555, 118, 75 and 75 times the fire's size
  expect_equal(10 * multiplier(c(9, 100, 3000, 18000)) / c(9, 100, 3000, 18000),
               c(5000 / 9, 117.9, 75, 75))
})
