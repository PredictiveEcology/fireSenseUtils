## Tests for multiplier() in makeBufferedFires.R

test_that("multiplier: respects minSize floor for small fires", {
  expect_equal(multiplier(1, minSize = 1000, baseMultiplier = 5), 1000)
  expect_equal(multiplier(10, minSize = 1000, baseMultiplier = 5), 1000)
})

test_that("multiplier: uses baseMultiplier when 14 - log(size) drops below it", {
  ## For large size, 14 - log(size) becomes small; baseMultiplier wins via pmax
  size <- 1e6  ## log(1e6) ~ 13.8; 14 - 13.8 = 0.2 < 5
  expect_equal(multiplier(size, minSize = 1, baseMultiplier = 5),
               round(5 * size, 0))
})

test_that("multiplier: uses (14 - log(size)) when above baseMultiplier", {
  ## At size = 100, log(100) ~ 4.6; 14 - 4.6 = 9.4 > 5 (baseMultiplier)
  size <- 100
  expected <- round(pmax(5, 14 - log(size)) * size, 0)
  expect_equal(multiplier(size, minSize = 1, baseMultiplier = 5), expected)
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
