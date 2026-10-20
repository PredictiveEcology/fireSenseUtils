## The year effect on expected counts: yearEffectMultiplier() and its estimator yearEffectSD()

simTotals <- function(sigma, mu, seed) {
  set.seed(seed)
  list(observed = stats::rpois(length(mu), mu * yearEffectMultiplier(sigma, stats::rnorm(length(mu)))),
       predicted = mu)
}

test_that("the multiplier is exp(sd * z - sd^2 / 2), exactly 1 at sd 0, and NA sd is 0", {
  expect_equal(yearEffectMultiplier(0.5, c(-2, 0, 2)), exp(0.5 * c(-2, 0, 2) - 0.125))
  expect_identical(yearEffectMultiplier(0, c(-2, 0, 2)), c(1, 1, 1))
  expect_identical(yearEffectMultiplier(NA, 1.7), 1)
  ## one z, an sd per pixel (several ELFs sharing a year)
  expect_equal(yearEffectMultiplier(c(0, 0.3, NA), 1), c(1, exp(0.3 - 0.045), 1))
})

test_that("one z moves every mean the same way: more for z > 0, less for z < 0", {
  expect_true(all(yearEffectMultiplier(c(0.2, 0.5, 0.9), 1.5) > 1))
  expect_true(all(yearEffectMultiplier(c(0.2, 0.5, 0.9), -1.5) < 1))
})

test_that("the mean is preserved over many draws", {
  set.seed(7)
  z <- stats::rnorm(2e5)
  for (sd in c(0.3, 0.7, 1.2))
    expect_equal(mean(yearEffectMultiplier(sd, z)), 1, tolerance = 0.02)
})

test_that("the estimator recovers a known sigma on simulated data", {
  mu <- rep(c(30, 80, 150, 300), length.out = 4000)   # yearly expected totals, as in an ELF
  for (sigma in c(0.3, 0.5, 0.9))
    expect_equal(yearEffectSD(simTotals(sigma, mu, seed = 1)$observed, mu), sigma, tolerance = 0.05)
})

test_that("with no year effect the estimate is about 0: the Poisson variance is removed", {
  mu <- rep(c(30, 80, 150, 300), length.out = 4000)
  s <- simTotals(0, mu, seed = 2)
  expect_lt(yearEffectSD(s$observed, s$predicted), 0.03)
  ## the plain sd of log(O / mu) would count the Poisson noise: about 0.1 here
  expect_gt(stats::sd(log(s$observed / mu)), 0.08)
})

test_that("small expected totals with many zero years work (escapes)", {
  mu <- rep(1.5, 4000)
  s <- simTotals(0.8, mu, seed = 3)
  expect_gt(sum(s$observed == 0), 1000)
  expect_equal(yearEffectSD(s$observed, s$predicted), 0.8, tolerance = 0.1)
})

test_that("a bias in the overall level is not counted as a year effect", {
  mu <- rep(c(30, 80, 150, 300), length.out = 4000)
  s <- simTotals(0.4, mu, seed = 4)
  expect_equal(yearEffectSD(s$observed, s$predicted / 2), yearEffectSD(s$observed, s$predicted))
})

test_that("less spread than Poisson is floored at 0; degenerate input gives 0; NA years are dropped", {
  expect_identical(yearEffectSD(c(100, 100, 100), c(100, 100, 100)), 0)
  expect_identical(yearEffectSD(5, 5), 0)                       # one year
  expect_identical(yearEffectSD(c(0, 0), c(1, 2)), 0)           # nothing observed
  expect_identical(yearEffectSD(c(1, 2), c(0, 0)), 0)           # nothing predicted
  expect_equal(yearEffectSD(c(10, NA, 40), c(20, 5, 20)), yearEffectSD(c(10, 40), c(20, 20)))
})
