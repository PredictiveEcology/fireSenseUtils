## The escape threshold, the annual-area term and the area-weighted size-distribution term of
## .objfunSpreadFit() (all off by default).

test_that("escapePixels() turns hectares into pixels from the landscape's resolution", {
  r240 <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 2400, ymin = 0, ymax = 2400, crs = "EPSG:3005")
  expect_equal(prod(terra::res(r240)) / 1e4, 5.76)
  expect_identical(fireSenseUtils:::escapePixels(50, r240), 9L)      # 50 / 5.76 = 8.7 -> 9
  expect_identical(fireSenseUtils:::escapePixels(NULL, r240), 2L)    # the original "more than 1 pixel"
  expect_identical(fireSenseUtils:::escapePixels(1, r240), 2L)       # never below 2
  r250 <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 2500, ymin = 0, ymax = 2500, crs = "EPSG:3005")
  expect_identical(fireSenseUtils:::escapePixels(50, r250), 8L)      # exactly 8 x 6.25 ha
})

test_that("areaWeightedCvM() is 0 for identical samples and matches a hand calculation", {
  awc <- fireSenseUtils:::areaWeightedCvM
  expect_equal(awc(c(1, 3, 10), c(1, 3, 10)), 0)
  expect_equal(awc(rep(c(1, 3, 10), 5), c(1, 3, 10)), 0)  # same shares of area
  ## obs 1 and 3: area shares at x = 1 and 3 are 0.25 and 1; sim all 1s: 1 and 1
  ## 2 * (0.25 * 0.75^2 + 0.75 * 0^2) = 0.28125
  expect_equal(awc(c(1, 1, 1, 1), c(1, 3)), 0.28125)
  ## too many large simulated fires scores worse than a close match
  obs <- c(2, 5, 20, 150, 900)
  expect_gt(awc(c(obs, 20000), obs), awc(c(obs, 1000), obs))
})

test_that("yearAreaNLL() equals the per-fire t likelihood on annual totals and responds to a shift", {
  set.seed(1)
  sims <- list(year2001 = round(rlnorm(50, log(400), 0.4)), year2002 = round(rlnorm(50, log(60), 0.4)))
  obs <- c(year2001 = 400, year2002 = 60)
  nll <- fireSenseUtils:::yearAreaNLL(obs, sims, sizeLik = "t", sizeLikDf = 5)
  byHand <- -log(fireSenseUtils:::sizeLikT(400, sims$year2001, 5)) - log(fireSenseUtils:::sizeLikT(60, sims$year2002, 5))
  expect_equal(nll, byHand)
  shifted <- fireSenseUtils:::yearAreaNLL(c(year2001 = 1600, year2002 = 60), sims, sizeLik = "t", sizeLikDf = 5)
  expect_gt(shifted, nll)
  ## years without simulated totals are skipped, not failed
  expect_equal(fireSenseUtils:::yearAreaNLL(obs, list(year2001 = sims$year2001, year2002 = NULL), "t", 5),
               -log(fireSenseUtils:::sizeLikT(400, sims$year2001, 5)))
  expect_true(is.finite(fireSenseUtils:::yearAreaNLL(obs, sims, sizeLik = "kde")))
})

test_that("the new arguments are off by default", {
  f <- formals(fireSenseUtils:::.objfunSpreadFit)
  expect_null(f$escapeSizeHa)
  expect_identical(f$yearAreaWeight, 0)
  expect_identical(f$areaDistWeight, 0)
  expect_false(f$returnTerms)
})
