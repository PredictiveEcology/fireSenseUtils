## adStatistic(), the two-sample Anderson-Darling statistic. The spread objective no longer has the
## pooled Anderson-Darling term: it ranked fires cut at the buffers among the observed ones, so it
## depended on the buffers. The statistic stays for diagnostics (scoreFireSizes()).

test_that("adStatistic() is the statistic kSamples::ad.test() reports, ties included", {
  skip_if_not_installed("kSamples")
  set.seed(2)
  for (shift in c(0, 1)) {
    obs <- round(rlnorm(300, 3, 2))
    sim <- round(rlnorm(300 * 25, 3 + shift, 2))
    expect_equal(adStatistic(sim, obs), kSamples::ad.test(sim, obs)$ad[1, 1], tolerance = 1e-4)
  }
})

test_that("adStatistic() errors on an empty sample, so the objective returns its fail value, not NaN", {
  ## every year can bail before simulating; DEoptim stops on a NaN objective
  expect_error(adStatistic(NULL, c(5, 10)))
  expect_error(adStatistic(c(5, 10), numeric(0)))
})

test_that("adStatistic() does not overflow on the sample sizes of a real fit", {
  ## 1800 fires x 50 replicates: N * M passes .Machine$integer.max
  set.seed(3)
  sim <- round(rlnorm(90000, 3, 2)); obs <- round(rlnorm(1800, 3, 2))
  expect_no_warning(ad <- adStatistic(sim, obs))
  expect_true(is.finite(ad))
})
