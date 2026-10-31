test_that("hindcastSetup runs 10 reps from the year, and no vegetation modules", {
  hc <- hindcastSetup(1991)
  expect_identical(hc$times, list(start = 1991, end = 1991 + 9.5 * 0.001))
  expect_false(any(grepl("Biomass|summary", hc$modules)))
  expect_identical(hc$params$climateYear$samplingRange, 1991)
})

test_that("hindcastSetup sets every module's own interval parameter to step", {
  p <- hindcastSetup(2001, step = 0.01)$params
  expect_identical(p$fireSense_dataPrepPredict$fireTimeStep, 0.01)
  for (m in c("fireSense_ignitionPredict", "fireSense_spreadPredict", "fireSense_burn"))
    expect_identical(p[[m]]$.runInterval, 0.01)
  expect_identical(p$fireSense_summary$hindcastInterval, 0.01)
  expect_identical(p$fireSense_summary$hindcastYear, 2001)
  expect_true(all(vapply(p[-1], function(x) is.na(x$.plots), logical(1))))
})

test_that("allSets runs 20 reps with a set per replicate; ignitions reach fireSense_burn only when given", {
  hc <- hindcastSetup(1986, allSets = TRUE, ignitions = "observedEscapes")
  expect_equal(hc$times$end, 1986 + 19.5 * 0.001)
  expect_true(hc$params$fireSense_spreadPredict$paramSetByReplicate)
  expect_identical(hc$params$fireSense_burn$ignitionSource, "observedEscapes")
  plain <- hindcastSetup(1986)$params
  expect_null(plain$fireSense_spreadPredict$paramSetByReplicate)
  expect_null(plain$fireSense_burn$ignitionSource)
})
