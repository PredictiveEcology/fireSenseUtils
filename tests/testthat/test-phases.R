test_that("fireSensePhases are all, inputs and fit", {
  expect_identical(fireSensePhases, c("all", "inputs", "fit"))
})

test_that("phase 'all' has no barriers, and is the default", {
  expect_null(phaseEvents("all"))
  expect_null(phaseEvents())
})

test_that("phase 'inputs' stops before the threshold calibration and the climate rasters", {
  expect_identical(phaseEvents("inputs"),
                   list(.stopBefore = list(fireSense_spreadFit = "estimateThreshold",
                                           fireSense_dataPrepPredict = "getClimateRasters")))
})

test_that("phase 'fit' stops after the fit, or after cross-validation when it validates", {
  fit <- phaseEvents("fit")
  expect_identical(fit$.stopAfter, list(fireSense_spreadFit = "run"))
  expect_identical(fit$.stopBefore, list(fireSense_dataPrepPredict = "getClimateRasters"))
  expect_identical(phaseEvents("fit", validates = TRUE)$.stopAfter,
                   list(fireSense_spreadFit = "crossValidate"))
})

test_that("phaseEvents() matches the phase partially and rejects unknown ones", {
  expect_identical(phaseEvents("inp"), phaseEvents("inputs"))
  expect_error(phaseEvents("predict"), "should be one of")
})
