test_that("fireSensePhases are all, inputs and fit", {
  expect_identical(fireSensePhases, c("all", "inputs", "fit"))
})

test_that("phase 'all' has no barriers, and is the default", {
  expect_null(phaseEvents("all"))
  expect_null(phaseEvents())
})

test_that("phase 'inputs' stops before the fit and the climate rasters", {
  ## fireSense_spreadFit has no "estimateThreshold" event any more; a barrier on an event that is never
  ## scheduled does not stop anything, so the job would go on to fit
  expect_identical(phaseEvents("inputs"),
                   list(.stopBefore = list(fireSense_spreadFit = "run",
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
