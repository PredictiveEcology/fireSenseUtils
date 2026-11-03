test_that("only an objective with the SNLL_FS term has an early stop", {
  expect_false(objectiveHasEarlyStop("crps"))
  expect_false(objectiveHasEarlyStop(c("adTest", "SNLL")))
  expect_true(objectiveHasEarlyStop(c("adTest", "SNLL_FS")))
  expect_true(objectiveHasEarlyStop("snll_fs"))
  expect_true(objectiveHasEarlyStop("Snll_Fs"))
})
