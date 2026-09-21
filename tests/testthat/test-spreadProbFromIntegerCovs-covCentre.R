test_that("covCentre is subtracted after rescaling, and only from the named columns", {
  mk <- function() data.table::data.table(pixelID = 1:3, a = c(0L, 500L, 1000L), b = c(0L, 200L, 400L))
  run <- function(covCentre, doAssertions = FALSE) {
    spreadProbFromIntegerCovs(
      shortAnnDTx1000 = mk(), yr = 2000, covMinMax = NULL, mutuallyExclusive = NULL,
      colsToUse = c("a", "b"), doAssertions = doAssertions, logisticPars = c(0.2, 1, 1),
      covPars = c(1, 1), maxFireSpread = 0.28, lowerSpreadProb = 0.1, covCentre = covCentre)
  }
  plain <- run(NULL)
  expect_equal(plain$a, c(0, 0.5, 1))

  centred <- run(list(a = 0.5))
  expect_equal(centred$a, c(-0.5, 0, 0.5))
  expect_equal(centred$b, plain$b)

  ## centred values are negative by design; the non-negative range assertion must not see them
  expect_no_error(run(list(a = 0.5), doAssertions = TRUE))
})
