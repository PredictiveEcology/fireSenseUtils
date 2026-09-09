## objFunInner() returns `ret` unconditionally, but only built it inside
## `if (isTRUE(doFitting))`. doFitting is `any(c(doSNLL_FSTest, doMADTest,
## doADTest))`, so with no test selected -- which is exactly how
## fireSense_SpreadFit's `debug` mode calls the chain, passing tests = "" -- the
## branch was skipped and the return failed with "object 'ret' not found". The
## comment at the far end of that branch ("Object ret doesn't exist") shows this
## had been noticed before.
##
## The work before the branch is stubbed out here: what is under test is the
## contract that the function returns a list whatever the flags say.

test_that("objFunInner returns a list when no test was asked for", {
  skip_if_not_installed("data.table")
  local_mocked_bindings(
    paramsSeparate = function(...) list(logisticPars = c(1, 1), covPars = 1),
    spreadProbFromIntegerCovs = function(...) data.table::data.table(pixelID = 1L, cov = 0),
    logisticAll = function(...) 0.2
  )
  out <- objFunInner(
    yr = "2004", annDTx1000 = NULL, par = 1, parsModel = NULL,
    annualFires = NULL, nonAnnualDTx1000 = NULL, annualFireBufferedDT = NULL,
    indexNonAnnual = NULL, colsToUse = "cov", covMinMax = NULL,
    mutuallyExclusive = NULL, doAssertions = FALSE, maxFireSpread = 0.28,
    lowerSpreadProb = 0.1, cells = integer(4), lanscape1stQuantileThresh = 0.1,
    weighted = TRUE, r = NULL, Nreps = 1L,
    doSNLL_FSTest = FALSE, doMADTest = FALSE, doADTest = FALSE,
    plot.it = FALSE, verbose = 0
  )
  expect_type(out, "list")
  expect_length(out, 0L)   ## nothing was asked for, so nothing is reported
})
