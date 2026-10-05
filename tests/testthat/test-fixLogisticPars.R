## hillSlope1 and inflectionPoint1 (the spread link's slope and Richards asymmetry exponent) are fixed
## at 1, not fitted by DEoptim. hillSlope1 is not identifiable: with the linear predictor
## x = covariates %*% beta it enters logistic3p()/logistic3pUpper() only as hillSlope1 * x, so scaling
## every coefficient by k and dividing hillSlope1 by k changes nothing. inflectionPoint1 is not a
## location but the exponent of the logistic, and in the 2026-10-04 fits it was bimodal and traded off
## against the coefficients. fixLogisticPars() reinserts both as the fixed 2nd and 3rd logistic
## parameters; .objfunSpreadFit() calls it on every `par` it receives (see test-yearSpreadSD.R for the
## integration test of that call).

test_that("fixLogisticPars() inserts hillSlope1 = 1 and inflectionPoint1 = 1 as elements 2 and 3, named or not", {
  expect_identical(fireSenseUtils:::fixLogisticPars(c(maxAsymptote = 0.27, cov = 2)),
                   c(maxAsymptote = 0.27, hillSlope1 = 1, inflectionPoint1 = 1, cov = 2))
  ## unnamed (as DEoptim's `par` arrives): still inserted after the 1st, by position
  expect_identical(fireSenseUtils:::fixLogisticPars(c(0.27, 4, 2)),
                   c(0.27, hillSlope1 = 1, inflectionPoint1 = 1, 4, 2))
  ## with upperTail1 (logistic3pUpper): maxAsymptote, hillSlope1, inflectionPoint1, upperTail1
  expect_identical(fireSenseUtils:::fixLogisticPars(c(maxAsymptote = 0.27, upperTail1 = -0.3, cov = 2)),
                   c(maxAsymptote = 0.27, hillSlope1 = 1, inflectionPoint1 = 1, upperTail1 = -0.3, cov = 2))
})

test_that("a stored full par (fitted inflectionPoint1) is not given a second one", {
  full <- c(maxAsymptote = 0.27, hillSlope1 = 1, inflectionPoint1 = 3, cov = 2)
  expect_identical(fireSenseUtils:::fixLogisticPars(full), full)
  ## the pre-hillSlope1-fix style of a fit that fitted inflectionPoint1 only: hillSlope1 goes in 2nd
  expect_identical(fireSenseUtils:::fixLogisticPars(c(maxAsymptote = 0.27, inflectionPoint1 = 3, cov = 2)),
                   c(maxAsymptote = 0.27, hillSlope1 = 1, inflectionPoint1 = 3, cov = 2))
})

test_that("the spread link on the short par equals the old link with both fixed values written out", {
  mat <- matrix(c(0.1, 0.4, 0.7, 1), ncol = 1, dimnames = list(NULL, "cov"))
  covPars <- c(cov = 2)
  for (named in c(TRUE, FALSE)) {
    parShort <- c(maxAsymptote = 0.27)                                              # DEoptim's par
    parOld <- c(maxAsymptote = 0.27, hillSlope1 = 1, inflectionPoint1 = 1)          # written out by hand
    if (!named) { parShort <- unname(parShort); parOld <- unname(parOld) }
    newPred <- fireSenseUtils::logisticAll(fireSenseUtils:::fixLogisticPars(parShort), mat, covPars,
                                           lowerSpreadProb = 0.13)
    oldPred <- fireSenseUtils::logisticAll(parOld, mat, covPars, lowerSpreadProb = 0.13)
    expect_identical(unname(newPred), unname(oldPred))
    ## and it is the plain logistic
    expect_equal(as.numeric(newPred), 0.13 + (0.27 - 0.13) * stats::plogis(as.numeric(mat %*% covPars)))
  }
  ## with upperTail1 (logistic3pUpper), named and unnamed (unnamed needs link =)
  for (named in c(TRUE, FALSE)) {
    parShort <- c(maxAsymptote = 0.27, upperTail1 = -0.3)
    parOld <- c(maxAsymptote = 0.27, hillSlope1 = 1, inflectionPoint1 = 1, upperTail1 = -0.3)
    if (!named) { parShort <- unname(parShort); parOld <- unname(parOld) }
    newPred <- fireSenseUtils::logisticAll(fireSenseUtils:::fixLogisticPars(parShort), mat, covPars,
                                           lowerSpreadProb = 0.13, link = "logistic3pUpper")
    oldPred <- fireSenseUtils::logisticAll(parOld, mat, covPars, lowerSpreadProb = 0.13,
                                           link = "logistic3pUpper")
    expect_identical(unname(newPred), unname(oldPred))
  }
  ## it is not ignoring the fixed values: a different inflectionPoint1 changes the prediction
  offPred <- fireSenseUtils::logisticAll(c(maxAsymptote = 0.27, hillSlope1 = 1, inflectionPoint1 = 0.4),
                                         mat, covPars, lowerSpreadProb = 0.13)
  expect_false(isTRUE(all.equal(unname(fireSenseUtils::logisticAll(
    fireSenseUtils:::fixLogisticPars(c(maxAsymptote = 0.27)), mat, covPars, lowerSpreadProb = 0.13)),
    unname(offPred))))
})
