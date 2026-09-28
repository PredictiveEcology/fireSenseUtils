## hillSlope1 (the spread link's slope) is fixed at 1, not fitted by DEoptim: with the linear
## predictor x = covariates %*% beta, hillSlope1 enters logistic3p()/logistic3pUpper() only as
## hillSlope1 * x, so scaling every covariate coefficient by k and dividing hillSlope1 by k leaves
## every prediction unchanged -- hillSlope1 is not identifiable. fixHillSlope1() reinserts it as the
## fixed 2nd logistic parameter; .objfunSpreadFit() calls it on every `par` it receives (see
## test-yearSpreadSD.R for the integration test of that call).

test_that("fixHillSlope1() inserts hillSlope1 = 1 as the 2nd element, named or not", {
  expect_identical(fireSenseUtils:::fixHillSlope1(c(maxAsymptote = 0.27, inflectionPoint1 = 4)),
                    c(maxAsymptote = 0.27, hillSlope1 = 1, inflectionPoint1 = 4))
  ## unnamed (as DEoptim's `par` arrives): still inserted 2nd, by position
  expect_identical(fireSenseUtils:::fixHillSlope1(c(0.27, 4, 2)),
                    c(0.27, hillSlope1 = 1, 4, 2))
  ## with upperTail1 (logistic3pUpper): maxAsymptote, hillSlope1, inflectionPoint1, upperTail1
  expect_identical(fireSenseUtils:::fixHillSlope1(c(maxAsymptote = 0.27, inflectionPoint1 = 4, upperTail1 = -0.3)),
                    c(maxAsymptote = 0.27, hillSlope1 = 1, inflectionPoint1 = 4, upperTail1 = -0.3))
})

test_that("the spread link on the new (short) par equals the old link with hillSlope1 = 1 by hand", {
  ## same toy covariate matrix and coefficient either way
  mat <- matrix(c(0.1, 0.4, 0.7, 1), ncol = 1, dimnames = list(NULL, "cov"))
  covPars <- c(cov = 2)
  parShort <- c(maxAsymptote = 0.27, inflectionPoint1 = 4)                # DEoptim's par: no hillSlope1
  parOld <- c(maxAsymptote = 0.27, hillSlope1 = 1, inflectionPoint1 = 4)  # the pre-fix style, hand-fixed at 1

  newPred <- fireSenseUtils::logisticAll(fireSenseUtils:::fixHillSlope1(parShort), mat, covPars,
                                         lowerSpreadProb = 0.13)
  oldPred <- fireSenseUtils::logisticAll(parOld, mat, covPars, lowerSpreadProb = 0.13)
  expect_identical(newPred, oldPred)
  ## and it is not simply ignoring the slope: a different hillSlope1 changes the prediction
  offPred <- fireSenseUtils::logisticAll(c(maxAsymptote = 0.27, hillSlope1 = 0.4, inflectionPoint1 = 4),
                                         mat, covPars, lowerSpreadProb = 0.13)
  expect_false(isTRUE(all.equal(newPred, offPred)))
})
