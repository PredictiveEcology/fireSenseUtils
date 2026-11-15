## Small real boosters, as the folds of an ignition fit
foldFit <- function() {
  set.seed(2)
  x <- matrix(rnorm(200), ncol = 2, dimnames = list(NULL, c("MDC", "youngAge")))
  y <- rpois(100, exp(0.3 * x[, 1]))
  b <- function(n) xgboost::xgboost(x, y, objective = "count:poisson", nrounds = n, nthreads = 1, max_depth = 2)
  list(fit = list(modelList = list(model = list(Fold1 = b(2), Fold2 = b(4), rocs = list(), yearTotals = NULL))),
       x = x)
}

test_that("foldMeanPrediction() is the mean of the fold models' predictions, floored at 0", {
  skip_if_not_installed("xgboost")
  f <- foldFit()
  expect_identical(names(foldModels(f$fit)), c("Fold1", "Fold2"))
  m <- f$fit$modelList$model
  expect_equal(foldMeanPrediction(f$fit, f$x), (predict(m$Fold1, f$x) + predict(m$Fold2, f$x)) / 2)
  expect_error(foldModels(list(modelList = list(model = list(rocs = list())))), "no fold models")
})

test_that("escapeProbPerIgnition() keeps an ignited pixel's expected escapes and caps at 1", {
  lambda <- c(0.01, 0.5, 2, 8)
  e <- 0.4
  q <- escapeProbPerIgnition(e, lambda)
  ## E[escapes] = lambda * q = e * P(ignitions > 0): e per ignited pixel
  expect_equal(lambda * q, e * (1 - exp(-lambda)))
  expect_identical(escapeProbPerIgnition(5, 0.01), 1)
  expect_identical(escapeProbPerIgnition(NA, 1), 0)
  expect_identical(escapeProbPerIgnition(c(0.3, 2), 0), c(0.3, 1))
  ## byte-identical to fireSense_burn's drawEscapes() before it called this
  e <- c(0.01, 0.37, 1.3); lam <- c(0.003, 0.7, 4.1)
  expect_identical(escapeProbPerIgnition(e, lam), pmin(1, e * (-expm1(-lam)) / lam))
})
