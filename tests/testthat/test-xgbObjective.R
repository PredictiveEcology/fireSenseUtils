## Real small boosters. xgboost >= 3 boosters have no `$params`; saved and read back they are
## external pointers, as the ignition-fit ledger holds them.
realBooster <- function(objective, reload = FALSE) {
  set.seed(1)
  x <- matrix(rnorm(200), ncol = 2, dimnames = list(NULL, c("MDC", "youngAge")))
  y <- rpois(100, exp(0.3 * x[, 1]))
  b <- xgboost::xgboost(x, y, objective = objective, nrounds = 3, nthreads = 1, max_depth = 2)
  if (reload) {
    f <- tempfile(fileext = ".rds")
    on.exit(unlink(f))
    saveRDS(b, f)
    b <- readRDS(f)
  }
  b
}

test_that("xgbObjective() reads the objective of a real booster, fresh or reloaded", {
  skip_if_not_installed("xgboost")
  for (reload in c(FALSE, TRUE)) {
    expect_identical(xgbObjective(realBooster("count:poisson", reload)), "count:poisson")
    expect_identical(xgbObjective(realBooster("reg:tweedie", reload)), "reg:tweedie")
  }
})

test_that("xgbObjective() of a list of boosters is one objective per booster, names kept", {
  skip_if_not_installed("xgboost")
  folds <- list(Fold1 = realBooster("count:poisson"), Fold2 = realBooster("reg:tweedie"))
  expect_identical(xgbObjective(folds), c(Fold1 = "count:poisson", Fold2 = "reg:tweedie"))
})

test_that("xgbObjective() falls back to attr(, 'params'), and gives NA when nothing is readable", {
  withParams <- structure(list(), params = list(objective = "count:poisson"))
  expect_identical(xgbObjective(withParams), "count:poisson")
  expect_identical(xgbObjective(structure(list(a = 1), class = "notABooster")), NA_character_)
  expect_identical(xgbObjective(list(Fold1 = withParams, Fold2 = NULL)), c(Fold1 = "count:poisson", Fold2 = NA_character_))
  ## a classed list (e.g. a model object that carries its params as an attribute) is one model, not a list of models
  expect_identical(xgbObjective(structure(list(fun = 1), params = list(objective = "reg:tweedie"), class = "toy")), "reg:tweedie")
  expect_identical(xgbObjective(NULL), NA_character_)
})

test_that("ignitionObjectiveTxt is the objective the ignition and escape fits use", {
  expect_identical(ignitionObjectiveTxt, "count:poisson")
})

test_that("the ignition ledger tag names the objective, so changing the objective needs a new tag", {
  family <- c("count:poisson" = "Poisson", "reg:tweedie" = "Tweedie")[[ignitionObjectiveTxt]]
  expect_match(ignitionFitFileTag, family, fixed = TRUE)
})
