covs <- function(cols) as.data.frame(setNames(rep(list(1), length(cols)), cols))

test_that("xgbFeatureNames() drops ids, coreFraction, yearChar and responses, and sorts", {
  expect_identical(xgbFeatureNames(c("pixelID", "year", "ignitions", "z", coreFractionTxt, "escapes",
                                     "yearChar2001", "a")), c("a", "z"))
})

test_that("storedFitProblem(): NULL when every feature is present; names the missing ones otherwise", {
  skip_if_not_installed("xgboost")
  fit <- storedFit()
  cur <- covs(c("pixelID", "year", "ignitions", "CMDsm", "Pice_mar", "youngAge"))
  expect_null(storedFitProblem(fit, cur))
  ## a covariate the model never used does not matter
  expect_null(storedFitProblem(fit, cbind(cur, extra = 1)))
  names(cur)[names(cur) == "Pice_mar"] <- "Pice_new"
  expect_match(storedFitProblem(fit, cur), "absent .*\\(Pice_mar\\).*not in the model: Pice_new")
})

test_that("storedFitProblem(): no model, or any fold with another objective, is a problem", {
  skip_if_not_installed("xgboost")
  cur <- covs(c("CMDsm", "Pice_mar", "youngAge"))
  expect_match(storedFitProblem(NULL, cur), "no stored model")
  expect_match(storedFitProblem(list(modelList = list(model = list(rocs = 1))), cur), "no stored model")
  expect_match(storedFitProblem(storedFit("reg:tweedie"), cur), "reg:tweedie, not count:poisson")
  mixed <- storedFit()
  mixed$modelList$model$Fold2 <- storedFit("reg:tweedie")$modelList$model$Fold1
  expect_identical(foreignObjectives(mixed), "reg:tweedie")
  expect_match(storedFitProblem(mixed, cur), "reg:tweedie")
})
