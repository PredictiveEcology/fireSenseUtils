## Small real boosters (test-storedFitProblem.R, test-latestIgnitionFits.R) in the shape of a stored ledger fit, trained on `features`.
storedFit <- function(objective = "count:poisson", features = c("CMDsm", "Pice_mar", "youngAge"), folds = 2) {
  set.seed(1)
  x <- matrix(runif(100 * length(features)), ncol = length(features), dimnames = list(NULL, features))
  y <- rpois(100, exp(x[, 1]))
  model <- lapply(seq_len(folds), function(i)
    xgboost::xgboost(x, y, objective = objective, nrounds = 3, max_depth = 2, nthreads = 1))
  names(model) <- paste0("Fold", seq_len(folds))
  list(modelList = list(model = model, rocs = list()))
}
