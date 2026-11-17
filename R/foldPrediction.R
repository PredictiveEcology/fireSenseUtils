#' The fold models of an ignition or escape fit, and their mean prediction
#'
#' `fireSense_ignitionFit` fits one xgboost model per cross-validation fold and stores them in
#' `fit$modelList$model`, named `Fold1`, `Fold2`, ... (beside non-model elements such as `rocs` and
#' `yearTotals`). `fireSense_ignitionPredict` predicts with their mean, and `fireSense_ignitionFit`
#' uses the same mean for the fitted year effects, so the two agree.
#'
#' @param fit A `fireSense_IgnitionFitted` or `fireSense_EscapeFitted`: a list with `modelList$model`.
#'
#' @return `foldModels()`: the list of fold models (the elements of `fit$modelList$model` whose name
#'   contains `"Fold"`). It stops if there are none.
#' @export
#' @rdname foldPrediction
foldModels <- function(fit) {
  mods <- fit$modelList$model
  modsOnly <- mods[grep("Fold", names(mods))]
  if (is.null(modsOnly) || !length(modsOnly))
    stop("the fit has no fold models (elements named 'Fold*' in modelList$model)")
  modsOnly
}

#' @param newdata The covariates to predict, prepared as for the fit (standardised with its
#'   `scaleData` when it was fitted with `rescaleVars`).
#'
#' @return `foldMeanPrediction()`: numeric, one per row of `newdata`: the mean over the fold models of
#'   `predict(model, newdata)`, each floored at 0 (an expected count).
#' @export
#' @importFrom stats predict
#' @rdname foldPrediction
foldMeanPrediction <- function(fit, newdata) {
  preds <- lapply(foldModels(fit), function(model) pmax(predict(model, newdata = newdata), 0))
  rowMeans(do.call(cbind, preds))
}

#' The chance that one ignition escapes, from the expected escapes of an ignited pixel
#'
#' The escape model is fitted on ignited pixel-years only, so its mean `escapeMean` is the expected escapes of
#' a pixel that has at least one ignition. With Poisson ignitions of mean `ignitionMean`, each ignition
#' escapes with probability `q = escapeMean * (1 - exp(-ignitionMean)) / ignitionMean`, capped at 1:
#' `(1 - exp(-ignitionMean)) / ignitionMean` is `1 / E[ignitions | ignitions > 0]`, so the expected escapes of
#' an ignited pixel are `escapeMean` and never exceed its ignitions. `fireSense_burn` draws a pixel's escapes
#' as `rbinom(ignitions, q)`, so the expected escapes of a pixel are `ignitionMean * q`.
#'
#' @param escapeMean Numeric: expected escapes of an ignited pixel; `NA` is 0.
#' @param ignitionMean Numeric: the mean the pixel's ignitions are drawn from. At 0, `q` is its limit,
#'   `min(1, escapeMean)`.
#'
#' @return Numeric, `q`, recycled to the longer argument.
#' @export
#' @examples
#' escapeProbPerIgnition(0.3, c(0.01, 1, 5))
escapeProbPerIgnition <- function(escapeMean, ignitionMean) {
  n <- max(length(escapeMean), length(ignitionMean))
  escapeMean <- rep_len(escapeMean, n)
  ignitionMean <- rep_len(ignitionMean, n)
  escapeMean[is.na(escapeMean)] <- 0
  pmin(1, ifelse(ignitionMean > 0, escapeMean * (-expm1(-ignitionMean)) / ignitionMean, escapeMean))
}

#' The features of an ignition or escape fit, and whether a stored fit can predict
#'
#' `fireSense_ignitionFit` gives xgboost every covariate column except `pixelID`, `year`,
#' [coreFractionTxt], the `yearChar*` columns and the responses ([ignitionsTxt], [escapesTxt]), in
#' alphabetical order (`xgbFeatureNames()`). A fit stored in the ledger records the features it was trained on
#' (its fold boosters' `variable.names()`); if the covariates have changed since (e.g. non-forest land-cover
#' groups merged into one), `predict()` fails on the missing columns. `fireSense_ignitionFit` checks this
#' before reusing a stored fit, against the covariates it would fit with; `fireSense_ignitionPredict` checks it
#' before predicting, against the covariates it predicts with, so the check holds when no fit data were
#' prepared.
#'
#' @param cols Character vector of the column names of the covariates.
#'
#' @return `xgbFeatureNames()`: character vector.
#' @export
#' @rdname storedFitProblem
xgbFeatureNames <- function(cols) {
  cols <- setdiff(cols, c("pixelID", "year", coreFractionTxt))
  cols <- grep("yearChar", cols, invert = TRUE, value = TRUE)
  sort(grep(paste0(ignitionsTxt, "|", escapesTxt), cols, value = TRUE, invert = TRUE))
}

#' @param fitted A stored `fireSense_IgnitionFitted` or `fireSense_EscapeFitted`: a list with
#'   `modelList$model`, or `NULL`.
#'
#' @return `foreignObjectives()`: the objectives of the fold boosters that differ from
#'   [ignitionObjectiveTxt] (empty when all match); `"unknown"` stands for a booster whose objective cannot
#'   be read ([xgbObjective()]).
#' @export
#' @rdname storedFitProblem
foreignObjectives <- function(fitted) {
  boosters <- fitted$modelList$model
  boosters <- boosters[grep("^Fold", names(boosters))]
  obj <- xgbObjective(boosters)
  obj[is.na(obj)] <- "unknown"
  unique(obj[obj != ignitionObjectiveTxt])
}

#' @param covariates The covariates the stored fit would predict (or be refitted) with: a `data.frame`, or
#'   anything with `colnames()`.
#'
#' @return `storedFitProblem()`: `NULL` when `fitted` has fold boosters, all fitted with
#'   [ignitionObjectiveTxt], whose features are all among `xgbFeatureNames(colnames(covariates))`; otherwise
#'   one string saying why not (no model; the objectives; the missing features, and the covariates the model
#'   does not use), for the caller's message.
#' @export
#' @rdname storedFitProblem
storedFitProblem <- function(fitted, covariates) {
  boosters <- fitted$modelList$model
  boosters <- boosters[grep("^Fold", names(boosters))]
  if (!length(boosters))
    return("there is no stored model")
  foreign <- foreignObjectives(fitted)
  if (length(foreign))
    return(paste0("it was fitted with objective ", paste(foreign, collapse = ", "), ", not ", ignitionObjectiveTxt))
  stored <- unique(unlist(lapply(boosters, stats::variable.names)))
  current <- xgbFeatureNames(colnames(covariates))
  missingNow <- setdiff(stored, current)
  if (length(missingNow))
    return(paste0("it uses features absent from the current covariates (", paste(missingNow, collapse = ", "),
                  "); current covariates not in the model: ", paste(setdiff(current, stored), collapse = ", ")))
  NULL
}
