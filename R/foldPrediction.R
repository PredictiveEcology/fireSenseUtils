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

#' @return `foldFeatureNames()`: the features the fold models were trained with (their `variable.names()`),
#'   without repeats.
#' @export
#' @rdname foldPrediction
foldFeatureNames <- function(fit) unique(unlist(lapply(foldModels(fit), stats::variable.names)))
