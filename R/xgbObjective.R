#' The objective an xgboost model was fitted with
#'
#' Reads `xgboost::xgb.config(model)$learner$objective$name`, then `attr(model, "params")$objective`.
#' The second is the only source for models that carry their parameters as an attribute; the first
#' is the only one for xgboost (>= 3) boosters, which have no `$params` and, once read back from a
#' file, are external pointers. Every fireSense module that needs the objective (to check a stored fit
#' against [ignitionObjectiveTxt], or to tell `fireSense_burn` how to draw counts) uses this function.
#'
#' @param model An `xgb.Booster`, or a plain (unclassed) list of them, e.g. the folds of a fit.
#'
#' @return For one booster a single string; for a list a character vector with one element per booster,
#'   keeping the list's names. `NA_character_` stands for a model whose objective cannot be read.
#'
#' @export
#' @examples
#' if (requireNamespace("xgboost", quietly = TRUE)) {
#'   x <- matrix(rnorm(100), ncol = 2)
#'   b <- xgboost::xgboost(x, rpois(50, 1), objective = "count:poisson", nrounds = 2, nthreads = 1)
#'   xgbObjective(b)
#'   xgbObjective(list(Fold1 = b, Fold2 = b))
#' }
xgbObjective <- function(model) {
  if (is.list(model) && is.null(oldClass(model)) && is.null(attr(model, "params")))
    ## a plain list is a list of models; a classed list, or one with `params`, is one model
    return(vapply(model, xgbObjective, character(1)))
  isOneString <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  o <- tryCatch(xgboost::xgb.config(model)$learner$objective$name, error = function(e) NULL)
  if (!isOneString(o)) o <- tryCatch(attr(model, "params")$objective, error = function(e) NULL)
  if (isOneString(o)) o else NA_character_
}
