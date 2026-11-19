#' Resolve the `doFit` parameter of the FireSense fit modules
#'
#' `fireSense_spreadFit` and `fireSense_ignitionFit` decide whether to fit from one three-valued
#' parameter, `doFit`: `FALSE` (the modules' default) never fits (use the ledger fit; stop when
#' there is none), `NA` fits only when the ledger has no usable fit for the polygon, `TRUE` always
#' fits (replacing a ledger fit). `refitExisting` and `stopIfNoPreRunFit` are deprecated aliases:
#' `refitExisting = TRUE` is `doFit = TRUE`; `stopIfNoPreRunFit = TRUE` is `doFit = FALSE`
#' (`stopIfNoPreRunFit = FALSE` only allowed a fit, which is what `NA` does); when both are set
#' `refitExisting` wins. A `doFit` of `NA` or `TRUE` is never overridden by an alias. A `doFit` of
#' `FALSE` is, because it is the modules' default and so cannot be told apart from an unset `doFit`.
#' A module that does not have an alias leaves it `NA`.
#'
#' @param doFit logical, length one: `NA`, `TRUE` or `FALSE`.
#' @param refitExisting,stopIfNoPreRunFit the deprecated aliases; `NA` means not set.
#' @param module character; module name used in the deprecation message.
#' @param message logical; send the deprecation message when an alias is set (once per call).
#' @return `NA`, `TRUE` or `FALSE`.
#' @export
resolveDoFit <- function(doFit = NA, refitExisting = NA, stopIfNoPreRunFit = NA,
                         module = "the fit module", message = TRUE) {
  aliasSet <- !isTRUE(is.na(refitExisting)) || !isTRUE(is.na(stopIfNoPreRunFit))
  if (aliasSet && isTRUE(message))
    base::message(module, ": parameters `refitExisting` and `stopIfNoPreRunFit` are deprecated; ",
                  "use `doFit` (NA: fit only if the ledger has no fit; TRUE: always fit; FALSE: never fit).")
  if (isTRUE(doFit) || (isFALSE(doFit) && !aliasSet)) return(isTRUE(doFit))
  if (isTRUE(refitExisting)) return(TRUE)
  if (isTRUE(stopIfNoPreRunFit)) return(FALSE)
  NA
}

#' Must this run fit, given `doFit` and the ledger?
#'
#' @param doFit `NA`, `TRUE` or `FALSE`, as returned by [resolveDoFit()].
#' @param hasLedgerFit logical: does the ledger hold a usable fit for this polygon? Evaluated
#'   lazily, only when `doFit` is not `TRUE`, so a forced fit never consults the ledger.
#' @param polygon the polygon identifier, for the error message.
#' @return `TRUE` to fit, `FALSE` to use the ledger's fit. Stops when `doFit` is `FALSE` and the
#'   ledger has no usable fit.
#' @export
mustFit <- function(doFit, hasLedgerFit, polygon = "") {
  if (isTRUE(doFit)) return(TRUE)
  hasFit <- isTRUE(hasLedgerFit)
  if (!hasFit && isFALSE(doFit))
    stop("There is no usable stored fit for polygon ", if (length(polygon)) polygon[1] else "<none>",
         " in the ledger",
         ", and `doFit` is FALSE (the default: never fit, only predict from a stored fit). ",
         "To fit it, set `doFit = NA` (fit only when there is no stored fit) or `doFit = TRUE` (always fit).",
         call. = FALSE)
  !hasFit
}

#' Will this job read the data prepared for a fit?
#'
#' The one rule `fireSense_dataPrepFit` and `fireSense_spreadFit` share for `skipUnusedFitPrep`: a job that only
#' predicts from the ledger's fit does not need the spread covariates. A job needs them when it holds out a
#' fold (`heldOutFold`), validates or profiles (`mode`), always fits (`doFit` TRUE), or fits when the ledger has
#' no fit (`doFit` NA and no ledger fit, or not known). `doFit` FALSE never fits, so it does not.
#'
#' @param doFit `NA`, `TRUE` or `FALSE`, as returned by [resolveDoFit()].
#' @param heldOutFold `NA` (default), `1L` or `2L`.
#' @param mode character, `fireSense_spreadFit`'s `mode`; `NULL` if that module is not in the project.
#' @param hasLedgerFit logical: does the ledger hold a fit for this polygon? `NA` (default): not known, so
#'   `doFit = NA` may fit.
#' @return logical, length one.
#' @export
fitDataNeeded <- function(doFit, heldOutFold = NA, mode = NULL, hasLedgerFit = NA) {
  if (!isTRUE(is.na(heldOutFold)) || any(c("validate", "profile") %in% mode)) return(TRUE)
  if (isFALSE(doFit)) return(FALSE)
  isTRUE(doFit) || !isTRUE(hasLedgerFit)
}
