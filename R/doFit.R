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
#' @param fitFileSuffix the module's `fitFileSuffix` parameter, i.e., which ledger was searched, for the
#'   error message. `NULL` (default) leaves it out of the message.
#' @return `TRUE` to fit, `FALSE` to use the ledger's fit. Stops when `doFit` is `FALSE` and the
#'   ledger has no usable fit.
#' @export
mustFit <- function(doFit, hasLedgerFit, polygon = "", fitFileSuffix = NULL) {
  if (isTRUE(doFit)) return(TRUE)
  hasFit <- isTRUE(hasLedgerFit)
  if (!hasFit && isFALSE(doFit))
    stop("There is no usable stored fit for polygon ", if (length(polygon)) polygon[1] else "<none>",
         " in the ledger", if (!is.null(fitFileSuffix)) paste0(" with fitFileSuffix \"", fitFileSuffix, "\""),
         ", and `doFit` is FALSE (the default: never fit, only predict from a stored fit). ",
         "To fit it, set `doFit = NA` (fit only when there is no stored fit) or `doFit = TRUE` (always fit).",
         call. = FALSE)
  !hasFit
}
