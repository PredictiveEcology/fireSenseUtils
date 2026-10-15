#' Resolve the `doFit` parameter of the FireSense fit modules
#'
#' `fireSense_spreadFit` and `fireSense_ignitionFit` decide whether to fit from one three-valued
#' parameter, `doFit`: `NA` (default) fits only when the ledger has no usable fit for the polygon,
#' `TRUE` always fits (replacing a ledger fit), `FALSE` never fits (use the ledger fit; stop when
#' there is none). `refitExisting` and `stopIfNoPreRunFit` are deprecated aliases: `refitExisting =
#' TRUE` is `doFit = TRUE`; `stopIfNoPreRunFit = TRUE` is `doFit = FALSE` (`stopIfNoPreRunFit =
#' FALSE` only allowed a fit, which is what `NA` does); when both are set `refitExisting` wins. An
#' explicit (non-`NA`) `doFit` is never overridden by an alias. A module that does not have an
#' alias leaves it `NA`.
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
  if (!isTRUE(is.na(doFit))) return(isTRUE(doFit))
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
    stop("There is no usable pre-run fit in the ledger for polygon ", polygon[1],
         " and `doFit` is FALSE. Set `doFit = NA` (or TRUE) to fit it.", call. = FALSE)
  !hasFit
}
