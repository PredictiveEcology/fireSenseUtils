#' Does a spread objective have an early stop?
#'
#' The early stop compares each fire's first-block value (the `"SNLL_FS"` term of
#' [.objfunSpreadFit()]) with a threshold. An objective without that term, such as
#' `tests = "crps"`, has no first-block value, so there is nothing to calibrate a threshold on.
#' This is the one place that rule is written; [.objfunSpreadFit()] uses it too.
#'
#' @param tests character; the `tests` of [.objfunSpreadFit()]. Case-insensitive.
#'
#' @return `TRUE` if `tests` includes `"SNLL_FS"`, else `FALSE`.
#'
#' @export
#' @examples
#' objectiveHasEarlyStop("crps")                    # FALSE
#' objectiveHasEarlyStop(c("adTest", "SNLL_FS"))    # TRUE
objectiveHasEarlyStop <- function(tests) {
  any(grepl("snll_fs", tolower(tests)))
}
