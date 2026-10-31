#' What a fireSense hindcast job changes in a `setupProject()` call
#'
#' A hindcast holds the vegetation of SCANFI year `year` static and draws that year's climate. It runs
#' 10 replicates as decimal years `year` to `year + 9 * step` (20, with `allSets`, one per spread parameter
#' set). It runs the climate and fireSense modules only: no vegetation or summary modules, so
#' `modules` replaces the project's own. The modules name their interval parameter differently
#' (`fireTimeStep`, `.runInterval`, `hindcastInterval`); this sets each of them to `step`.
#' `.plots = NA` keeps the modules that run every replicate from plotting. The data modules keep the
#' fit's `.plots`, which is in their cache keys, so a hindcast reuses the fit's data preparation.
#'
#' @param year The SCANFI year whose vegetation is held static.
#' @param step Years between replicates.
#' @param allSets If `TRUE`, one job runs every fitted spread parameter set (20 replicates, 5 sets of 4).
#' @param ignitions `NA` for the `fireSense_burn` default (`"predicted"`); `"observedEscapes"` spreads
#'   each year's observed escapes (spread alone).
#'
#' @return A list of `times`, `modules` and `params`. Use `times` and `modules` instead of the
#'   project's, and merge `params` into the project's with `Require::modifyList2()`.
#'
#' @export
hindcastSetup <- function(year, step = 0.001, allSets = FALSE, ignitions = NA_character_) {
  stopifnot(length(year) == 1L, !is.na(year), step > 0)
  interval <- function(arg) c(stats::setNames(list(step), arg), list(.plots = NA))
  list(
    # end at 9.5 steps, as nine additions of `step` land just below 9 steps (19.5: nineteen)
    times = list(start = year, end = year + (if (isTRUE(allSets)) 19.5 else 9.5) * step),
    modules = c("PredictiveEcology/canClimateData@modsForFireSense",
                "PredictiveEcology/climateYear@modsForFireSense",
                "PredictiveEcology/fireSense@development"),
    params = list(
      climateYear = list(samplingRange = year),
      fireSense_dataPrepPredict = interval("fireTimeStep"),
      fireSense_ignitionPredict = interval(".runInterval"),
      fireSense_spreadPredict = c(interval(".runInterval"),
                                  if (isTRUE(allSets)) list(paramSetByReplicate = TRUE)),
      fireSense_burn = c(interval(".runInterval"),
                         if (!is.na(ignitions)) list(ignitionSource = ignitions)),
      fireSense_summary = c(interval("hindcastInterval"), list(hindcastYear = year))))
}
