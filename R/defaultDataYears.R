#' SCANFI V2 vegetation years
#'
#' The years that SCANFI V2 has, every five years from 1985 to 2025. This repeats
#' `LandR:::.scanfi_v2_years`, which LandR does not export: keep the two in step.
#'
#' @export
scanfiV2Years <- seq(1985L, 2025L, 5L)

#' Default vegetation data years for the fireSense modules
#'
#' The vegetation years that `fireSense_dataPrepFit` builds by default ("all years that
#' exist"): every SCANFI V2 year ([scanfiV2Years]) strictly before the last of
#' [defaultFireYears()]. A data year needs a fire year after it to be fitted against, so with
#' historical climate to 2024 the years are 1985, 1990, ..., 2020.
#'
#' @return An integer vector of years.
#'
#' @export
defaultDataYears <- function() {
  scanfiV2Years[scanfiV2Years < max(defaultFireYears())]
}
