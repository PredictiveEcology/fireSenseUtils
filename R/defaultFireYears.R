#' Default fire years for the fireSense modules
#'
#' The years of fire records that `fireSense_dataPrepFit` fits to, and over which
#' `fireSense_ELFs` counts fires to merge ELFs with too few of them. Both modules take
#' their `fireYears` default from here so that the ELFs and the fit always use the
#' same years.
#'
#' The years run from 1985, the first SCANFI V2 year, to the latest year with historical
#' climate for every tile ([climateData::latestHistoricalYear()]); climate is the last of
#' the inputs to reach a year.
#'
#' @return An integer vector of years.
#'
#' @export
defaultFireYears <- function() {
  1985L:climateData::latestHistoricalYear()
}
