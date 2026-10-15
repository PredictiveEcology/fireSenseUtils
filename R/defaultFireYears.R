#' Default fire years for the fireSense modules
#'
#' The years of fire records that `fireSense_dataPrepFit` fits to, and over which
#' `fireSense_ELFs` counts fires to merge ELFs with too few of them. Both modules take
#' their `fireYears` default from here so that the ELFs and the fit always use the
#' same years.
#'
#' The years run from 1986, the year after the first SCANFI V2 year (1985), to the latest
#' year with historical climate for every tile ([climateData::latestHistoricalYear()]);
#' climate is the last of the inputs to reach a year. A fire year is fitted against the latest
#' vegetation year strictly before it (the vegetation of year Y already shows year Y's burns),
#' so 1985 has no vegetation to pair with.
#'
#' @return An integer vector of years.
#'
#' @export
defaultFireYears <- function() {
  1986L:climateData::latestHistoricalYear()
}
