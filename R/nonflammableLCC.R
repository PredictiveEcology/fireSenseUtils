#' Non-flammable land-cover codes for the FireSense land cover
#'
#' The single source of truth for which codes in the land cover built by
#' [makeFireSenseLCC()] are non-flammable. `fireSense_dataPrepFit` and
#' `fireSense_dataPrepPredict` take their `nonflammableLCC` parameter default from this
#' constant, and [makeFireSenseLCC()] and [ELFflammableArea()] take theirs from it too,
#' so the land cover a fit/prediction is built from and the land cover it is run on
#' cannot silently disagree on what counts as flammable.
#'
#' Every code below comes from [makeFireSenseLCC()]'s two land-cover sources,
#' [LandR::prepInputs_NTEMS_LCC_FAO()] and [LandR::prepInputs_SCANFI_LCC_FAO()]
#' (documented in the "Data codes" comment in `LandR/R/prepInputs_NTEMS.R` and
#' `LandR/R/maps.R`), plus SCANFI's raw-class recode in
#' [LandR::convert_SCANFI_LCC_codes()]:
#'
#' - `0`: no data / unclassified -- non-flammable.
#' - `20`: water -- non-flammable.
#' - `30`: rock/exposed land, SCANFI's combined rock+exposed class
#'   ([LandR::convert_SCANFI_LCC_codes()]) -- non-flammable. This was missing from the
#'   modules' old default, so SCANFI-derived rock entered fits as flammable non-forest.
#' - `31`: snow/ice (NTEMS) -- non-flammable.
#' - `32`: rock/rubble (NTEMS) -- non-flammable.
#' - `33`: exposed/barren land (NTEMS) -- non-flammable.
#' - `40`: bryoids -- flammable ground cover, not included.
#' - `50`: shrubs -- flammable, not included.
#' - `80`, `81`: wetland, wetland-treed -- flammable, not included.
#' - `100`: herbs -- flammable, not included.
#' - `210`, `220`, `230`: coniferous, broadleaf, mixedwood forest -- flammable, not
#'   included.
#' - `240`: FAO-forest disturbed code (`disturbedCode` in
#'   [LandR::prepInputs_SCANFI_LCC_FAO()]/[LandR::prepInputs_NTEMS_LCC_FAO()]) --
#'   flammable forest, not included.
#'
#' @format `numeric` vector, `c(0, 20, 30, 31, 32, 33)`.
#'
#' @export
fireSenseNonflammableLCC <- c(0, 20, 30, 31, 32, 33)
