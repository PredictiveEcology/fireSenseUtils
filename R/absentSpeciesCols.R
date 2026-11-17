#' Species covariate columns a run's vegetation may lack
#'
#' A fit trained with a species (fuel class) covariate, e.g. `Pinu_con`, needs that column to predict. A run
#' whose vegetation lacks the species has no such column, but the species has no biomass there: the column
#' is 0 (or, for fuel on the log scale, `logMinB(0)`; see `value`). Any other missing covariate
#' (climate, a non-forest group, `youngAge`) cannot be filled in this way.
#'
#' `speciesCovCols()` names the species covariates among `covNames`: the fuel columns
#' ([spreadFuelCols()]: everything that is not an indicator) that are not climate variables
#' ([climateNAnames()] names, with or without underscores), lightning (`lightning*`), interaction terms (`a:b`),
#' `pixelID` or `year`. Drought indices (`droughtIndexPattern`: `CMD`, `MDC`, also in `cumMDC`) are climate too, whether or not
#' `climateData` lists them.
#' A column that is none of these is taken to be a species, so a covariate of another kind must be named
#' here before a model uses it.
#'
#' @param covNames Character vector of column names.
#' @return `speciesCovCols()`: the names among `covNames` that are species covariates.
#' @export
#' @rdname absentSpeciesCols
speciesCovCols <- function(covNames) {
  climate <- getFromNamespace(".allowedClimateVars", ns = asNamespace("climateData"))
  isClimate <- gsub("_", "", covNames) %in% gsub("_", "", climate)
  fuel <- spreadFuelCols(covNames[!isClimate])
  setdiff(fuel, c("pixelID", yearTxt, spreadInteractionCols(fuel), grep("^lightning", fuel, value = TRUE),
                  grep(droughtIndexPattern, fuel, value = TRUE)))
}

#' @param features Character vector: the columns a model needs (see [foldFeatureNames()]).
#' @return `absentSpeciesCols()`: the species covariates among `features` that are not in `covNames`, but
#'   none (`character(0)`) when `covNames` has a species covariate that is not in `features`: that is a renamed or
#'   merged class, not an absent species, and the model does not match the covariates.
#' @export
#' @rdname absentSpeciesCols
absentSpeciesCols <- function(features, covNames) {
  if (length(setdiff(speciesCovCols(covNames), features))) return(character())
  intersect(setdiff(features, covNames), speciesCovCols(features))
}

#' @param dt A `data.table` of covariates.
#' @param value What an absent species is, in the units of `dt` (and before any rescaling or centring, so
#'   that the model's own rescaling treats it as it treats the fit's data): 0, the default, for a
#'   biomass or cover; `logMinB(0)` for fuel stored as log biomass.
#' @return `addAbsentSpeciesCols()`: `dt`, or a copy with a column of `value` for each of
#'   `absentSpeciesCols(features, names(dt))`. A message names the added columns. `dt` is not changed.
#' @export
#' @rdname absentSpeciesCols
addAbsentSpeciesCols <- function(dt, features, value = 0) {
  absent <- absentSpeciesCols(features, names(dt))
  if (!length(absent)) return(dt)
  message("Species covariate(s) ", paste(absent, collapse = ", "), " are absent from this run's vegetation; ",
          "setting them to ", value, ".")
  dt <- data.table::copy(dt)
  for (cn in absent) data.table::set(dt, NULL, cn, value)
  dt
}
