## Which climate variables the fire models use, and the climate layers that requires. Moved from
## fireSense_dataPrepFit, so fireSense_dataPrepPredict can ask for the climate of stored fits in a run
## without fireSense_dataPrepFit.
##
## `climateVariablesForFire` names the variables for each process: IgnitionFit (xgboost) takes all of
## `ignition`, as it copes with correlated covariates; SpreadFit takes `spread`. Unless a user supplies
## them, fireSense_dataPrepFit sets both it and `climateVariables` (what canClimateData prepares) from the
## default below, so the two always agree and a user need not choose. canClimateData leaves
## `climateVariables` to the fireSense module that lists it as an output, and keeps its own default
## otherwise.

#' Default climate variables for the fire models
#'
#' ClimateNA names; `"cumMDC"` is derived by `climateData::calcCumMDC()`. `spread = "auto"` picks,
#' per study area, the `ignition` variable that best tracks annual area burned
#' (`fireSense_dataPrepFit`'s `selectSpreadClimateVariable()`).
#' @export
defaultClimateVariablesForFire <- list(ignition = c("CMD", "cumMDC", "CMD_sm", "CMD_sp"),
                                       spread = "auto")

#' ClimateNA names of climate variables given with or without underscores
#'
#' The fire modules name climate layers without underscores (`"CMDsm"`), ClimateNA with them
#' (`"CMD_sm"`). Both forms map to the ClimateNA name.
#'
#' @param x Character vector of variable names.
#' @return Character vector of ClimateNA names, same length as `x`.
#' @export
climateNAnames <- function(x) {
  allowed <- getFromNamespace(".allowedClimateVars", ns = asNamespace("climateData"))
  out <- allowed[match(gsub("_", "", x), gsub("_", "", allowed))]
  if (anyNA(out))
    stop("Unknown climate variable(s): ", paste(x[is.na(out)], collapse = ", "),
         ". Use ClimateNA names (e.g., CMD_sm) or climateData's derived variables (e.g., cumMDC).")
  out
}

#' Climate layers needed by the fire models
#'
#' @param climateVariablesForFire List with `ignition` and `spread` character vectors.
#' @param historicalYears Years of historical climate: the fire years for a fit, the simulated years
#'   with historical climate for prediction. None (`NULL` or empty) requests no historical climate.
#' @param projectedYears Projected years, from `canClimateProjectedYears()`. Projected climate is
#'   requested for these years only; none (`NULL` or empty) requests no projected climate.
#' @return A list for `sim$climateVariables`, as built by `climateData::climateLayers()`.
#' @export
fireClimateLayers <- function(climateVariablesForFire, historicalYears, projectedYears = NULL) {
  vars <- unlist(climateVariablesForFire[c("ignition", "spread")], use.names = FALSE)
  vars <- unique(climateNAnames(setdiff(vars, "auto")))    # "auto" chooses among the ignition variables
  climateData::climateLayers(vars, fun = quote(calcAsIs), historical = length(historicalYears) > 0,
                             historicalYears = historicalYears,
                             projected = length(projectedYears) > 0, projectedYears = projectedYears)
}

#' Add the climate variables of existing fits to what is being prepared
#'
#' For prediction from existing fits, possibly of several ELFs that chose different spread variables
#' (`spread = "auto"`): every variable a fit uses must be prepared. Variables already in
#' `climateVariables` keep their definition (and years); only missing ones are added.
#'
#' @param climateVariables The list for `sim$climateVariables` (names like `historical_CMDsm`).
#' @param fitClimVars Climate variables the fits use, ClimateNA names.
#' @inheritParams fireClimateLayers
#' @return `climateVariables`, with any missing variable added.
#' @export
addFitClimateVariables <- function(climateVariables, fitClimVars, historicalYears, projectedYears = NULL) {
  prepared <- unique(sub("^[^_]+_", "", names(climateVariables)))
  toAdd <- fitClimVars[!gsub("_", "", fitClimVars) %in% prepared]
  if (!length(toAdd)) return(climateVariables)
  c(climateVariables, fireClimateLayers(list(ignition = toAdd), historicalYears = historicalYears,
                                        projectedYears = projectedYears))
}

#' The years of projected climate that canClimateData prepares
#'
#' canClimateData's `projectedClimateYears` is the one place these are set (its default is
#' 2011:2100). Empty, or without canClimateData in the run, means no projected climate: runs that
#' sample past climate (NRV) set it empty.
#'
#' @param sim A `simList`.
#' @return Integer/numeric years; `NULL` if canClimateData is not in the run.
#' @export
canClimateProjectedYears <- function(sim) {
  tryCatch(SpaDES.core::P(sim, module = "canClimateData")$projectedClimateYears, error = function(e) NULL)
}

#' The fire models' climate variables with those of stored fits
#'
#' `spread` becomes the fits' variables, all of them (prediction blends neighbour ELFs, which may have
#' chosen different ones); `ignition` (xgboost) keeps its variables and gains the fits'. Names without
#' underscores, as the fire modules use them (`"CMDsm"`).
#'
#' @param climateVariablesForFire List with `ignition` and `spread`.
#' @param fitClimVars Climate variables the fits use, from [spreadFitClimateVariables()].
#' @return `climateVariablesForFire`, updated.
#' @export
addFitClimateVariablesForFire <- function(climateVariablesForFire, fitClimVars) {
  fitVars <- gsub("_", "", fitClimVars)
  climateVariablesForFire[["spread"]] <- fitVars
  climateVariablesForFire[["ignition"]] <- sort(unique(c(gsub("_", "", climateVariablesForFire[["ignition"]]), fitVars)))
  climateVariablesForFire
}

#' The climate variables of stored spread fits
#'
#' Every climate variable any of the ledger rows' fits uses (each ELF chose its own with
#' `spread = "auto"`), for preparing the climate of prediction, which blends neighbour ELFs.
#'
#' @param ledgerInputs As returned by [spreadFitLedgerInputs()].
#' @return Character vector of ClimateNA names (e.g. `"CMD_sm"`), each once.
#' @export
spreadFitClimateVariables <- function(ledgerInputs) {
  unique(unlist(lapply(ledgerInputs, `[[`, "theseClimVars"), use.names = FALSE))
}
