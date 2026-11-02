#' Per-ELF species, fuel and climate inputs from SpreadFit ledger rows
#'
#' Each row of the SpreadFit ledger is one ELF's fit, with the species table (`sppEquiv`), non-forest
#' land-cover groups and missing-land-cover group it was fitted with, and its parameters.
#' `fireSense_dataPrepFit` builds the per-ELF inputs from the ledger it reads before the fit, and
#' `fireSense_spreadFit` rebuilds them from the ledger after it writes its own row, so prediction in the
#' same run has the new ELF.
#'
#' `spreadFitLedgerInputs()` harmonizes each row's `sppEquiv` with [LandR::sppHarmonize()] (species
#' without a name in `sppEquivCol` are dropped and a colour vector is made), and finds the climate
#' variables of its parameters: the columns of `params` that are neither logistic parameters
#' ([logisticParamNames]), `youngAge`, fuel classes nor non-forest groups, and are a climate variable of
#' `climateData`.
#'
#' @param ledger `data.frame` or `sf`, SpreadFit ledger rows, with columns `polygonID`
#'   ([polygonIDTxt]), `params` and the list-columns in `colNames`.
#' @param sppEquivCol character, the column of `sppEquiv` naming the species.
#' @param fuelClassCol character, the column of `sppEquiv` naming the fuel classes.
#' @param studyArea passed to [LandR::sppHarmonize()]; used only when a row has no `sppEquiv`.
#' @param colNames character, the list-columns to read (those holding `"param"` are left out). Columns
#'   the ledger lacks (rows written before they existed) are skipped.
#'
#' @return `spreadFitLedgerInputs()`: a list with one element per row, named by `polygonID`, each a list
#'   of `missingLCCgroup` and `nonForestedLCCGroups` (when the ledger has them), `sppColorVect`,
#'   `sppEquiv` and `sppNameVector` (harmonized), `theseClimVars` and `theseClimVarsNoUnderscore`.
#' @export
#' @rdname spreadFitLedgerInputs
spreadFitLedgerInputs <- function(ledger, sppEquivCol, fuelClassCol, studyArea,
                                  colNames = spreadFitAdditionalColNamesTxt) {
  colNames <- setdiff(colNames, grep(value = TRUE, "param", colNames))
  colNames <- intersect(colNames, names(ledger))
  allVars <- unique(utils::getFromNamespace(".allowedClimateVars", ns = "climateData")) # it lists CMD twice
  allVarsNoUnderscore <- gsub("_", "", allVars)
  lpn <- logisticParamNames
  ids <- as.character(ledger[[polygonIDTxt]])
  out <- lapply(seq_len(NROW(ledger)), function(ind) {
    row <- lapply(stats::setNames(colNames, colNames), function(nm) ledger[[nm]][[ind]])
    sppEquiv <- row[["sppEquiv"]]
    sppOuts <- LandR::sppHarmonize(sppEquiv, unique(sppEquiv[[sppEquivCol]]), sppEquivCol, sppColorVect = NULL,
                                   vegLeadingProportion = NULL, studyArea = studyArea)

    pars <- ledger[["params"]][[ind]]
    allMatched <- sapply(lpn, function(lpn) all(lpn %in% colnames(pars)))
    numMatches <- sapply(lpn, function(lpn) sum(colnames(pars) %in% lpn))
    whLogistic <- which(allMatched & numMatches == max(numMatches))
    FuelNames <- c(sppOuts$sppEquiv[[fuelClassCol]], names(row[["nonForestedLCCGroups"]]))
    if (youngAgeTxt %in% colnames(pars))
      FuelNames <- c(youngAgeTxt, FuelNames)
    ClimateNames <- setdiff(setdiff(colnames(pars), lpn[[whLogistic]]), FuelNames)
    whClimateVar <- which(allVarsNoUnderscore %in% ClimateNames)

    c(row[intersect(c("missingLCCgroup", "nonForestedLCCGroups"), colNames)],
      sppOuts[c("sppColorVect", "sppEquiv", "sppNameVector")],
      list(theseClimVars = allVars[whClimateVar], theseClimVarsNoUnderscore = allVarsNoUnderscore[whClimateVar]))
  })
  stats::setNames(out, ids)
}

#' @description
#' `spreadFitLedgerInputLists()` regroups them into the per-ELF lists that `fireSense_dataPrepPredict`
#' reads.
#'
#' @param ledgerInputs the result of `spreadFitLedgerInputs()`.
#'
#' @return `spreadFitLedgerInputLists()`: a named list of five lists, each with one element per ledger
#'   row (named by `polygonID`): `sppEquivs`, `sppColorVects`, `sppNameVectors`,
#'   `nonForestedLCCGroupsList` and `missingLCCgroupList`.
#' @export
#' @rdname spreadFitLedgerInputs
spreadFitLedgerInputLists <- function(ledgerInputs) {
  elements <- c(sppEquivs = "sppEquiv", sppColorVects = "sppColorVect", sppNameVectors = "sppNameVector",
                nonForestedLCCGroupsList = "nonForestedLCCGroups", missingLCCgroupList = "missingLCCgroup")
  lapply(elements, function(el) lapply(ledgerInputs, function(x) x[[el]]))
}
