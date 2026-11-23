#' The fuel terms of a spread model
#'
#' @param covNames Character; covariate names, e.g. [spreadSourceCols()] of the formula or the names of a fit's parameters.
#' @return The names that are fuel biomass: the dominant and secondary fuel classes ([domAgbTxt], [secAgbTxt]) and
#'   [treedWetlandAgbTxt].
#' @export
spreadFuelTerms <- function(covNames)
  covNames[startsWith(covNames, domAgbTxt) | startsWith(covNames, secAgbTxt) | covNames == treedWetlandAgbTxt]

#' The change in spread probability that a fuel makes
#'
#' The coefficients of the fuel terms are not comparable between fuels: a fuel covariate is the class's biomass over
#' `1e4` (0 where the class is absent), classes differ hugely in how much biomass they have, the link saturates and
#' correlated classes trade off. This puts every fuel on one scale. For each fuel term of a ledger row it is the change in
#' spread probability from 0 to the term's median non-zero value in the fit's data (`covMedianNonZero_spread`,
#' [spreadFuelSummary()]), with every other covariate at its centre (`covCentre_spread`), through the fitted link, averaged over the
#' row's parameter sets. The link is evaluated by [spreadProbFromCovs()], the function the fit and `fireSense_spreadPredict` use,
#' with the link's parameters read from the row as `fireSense_spreadPredict` reads them, so a fit near its ceiling shows
#' an effect near 0 whatever its coefficient. With climate outside the ceiling or in it (`climateRef_spread`), climate is at its
#' reference; with the activation, at the reference too.
#'
#' @param row A one-row ledger row, as [latestSpreadFits()] gives.
#' @param lowerSpreadProb The spread probability floor, `fireSense_spreadPredict`'s `lowerSpreadProb`.
#' @return A named numeric vector, one change in probability per fuel term of the row's parameters (empty without
#'   fuel terms); `NA` for all of them when the row has no `covMedianNonZero_spread` (a row written before it), and for
#'   a fuel never non-zero.
#' @export
spreadFuelEffect <- function(row, lowerSpreadProb = spreadProbFloor) {
  one <- function(col) if (col %in% names(row) && length(row[[col]])) row[[col]][[1L]]
  actTxt <- activationInterceptTxt
  params <- as.data.frame(one("params"))
  fuelNames <- spreadFuelTerms(names(params))
  na <- stats::setNames(rep(NA_real_, length(fuelNames)), fuelNames)
  medians <- one(spreadFitCovMedianNonZeroTxt)
  if (!length(fuelNames) || !is.list(medians)) return(na)

  ## the row's link, climate and centres, read as fireSense_spreadPredict reads them (ledger*(), R/ledgerRow.R)
  linkArgs <- ledgerLinkArgs(row, 1L, floor = lowerSpreadProb)
  link <- linkArgs$link
  mode <- ledgerSpreadClimate(row, 1L)
  climateRef <- ledgerClimateRef(row, 1L)
  if (!is.null(link) && identical(mode, "predictor")) climateRef <- NULL
  ## the covariate values of "everything at its centre": 0 centred, and with climate in or outside the ceiling, the reference
  centre <- ledgerCovCentre(row, 1L)
  centreOf <- function(cn) if (is.null(centre[[cn]])) 0 else centre[[cn]]
  design <- setdiff(names(params), c(unlist(logisticParamNames), spreadRandomEffectTxt, actTxt))
  sources <- spreadSourceCols(design)
  atCentre <- stats::setNames(numeric(length(sources)), sources)
  if (length(climateRef)) atCentre[names(climateRef)] <- unlist(climateRef)
  ## two rows per fuel: it absent (0) and at its median non-zero, as the objective's design matrix has them
  at <- function(fuel, value) {
    v <- atCentre
    v[fuel] <- value - centreOf(fuel)
    vapply(design, function(cn) if (cn == spreadInterceptTxt) 1 else prod(v[strsplit(cn, ":", fixed = TRUE)[[1L]]]), numeric(1))
  }
  fuels <- intersect(fuelNames, names(medians)[!vapply(medians, is.na, logical(1))])
  if (!length(fuels)) return(na)
  mat <- do.call(rbind, lapply(fuels, function(f) rbind(at(f, 0), at(f, medians[[f]]))))
  colnames(mat) <- design
  isHigh <- rep(c(FALSE, TRUE), length(fuels))
  powerLawArgs <- c(list(climateRef = climateRef, powerLawPc50 = if (!is.null(climateRef)) ledgerPowerLawPc50(row, 1L)),
                    linkArgs[intersect(names(linkArgs), c("link", "powerLawGamma", "powerLawShape", "lowerAsymptote", "richardsNu"))])
  change <- vapply(seq_len(NROW(params)), function(k) {
    par <- unlist(params[k, , drop = FALSE])
    p <- as.vector(do.call(spreadProbFromCovs, c(list(mat, logisticPars = par[names(par) %in% unlist(logisticParamNames)],
                                                       covPars = par[design], lowerSpreadProb = lowerSpreadProb,
                                                       activationIntercept = if (identical(mode, "activation")) params[[actTxt]][k]),
                                                  powerLawArgs)))
    p[isHigh] - p[!isHigh]
  }, numeric(length(fuels)))
  na[fuels] <- rowMeans(matrix(change, nrow = length(fuels)))
  na
}
