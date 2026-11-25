## The ledger row of a spread fit, read as `fireSense_spreadPredict` reads it: the covariate ranges and centres, the link and its
## parameters, and where climate enters. `fireSense_spreadPredict` (prediction) and `spreadFuelEffect()` call these, so a row is
## read in one place. Each takes `sa`, the rows of a ledger (`sim$studyAreaWithSpreadParams`, or [latestSpreadFits()]), and the row `i`.
## A row written before a column existed has no such column (or `NA`) and is read as it always was.

#' The covariate ranges of one ledger row
#'
#' A fit rescales its covariates with fixed ranges and stores them in its ledger row (`covMinMax_spread`).
#' Prediction rescales with the same values, so it never depends on `sim$covMinMax_spread`, which
#' `fireSense_spreadFit` rebuilds in every run from that run's own covariate tables (in a hindcast, from a
#' formula rebuilt for another year, which can lack a fitted term).
#'
#' @param sa `sim$studyAreaWithSpreadParams`, one row per ELF.
#' @param i integer; the row.
#' @param default what to return for a row without ranges (an older ledger row): `sim$covMinMax_spread`.
#' @return The row's `covMinMax_spread` (a `data.table`, 2 rows: min, max), else `default`.
#' @export
ledgerCovMinMax <- function(sa, i, default = NULL) {
  if (!"covMinMax_spread" %in% names(sa)) return(default)
  cmm <- sa[["covMinMax_spread"]][[i]]
  if (is.null(cmm) || !NROW(cmm)) return(default)
  cmm
}

#' The climate reference of one ledger row
#'
#' A fit made with climate outside the `maxAsymptote` ceiling (`fireSense_spreadFit`'s `spreadLink = "climateOutside"`)
#' stores the reference value of each climate covariate in the ledger row (`spreadFitClimateRefTxt`),
#' with `spreadLink` beside it. Prediction applies the same link with the same reference
#' (`spreadProbFromCovs()`).
#'
#' @param sa `sim$studyAreaWithSpreadParams`, one row per ELF.
#' @param i integer; the row.
#' @return A named list, or `NULL` for a row without a reference: a fit with climate under the ceiling, or a row
#'   written before the link existed. Such a row predicts exactly as it always did.
#' @export
ledgerClimateRef <- function(sa, i) {
  col <- spreadFitClimateRefTxt
  if (!col %in% names(sa)) return(NULL)
  ref <- sa[[col]][[i]]
  if (is.list(ref) && length(ref)) ref else NULL
}

#' The climate hinge knots of one ledger row
#'
#' A fit made with `fireSense_dataPrepFit`'s `spreadClimateHinge` has a hinge term (`climateHinge`) for each climate
#' variable, at the ELF's median of the variable over the fit's burnable pixels. The fit stores those knots in the ledger
#' row (`spreadFitClimateHingeKnotsTxt`). Prediction builds the terms from them (`addClimateHinge()`).
#'
#' @param sa `sim$studyAreaWithSpreadParams`, one row per ELF.
#' @param i integer; the row.
#' @return A named list (variable to knot), or `NULL` for a row without knots: a fit without hinge terms, or a row
#'   written before they existed. Such a row predicts exactly as it always did.
#' @export
ledgerClimateHingeKnots <- function(sa, i) {
  col <- spreadFitClimateHingeKnotsTxt
  if (!col %in% names(sa)) return(NULL)
  knots <- sa[[col]][[i]]
  if (is.list(knots) && length(knots)) knots else NULL
}

#' The runaway threshold of one ledger row
#'
#' A fit made with `fireSense_spreadFit`'s `spreadClimate = "ceiling"` moves `maxAsymptote` toward the landscape's
#' runaway threshold `p50`, stored in the ledger row (`spreadFitPowerLawPc50Txt`) as for every power-law
#' fit. Only a row with `climateRef_spread` (`ledgerClimateRef()`) uses it.
#'
#' @param sa `sim$studyAreaWithSpreadParams`, one row per ELF.
#' @param i integer; the row.
#' @return `p50`, or `NULL` for a row without the column or with `NA`. A row with `climateRef_spread` and a power-law
#'   link but no finite `p50` stops, naming the row.
#' @export
ledgerPowerLawPc50 <- function(sa, i) {
  col <- spreadFitPowerLawPc50Txt
  p50 <- if (col %in% names(sa)) sa[[col]][[i]]
  if (!is.null(ledgerClimateRef(sa, i)) && any(ledgerSpreadLink(sa, i) %in% spreadPowerLawLinks) &&
      !(is.numeric(p50) && length(p50) == 1L && is.finite(p50)))
    stop("spread ledger row ", i, " has climate in the ceiling ('", spreadFitClimateRefTxt,
         "' and spreadLink '", ledgerSpreadLink(sa, i), "') but no finite '", col, "' (", format(p50), ")", call. = FALSE)
  if (is.null(p50) || is.na(p50)) NULL else p50
}

#' The power-law exponent of one ledger row
#'
#' A fit made with `fireSense_spreadFit`'s `spreadLink = "powerLaw"` or `"powerLawSmooth"` stores its `gamma` in the ledger row
#' (`spreadFitPowerLawGammaTxt`), with the link in `spreadLinkTxt`. Prediction
#' applies the same link with the same `gamma` (`spreadProbFromCovs()`).
#'
#' @param sa `sim$studyAreaWithSpreadParams`, one row per ELF.
#' @param i integer; the row.
#' @return `gamma`, or `NULL` for a row of any other link (`"richards"` has none), such as every row written before the
#'   link existed. Such a row predicts exactly as it always did. A `"powerLaw"` or `"powerLawSmooth"` row without a valid
#'   `gamma` stops.
#' @export
ledgerPowerLawGamma <- function(sa, i) {
  link <- ledgerSpreadLink(sa, i)
  if (!any(link %in% spreadPowerLawLinks)) return(NULL)
  gamma <- sa[[spreadFitPowerLawGammaTxt]][[i]]
  if (!(is.numeric(gamma) && length(gamma) == 1L && is.finite(gamma) && gamma > 0))
    stop("spread ledger row ", i, " has spreadLink '", link, "' but no valid '",
         spreadFitPowerLawGammaTxt, "' (", format(gamma), ")", call. = FALSE)
  gamma
}

#' The spread link of one ledger row
#'
#' @param sa `sim$studyAreaWithSpreadParams`, one row per ELF.
#' @param i integer; the row.
#' @return The row's `spreadLinkTxt`, or `NULL` for a row without one (the default link, and every row
#'   written before the column existed).
#' @export
ledgerSpreadLink <- function(sa, i) {
  linkCol <- spreadLinkTxt
  if (linkCol %in% names(sa)) as.character(sa[[linkCol]][[i]])
}

#' A number of one ledger row, for the links that record it
#'
#' The row's `lowerAsymptote_spread` and `richardsNu_spread` are checked together by
#' `checkLinkParams()`, the check `fireSense_spreadFit` makes of its parameters; an invalid or missing
#' value stops, naming the row.
#'
#' @param sa `sim$studyAreaWithSpreadParams`, one row per ELF.
#' @param i integer; the row.
#' @param col the ledger column.
#' @param links the links that write `col`; any other row gives `NULL`.
#' @param floor the run's `lowerSpreadProb`, which an activation row's lower asymptote must equal.
#' @return The value, or `NULL` for a row of another link.
#' @export
ledgerLinkNumber <- function(sa, i, col, links, floor = spreadProbFloor) {
  link <- ledgerSpreadLink(sa, i)
  if (!any(link %in% links)) return(NULL)
  num <- function(cn) if (cn %in% names(sa)) sa[[cn]][[i]]
  tryCatch(checkLinkParams(link, num(spreadFitLowerAsymptoteTxt),
                                           num(spreadFitRichardsNuTxt),
                                           activation = identical(ledgerSpreadClimate(sa, i), "activation"), floor = floor),
           error = function(e) stop("spread ledger row ", i, " has spreadLink '", link, "': ",
                                    conditionMessage(e), " (columns '", spreadFitLowerAsymptoteTxt, "', '",
                                    spreadFitRichardsNuTxt, "')", call. = FALSE))
  num(col)
}

#' The lower asymptote of one ledger row
#'
#' `fireSense_spreadFit`'s `lowerAsymptote` of a `"powerLawSmooth"` or `"richards"` fit, in
#' `spreadFitLowerAsymptoteTxt`.
#'
#' @inheritParams ledgerLinkNumber
#' @return A number from 0 to `spreadProbFloor`; `NULL` for a row of another link.
#' @export
ledgerLowerAsymptote <- function(sa, i, floor = spreadProbFloor)
  ledgerLinkNumber(sa, i, spreadFitLowerAsymptoteTxt, spreadAsymptoteLinks, floor)

#' The `nu` of one ledger row
#'
#' `fireSense_spreadFit`'s `richardsNu` of a `"richards"` fit, in `spreadFitRichardsNuTxt`.
#'
#' @inheritParams ledgerLinkNumber
#' @return A number, 0 or more; `NULL` for a row of another link.
#' @export
ledgerRichardsNu <- function(sa, i, floor = spreadProbFloor)
  ledgerLinkNumber(sa, i, spreadFitRichardsNuTxt, "richards", floor)

#' The link of one ledger row and its parameters, as the arguments of `spreadProbOneELF()`
#'
#' @inheritParams ledgerLinkNumber
#' @return A list of `link` (only for a link of `spreadLinearLinks`), `powerLawGamma` (not `"richards"`),
#'   `powerLawShape`, `lowerAsymptote` and `richardsNu` (the links that have them), and `spreadClimate` when the row
#'   records it; the others are left out, so they take `spreadProbOneELF()`'s defaults. A row of any other link,
#'   or without one, has `powerLawShape` only.
#' @export
ledgerLinkArgs <- function(sa, i, floor = spreadProbFloor) {
  link <- ledgerSpreadLink(sa, i)
  if (!any(link %in% spreadLinearLinks)) link <- NULL
  col <- spreadClimateTxt
  recorded <- col %in% names(sa) && !is.na(sa[[col]][[i]])
  Filter(Negate(is.null), list(link = link, powerLawGamma = ledgerPowerLawGamma(sa, i),
                               powerLawShape = ledgerPowerLawShape(sa, i),
                               lowerAsymptote = ledgerLowerAsymptote(sa, i, floor), richardsNu = ledgerRichardsNu(sa, i, floor),
                               spreadClimate = if (recorded) ledgerSpreadClimate(sa, i)))
}

#' Where the fit of one ledger row put climate
#'
#' `fireSense_spreadFit` records its `spreadClimate` in the row (`spreadClimateTxt`). A row without it
#' (written before it existed) is read as the fit made it: `"ceiling"` when it has `climateRef_spread`
#' (`ledgerClimateRef()`) and a power-law link, else `"predictor"`.
#'
#' @inheritParams ledgerLinkNumber
#' @return One of `spreadClimateOptions`; an unknown value stops, naming the row.
#' @export
ledgerSpreadClimate <- function(sa, i) {
  col <- spreadClimateTxt
  mode <- if (col %in% names(sa)) as.character(sa[[col]][[i]])
  if (!length(mode) || is.na(mode))
    return(if (!is.null(ledgerClimateRef(sa, i)) && any(ledgerSpreadLink(sa, i) %in% spreadPowerLawLinks))
      "ceiling" else "predictor")
  if (!mode %in% spreadClimateOptions)
    stop("spread ledger row ", i, " has an unknown '", col, "' (", mode, "); one of ",
         paste(spreadClimateOptions, collapse = ", "), " expected", call. = FALSE)
  mode
}

#' The power-law curvature exponent of one ledger row
#'
#' `fireSense_spreadFit`'s `powerLawShape` is stored in the ledger row (`spreadFitPowerLawShapeTxt`).
#' A row without the column (every row written before it existed) or with `NA` has shape 1, the original link.
#'
#' @param sa `sim$studyAreaWithSpreadParams`, one row per ELF.
#' @param i integer; the row.
#' @return `shape`, a positive number; 1 when the row has none.
#' @export
ledgerPowerLawShape <- function(sa, i) {
  col <- spreadFitPowerLawShapeTxt
  shape <- if (col %in% names(sa)) sa[[col]][[i]]
  if (is.null(shape) || is.na(shape)) return(spreadPowerLawShape)
  if (!(is.numeric(shape) && length(shape) == 1L && is.finite(shape) && shape > 0))
    stop("spread ledger row ", i, " has an invalid '", col, "' (", format(shape), ")", call. = FALSE)
  shape
}

#' The covariate centres of one ledger row
#'
#' A fit with an intercept centres its covariates and stores the centres in the ledger row
#' (`spreadFitCovCentreTxt`). Prediction subtracts the same values.
#'
#' @param sa `sim$studyAreaWithSpreadParams`, one row per ELF.
#' @param i integer; the row.
#' @return A named list, or `NULL` for a row without centres: a fit without an intercept, or a row written
#'   before they were stored. Such a row predicts exactly as it always did.
#' @export
ledgerCovCentre <- function(sa, i) {
  col <- spreadFitCovCentreTxt
  if (!col %in% names(sa)) return(NULL)
  centre <- sa[[col]][[i]]
  if (is.list(centre) && length(centre)) centre else NULL
}
