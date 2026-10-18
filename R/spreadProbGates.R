## The spreadProb gates of `.objfunSpreadFit()`, as functions of their own, so that the objective
## and the threshold calibration (fireSense_spreadFit) apply the same test.

#' The objective's spread-probability gates, on one year's spreadProb
#'
#' `.objfunSpreadFit()` refuses to simulate a year (it scores the `minLik` floor, and the parameter set
#' is "bailed") unless its spreadProb passes these checks. The objective calls this function, so the
#' checks cannot differ from [spreadProbGates()].
#'
#' @param spreadProb numeric; the spreadProb of each pixel of the year.
#' @param ceiling numeric; the logistic's asymptote (its first parameter). Pixels at or above 99% of
#'   it (or within 2.5% of `lowerSpreadProb`) are "edge" values and are left out of the spread checks.
#' @param lowerSpreadProb,maxFireSpread numeric; the median must be within them.
#' @param lanscape1stQuantileThresh numeric; the first quartile of the non-edge values must be below this.
#' @return a list: `pass`; `burny` (first quartile too high, "Too burny a landscape"); `notSpread`
#'   ((q90 - q10) / median <= 0.025, "Not spread out enough"); `medianOK`; and the `nonEdgeValues` and
#'   their `summ`ary that the objective logs.
#' @export
spreadProbGateTest <- function(spreadProb, ceiling, lowerSpreadProb, maxFireSpread,
                               lanscape1stQuantileThresh) {
  medSP <- median(spreadProb, na.rm = TRUE)
  ## Taken from the spreadProb column, not by scanning `cells`. `cells` is landscape-length and
  ## zero everywhere except this year's pixels (6.0M cells vs a median 41k pixels on ELF 5.3.1),
  ## so `cells[cells > a | cells > b]` made four passes over the landscape to recover values that
  ## are already here. Same values: pixelID is unique within a year, the zeros never pass a
  ## non-negative threshold, and `x > a | x > b` is `x > min(a, b)`. Only quantile() and summary()
  ## read it, so order does not matter. Measured on ELF 5.3.1 with identical seeds: identical
  ## objective values, and an evaluation 1.1-1.4x faster.
  nonEdgeValues <- spreadProb[spreadProb > min(lowerSpreadProb * 1.025, ceiling * 0.99)]
  sdSP <- diff(quantile(nonEdgeValues, c(0.1, 0.9)))
  if (is.na(sdSP)) sdSP <- 0
  medianOK <- medSP <= maxFireSpread & medSP >= lowerSpreadProb
  spreadOutEnough <- sdSP / medSP > 0.025
  summ <- summary(nonEdgeValues)
  lowSPLowEnough <- summ[2] < lanscape1stQuantileThresh
  list(pass = isTRUE(medianOK && spreadOutEnough && lowSPLowEnough),
       burny = isTRUE(!lowSPLowEnough), notSpread = isTRUE(!spreadOutEnough),
       medianOK = isTRUE(medianOK), nonEdgeValues = nonEdgeValues, summ = summ)
}

#' A year's covariates, rescaled as the objective does, before the logistic
#' @keywords internal
spreadProbCovsOfYear <- function(yr, annDTx1000, nonAnnualDTx1000, indexNonAnnual, covMinMax,
                                 mutuallyExclusive, colsToUse, doAssertions, logisticPars, covPars,
                                 maxFireSpread, lowerSpreadProb, covCentre = NULL) {
  spreadProbFromIntegerCovs(
    shortAnnDTx1000 = NULL, annDTx1000, nonAnnualDTx1000,
    indexNonAnnual, yr, covMinMax, mutuallyExclusive, colsToUse,
    doAssertions, logisticPars, covPars, maxFireSpread, lowerSpreadProb, covCentre = covCentre
  )
}

#' Which non-annual table (`nonAnnualDTx1000`) belongs to which years
#'
#' @param nonAnnualDTx1000 list of `data.table`s named by the years they cover (`"1985_1986"`).
#' @return a `data.table`: `ind`, the position in the list, and `date`, the years it covers.
#' @keywords internal
nonAnnualIndex <- function(nonAnnualDTx1000) {
  yearSplit <- strsplit(names(nonAnnualDTx1000), "_")
  data.table::rbindlist(Map(
    ind = seq_along(nonAnnualDTx1000), date = yearSplit,
    function(ind, date) data.table::data.table(ind = ind, date = date)))
}

#' The centre of each covariate of a spread fit: its mean over the data the objective uses
#'
#' With an intercept in the model the covariates are centred (`covCentre` of [.objfunSpreadFit()]), so the
#' coefficients describe variation and the intercept the level. The centre is the mean of each covariate
#' of the formula, over every pixel-year in `annualDTx1000`, after the rescaling and the mutual
#' exclusivity the objective applies and before it centres: the very values it goes on to centre.
#'
#' @param annualDTx1000,nonAnnualDTx1000,covMinMax,mutuallyExclusive As in [.objfunSpreadFit()].
#' @param formulaToFit The spread formula; the intercept is not a covariate and has no centre.
#' @return `NULL` when the formula has no intercept (nothing is centred), else a named list, one mean per
#'   covariate, to give as `covCentre`.
#' @export
spreadCovCentre <- function(annualDTx1000, nonAnnualDTx1000, formulaToFit, covMinMax = NULL,
                            mutuallyExclusive = list("youngAge" = c("class", "nf"))) {
  colsToUse <- spreadDesignCols(formulaToFit)
  if (!spreadInterceptTxt %in% colsToUse) return(NULL)
  spreadCovMeans(annualDTx1000, nonAnnualDTx1000, spreadCovCols(colsToUse), covMinMax, mutuallyExclusive)
}

## The mean of each of `covCols` over every pixel-year of `annualDTx1000`, on the scale the objective's
## design matrix has them: rescaled, made mutually exclusive and, with `covCentre`, centred. Used by
## spreadCovCentre() (no centre: the values it goes on to centre) and spreadClimateRef() (the fit's centre, so
## the reference is on the design scale).
spreadCovMeans <- function(annualDTx1000, nonAnnualDTx1000, covCols, covMinMax = NULL, mutuallyExclusive = NULL,
                           covCentre = NULL) {
  lapply(nonAnnualDTx1000, data.table::setDT)
  indexNonAnnual <- nonAnnualIndex(nonAnnualDTx1000)
  sums <- numeric(length(covCols))
  n <- numeric(length(covCols))
  for (yr in names(annualDTx1000)) {
    covs <- spreadProbCovsOfYear(yr, annualDTx1000[[yr]], nonAnnualDTx1000, indexNonAnnual, covMinMax,
                                 mutuallyExclusive, covCols, FALSE, NULL, NULL, NULL, NULL, covCentre = covCentre)
    sums <- sums + vapply(covCols, function(cn) sum(covs[[cn]], na.rm = TRUE), numeric(1))
    n <- n + vapply(covCols, function(cn) sum(!is.na(covs[[cn]])), numeric(1))
  }
  as.list(stats::setNames(sums / n, covCols))
}

#' The climate covariates of a spread formula
#'
#' The covariates of the formula that are annual (columns of `annualDTx1000`) and not 0/1 indicators
#' (`youngAge`, `nfLCC_*`, `treedWetland`; see [spreadIndicatorRanges()]): the rule `fireSense_spreadFit`
#' uses to require a fixed range for a climate covariate.
#'
#' @param formulaToFit The spread formula, a character string or a formula.
#' @param annualDTx1000 The annual covariate tables, as in [.objfunSpreadFit()]; only the column names of
#'   the first are used.
#' @return A character vector; empty when the formula has no climate covariate.
#' @export
spreadClimateCols <- function(formulaToFit, annualDTx1000) {
  covCols <- spreadCovCols(spreadDesignCols(formulaToFit))
  annual <- intersect(covCols, names(annualDTx1000[[1]]))
  setdiff(annual, spreadIndicatorCols(annual))
}

#' The climate reference of a spread fit with climate outside the ceiling
#'
#' With climate outside the ceiling (`climateRef` of [spreadProbFromCovs()]; `spreadLink = "climateOutside"`
#' in `fireSense_spreadFit`) the climate term is `beta_clim * (climate - climateRef)`, added on the logit of the
#' spread probability after the fuel ceiling. `climateRef` is the fit's typical climate: the mean of each
#' climate covariate ([spreadClimateCols()]) over every pixel-year in `annualDTx1000`, on the scale of the
#' objective's design matrix (rescaled with `covMinMax`, and centred with `covCentre` when the fit has an
#' intercept, so then it is 0 up to rounding). It is found once, from the data the objective uses, stored in
#' the ledger row (`spreadFitClimateRefTxt`), and `fireSense_spreadPredict` uses the same value, so the fit
#' and the prediction cannot drift apart.
#'
#' @inheritParams spreadCovCentre
#' @param covCentre The fit's covariate centres ([spreadCovCentre()]); `NULL` without an intercept.
#' @return A named list, one mean per climate covariate, to give as `climateRef`; `NULL` when the formula has
#'   no climate covariate.
#' @export
spreadClimateRef <- function(annualDTx1000, nonAnnualDTx1000, formulaToFit, covMinMax = NULL,
                             mutuallyExclusive = list("youngAge" = c("class", "nf")), covCentre = NULL) {
  climCols <- spreadClimateCols(formulaToFit, annualDTx1000)
  if (!length(climCols)) return(NULL)
  spreadCovMeans(annualDTx1000, nonAnnualDTx1000, climCols, covMinMax, mutuallyExclusive, covCentre = covCentre)
}

#' The spread probability of each pixel from its covariates: the one place the spread link is computed
#'
#' Every spread probability the model gives -- in the objective ([.objfunSpreadFit()]), its gates
#' ([spreadProbGates()]), the validation ([spreadFitValidationData()], [plotSpreadFitResponse()]) and
#' the prediction (`fireSense_spreadPredict`) -- comes from this function, so the fit and the prediction
#' use one formula.
#'
#' Two links, chosen by `climateRef`:
#'
#' * `climateRef = NULL` (the default; `spreadLink = "ceiling"` in `fireSense_spreadFit`): the logistic of
#'   the whole linear predictor, `p = f + (maxAsymptote - f) * plogis(x)` with `x = mat %*% covPars`
#'   ([logisticAll()]; `f` is `lowerSpreadProb`). Climate is one of the covariates in `x`, so it acts under
#'   the `maxAsymptote` ceiling: where fuel alone puts a pixel near the ceiling a dry year adds almost
#'   nothing (in the 2026-10 hindcasts of three ELFs, 59-77% of pixels sat within 0.01 of the ceiling every
#'   year, and spread barely changed between years).
#' * `climateRef` a named numeric (`spreadLink = "climateOutside"`): the climate covariates (its names,
#'   columns of `mat`) are taken out of the ceiling. `p0 = f + (maxAsymptote - f) * plogis(x_fuel)`, with
#'   `x_fuel` the linear predictor of every other column (intercept, fuels, youngAge, non-forest, ...), is
#'   the spread probability at the reference climate; then
#'   `p = plogis(qlogis(p0) + sum(beta_clim * (climate - climateRef)))`, uncapped. The year effect
#'   `eps_Y` (`yearSpreadSD` in the objective, `fireSense_SpreadSD` in `fireSense_burn`) is added on the
#'   same logit afterwards, so climate and the year effect act together after the fuel ceiling. At
#'   `climate = climateRef` the two links give the same `p` for the same `p0`. The sign constraint on the
#'   climate coefficient (`>= 0`, from `fireSense_spreadFit`'s bounds) is unchanged.
#'
#' @param mat Numeric matrix, one row per pixel, one column per design column (with `colnames` when
#'   `climateRef` is given), as `spreadProbCovsOfYear()` gives them: `as.matrix(covs[, ..colsToUse])`.
#' @param logisticPars,covPars,lowerSpreadProb,link As for [logisticAll()]; `covPars` in the order of the
#'   columns of `mat`.
#' @param climateRef `NULL`, or a named numeric or list: the reference value of each climate column, on the
#'   scale of `mat` ([spreadClimateRef()]). Every name must be a column of `mat`.
#' @param climateAtRef Logical. With `climateRef`, `TRUE` returns `p0`, the spread probability with climate
#'   held at its reference: what the ceiling applies to. The objective's gates and its saturation summary
#'   look at `p0`, because under this link a dry year legitimately takes the landscape past `maxFireSpread`
#'   and a gate on `p` would refuse the years that carry the climate signal. Ignored without `climateRef`.
#' @return A one-column numeric matrix, one value per row of `mat` (as [logisticAll()]).
#' @export
spreadProbFromCovs <- function(mat, logisticPars, covPars, lowerSpreadProb, link = NULL, climateRef = NULL,
                               climateAtRef = FALSE) {
  if (is.null(climateRef) || !length(climateRef))
    return(logisticAll(logisticPars, mat = mat, covPars, lowerSpreadProb, link = link))
  climCols <- names(climateRef)
  if (is.null(climCols) || !all(climCols %in% colnames(mat)))
    stop("spreadProbFromCovs(): every name of `climateRef` must be a column of `mat`; climateRef: ",
         paste(climCols, collapse = ", "), "; columns: ", paste(colnames(mat), collapse = ", "))
  isClim <- colnames(mat) %in% climCols
  p0 <- logisticAll(logisticPars, mat = mat[, !isClim, drop = FALSE], covPars[!isClim], lowerSpreadProb,
                    link = link)
  if (isTRUE(climateAtRef)) return(p0)
  ref <- unlist(climateRef)[colnames(mat)[isClim]]
  shift <- (mat[, isClim, drop = FALSE] - rep(ref, each = nrow(mat))) %*% covPars[isClim]
  stats::plogis(stats::qlogis(p0) + shift)
}

#' Split `par` as the objective does: drop the trailing `yearSpreadSD` and insert the fixed `hillSlope1` and `inflectionPoint1`
#' @keywords internal
splitSpreadPar <- function(par, fitYearSpreadSD = NULL) {
  yearSpreadSD <- 0
  if (is.null(fitYearSpreadSD)) fitYearSpreadSD <- yearSpreadSDTxt %in% names(par)
  if (isTRUE(fitYearSpreadSD)) {
    if (!is.null(names(par)) && !identical(names(par)[length(par)], yearSpreadSDTxt))
      stop("`", yearSpreadSDTxt, "` must be the last parameter")
    yearSpreadSD <- unname(par[length(par)])
    par <- par[-length(par)]
  }
  ## hillSlope1 and inflectionPoint1 are fixed at 1, not fitted -- see fixLogisticPars() for why.
  list(par = fixLogisticPars(par), yearSpreadSD = yearSpreadSD)
}

#' The observed fires of each year, at or above `minFireSize`, one row per cell; years left with none
#' are dropped. Used by the objective and by [firstBlockYears()].
#' @keywords internal
firesAboveMinSize <- function(historicalFires, minFireSize) {
  x <- lapply(historicalFires, function(x) {
    x <- x[x$size >= minFireSize, ]
    x[!duplicated(x$cells), ]
  })
  ## can't fit fires for years with no data; drop these years
  omitYears <- names(x[which(lapply(x, nrow) == 0)])
  if (length(omitYears) > 0) x[omitYears] <- NULL
  x
}

#' The objective's first block of years: the two with the most burned area
#' @param historicalFiresAboveMin from `firesAboveMinSize()`
#' @keywords internal
firstBlockYears <- function(historicalFiresAboveMin) {
  fireSizesByYear <- unlist(lapply(historicalFiresAboveMin, function(x) sum(x$size)))
  names(head(sort(fireSizesByYear, decreasing = TRUE), 2))
}

## an escaped fire has at least `escapeMinPx` pixels, so the objective fits only fires that large
effectiveMinFireSize <- function(minFireSize, escapeMinPx) {
  if (is.null(escapeMinPx)) minFireSize else max(minFireSize, escapeMinPx)
}

#' Would parameter sets pass the objective's spreadProb gates in its first block of years?
#'
#' Computes each parameter set's spreadProb on the covariates of the first block (the two years with
#' the most burned area) exactly as [.objfunSpreadFit()] does, and applies [spreadProbGateTest()],
#' which the objective itself calls, without simulating any fire. A parameter set that fails is one
#' the objective would bail on in the first block. It is cheap enough to screen thousands of draws.
#'
#' @param par a named numeric vector as `.objfunSpreadFit()` takes (logistic parameters, then the
#'   covariate coefficients, then `yearSpreadSD` if fitted), or a list of them.
#' @param annualDTx1000,nonAnnualDTx1000,historicalFires,covMinMax,mutuallyExclusive,covCentre,link,climateRef,maxFireSpread,lowerSpreadProb,lanscape1stQuantileThresh,minFireSize,escapeSizeHa,landscape
#'   As in [.objfunSpreadFit()] (`landscape` is needed only with `escapeSizeHa`). With `climateRef` the gates
#'   test the spread probability at the reference climate, as the objective does ([spreadProbFromCovs()]).
#' @param formulaToFit character; the spread formula, as in the objective.
#' @param years character; the years to test. Default: the objective's first block.
#' @param doAssertions logical; as in the objective.
#' @return for one `par`, a list of logicals `pass` (every tested year passes), `burny` and `notSpread`
#'   (any tested year fails that check), `medianOK` (every year's median is in range) and the tested
#'   `years`; for a list of parameter sets, a `data.frame` with one row of the same per set.
#' @export
spreadProbGates <- function(par, annualDTx1000, nonAnnualDTx1000, historicalFires, formulaToFit,
                            covMinMax = NULL, mutuallyExclusive = list("youngAge" = c("class", "nf")),
                            covCentre = NULL, link = NULL, maxFireSpread = spreadProbCeiling, lowerSpreadProb = spreadProbFloor,
                            lanscape1stQuantileThresh = 0.265, minFireSize = 2, escapeSizeHa = NULL,
                            landscape = NULL, years = NULL, doAssertions = FALSE, climateRef = NULL) {
  data.table::setDTthreads(1)
  formulaToFit <- as.formula(formulaToFit, env = .GlobalEnv)
  colsToUse <- spreadDesignCols(formulaToFit)
  parsModel <- length(colsToUse)
  if (is.null(years)) {
    escapeMinPx <- if (!is.null(escapeSizeHa)) escapeSizePixels(escapeSizeHa, landscape)
    years <- firstBlockYears(firesAboveMinSize(historicalFires,
                                               effectiveMinFireSize(minFireSize, escapeMinPx)))
  }
  years <- as.character(years)
  lapply(nonAnnualDTx1000, data.table::setDT)
  indexNonAnnual <- nonAnnualIndex(nonAnnualDTx1000)
  single <- !is.list(par)
  pars <- if (single) list(par) else par
  pars <- lapply(pars, function(p) splitSpreadPar(p)$par)
  firstLogistic <- paramsSeparate(pars[[1]], parsModel)
  ## covariates do not depend on par: prepare each year once
  covs <- lapply(years, function(yr) as.matrix(
    spreadProbCovsOfYear(yr, annualDTx1000[[yr]], nonAnnualDTx1000, indexNonAnnual, covMinMax,
                         mutuallyExclusive, colsToUse, doAssertions, firstLogistic$logisticPars,
                         firstLogistic$covPars, maxFireSpread, lowerSpreadProb, covCentre = covCentre)[, ..colsToUse]))
  res <- lapply(pars, function(p) {
    ps <- paramsSeparate(p, parsModel)
    g <- lapply(covs, function(cv)
      spreadProbGateTest(spreadProbFromCovs(cv, ps$logisticPars, ps$covPars, lowerSpreadProb, link,
                                            climateRef = climateRef, climateAtRef = TRUE),
                         ps$logisticPars[1], lowerSpreadProb, maxFireSpread, lanscape1stQuantileThresh))
    list(pass = all(vapply(g, `[[`, logical(1), "pass")),
         burny = any(vapply(g, `[[`, logical(1), "burny")),
         notSpread = any(vapply(g, `[[`, logical(1), "notSpread")),
         medianOK = all(vapply(g, `[[`, logical(1), "medianOK")))
  })
  if (single) return(c(res[[1]], list(years = years)))
  data.frame(pass = vapply(res, `[[`, logical(1), "pass"), burny = vapply(res, `[[`, logical(1), "burny"),
             notSpread = vapply(res, `[[`, logical(1), "notSpread"),
             medianOK = vapply(res, `[[`, logical(1), "medianOK"))
}
