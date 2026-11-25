#' The year-level climate term of the expected ignitions
#'
#' `fireSense_ignitionFit` can scale each year's expected ignitions by `exp(beta * z_y) / norm`, with `beta >= 0`
#' and `z_y` the standardised yearly mean of one climate covariate over the ELF (by default summer CMD,
#' `"CMDsm"`). The fit stores `list(variable, center, scale, beta, norm)` (`yearClimate`): `center` and `scale` are
#' the mean and sd of the yearly means over the fit years, and `norm` is the mean of `exp(beta * z_y)` over the fit
#' years, so the multiplier averages 1 over them. `fireSense_ignitionPredict` applies it. These functions are the
#' one place the definition lives: the yearly mean is the same unweighted mean over the same rows (core pixels, with
#' the variable) in the fit and in the prediction, provided both covariate tables carry `coreFraction` from [addCoreFraction()].
#'
#' @param variable Character: a covariate name, with or without underscores (`"CMD_sm"` is `"CMDsm"`).
#' @param columns Character: the column names to search.
#' @param covariates A `data.frame` of covariates, one row per pixel (and year): the raw (not standardised)
#'   covariates of the fit, or of the pixels to predict, both with a `coreFraction` column ([coreFractionTxt],
#'   made by [addCoreFraction()]). A row counts if `variable` and `coreFraction` are not `NA` and `coreFraction` is above 0;
#'   no other column decides. Without a `coreFraction` column every row with `variable` counts. A `year`
#'   column ([yearTxt]) gives the years; without one all rows are one year.
#' @param yearClimate `list(variable, center, scale, beta, norm)`, or `NULL` (no term). `norm` may be absent (1).
#' @param yearMean The yearly mean(s) of `yearClimate$variable`, on the covariate's raw scale.
#'
#' @return `yearClimateColumn()`: the element of `columns` that is `variable`; an error if there is none.
#'   `yearClimateMeans()`: a `data.frame` with `year` and `mean`, one row per year that has rows left.
#'   `yearClimateMultiplier()`: `exp(beta * (yearMean - center) / scale) / norm`, 1 when `yearClimate` is `NULL`.
#' @export
#' @rdname yearClimate
#' @examples
#' yearClimateColumn("CMD_sm", c("CMDsm", "youngAge"))
#' yearClimateMeans(data.frame(year = c(1, 1, 2, 2), CMDsm = c(1, 3, 10, NA)), "CMDsm")
#' yc <- list(variable = "CMDsm", center = 5, scale = 2, beta = 0.5, norm = 1)
#' yearClimateMultiplier(yc, c(3, 5, 7))
#' yearClimateMultiplier(NULL, 7)
yearClimateColumn <- function(variable, columns) {
  out <- matchIgnoringUnderscores(variable, columns)
  if (is.na(out))
    stop("The year climate variable '", variable, "' is not among the covariates (", paste(columns, collapse = ", "), ").",
         call. = FALSE)
  out
}

#' @export
#' @rdname yearClimate
yearClimateMeans <- function(covariates, variable) {
  x <- as.data.frame(covariates)
  col <- yearClimateColumn(variable, names(x))
  keep <- !is.na(x[[col]])
  if (!is.null(x[[coreFractionTxt]])) keep <- keep & !is.na(x[[coreFractionTxt]]) & x[[coreFractionTxt]] > 0
  x <- x[keep, , drop = FALSE]
  year <- if (is.null(x[[yearTxt]])) rep(1, nrow(x)) else x[[yearTxt]]
  if (!nrow(x)) return(data.frame(year = numeric(0), mean = numeric(0)))
  tot <- rowsum(cbind(x = x[[col]], n = 1), group = year)
  data.frame(year = as.numeric(rownames(tot)), mean = unname(tot[, "x"] / tot[, "n"]), row.names = NULL)
}

#' @export
#' @rdname yearClimate
yearClimateMultiplier <- function(yearClimate, yearMean) {
  if (is.null(yearClimate)) return(rep(1, length(yearMean)))
  norm <- if (is.null(yearClimate$norm)) 1 else yearClimate$norm
  exp(yearClimate$beta * (yearMean - yearClimate$center) / yearClimate$scale) / norm
}

#' The standardised climate of each fit year of a spread fit
#'
#' The yearly mean of the formula's first annual climate covariate ([spreadClimateCols()]) over the pixels of
#' `annualDTx1000` ([yearClimateMeans()]), standardised over the fit years: `(mean - mean(means)) / sd(means)`.
#' The same year-level climate as the ignition year-climate term, for the spread fit's own climate covariate.
#'
#' @param annualDTx1000 Named list of the annual covariate tables, one per year ([.objfunSpreadFit()]).
#' @param formulaToFit The spread formula ([.objfunSpreadFit()]).
#' @param years The fit years (character); default all of `annualDTx1000`.
#' @return A named numeric (by year), with mean 0 and sd 1 (all 0 if the years do not differ); `NULL` if the
#'   formula has no annual climate covariate.
#' @export
spreadYearClimate <- function(annualDTx1000, formulaToFit, years = names(annualDTx1000)) {
  climCols <- spreadClimateCols(formulaToFit, annualDTx1000)
  if (!length(climCols)) return(NULL)
  ## years are labelled by position, so any year names work (yearClimateMeans() reads the year as a number)
  yrDT <- data.table::rbindlist(stats::setNames(annualDTx1000[years], seq_along(years)), idcol = yearTxt,
                                use.names = TRUE, fill = TRUE)
  yrDT <- yrDT[, c(yearTxt, climCols[[1]]), with = FALSE]
  m <- yearClimateMeans(yrDT, climCols[[1]])
  mean <- stats::setNames(m$mean, years[m$year])[years]
  sdm <- stats::sd(mean)
  stats::setNames(if (is.finite(sdm) && sdm > 0) (mean - base::mean(mean)) / sdm else rep(0, length(mean)), years)
}
