#' The year-level climate term of the expected ignitions and escapes
#'
#' `fireSense_ignitionFit` can scale each year's expected ignitions, and separately its expected escapes, by `exp(beta * z_y) / norm`, with `beta >= 0`
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
#'   With `shape = "logistic"` (absent: the log-linear term above, every fit made before) it is
#'   `list(variable, center, scale, shape, k, ratio, norm)`: the multiplier is
#'   `(1 + (ratio - 1) * plogis(k * z)) / norm`, an S-curve through its midpoint at `z = 0` (the ELF's average
#'   year) from a floor of 1 to an upper level of `ratio` (`ratio >= 1`, `k >= 0`), which saturates beyond the fit range.
#' @param yearMean The yearly mean(s) of `yearClimate$variable`, on the covariate's raw scale.
#'
#' @return `yearClimateColumn()`: the element of `columns` that is `variable`; an error if there is none.
#'   `yearClimateMeans()`: a `data.frame` with `year` and `mean`, one row per year that has rows left.
#'   `yearClimateMultiplier()`: `exp(beta * z) / norm` (log-linear) or `(1 + (ratio - 1) * plogis(k * z)) / norm`
#'   (logistic) with `z = (yearMean - center) / scale`, 1 when `yearClimate` is `NULL`.
#' @export
#' @rdname yearClimate
#' @examples
#' yearClimateColumn("CMD_sm", c("CMDsm", "youngAge"))
#' yearClimateMeans(data.frame(year = c(1, 1, 2, 2), CMDsm = c(1, 3, 10, NA)), "CMDsm")
#' yc <- list(variable = "CMDsm", center = 5, scale = 2, beta = 0.5, norm = 1)
#' yearClimateMultiplier(yc, c(3, 5, 7))
#' yearClimateMultiplier(NULL, 7)
#' yearClimateMultiplier(list(variable = "CMDsm", center = 5, scale = 2, shape = "logistic", k = 2, ratio = 4,
#'                             norm = 1), c(3, 5, 7, 50))  # saturates
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
  x <- yearMean - yearClimate$center
  if (identical(yearClimate$shape, "logistic"))
    return((1 + (yearClimate$ratio - 1) * stats::plogis(yearClimate$k * x / yearClimate$scale)) / norm)
  exp(yearClimate$beta * x / yearClimate$scale) / norm   # a term without `shape`: the arithmetic of old, to the last bit
}
