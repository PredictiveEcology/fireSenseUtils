#' The year-level climate term of the expected ignitions
#'
#' `fireSense_ignitionFit` can scale each year's expected ignitions by `exp(beta * z_y) / norm`, with `beta >= 0`
#' and `z_y` the standardised yearly mean of one climate covariate over the ELF (by default summer CMD,
#' `"CMDsm"`). The fit stores `list(variable, center, scale, beta, norm)` (`yearClimate`): `center` and `scale` are
#' the mean and sd of the yearly means over the fit years, and `norm` is the mean of `exp(beta * z_y)` over the fit
#' years, so the multiplier averages 1 over them. `fireSense_ignitionPredict` applies it. These functions are the
#' one place the definition lives: the yearly mean is the same unweighted mean over the same rows in the fit and in
#' the prediction.
#'
#' @param variable Character: a covariate name, with or without underscores (`"CMD_sm"` is `"CMDsm"`).
#' @param columns Character: the column names to search.
#' @param covariates A `data.frame` of covariates, one row per pixel (and year): the raw (not standardised)
#'   covariates of the fit, or of the pixels to predict. Rows with an `NA` in any column are dropped, and, when
#'   there is a `coreFraction` column ([coreFractionTxt]), rows without core (`coreFraction` not above 0). A `year`
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
  x <- x[stats::complete.cases(x), , drop = FALSE]
  if (!is.null(x[[coreFractionTxt]])) x <- x[x[[coreFractionTxt]] > 0, , drop = FALSE]
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
