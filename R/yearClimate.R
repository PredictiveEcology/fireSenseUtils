#' The year-level climate term of the expected ignitions
#'
#' `fireSense_ignitionFit` can scale each year's expected ignitions by `exp(beta * z_y)`, with `beta >= 0`
#' and `z_y` the standardised yearly mean of one climate covariate over the ELF (by default summer CMD,
#' `"CMDsm"`). The fit stores `list(variable, center, scale, beta)` (`yearClimate`), with
#' `center` and `scale` the mean and sd of the yearly means over the fit years, and
#' `fireSense_ignitionPredict` applies it. These functions are the one place the formula lives.
#'
#' @param variable Character: a covariate name, with or without underscores (`"CMD_sm"` is `"CMDsm"`).
#' @param columns Character: the column names to search.
#' @param x Numeric: the covariate, one element per row.
#' @param year The year of each row.
#' @param weight Numeric: the weight of each row (the share of the pixel in the core,
#'   `fireSense_ignitionFit`'s `exposureOf()`); 1 for a plain mean.
#' @param yearClimate `list(variable, center, scale, beta)`, or `NULL` (no term).
#' @param yearMean The yearly mean(s) of `yearClimate$variable`, on the covariate's raw scale.
#'
#' @return `yearClimateColumn()`: the element of `columns` that is `variable`; an error if there is none.
#'   `yearClimateMeans()`: a `data.frame` with `year` and `mean`, one row per year, the weighted mean of `x`.
#'   `yearClimateMultiplier()`: `exp(beta * (yearMean - center) / scale)`, 1 when `yearClimate` is `NULL`.
#' @export
#' @rdname yearClimate
#' @examples
#' yearClimateColumn("CMD_sm", c("CMDsm", "youngAge"))
#' yearClimateMeans(c(1, 3, 10, 20), year = c(1, 1, 2, 2), weight = c(1, 1, 1, 0))
#' yearClimateMultiplier(list(variable = "CMDsm", center = 5, scale = 2, beta = 0.5), c(3, 5, 7))
#' yearClimateMultiplier(NULL, 7)
yearClimateColumn <- function(variable, columns) {
  out <- columns[match(gsub("_", "", variable), gsub("_", "", columns))]
  if (is.na(out))
    stop("The year climate variable '", variable, "' is not among the covariates (", paste(columns, collapse = ", "), ").",
         call. = FALSE)
  out
}

#' @export
#' @rdname yearClimate
yearClimateMeans <- function(x, year, weight = 1) {
  weight <- rep_len(weight, length(x))
  tot <- rowsum(cbind(wx = weight * x, w = weight), group = year)
  data.frame(year = as.numeric(rownames(tot)), mean = unname(tot[, "wx"] / tot[, "w"]), row.names = NULL)
}

#' @export
#' @rdname yearClimate
yearClimateMultiplier <- function(yearClimate, yearMean) {
  if (is.null(yearClimate)) return(rep(1, length(yearMean)))
  exp(yearClimate$beta * (yearMean - yearClimate$center) / yearClimate$scale)
}
