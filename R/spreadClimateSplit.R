## The spread climate covariate split into a between-year and a within-year term
## (`spreadClimateSplit` of fireSense_dataPrepFit). One place for the names, the ELF median, the split and the
## ranges, used by the fit (fireSense_dataPrepFit, fireSense_spreadFit, spreadFitLedgerInputs()) and by the
## prediction (fireSense_spreadPredict), so a pixel's two terms are the same numbers in both.

#' The spread climate covariate split into a between-year and a within-year term
#'
#' With `spreadClimateSplit` the spread formula has two terms for a climate variable `C` instead of one:
#' `<C>_yr`, the year's ELF median of `C` over the ELF's burnable pixels (one value for every pixel of
#' the year: it carries the year-to-year swings), and `<C>_dev`, `C` minus that median (each pixel's
#' departure from its year's median: it carries the spatial gradient within a year). Each has its own coefficient.
#'
#' The median is `stats::median(na.rm = TRUE)` of the climate layer's values over the ELF's burnable pixels,
#' before the x1000 integer storage. The fit and the prediction both use `splitSpreadClimate()` with the climate
#' of those pixels as `ref`, so the terms agree.
#'
#' The terms keep the sign rule of their variable, which `fireSense_spreadFit` applies by name: `CMD_yr` and `CMD_dev`
#' (names with `CMD` or `MDC`) are bounded below by 0; the others, e.g. `Tmaxsm_yr`, are not constrained.
#' Both are rescaled with fixed ranges ([climateSplitRanges()]).
#'
#' @param vars character; climate variables, as named in the annual tables (`"CMD"`, `"Tmaxsm"`).
#' @return `climateSplitCols()`: the two terms of each variable, `<var>_yr` then `<var>_dev`, variable by variable.
#' @export
#' @rdname spreadClimateSplit
climateSplitCols <- function(vars) as.vector(rbind(paste0(vars, climateSplitYrSuffix),
                                                   paste0(vars, climateSplitDevSuffix)))

#' @export
#' @rdname spreadClimateSplit
climateSplitYrSuffix <- "_yr"

#' @export
#' @rdname spreadClimateSplit
climateSplitDevSuffix <- "_dev"

#' @param x character; column or parameter names.
#' @return `isClimateSplitCol()`: logical, `TRUE` for the names that end in `_yr` or `_dev`.
#'   `climateSplitBase()`: `x` with those endings removed, i.e. the climate variable of a term.
#' @export
#' @rdname spreadClimateSplit
isClimateSplitCol <- function(x) grepl(paste0("(", climateSplitYrSuffix, "|", climateSplitDevSuffix, ")$"), x)

#' @export
#' @rdname spreadClimateSplit
climateSplitBase <- function(x) sub(paste0("(", climateSplitYrSuffix, "|", climateSplitDevSuffix, ")$"), "", x)

#' @param dt `data.table` with the climate variables `vars` and the rows to give the two terms to.
#'   Changed by reference: each variable is replaced by its `_yr` and `_dev` columns.
#' @param ref `data.table` with the climate variables `vars` (and `by`) of the ELF's burnable pixels, the rows
#'   the median is taken over.
#' @param by character; the columns of `dt` and `ref` that define a year (`yearTxt` in the fit's tables of
#'   all years); `NULL` when both are one year (a prediction).
#' @return `splitSpreadClimate()`: `dt`, invisibly. A year with no value in `ref` gives `NA`.
#' @export
#' @rdname spreadClimateSplit
splitSpreadClimate <- function(dt, vars, ref, by = NULL) {
  yr <- paste0(vars, climateSplitYrSuffix)
  dev <- paste0(vars, climateSplitDevSuffix)
  med <- ref[, lapply(.SD, stats::median, na.rm = TRUE), by = by, .SDcols = vars]
  data.table::setnames(med, vars, yr)
  if (length(by)) {
    dt[med, (yr) := mget(paste0("i.", yr)), on = by]
  } else {
    for (i in seq_along(vars)) data.table::set(dt, NULL, yr[i], med[[yr[i]]])
  }
  for (i in seq_along(vars)) {
    data.table::set(dt, NULL, dev[i], dt[[vars[i]]] - dt[[yr[i]]])
    data.table::set(dt, NULL, vars[i], NULL)
  }
  invisible(dt)
}

#' Fixed ranges of the split climate terms
#'
#' `<var>_yr` is rescaled with the range of `<var>`, `<var>_dev` with a range of the same width centred on 0
#' (`c(-w / 2, w / 2)` for a range of width `w`), so a unit of either term is the same amount of the climate
#' variable and the two coefficients compare. A range already in `ranges` for a split name is kept.
#'
#' @param ranges named list of `c(min, max)`, as `climateCovRanges`.
#' @return `ranges` with the two split terms of each of its variables added.
#' @export
climateSplitRanges <- function(ranges) {
  if (!length(ranges)) return(ranges)
  derived <- c(stats::setNames(ranges, paste0(names(ranges), climateSplitYrSuffix)),
               stats::setNames(lapply(ranges, function(r) c(-1, 1) * diff(r) / 2),
                               paste0(names(ranges), climateSplitDevSuffix)))
  c(ranges, derived[setdiff(names(derived), names(ranges))])
}
