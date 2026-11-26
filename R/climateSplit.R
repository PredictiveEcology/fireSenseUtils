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
#' @rdname climateSplit
climateSplitCols <- function(vars) climateTermCols(vars, c(climateSplitYrSuffix, climateSplitDevSuffix))

#' @export
#' @rdname climateSplit
climateSplitYrSuffix <- "_yr"

#' @export
#' @rdname climateSplit
climateSplitDevSuffix <- "_dev"

#' @param x character; column or parameter names.
#' @return `isClimateSplitCol()`: logical, `TRUE` for the names that end in `_yr` or `_dev`.
#'   `climateSplitBase()`: `x` with those endings removed, i.e. the climate variable of a term.
#' @export
#' @rdname climateSplit
isClimateSplitCol <- function(x) isClimateTermCol(x, c(climateSplitYrSuffix, climateSplitDevSuffix))

#' @export
#' @rdname climateSplit
climateSplitBase <- function(x) climateTermBase(x, c(climateSplitYrSuffix, climateSplitDevSuffix))

#' The climate terms made from a climate variable by a suffix
#'
#' The one set of helpers for the terms of a climate variable `C` that are `C` and a suffix: `C_yr` and `C_dev`
#' ([climateSplitCols()]), `C_curve` ([climateCurveSuffix]). The split and curve helpers call them.
#'
#' @param vars character; climate variables.
#' @param suffixes character; the suffixes of the terms, e.g. `climateTermSuffixes` (every term there is).
#' @return `climateTermCols()`: the terms of each variable, suffix by suffix within a variable, variable by variable.
#' @export
#' @rdname climateTerm
climateTermCols <- function(vars, suffixes) as.vector(do.call(rbind, lapply(suffixes, function(sfx) paste0(vars, sfx))))

#' @param x character; column or parameter names.
#' @return `isClimateTermCol()`: logical, `TRUE` for the names that end in one of `suffixes`. `climateTermBase()`: `x`
#'   with that ending removed, i.e. the climate variable of a term.
#' @export
#' @rdname climateTerm
isClimateTermCol <- function(x, suffixes) grepl(paste0("(", paste(suffixes, collapse = "|"), ")$"), x)

#' @export
#' @rdname climateTerm
climateTermBase <- function(x, suffixes) sub(paste0("(", paste(suffixes, collapse = "|"), ")$"), "", x)

#' @format `climateTermSuffixes`: every suffix of a climate term, `_yr`, `_dev` and `_curve`.
#' @export
#' @rdname climateTerm
climateTermSuffixes <- c(climateSplitYrSuffix, climateSplitDevSuffix, climateCurveSuffix)

## `ranges` with the derived ranges of each term added: `terms` is a named list, the suffix of a term to the function that
## makes its range from the variable's `c(min, max)`. A range already in `ranges` for a term's name is kept.
climateTermRanges <- function(ranges, terms) {
  if (!length(ranges)) return(ranges)
  derived <- unlist(lapply(names(terms), function(sfx)
    stats::setNames(lapply(ranges, terms[[sfx]]), paste0(names(ranges), sfx))), recursive = FALSE)
  c(ranges, derived[setdiff(names(derived), names(ranges))])
}

#' @param dt `data.table` with the climate variables `vars` and the rows to give the two terms to.
#'   Changed by reference: each variable is replaced by its `_yr` and `_dev` columns.
#' @param ref `data.table` with the climate variables `vars` (and `by`) of the ELF's burnable pixels, the rows
#'   the median is taken over.
#' @param by character; the columns of `dt` and `ref` that define a year (`yearTxt` in the fit's tables of
#'   all years); `NULL` when both are one year (a prediction).
#' @return `splitSpreadClimate()`: `dt`, invisibly. A year with no value in `ref` gives `NA`.
#' @export
#' @rdname climateSplit
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
climateSplitRanges <- function(ranges)
  climateTermRanges(ranges, stats::setNames(list(identity, function(r) c(-1, 1) * diff(r) / 2),
                                            c(climateSplitYrSuffix, climateSplitDevSuffix)))

#' The pixels a year's climate median is taken over
#'
#' The one choice of reference pixels for [splitSpreadClimate()], in the fit and in the prediction: the ELF's
#' burnable (flammable) pixels. The fit gives its `burnablePixels` (from the data year's `landcoverDT`, already
#' flammable); the prediction gives the rows of its covariate table with `flammableRTM`, so a non-flammable row is
#' not in the median, and `elfPixels`, the ELF's own cells, so neither is a neighbour's blend-zone row.
#'
#' @param pixelID integer; candidate pixels (cells).
#' @param flammableRTM optional `SpatRaster`; only its cells equal to 1 are kept.
#' @param elfPixels optional integer; only these cells are kept (`NULL`: all).
#' @return the kept `pixelID`s.
#' @export
climateSplitRefPixels <- function(pixelID, flammableRTM = NULL, elfPixels = NULL) {
  if (!is.null(flammableRTM))
    pixelID <- pixelID[terra::values(flammableRTM, mat = FALSE)[pixelID] %in% 1]
  if (!is.null(elfPixels)) pixelID <- pixelID[pixelID %in% elfPixels]
  pixelID
}

#' The ELF of each cell of `rasterToMatchLargeELF`
#'
#' The one lookup of the ELF labels: a categorical raster gives its category labels (the second column of its
#' levels), any other raster its values as text. `unique()` of the result is the set of ELFs a raster labels.
#'
#' @param r a one-layer `SpatRaster`.
#' @return character, one label per cell in cell order (`NA` where the cell is).
#' @export
elfCellLabels <- function(r) {
  v <- terra::values(r, mat = FALSE)
  lv <- terra::levels(r)[[1]]
  if (is.data.frame(lv) && NCOL(lv) >= 2) as.character(lv[[2]][match(v, lv[[1]])]) else as.character(v)
}
