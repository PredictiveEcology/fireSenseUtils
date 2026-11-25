## The spread climate covariate with a smooth curve term (`spreadClimateCurve` of fireSense_dataPrepFit). One place for
## the names, the centre and scale, the term and its range, used by the fit (fireSense_dataPrepFit, fireSense_spreadFit) and by
## the prediction (fireSense_spreadPredict), so a pixel's curve term is the same number in both. The centre and scale are
## stored in the ledger row.

#' The softplus function
#'
#' `softplus(x, s) = s * log(1 + exp(x / s))`: a smooth version of `pmax(0, x)` with no corner. It is near 0 well below 0
#' and near `x` well above it, and bends over a distance of about `s` around 0 (`softplus(0, s) = s * log(2)`). Computed
#' as `pmax(x, 0) + s * log1p(exp(-abs(x) / s))`, so a large `x / s` does not overflow.
#'
#' @param x numeric.
#' @param s numeric, the scale (> 0) of the bend; smaller is closer to `pmax(0, x)`.
#' @return numeric, the length of `x`; `NA` stays `NA`.
#' @export
softplus <- function(x, s = 1) pmax(x, 0) + s * log1p(exp(-abs(x) / s))

#' A smooth curve term for the spread climate covariate
#'
#' With `spreadClimateCurve` the spread formula has an extra term for each climate variable `C`: `<C>_curve =
#' softplus(C - k, s)` ([softplus()]), which is near 0 below the centre `k` and near `C - k` above it, bending gradually
#' with no corner, with its own coefficient. `C` keeps its own term, so the climate response is linear well below `k`
#' and steeper above it when the coefficient is positive. The coefficient's sign is free; a negative value flattens the
#' dry side instead.
#'
#' `k` is the `stats::median(na.rm = TRUE)` and `s` the `stats::sd(na.rm = TRUE)` of `C` over the burnable pixels of the
#' ELF in all the fit's years (`climateCurvePars()`), on the scale of the climate layer before the x1000 integer storage.
#' They are one pair of numbers per variable for the ELF, stored in the ledger row (`spreadFitClimateCurveParsTxt`;
#' `ledgerClimateCurvePars()`), and the prediction builds the term from them with `addClimateCurve()`, so the terms agree.
#'
#' The term keeps the sign rule of its variable, which `fireSense_spreadFit` applies by name: with
#' `coefSigns = "constrained"` `CMD_curve` (a name with `CMD` or `MDC`) is bounded below by 0. It is rescaled with a
#' fixed range ([climateCurveRanges()]).
#'
#' @param vars character; climate variables, as named in the annual tables (`"CMD"`, `"Tmaxsm"`).
#' @return `climateCurveCols()`: the curve term of each variable, `<var>_curve`.
#' @export
#' @rdname climateCurve
climateCurveCols <- function(vars) paste0(vars, climateCurveSuffix)

#' @export
#' @rdname climateCurve
climateCurveSuffix <- "_curve"

#' @param x character; column or parameter names.
#' @return `isClimateCurveCol()`: logical, `TRUE` for the names that end in `_curve`. `climateCurveBase()`: `x`
#'   with that ending removed, i.e. the climate variable of a term.
#' @export
#' @rdname climateCurve
isClimateCurveCol <- function(x) grepl(paste0(climateCurveSuffix, "$"), x)

#' @export
#' @rdname climateCurve
climateCurveBase <- function(x) sub(paste0(climateCurveSuffix, "$"), "", x)

#' @param ref `data.table` with the climate variables `vars`: the rows (the ELF's burnable pixels in all the fit's
#'   years) the centre and scale are taken over.
#' @return `climateCurvePars()`: a named list, one element per variable of `vars`: `c(k = , s = )`, the median and the
#'   sd. Stops for a variable whose sd is not above 0 (a constant or missing climate has no curve).
#' @export
#' @rdname climateCurve
climateCurvePars <- function(ref, vars) {
  pars <- lapply(stats::setNames(vars, vars), function(v)
    c(k = stats::median(ref[[v]], na.rm = TRUE), s = stats::sd(ref[[v]], na.rm = TRUE)))
  checkClimateCurvePars(pars)
  pars
}

#' @param dt `data.table` with the climate variables named in `pars`, in the unit `pars` were found in. Changed
#'   by reference: each `<var>_curve` column is added; the variable stays.
#' @param pars named list: `c(k = , s = )` of each climate variable, as `climateCurvePars()` or
#'   `ledgerClimateCurvePars()` give them.
#' @return `addClimateCurve()`: `dt`, invisibly. `NA` climate gives `NA`.
#' @export
#' @rdname climateCurve
addClimateCurve <- function(dt, pars) {
  checkClimateCurvePars(pars)
  for (v in names(pars))
    data.table::set(dt, NULL, climateCurveCols(v), softplus(dt[[v]] - pars[[v]][["k"]], pars[[v]][["s"]]))
  invisible(dt)
}

checkClimateCurvePars <- function(pars) {
  ok <- vapply(pars, function(p) all(c("k", "s") %in% names(p)) && is.finite(p[["k"]]) && is.finite(p[["s"]]) && p[["s"]] > 0,
               logical(1))
  if (!all(ok))
    stop("the climate curve needs a finite centre `k` and a scale `s` above 0 for ", paste(names(pars)[!ok], collapse = ", "),
         call. = FALSE)
  invisible(pars)
}

#' Fixed ranges of the curve climate terms
#'
#' `<var>_curve` is rescaled with `c(0, w / 2)` for a variable of range width `w`: the term is at most about the distance
#' from the centre to the top of the range, and the terms of the same variable (`_yr`, `_dev`) use that width, so a unit
#' of each is the same amount of the climate variable. A range already in `ranges` for a curve name is kept.
#'
#' @param ranges named list of `c(min, max)`, as `climateCovRanges`.
#' @return `ranges` with the curve term of each of its variables added.
#' @export
climateCurveRanges <- function(ranges) {
  if (!length(ranges)) return(ranges)
  derived <- stats::setNames(lapply(ranges, function(r) c(0, diff(r) / 2)), climateCurveCols(names(ranges)))
  c(ranges, derived[setdiff(names(derived), names(ranges))])
}
