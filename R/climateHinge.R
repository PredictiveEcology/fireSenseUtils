## The spread climate covariate with a hinge term (`spreadClimateHinge` of fireSense_dataPrepFit). One place for the names,
## the knot, the term and its range, used by the fit (fireSense_dataPrepFit, fireSense_spreadFit) and by the prediction
## (fireSense_spreadPredict), so a pixel's hinge is the same number in both. The knot is stored in the ledger row.

#' A hinge term for the spread climate covariate
#'
#' With `spreadClimateHinge` the spread formula has an extra term for each climate variable `C`: `<C>_hinge =
#' pmax(0, C - k)`, which is 0 below the knot `k` and `C - k` above it, with its own coefficient. `C` keeps its
#' own term, so the climate response is linear below `k` and steeper (or flatter) above it.
#'
#' The knot is the `stats::median(na.rm = TRUE)` of `C` over the burnable pixels of the ELF in all the fit's years
#' (`climateHingeKnots()`), on the scale of the climate layer before the x1000 integer storage. It is one number
#' per variable for the ELF, stored in the ledger row (`spreadFitClimateHingeKnotsTxt`; `ledgerClimateHingeKnots()`),
#' and the prediction builds the term from it with `addClimateHinge()`, so the terms agree.
#'
#' The term keeps the sign rule of its variable, which `fireSense_spreadFit` applies by name: with
#' `coefSigns = "constrained"` `CMD_hinge` (a name with `CMD` or `MDC`) is bounded below by 0. It is rescaled with a
#' fixed range ([climateHingeRanges()]).
#'
#' @param vars character; climate variables, as named in the annual tables (`"CMD"`, `"Tmaxsm"`).
#' @return `climateHingeCols()`: the hinge term of each variable, `<var>_hinge`.
#' @export
#' @rdname climateHinge
climateHingeCols <- function(vars) paste0(vars, climateHingeSuffix)

#' @export
#' @rdname climateHinge
climateHingeSuffix <- "_hinge"

#' @param x character; column or parameter names.
#' @return `isClimateHingeCol()`: logical, `TRUE` for the names that end in `_hinge`. `climateHingeBase()`: `x`
#'   with that ending removed, i.e. the climate variable of a term.
#' @export
#' @rdname climateHinge
isClimateHingeCol <- function(x) grepl(paste0(climateHingeSuffix, "$"), x)

#' @export
#' @rdname climateHinge
climateHingeBase <- function(x) sub(paste0(climateHingeSuffix, "$"), "", x)

#' @param ref `data.table` with the climate variables `vars`: the rows (the ELF's burnable pixels in all the fit's
#'   years) the knot is the median over.
#' @return `climateHingeKnots()`: a named numeric, the knot of each of `vars`.
#' @export
#' @rdname climateHinge
climateHingeKnots <- function(ref, vars) {
  vapply(stats::setNames(vars, vars), function(v) stats::median(ref[[v]], na.rm = TRUE), numeric(1))
}

#' @param dt `data.table` with the climate variables named in `knots`, in the unit the knots were found in. Changed
#'   by reference: each `<var>_hinge` column is added; the variable stays.
#' @param knots named numeric (or list): the knot of each climate variable, as `climateHingeKnots()` or
#'   `ledgerClimateHingeKnots()` give them.
#' @return `addClimateHinge()`: `dt`, invisibly. `NA` climate gives `NA`.
#' @export
#' @rdname climateHinge
addClimateHinge <- function(dt, knots) {
  knots <- unlist(knots)
  for (v in names(knots))
    data.table::set(dt, NULL, climateHingeCols(v), pmax(0, dt[[v]] - knots[[v]]))
  invisible(dt)
}

#' Fixed ranges of the hinge climate terms
#'
#' `<var>_hinge` is rescaled with `c(0, w / 2)` for a variable of range width `w`: the hinge is at most the distance
#' from the median to the top of the range, and the terms of the same variable (`_yr`, `_dev`) use that width, so a unit
#' of each is the same amount of the climate variable. A range already in `ranges` for a hinge name is kept.
#'
#' @param ranges named list of `c(min, max)`, as `climateCovRanges`.
#' @return `ranges` with the hinge term of each of its variables added.
#' @export
climateHingeRanges <- function(ranges) {
  if (!length(ranges)) return(ranges)
  derived <- stats::setNames(lapply(ranges, function(r) c(0, diff(r) / 2)), climateHingeCols(names(ranges)))
  c(ranges, derived[setdiff(names(derived), names(ranges))])
}
