#' Fixed ranges of the spread covariates
#'
#' `climateCovRanges` is the one table of fixed `c(min, max)` ranges, one per climate variable of the
#' spread model, that `fireSense_spreadFit` uses by default (its `covFixedRange`) instead of the range of
#' each ELF's data, so a coefficient means the same thing in every ELF. The variables and their units:
#' `CMD`, `CMDsm` and `CMDsp` are climate moisture deficit in mm (annual, summer-to-date and spring);
#' `cumMDC` is a cumulative Monthly Drought Code, an index with no unit (it is not mm). The values are
#' provisional: all four are 0-100 for now, to be replaced after the `cumMDC` values have been checked.
#' A climate covariate selected for the spread fit that is not in this table stops the fit; add its range
#' here.
#'
#' @format A named list of `c(min, max)`.
#' @export
climateCovRanges <- list(CMDsm = c(0, 100), CMD = c(0, 100), CMDsp = c(0, 100), cumMDC = c(0, 100))

#' Indicator covariates of the spread model
#'
#' The covariates that are 0 or 1: `youngAge`, the non-forest land-cover groups (`nfLCC_*`) and the
#' species-mode `treedWetland` indicator (not the biomass `treedWetland_agb`). Their range is always
#' `c(0, 1)` ([spreadIndicatorRanges()]), not the range of the data: in an ELF where one is constant
#' the data's range is zero wide and rescaling divides by zero.
#'
#' @param covNames character vector of covariate column names.
#' @return `spreadIndicatorCols()`: the names among `covNames` that are indicators.
#'   `spreadIndicatorRanges()`: a named list, `c(0, 1)` for each of them.
#' @export
#' @rdname spreadIndicators
spreadIndicatorCols <- function(covNames) {
  c(grep("^nfLCC_", covNames, value = TRUE), intersect(c(treedWetlandTxt, youngAgeTxt), covNames))
}

#' @export
#' @rdname spreadIndicators
spreadIndicatorRanges <- function(covNames) {
  ind <- spreadIndicatorCols(covNames)
  stats::setNames(rep(list(c(0, 1)), length(ind)), ind)
}

#' Spread covariates that carry no information
#'
#' A covariate is empty when it is zero everywhere, or, for fuel biomass, when it is on the
#' [logMinB()] floor everywhere: fuel is stored as log biomass, so a class with no biomass in the
#' buffers is 3.6 everywhere, not 0. The floor is recognised by its value, so the test does not
#' depend on how a column is named. An indicator with any 1, or a fuel class with any biomass, is not empty.
#' `NA`s are ignored, as `colSums(na.rm = TRUE)` did.
#'
#' @param dt `data.table` of covariates.
#' @param cols character vector: the columns of `dt` to test.
#' @return the names among `cols` that are empty.
#' @export
emptySpreadCovariates <- function(dt, cols) {
  isEmpty <- vapply(cols, function(cn) {
    x <- dt[[cn]]
    all(x == 0, na.rm = TRUE) || all(abs(x - .logMinBFloor) <= .logMinBTol, na.rm = TRUE)
  }, logical(1))
  cols[isEmpty]
}
