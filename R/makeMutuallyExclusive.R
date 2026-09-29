#' guarantees mutually exclusive values in a data table
#'
#' @param dt a data.table with columns that should be mutually exclusive
#' @param mutuallyExclusiveCols A named list, where the name of the list element must be a single
#'   covariate column name in `dt`. The list
#'   content should be a "grep" pattern with which to match column names, e.g., `"vegPC"`.
#'   The values of all column names that match the grep value will be set to `0`, whenever
#'   the name of that list element is non-zero. Default is `list("youngAge" = list("vegPC"))`,
#'   meaning that all columns with `vegPC` in their name will be set to zero wherever `youngAge`
#'   is non-zero. Patterns are matched anchored to the start of the column name, and the key
#'   column itself (`cov1`) is never zeroed, even if it matches its own pattern.
#'
#' @return a data.table with relevant columns made mutually exclusive
#'
#' @export
#' @importFrom data.table set
#'
makeMutuallyExclusive <- function(dt, mutuallyExclusiveCols = list("youngAge" = c("vegPC"))) {
  for (cov1 in names(mutuallyExclusiveCols)) {
    ## rows are fixed before any zeroing -- otherwise zeroing cov1 (or an earlier match) changes
    ## which rows later patterns see as non-zero
    whToZero <- which(dt[[cov1]] != 0)
    if (length(whToZero)) {
      for (grepVal in mutuallyExclusiveCols[[cov1]]) {
        cns <- grep(paste0("^", grepVal), colnames(dt), value = TRUE)
        cns <- setdiff(cns, cov1) ## the key column is never zeroed by its own pattern
        for (cn in cns) {
          set(dt, whToZero, cn, 0)
        }
      }
    }
  }
  dt
}

#' Columns `youngAge` is mutually exclusive with
#'
#' `youngAge` is mutually exclusive with every other non-climate spread covariate: wherever it is
#' non-zero, fuel biomass, non-forest land cover (always named `nfLCC_*`, see `fuelClassPrep()`)
#' and `treedWetland` are all zero. `fireSense_spreadFit` derives this from the covariates it holds
#' apart from climate (its non-annual covariate table); `fireSense_spreadPredict` has no such
#' table, so it identifies fuel columns itself (their `covMinMax` range) and calls this for the
#' rest, so both modules end up applying the identical rule.
#'
#' @param covNames character vector of all covariate column names present (excluding `pixelID`).
#' @param fuelCols character vector, the fuel biomass column names among `covNames`.
#' @return a named list, `list(youngAge = <matched column names>)`, suitable as
#'   `makeMutuallyExclusive()`'s `mutuallyExclusiveCols` argument.
#' @export
#' @importFrom stats setNames
youngAgeExclusiveCols <- function(covNames, fuelCols = character()) {
  lcc <- grep("^nfLCC_", covNames, value = TRUE)
  cols <- setdiff(unique(c(fuelCols, lcc, intersect(treedWetlandTxt, covNames))), youngAgeTxt)
  stats::setNames(list(cols), youngAgeTxt)
}
