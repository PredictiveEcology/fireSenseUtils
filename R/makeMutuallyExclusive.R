#' guarantees mutually exclusive values in a data table
#'
#' @param dt a data.table with columns that should be mutually exclusive
#' @param mutuallyExclusiveCols A named list, where the name of the list element must be a single
#'   covariate column name in `dt`. The list
#'   content should be a "grep" pattern with which to match column names, e.g., `"vegPC"`.
#'   The values of all column names that match the grep value will be set to `0`, whenever
#'   the name of that list element is non-zero. Default is `list("youngAge" = list("vegPC"))`,
#'   meaning that all columns with `vegPC` in their name will be set to zero wherever `youngAge`
#'   is non-zero.
#'
#' @return a data.table with relevant columns made mutually exclusive
#'
#' @export
#' @importFrom data.table set
#'
makeMutuallyExclusive <- function(dt, mutuallyExclusiveCols = list("youngAge" = c("vegPC"))) {
  for (cov1 in names(mutuallyExclusiveCols)) {
    for (grepVal in mutuallyExclusiveCols[[cov1]]) {
      cns <- grep(grepVal, colnames(dt), value = TRUE)
      whToZero <- which(dt[[cov1]] != 0)
      if (length(whToZero)) {
        for (cn in cns) {
          set(dt, whToZero, cn, 0)
        }
      }
    }
  }
  dt
}
