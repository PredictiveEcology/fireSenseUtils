#' @param mutuallyExclusive If there are any covariates, e.g,. youngAge, that should be
#'   considered mutually exclusive, i.e., "if youngAge is non-zero, should vegPC2 be set to zero", then
#'   this can be done here. A named list, where the name of the list element must be a single
#'   covariate column name in either \code{annualDTx1000} or \code{nonAnnualDTx1000}. The list
#'   content should be a "grep" pattern with which to match column names, e.g., \code{"vegPC"}.
#'   The values of all column names that match the grep value will be set to \code{0}, whenever
#'   the name of that list element is non-zero. Default is \code{list("youngAge" = list("vegPC"))},
#'   meaning that all columns with \code{vegPC} in their name will be set to zero wherever \code{youngAge}
#'   is non-zero.
