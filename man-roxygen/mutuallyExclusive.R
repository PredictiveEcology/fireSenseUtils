#' @param mutuallyExclusive If there are any covariates, e.g,. `youngAge`, that should be
#'   considered mutually exclusive, then this can be done here
#'   (i.e., "if `youngAge` is non-zero, should `vegPC2` be set to zero").
#'   A named list, where the name of the list element must be a single covariate column name
#'   in either `annualDTx1000` or `nonAnnualDTx1000`.
#'   The list content should be a "grep" pattern with which to match column names, e.g., `"vegPC"`.
#'   The values of all column names that match the grep value will be set to `0`, whenever the
#'   name of that list element is non-zero. Default is `list("youngAge" = list("vegPC"))`,
#'   meaning that all columns with `vegPC` in their name will be set to zero wherever `youngAge`
#'   is non-zero.
