#' `fireSenseUtils` package
#'
#' Utilities for working with the 'fireSense' group of 'SpaDES' modules.
#'
#' @import methods
#' @importFrom grDevices rgb
#' @importFrom stats terms.formula
#' @importFrom terra cats classify cover crs freq mask mosaic sprc
#' @importFrom reproducible .robustDigest Copy Filenames asPath cacheId
#' @importFrom SpaDES.tools rastFromDF
#' @importFrom data.table setattr
#'
"_PACKAGE"

## usethis namespace: start
#' @importFrom withr deferred_run
#' @importFrom withr local_package
## usethis namespace: end
NULL

## Data.table NSE column names, `..varname` selectors, and other symbols that
## codetools cannot resolve via static analysis. None of these are package-
## level objects; they are column names referenced inside `[ ]` calls on
## data.tables, or names looked up via `mget()`-style indirection.
utils::globalVariables(c(
  "..SDcols",
  "YA_NF",
  "as.int",
  "coverSums",
  "inRange",
  "isNonForest",
  "pixelId",
  "rowcheck",
  "rows",
  "run",
  "targetFile",
  ## Anonymous lapply over the package's own internal Filenames symbol; kept
  ## so codetools doesn't re-flag mosaic now that it is qualified above.
  "mosaic"
))
