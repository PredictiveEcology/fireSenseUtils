#' Data checks and assertions for `spreadFitRun`
#'
#' `chk_duplicatedStartPixels` enforces the invariant that no more than one
#' fire can be ignited in a given pixel within a single time interval. When
#' duplicate `cells` are detected it issues a warning and keeps only the
#' largest fire (by `size`) for each duplicated pixel.
#'
#' @param cells Integer vector of raster cell indices (pixel IDs) at which
#'   fires are scheduled to start.
#' @param size Integer vector the same length as `cells`, giving the size of
#'   each fire (number of pixels). Used to break ties when duplicates exist.
#'
#' @return A list with two elements: `loci` (deduplicated `cells`) and
#'   `sizes` (corresponding `size` values).
#'
#' @export
#' @rdname assertions-spreadFitRun
chk_duplicatedStartPixels <- function(cells, size) {
  if (anyDuplicated(cells)) {
    warning("> No more than one fire can start in a given pixel during",
      " the same time interval, keeping the largest fire.",
      immediate. = TRUE
    )

    to_rm <- unlist(
      lapply(
        unique(cells[duplicated(cells)]),
        function(locus) {
          wh <- which(cells == locus)
          sizes <- size[wh]
          wh[-base::which.max(sizes)]
        }
      )
    )

    list(loci = cells[-to_rm], sizes = size[-to_rm])
  } else {
    list(loci = cells, sizes = size)
  }
}

#' @param moduleName Character. Name of the calling SpaDES module, used as a
#'   prefix in any error messages thrown.
#' @param envir An environment (typically a SpaDES module's `sim` object) that
#'   is expected to contain the data object named `"fireAttributesFireSense_SpreadFit"`.
#' @param attribs Character. Name of the data object being checked, included
#'   in error messages so the caller can identify which input failed.
#' @param fml A model formula passed through to [stats::is.empty.model()].
#'
#' @return Called for its side effects. Throws an informative error if any
#'   check fails; returns invisibly otherwise.
#'
#' @export
#' @importFrom stats is.empty.model
#' @rdname assertions-spreadFitRun
.doDataChecks <- function(moduleName, envir, attribs, fml) {
  if (is.null(envir[["fireAttributesFireSense_SpreadFit"]])) {
    stop(moduleName, "> '", attribs, "' not found in data objects or NULL.")
  }

  if (!is(envir[["fireAttributesFireSense_SpreadFit"]], "SpatialPointsDataFrame")) {
    stop(moduleName, "> '", attribs, "' is not a SpatialPointsDataFrame.")
  }

  if (is.null(envir[["fireAttributesFireSense_SpreadFit"]][["size"]])) {
    stop(moduleName, "> The SpatialPointsDataFrame '", attribs, "' must have a 'size' column.")
  }

  if (is.empty.model(fml)) {
    stop(moduleName, "> The formula describes an empty model.")
  }
}

#' Assert that rescaled spread covariates are usable
#'
#' Covariates were historically required to lie in `[0, 1]`. That upper bound was a
#' convention of the log-biomass parameterisation, not a requirement of the objective
#' function: covariates enter only through `exp(mat %*% covPars)` inside [logistic3p()],
#' which saturates gracefully (`exp(Inf)^(-b)` is `0`, giving `maxAsymptote`) and never
#' produces `NaN`. Fuel covariates expressed as `biomass / 1e4` legitimately exceed 1,
#' so only non-negativity and finiteness are enforced.
#'
#' @param dt A `data.table` of rescaled covariates.
#' @param colsToUse Character vector of covariate column names to check.
#'
#' @return `invisible(TRUE)`, or an error.
#' @export
assertCovariateRange <- function(dt, colsToUse) {
  vals <- dt[, ..colsToUse]
  bad <- vapply(vals, function(x) any(!is.finite(x)) || any(round(x, 3) < 0), logical(1))
  if (any(bad)) {
    stop("Covariates must be non-negative and finite; these are not: ",
         paste(colsToUse[bad], collapse = ", "))
  }
  invisible(TRUE)
}
