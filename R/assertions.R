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
