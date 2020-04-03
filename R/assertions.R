#' Data checks and assertions for \code{spreadFitRun}
#'
#' @param cells DESCRIPTION NEEDED
#' @param size DESCRIPTION NEEDED
#'
#' @return DESCRIPTION NEEDED
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

#' @param moduleName DESCRIPTION NEEDED
#' @param envir DESCRIPTION NEEDED
#' @param attribs DESCRIPTION NEEDED
#' @param fml DESCRIPTION NEEDED
#'
#' @return DESCRIPTION NEEDED
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
