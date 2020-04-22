#' Convert a list of SpatialPointsDataFrame object to a list of data.table objects
#'
#' Must supply a raster so that points can be converted to the cells on a raster.
#' It is assumed that the \code{sizeCol} is accurate. If not, it should be reclaculated
#' before this function call
#'
#'
#' @param ras A raster that will be the template for cells (pixel ids)
#' @param pts A list of SpatialPointsDataFrame objects
#' @param useCentroids Logical
#'
#' @return
#' A list of
#'
#' @export
#' @importFrom data.table rbindlist as.data.table set
#' @importFrom purrr map
#' @importFrom raster extract
makeLociList <- function(ras, pts, idCol = "NFIREID", dateCol = "YEAR", sizeCol = "POLY_HA",
                         returnCols = c("size", "date", "cells", "ids")) {
  finalCols <- c("size", "date", "cells")
  #if (isTRUE(P(sim)$useCentroids)) {
    keepCols <- c(sizeCol, dateCol, idCol)
    # keepCols <- c("POLY_HA", "YEAR", "NFIREID")
    lociDF <- purrr::map(pts, ras = ras,
                         function(.x, ras) {
                           raster::extract(x = ras,
                                           y = spTransform(.x[, keepCols], crs(ras)),
                                           cellnumbers = TRUE,
                                           sp = TRUE,
                                           df = TRUE) %>%
                             as.data.table()
                         }) %>%
      rbindlist()
    set(lociDF, NULL, "size", round(lociDF[[sizeCol]] / (prod(res(ras))/1e4)))
    browser()
    set(lociDF, NULL, setdiff(colnames(lociDF), c("size", "YEAR", "cells", "NFIREID")), NULL)
    setnames(lociDF, c("YEAR", "NFIREID"), c("date", "ids"))


  #} else {
  #   lociDF <- raster::extract(x = ras,
  #                             y = sim[["fireAttributesFireSense_SpreadFit"]],
  #                             cellnumbers = TRUE,
  #                             df = TRUE,
  #                             sp = TRUE) %>%
  #     as.data.table() %>%
  #     set(NULL, setdiff(colnames(.), finalCols), NULL)
  # }
  lociList <- split(lociDF, f = lociDF$date, keep.by = FALSE)

  lociList
}
