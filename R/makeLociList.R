#' Convert a list of \code{SpatialPointsDataFrame} object to a list of \code{data.table} objects
#'
#' Must supply a raster so that points can be converted to the cells on a raster.
#' It is assumed that the \code{sizeCol} is accurate.
#' If not, it should be recalculated before this function call.
#'
#' @param ras A raster that will be the template for cells (pixel ids)
#' @param pts A list of \code{SpatialPointsDataFrame} objects
#' @param idsCol Character string identifying column name in \code{pts} that has unique
#'   id per event (i.e., fire)
#' @param dateCol Character string identifying column name in \code{pts} that has year
#' @param sizeCol Character string identifying column name in \code{pts} that has size of
#'   individual event. Can be in hectares or metres squared. Should set \code{sizeColUnits}
#' @param sizeColUnits Character string. Either \code{"ha"} or \code{"m2"}.
#'
#' @return
#' A list of data.table objects, each with 4 columns, \code{"size"} (in pixels), \code{"date"},
#' \code{"ids"} from \code{idsCol}, and \code{"cells"}, which are the pixel indices of the
#' \code{pts} points.
#'
#' @export
#' @importFrom data.table rbindlist as.data.table set
#' @importFrom magrittr %>%
#' @importFrom purrr map
#' @importFrom raster extract
makeLociList <- function(ras, pts, idsCol = "FIRE_ID", dateCol = "YEAR", sizeCol = "POLY_HA",
                         sizeColUnits = "ha") {
  returnCols = c("size", "date", "ids", "cells")
  keepCols <- c(sizeCol, dateCol, idsCol)
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
  dtColNames <- data.table(old = c(sizeCol, dateCol, idsCol, "cells"),
                           new = returnCols)
  ma <- match(dtColNames$old, colnames(lociDF))
  setnames(lociDF, old = dtColNames$old[ma], new = dtColNames$new[ma])
  # setnames(lociDF, c(dateCol, idsCol), c("date", "ids"))
  set(lociDF, NULL, setdiff(colnames(lociDF), returnCols), NULL)
  divisor <- switch(sizeColUnits,
                    "ha" = 1e4,
                    "m2" = 1,
                    stop("Must provide siceColUnits either ha or m2"))
  set(lociDF, NULL, "size", round(lociDF$size / (prod(res(ras))/divisor), 0))

  for (index in colnames(lociDF)) {
    if (is.numeric(lociDF[[index]]) && max(lociDF[[index]] < 1e9) && !is.integer(lociDF[[index]]))
      set(lociDF, NULL, index, as.integer(lociDF[[index]]))
  }

  split(lociDF, f = lociDF$date, keep.by = FALSE)
}
