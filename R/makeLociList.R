globalVariables(c(
"cell"
))
#' Convert a list of `SpatialPointsDataFrame` object to a list of `data.table` objects
#'
#' Must supply a raster so that points can be converted to the cells on a raster.
#' It is assumed that the `sizeCol` is accurate.
#' If not, it should be recalculated before this function call.
#'
#' @param ras A raster that will be the template for cells (pixel ids)
#' @param pts A list of `sf` point objects
#' @param idsCol Character string identifying column name in `pts` that has unique
#'   id per event (i.e., fire)
#' @param dateCol Character string identifying column name in `pts` that has year
#' @param sizeCol Character string identifying column name in `pts` that has size of
#'   individual event. Can be in hectares or metres squared. Should set `sizeColUnits`
#' @param sizeColUnits Character string. Either `"ha"` or `"m2"`.
#'
#' @return
#' A list of data.table objects, each with 4 columns, `"size"` (in pixels), `"date"`,
#' `"ids"` from `idsCol`, and `"cells"`, which are the pixel indices of the
#' `pts` points.
#'
#' @export
#' @importFrom data.table rbindlist as.data.table set
#' @importFrom purrr map
#' @importFrom terra extract res
#' @importFrom sf %>% st_crs st_transform
makeLociList <- function(ras, pts, idsCol = "FIRE_ID", dateCol = "YEAR", sizeCol = "POLY_HA",
                         sizeColUnits = "ha") {
  returnCols <- c("size", "date", "ids", "cells")
  keepCols <- c(sizeCol, dateCol, idsCol)
  #pts shoudl already be projected but no harm in forcing..
  pts <- lapply(pts, st_transform, crs = st_crs(ras))

  lociDF <- purrr::map(pts,
    ras = ras,
    function(x, ras) {
    extract(x = ras, y = x, bind = TRUE, cells = TRUE) %>%
        as.data.table()
    }
  ) %>%
    rbindlist()
  dtColNames <- data.table(
    old = c(sizeCol, dateCol, idsCol, "cell"),
    new = returnCols
  )
  ma <- match(dtColNames$old, colnames(lociDF))
  lociDF <- lociDF[, .SD, .SDcol = dtColNames$old]
  setnames(lociDF, old = dtColNames$old, new = dtColNames$new)
  divisor <- switch(sizeColUnits,
    "ha" = 1e4,
    "m2" = 1,
    stop("Must provide sizeColUnits either ha or m2")
  )
  set(lociDF, NULL, "size", round(lociDF$size / (prod(res(ras)) / divisor), 0))

  for (index in colnames(lociDF)) {
    if (is.numeric(lociDF[[index]]) &&
      max(as.numeric(lociDF[[index]]) < 1e9) && ## ~ 2^31 / 2 (half of signed integer bits)
      !is.integer(lociDF[[index]])) {
      set(lociDF, NULL, index, as.integer(lociDF[[index]]))
    }
  }
  lociDF[, date := paste0("year", date)] # fires now have year in front
  split(lociDF, f = lociDF$date, keep.by = FALSE)
}
