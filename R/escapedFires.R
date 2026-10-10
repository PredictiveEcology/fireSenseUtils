## An escaped fire is one that reached `escapeSizeHa` (the same threshold in the escape model, the spread fit
## and fireSense's burning). Fires below it, though recorded, are not spread; their sizes are kept so a
## forecast can give its non-escaped ignitions a size.

#' Fires that escaped
#'
#' @param fires Fire points or polygons (`sf`, `SpatVector` or `data.frame`) with a size column in ha.
#' @param escapeSizeHa Size (ha) a fire must reach to count as escaped.
#' @param pixSizeHa Pixel area (ha). A fire must also be larger than one pixel.
#' @param sizeCol The size column.
#' @return The rows of `fires` that escaped.
#' @export
escapedFires <- function(fires, escapeSizeHa, pixSizeHa, sizeCol = "SIZE_HA") {
  s <- as.data.frame(fires)[[sizeCol]]
  fires[which(s > pixSizeHa & s >= escapeSizeHa), ]
}
