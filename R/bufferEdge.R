#' The edge ring of each fire's buffer
#'
#' A simulated fire that reaches the outer edge of its own buffer has run away: the spread
#' objective ([.objfunSpreadFit()]) censors it. `bufferEdge()` marks the pixels that make that
#' edge: a pixel of a fire's buffer (burned or not) with at least one of its eight (queen)
#' neighbours outside that fire's own buffer. "Outside" is anything that is not a row of the same
#' fire in the table: unburned land of another fire's buffer, pixels removed as non-flammable (see
#' [removeBufferedFiresOutsideRTM()]), `NA` cells, and the raster's boundary.
#'
#' @param fireBufferedDT A table of one fire year with columns `ids` (the fire) and `pixelID`
#'   (cell number in `r`), as in `fireBufferedListDT`.
#' @param r A raster with the grid `pixelID` refers to; only its size is used.
#'
#' @return `bufferEdge()`: a logical vector, one element per row of `fireBufferedDT`.
#'   `addBufferEdge()`: `fireBufferedListDT` with that vector as a column `edge` in every table
#'   that does not have one already (the tables are copied, not modified).
#'
#' @export
#' @rdname bufferEdge
bufferEdge <- function(fireBufferedDT, r) {
  nr <- terra::nrow(r)
  nc <- terra::ncol(r)
  n <- as.numeric(nr) * nc
  pix <- as.numeric(fireBufferedDT$pixelID)
  fire <- match(fireBufferedDT$ids, unique(fireBufferedDT$ids)) - 1
  row <- (pix - 1) %/% nc
  col <- (pix - 1) %% nc
  key <- fire * n + pix
  edge <- logical(length(pix))
  for (dr in -1:1) {
    for (dc in -1:1) {
      if (dr == 0 && dc == 0) next
      rr <- row + dr
      cc <- col + dc
      inside <- rr >= 0 & rr < nr & cc >= 0 & cc < nc
      edge <- edge | !inside | !((fire * n + rr * nc + cc + 1) %in% key)
    }
  }
  edge
}

#' @param fireBufferedListDT A list of tables as `fireBufferedDT`, one per fire year.
#' @param landscape As `r`.
#' @export
#' @rdname bufferEdge
addBufferEdge <- function(fireBufferedListDT, landscape) {
  if (is.null(fireBufferedListDT)) return(NULL)
  lapply(fireBufferedListDT, function(x) {
    if (is.null(x) || "edge" %in% names(x)) return(x)
    x <- data.table::copy(data.table::as.data.table(x))
    data.table::set(x, NULL, "edge", bufferEdge(x, landscape))
    x
  })
}

## Once per session: `capSizes` and `penaliseCapHits` no longer do anything.
.capArgsWarned <- new.env(parent = emptyenv())
deprecatedCapArgs <- function(argNames) {
  old <- intersect(argNames, c("capSizes", "penaliseCapHits"))
  if (length(old) && is.null(.capArgsWarned$done)) {
    .capArgsWarned$done <- TRUE
    warning("`", paste(old, collapse = "`, `"), "` is deprecated and ignored: simulated fires are no ",
            "longer capped at a size. A fire that reaches the edge of its own buffer is penalised ",
            "instead; see `penaliseRunaways`.", call. = FALSE)
  }
  invisible()
}
