#' The edge ring of each fire's buffer, and its open part
#'
#' A simulated fire may spread only inside the year's buffers. `bufferEdge()` marks the pixels of
#' each fire's buffer (burned or not) with at least one of its eight (queen) neighbours outside that
#' fire's own buffer. "Outside" is anything that is not a row of the same fire in the table:
#' unburned land of another fire's buffer, pixels removed as non-flammable (see
#' [removeBufferedFiresOutsideRTM()]), `NA` cells, and the raster's boundary.
#'
#' `bufferOpen()` marks the edge pixels where a fire is actually stopped by the buffer: those with at
#' least one neighbour that is flammable (not `NA` in `r`) and outside every buffer of the year. An
#' edge pixel next to non-flammable land or next to another fire's buffer is not open: a fire there
#' could not have kept burning, or keeps burning in the neighbour's buffer. A simulated fire that burns
#' an open pixel of any buffer of its year is cut ("at least this big") by the spread objective
#' ([.objfunSpreadFit()]).
#'
#' @param fireBufferedDT A table of one fire year with columns `ids` (the fire) and `pixelID`
#'   (cell number in `r`), as in `fireBufferedListDT`. `bufferOpen()` uses its `edge` column if it has
#'   one.
#' @param r A raster with the grid `pixelID` refers to. `bufferEdge()` uses only its size;
#'   `bufferOpen()` also its values, `NA` being non-flammable. A raster without values is all flammable.
#'
#' @return `bufferEdge()`, `bufferOpen()`: a logical vector, one element per row of `fireBufferedDT`.
#'   `addBufferEdge()`: `fireBufferedListDT` with these vectors as columns `edge` and `open` in every
#'   table that does not have both already (the tables are copied, not modified).
#'
#' @export
#' @rdname bufferEdge
bufferEdge <- function(fireBufferedDT, r) {
  fire <- match(fireBufferedDT$ids, unique(fireBufferedDT$ids)) - 1
  pix <- as.numeric(fireBufferedDT$pixelID)
  n <- as.numeric(terra::ncell(r))
  key <- fire * n + pix
  edge <- logical(length(pix))
  forNeighbours(pix, r, function(nb, inside) {
    edge <<- edge | !inside | !((fire * n + nb) %in% key)
  })
  edge
}

#' @param flammable Internal: `!is.na()` of the values of `r`, when the caller has it already.
#' @export
#' @rdname bufferEdge
bufferOpen <- function(fireBufferedDT, r, flammable = cellsFlammable(r)) {
  edge <- if ("edge" %in% names(fireBufferedDT)) fireBufferedDT$edge else bufferEdge(fireBufferedDT, r)
  inAnyBuffer <- unique(as.numeric(fireBufferedDT$pixelID))
  e <- which(edge)
  open <- logical(length(edge))
  forNeighbours(as.numeric(fireBufferedDT$pixelID[e]), r, function(nb, inside) {
    free <- inside & !(nb %in% inAnyBuffer)
    if (!is.null(flammable)) free[free] <- flammable[nb[free]]
    open[e] <<- open[e] | free
  })
  open
}

## !is.na() of the raster's values; NULL (all flammable) for a raster without values
cellsFlammable <- function(r) {
  if (terra::hasValues(r)) !is.na(terra::values(r, mat = FALSE))
}

## calls fn(nb, inside) for each of the eight queen neighbours of the cells `pix` of `r`: `nb` the
## neighbour's cell number, `inside` whether it is on the raster
forNeighbours <- function(pix, r, fn) {
  nr <- terra::nrow(r)
  nc <- terra::ncol(r)
  row <- (pix - 1) %/% nc
  col <- (pix - 1) %% nc
  for (dr in -1:1) {
    for (dc in -1:1) {
      if (dr == 0 && dc == 0) next
      rr <- row + dr
      cc <- col + dc
      fn(rr * nc + cc + 1, rr >= 0 & rr < nr & cc >= 0 & cc < nc)
    }
  }
  invisible()
}

#' @param fireBufferedListDT A list of tables as `fireBufferedDT`, one per fire year.
#' @param landscape As `r`.
#' @export
#' @rdname bufferEdge
addBufferEdge <- function(fireBufferedListDT, landscape) {
  if (is.null(fireBufferedListDT)) return(NULL)
  done <- vapply(fireBufferedListDT, function(x) is.null(x) || all(c("edge", "open") %in% names(x)), logical(1))
  if (all(done)) return(fireBufferedListDT)
  flammable <- cellsFlammable(landscape) # once for all years
  fireBufferedListDT[!done] <- lapply(fireBufferedListDT[!done], function(x) {
    x <- data.table::copy(data.table::as.data.table(x))
    if (!"edge" %in% names(x)) data.table::set(x, NULL, "edge", bufferEdge(x, landscape))
    data.table::set(x, NULL, "open", bufferOpen(x, landscape, flammable))
    x
  })
  fireBufferedListDT
}

## Arguments of the spread objective that no longer do anything, and why. Each warns once per session.
retiredObjectiveArgs <- local({
  cap <- paste("simulated fires are no longer capped at a size. A fire that reaches an open edge of the",
               "year's buffers is cut instead; see `penaliseRunaways`.")
  noTerm <- function(what) paste0("the spread objective no longer has the ", what, " term; it is the per-fire ",
                                  "kernel log score only, which fires cut at the buffers cannot bias.")
  cut <- "a replicate cut at the edge of the buffers has no size in the per-fire kernel log score."
  c(capSizes = cap, penaliseCapHits = cap,
    yearAreaWeight = noTerm("annual-area"), areaDistWeight = noTerm("area-weighted size-distribution"),
    adTest = noTerm("pooled Anderson-Darling (`tests = \"adTest\"`)"),
    runawaySize = cut, runawayBufferMultiple = cut)
})
.deprecatedArgsWarned <- new.env(parent = emptyenv())
## `argNames`: the retired arguments a caller used (others are ignored here)
deprecatedArgs <- function(argNames) {
  new <- setdiff(intersect(argNames, names(retiredObjectiveArgs)), ls(.deprecatedArgsWarned))
  if (!length(new)) return(invisible())
  for (nm in new) assign(nm, TRUE, envir = .deprecatedArgsWarned)
  warning(paste0("`", new, "` is deprecated and ignored: ", retiredObjectiveArgs[new], collapse = "\n"),
          call. = FALSE)
  invisible()
}
