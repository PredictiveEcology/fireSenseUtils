## Fixtures shared by the ELF merge tests.

## ELFs as vertical bands of core cells on one grid (1000 m cells), optionally with a buffer band.
bandELFs <- function(bands, nrow = 4L, buffers = list()) {
  ncol <- max(unlist(c(bands, buffers)))
  r <- terra::rast(nrows = nrow, ncols = ncol, xmin = 0, xmax = ncol * 1000,
                   ymin = 0, ymax = nrow * 1000, crs = "EPSG:3978")
  layers <- Map(cols = bands, nam = names(bands), function(cols, nam) {
    v <- matrix(0L, nrow, ncol)
    if (!is.null(buffers[[nam]])) v[, buffers[[nam]]] <- 1L
    v[, cols] <- 2L
    terra::setValues(r, as.vector(t(v)))
  })
  out <- terra::rast(unname(layers))
  names(out) <- names(bands)
  out
}

neighboursOf <- function(ELF1, ELF2, sharedLength) {
  data.table::data.table(ELF1 = ELF1, ELF2 = ELF2, sharedEdges = sharedLength / 1000,
                         sharedLength = sharedLength)
}
