#' Create buffers around polygons based on area target for buffer
#'
#' @param verb Logical or numeric related to how much verbosity is printed. \code{FALSE} or
#'   \code{0} is none. \code{TRUE} or \code{1} is some. \code{2} is much more.
#' @export
#' @importFrom data.table data.table
#' @importFrom fasterize fasterize
#' @importFrom purrr map2
#' @importFrom sf st_as_sf st_crs st_transform
#' @param poly SpatialPolygons or a list of SpatialPolygons, containing polygons to buffer
#' @param rasterToMatch A RasterLayer with res, origin, extent, crs of desired outputted
#'   pixelID values
#' @param areaMultiplier Either a scalar that will buffer areaMultiplier * fireSize or
#'   a function of fireSize. Default is 1. See \code{\link{multiplier}} for an example.
#' @param polyName Optional character string of the polygon layer name (not the individual polygons
#'   on a SpatialPolygons object)
#' @param field Passed to \code{fasterize::fasterize}. If this is unique (such as polygon id),
#'   then each polygon will have its buffer calculated independently for each unique value
#'   in \code{field}
#' @param ... passed to \code{fasterize::fasterize}
#' @export
#' @rdname bufferToArea
#' @return
#' A data.table (or list of data.tables if \code{poly} was a list) with 2 columns
#' \code{buffer} and \code{pixelID}. \code{buffer} is either \code{1} (the
#' original polygon) or \code{0} (in the buffer)
bufferToArea <- function(poly, rasterToMatch, areaMultiplier,
                         verb = FALSE, polyName = NULL, field = NULL,
                         ...) {
  UseMethod("bufferToArea")
}

#' @export
#' @rdname bufferToArea
bufferToArea.list <- function(poly, rasterToMatch, areaMultiplier = 1,
                              verb = FALSE, polyName = NULL, field = NULL, ...) {
  if (is.null(polyName)) polyName <- names(poly)
  out <-  purrr::pmap(
    .l = list(poly = poly, polyName = polyName),
    rasterToMatch = rasterToMatch, verb = verb,
    areaMultiplier = areaMultiplier, field = field, ...,
    .f = bufferToArea)
  out
}

#' @export
#' @rdname bufferToArea
bufferToArea.SpatialPolygons <- function(poly, rasterToMatch, areaMultiplier = 1,
                                         verb = FALSE, polyName = NULL, field = NULL,
                                         ...) {
  if (is.null(polyName)) polyName <- "Layer 1"
  if  (as.integer(verb) >= 1) print(paste("Buffering polygons on",polyName))
  r <- fasterize::fasterize(
    sf::st_transform(sf::st_as_sf(poly), sf::st_crs(rasterToMatch)),
    raster = rasterToMatch, field = field, ...)

  loci <- which(!is.na(r[]))
  ids <- r[loci]
  initialDf <- data.table(loci, ids, id = seq(ids))
  am <- if (is(areaMultiplier, "function")) {
    areaMultiplier
  } else {
    function(x) areaMultiplier
  }
  fireSize <- initialDf[, list(actualSize = .N,
                               # simSize = .N,# needed for numIters
                               goalSize = areaMultiplier(.N)), by = "ids"]

  out <- list()
  simSizes <- initialDf[, list(simSize = .N), by = "ids"]
  simSizes <- fireSize[simSizes, on = "ids"]

  # if (!is.null(allowCells)) {
  #   spreadProb <- rep(NA, ncell(r))
  #   spreadProb[allowCells] <- 1
  # } else {
     spreadProb <- 1
  # }
  while(length(loci) > 0) {
    df <- data.table(loci, ids, id = seq(ids))
    r1 <- spread(r, loci = loci, iterations = 1,
                 spreadProb = spreadProb, quick = TRUE, returnIndices = TRUE)
    df <- df[r1, on = "id"]
    simSizes <- df[, list(simSize = .N), by = "ids"]
    simSizes <- fireSize[simSizes, on = "ids"]
    bigger <- simSizes$simSize > simSizes$goalSize

    if (any(bigger)) {
      idsBigger <- simSizes$ids[bigger]
      names(idsBigger) <- idsBigger
      out1 <- lapply(idsBigger, function(idBig) {
        browser(expr = idBig == "423" || idBig == 423)
        wh <- which(df$ids %in% idBig)
        if (as.integer(verb) >= 2) print(paste("  Fire id:,",idBig, "finished. Num pixels in buffer:",
                                               simSizes[ids == idBig]$goalSize - simSizes[ids == idBig]$actualSize,
                                               ", in fire:", simSizes[ids == idBig]$actualSize))
        lastIters <- !df[wh]$active
        needMore <- simSizes[ids == idBig]$goalSize - sum(lastIters)
        dt <- rbindlist(list(df[wh][lastIters],
                             df[wh][sample(which(df[wh]$active), needMore)]))
        dtOut <- dt[, list(buffer = 0, pixelID = indices, ids)]
        #browser(expr = exists("jj"))
        #rm(jj, envir = .GlobalEnv)

        dtOut[dtOut$pixelID %in% initialDf$loci[initialDf$ids %in% idBig], buffer := 1]
        dtOut
      })
      out <- append(out, out1)
    }
    if (any(!bigger)) {
      if (!all(!bigger)) {
        simSizes <- simSizes[simSize <= goalSize]
        df <- df[df$ids %in% simSizes$ids]
      }
      loci <- df$indices
      ids <- df$ids
    } else {
      loci <- integer()
    }
  }
  rbindlist(out)
}

#' @export
multiplier <- function(size) {
  round(pmax(5e2, pmax(2, 14 - log(size)) * size), 0)
}
