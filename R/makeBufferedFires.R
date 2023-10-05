utils::globalVariables(c(
  ".N", "goalSize", "simSize", "pixels"
))

#' Create buffers around polygons based on area target for buffer
#'
#' @param poly `sf` polygons or a list of `sf` containing polygons to buffer.
#' @param rasterToMatch A `SpatRaster` with `res`, `origin`, `extent`,
#'   `crs` of desired outputted `pixelID` values.
#' @param areaMultiplier Either a scalar that will buffer `areaMultiplier * fireSize` or
#'   a function of `fireSize.` Default is 1. See [multiplier()] for an example.
#' @param verb Logical or numeric related to how much verbosity is printed. `FALSE` or
#'   `0` is none. `TRUE` or `1` is some. `2` is much more.
#' @param polyName Optional character string of the polygon layer name (not the individual polygons
#'   on a `sf` polygon object)
#' @param field Passed to `fasterize::fasterize`. If this is unique (such as polygon id),
#'   then each polygon will have its buffer calculated independently for each unique value
#'   in `field`
#' @param minSize The absolute minimum size of the buffer & non-buffer together. This will
#'   be imposed after `areaMultiplier`.
#' @param cores number of processor cores to use
#' @param ... passed to `fasterize::fasterize`
#'
#' @return
#' A `data.table` (or list of `data.table`s if `poly` was a list) with 2 columns:
#' `buffer` and `pixelID`. `buffer` is either `1` (the original polygon) or
#' `0` (in the buffer).
#'
#' @export
#' @rdname bufferToArea
bufferToArea <- function(poly, rasterToMatch, areaMultiplier,
                         verb = FALSE, polyName = NULL, field = NULL,
                         minSize = 500, cores = 1, ...) {
  UseMethod("bufferToArea")
}

#' @export
#' @importFrom parallelly availableCores
#' @importFrom purrr pmap
#'
#' @rdname bufferToArea
bufferToArea.list <- function(poly, rasterToMatch, areaMultiplier = 10,
                              verb = FALSE, polyName = NULL, field = NULL,
                              minSize = 500, cores = 1, ...) {
  if (is.null(polyName)) {
    polyName <- names(poly)
  }
  maxCores <- parallelly::availableCores(constraints = "connections", omit = 1)
  cores <- min(min(length(poly), cores), maxCores)
  if (cores > 1) {
    out <- parallel::mcMap(
      mc.cores = cores,
      poly = poly,
      polyName = polyName,
      MoreArgs = list(
        rasterToMatch = rasterToMatch, verb = verb,
        areaMultiplier = areaMultiplier, field = field, minSize = minSize,
        cores = 1,
        ...
      ),
      bufferToArea
    )
  } else {
    out <- purrr::pmap(
      .l = list(poly = poly, polyName = polyName),
      rasterToMatch = rasterToMatch, verb = verb,
      areaMultiplier = areaMultiplier, field = field, minSize = minSize,
      cores = 1,
      ...,
      .f = bufferToArea
    )
  }
  names(out) <- names(poly)
  out
}

#' @export
#' @importFrom sf st_as_sf
#' @rdname bufferToArea
bufferToArea.SpatialPolygons <- function(poly, rasterToMatch, areaMultiplier = 10,
                                         verb = FALSE, polyName = NULL, field = NULL,
                                         minSize = 500, cores = 1, ...) {
  bufferToArea.sf(sf::st_as_sf(poly), rasterToMatch, areaMultiplier = areaMultiplier,
                  verb = verb, polyName = polyName, field = field,
                  minSize = minSize, cores = cores, ...)
}

#' @export
#' @importFrom data.table data.table rbindlist setorderv
#' @importFrom terra rasterize values
#' @importFrom sf st_crs st_transform
#' @importFrom SpaDES.tools spread2
#' @rdname bufferToArea
bufferToArea.sf <- function(poly, rasterToMatch, areaMultiplier = 10,
                            verb = FALSE, polyName = NULL, field = NULL,
                            minSize = 500, cores = 1, ...) {
  if (is.null(polyName)) polyName <- "Layer 1"
  if (as.integer(verb) >= 1) print(paste("Buffering polygons on", polyName))
  r <- rasterize(x = sf::st_transform(poly, sf::st_crs(rasterToMatch)),
                 y = rasterToMatch, field = field, ...
  )

  emptyDT <- data.table(pixelID = integer(0), buffer = integer(0), ids = integer(0))

  if (all(is.na(r[]))) {
    return(emptyDT)
  }
  rvals <- values(r, mat = FALSE)
  loci <- which(!is.na(rvals))
  ids <- as.integer(rvals[loci])

  initialDf <- data.table(loci, ids, id = seq(ids))
  am <- if (is(areaMultiplier, "function")) {
    areaMultiplier
  } else {
    function(x) areaMultiplier * x
  }
  fireSize <- initialDf[, list(
    actualSize = .N,
    # simSize = .N,# needed for numIters
    goalSize = asInteger(pmax(minSize, am(.N)))
  ), by = "ids"]

  out <- list()
  simSizes <- initialDf[, list(simSize = .N), by = "ids"]
  simSizes <- fireSize[simSizes, on = "ids"]

  # if (!is.null(allowCells)) {
  #   spreadProb <- rep(NA, ncell(r))
  #   spreadProb[allowCells] <- 1
  # } else {
  spreadProb <- 1
  # }
  it <- 1L
  maxIts <- 100L ## TODO: what's reasonable here???
  while ((length(loci) > 0) & (it <= maxIts)) {
    dups <- duplicated(loci)
    df <- data.table(loci = loci[!dups], ids = ids[!dups], id = seq_along(ids[!dups]))
    r1 <- SpaDES.tools::spread2(landscape = r, start = df$loci, iterations = 1,
                                spreadProb = spreadProb, asRaster = FALSE)
    df <- df[r1, on = c("loci" = "initialPixels")] #TODO: confirm this
    simSizes <- df[, list(simSize = .N), by = "ids"]
    simSizes <- fireSize[simSizes, on = "ids"]
    bigger <- simSizes$simSize > simSizes$goalSize

    if (any(bigger)) {
      idsBigger <- simSizes$ids[bigger]
      names(idsBigger) <- idsBigger
      out1 <- lapply(idsBigger, function(idBig) {
        wh <- which(df$ids %in% idBig)
        if (as.integer(verb) >= 2) {df}
        lastIters <- !df[wh]$state == "activeSource"
        needMore <- simSizes[ids == idBig]$goalSize - sum(lastIters)
        if (needMore > 0) {
          dt <- try(rbindlist(list(
            df[wh][lastIters],
            df[wh][sample(which(df[wh]$state == "activeSource"), needMore)]
          )))
        } else {
          dt <- df[wh][lastIters][sample(sum(lastIters), simSizes[ids == idBig]$goalSize)]
        }
        if (is(dt, "try-error")) browser() # stop("try error here")
        dtOut <- dt[, list(buffer = 0L, pixelID = pixels, ids)]

        dtOut[dtOut$pixelID %in% initialDf$loci[initialDf$ids %in% idBig], buffer := 1L]
        dtOut
      })
      out <- append(out, out1)
    }
    if (any(!bigger)) {
      if (!all(!bigger)) {
        simSizes <- simSizes[simSize <= goalSize]
        df <- df[df$ids %in% simSizes$ids]
      }
      loci <- df$pixels
      ids <- df$ids
      it <- it + 1L
    } else {
      loci <- integer(0)
    }
  }
  out3 <- if (length(out) > 0) {
    rbindlist(out)
  } else {
    emptyDT
  }
  setorderv(out3, "buffer", order = -1L)
  out4 <- out3[, list(buffer = buffer[1], ids = ids[1]), by = "pixelID"]
}

#' multiplier
#'
#' DESCRIPTION NEEDED
#'
#' @param size DESCRIPTION NEEDED
#' @param minSize DESCRIPTION NEEDED
#' @param baseMultiplier DESCRIPTION NEEDED
#'
#' @export
multiplier <- function(size, minSize = 1000, baseMultiplier = 5) {
  pmax(minSize, round(pmax(baseMultiplier, 14 - log(size)) * size, 0))
}

#' buffer ignition points to create non-ignitions for model
#'
#'
#' @param ignitionPoints SpatialPolygonsDataFrame with year of ignition
#' @param rtm a template raster
#' @param bufferSize the size of the buffers
#'
#' @return a list of data.tables containing indices inside buffered area of each year's ignitions
#'
#' @importFrom terra buffer rast set.values extract rasterize
#' @importFrom sf st_as_sf
#' @export
#' @rdname bufferIgnitionPoints
bufferIgnitionPoints <- function(ignitionPoints, rtm, bufferSize) {
  years <- sort(unique(ignitionPoints$YEAR))
  rtm <- setValues(rtm, 1:ncell(rtm))
  ignitionDT <- lapply(years, FUN = function(year, rtmIndices = rtm,
                                             ignitions = ignitionPoints) {
    # subset ignitions to year
    annIg <- ignitions[ignitions$YEAR == year, ]

    # buffer them
    buffed <- buffer(annIg, bufferSize)

    # get rtm indices of ignitions
    ignitionPix <- extract(rtmIndices, annIg)

    # get rtm indices of buffer
    buffed <- sf::st_as_sf(buffed)
    buffed <- rasterize(buffed, y = rtmIndices)
    annIndices <- rtmIndices[!is.na(buffed[])]

    # make data.table
    indices <- data.table(pixelID = annIndices, ignited = annIndices %in% ignitionPix)
    indices <- indices[!duplicated(pixelID), ]

    return(indices)
  })
  # this does not list which cell ignited
  names(ignitionDT) <- paste0("year", years)
  return(ignitionDT)
}
