## The bufferToArea.sf() of fireSenseUtils 0.2.3.9086, which re-spread every cell of every active fire
## each iteration with SpaDES.tools::spread2(). Kept as the reference the faster version is tested against.
bufferToAreaSpread2 <- function(poly, rasterToMatch, areaMultiplier = 10,
                            verb = FALSE, polyName = NULL, field = NULL,
                            minSize = 500, cores = 1, ...) {
  if (is.null(polyName)) polyName <- "Layer 1"
  if (as.integer(verb) >= 1) print(paste("Buffering polygons on", polyName))
  r <- rasterize(
    x = sf::st_transform(poly, sf::st_crs(rasterToMatch)),
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
  ## A buffer grows one cell per iteration, so it can take at most this many to cover the raster.
  maxIts <- max(terra::nrow(r), terra::ncol(r)) + 1L
  ## Every cell of every pixel of `idsToEmit` reached so far, as buffer rows.
  emitAll <- function(df, idsToEmit) {
    lapply(idsToEmit, function(idAll) {
      dtOut <- df[df$ids %in% idAll, list(buffer = 0L, pixelID = pixels, ids)]
      dtOut[dtOut$pixelID %in% initialDf$loci[initialDf$ids %in% idAll], buffer := 1L]
      dtOut
    })
  }
  prevSize <- setNames(fireSize$actualSize, fireSize$ids)
  while ((length(loci) > 0) & (it <= maxIts)) {
    dups <- duplicated(loci)
    df <- data.table(loci = loci[!dups], ids = ids[!dups], id = seq_along(ids[!dups]))
    r1 <- SpaDES.tools::spread2(
      landscape = r, start = df$loci, iterations = 1,
      spreadProb = spreadProb, asRaster = FALSE
    )
    df <- df[r1, on = c("loci" = "initialPixels")] # TODO: confirm this
    simSizes <- df[, list(simSize = .N), by = "ids"]
    simSizes <- fireSize[simSizes, on = "ids"]
    bigger <- simSizes$simSize > simSizes$goalSize
    ## A buffer that did not grow has filled the landscape (raster edge) before reaching its target:
    ## keep it, at the largest size the landscape allows, rather than lose the fire.
    stalled <- !bigger & simSizes$simSize <= prevSize[as.character(simSizes$ids)]
    prevSize[as.character(simSizes$ids)] <- simSizes$simSize
    if (any(stalled)) {
      idsStalled <- simSizes$ids[stalled]
      names(idsStalled) <- idsStalled
      out <- append(out, emitAll(df, idsStalled))
    }

    if (any(bigger)) {
      idsBigger <- simSizes$ids[bigger]
      names(idsBigger) <- idsBigger
      out1 <- lapply(idsBigger, function(idBig) {
        wh <- which(df$ids %in% idBig)
        if (as.integer(verb) >= 2) {
          df
        }
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
        if (is(dt, "try-error"))
          stop("could not sample pixels for fire ", idBig, ": ", conditionMessage(attr(dt, "condition")))
        dtOut <- dt[, list(buffer = 0L, pixelID = pixels, ids)]

        dtOut[dtOut$pixelID %in% initialDf$loci[initialDf$ids %in% idBig], buffer := 1L]
        dtOut
      })
      out <- append(out, out1)
    }
    done <- bigger | stalled
    if (any(!done)) {
      if (any(done)) {
        simSizes <- simSizes[!done]
        df <- df[df$ids %in% simSizes$ids]
      }
      loci <- df$pixels
      ids <- df$ids
      it <- it + 1L
    } else {
      loci <- integer(0)
    }
  }
  ## fires still growing when the iterations ran out
  if (length(loci) > 0) {
    idsLeft <- unique(df$ids)
    names(idsLeft) <- idsLeft
    out <- append(out, emitAll(df, idsLeft))
  }
  out3 <- if (length(out) > 0) {
    rbindlist(out)
  } else {
    emptyDT
  }
  setorderv(out3, "buffer", order = -1L)
  out4 <- out3[, list(buffer = buffer[1], ids = ids[1]), by = "pixelID"]
}
environment(bufferToAreaSpread2) <- asNamespace("fireSenseUtils")
