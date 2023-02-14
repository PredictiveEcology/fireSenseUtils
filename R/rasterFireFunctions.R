utils::globalVariables(c(
  "pixelID", "fireID", "ids", "id"
))

#' prepare fire
#' @param fireRaster a \code{SpatRaster} with values representing fire years
#' @param year numeric fire year
#' @template flammableRTM
#' @param bufferForFireRaster buffer size used to group discrete patches of burned pixels as
#' belonging to the same fire
#' @param areaMultiplier A scalar that will buffer \code{areaMultiplier * fireSize}
#'  @param verb Logical or numeric related to how much verbosity is printed. \code{FALSE} or
#'   \code{0} is none. \code{TRUE} or \code{1} is some. \code{2} is much more.
#' @param minSize The absolute minimum size of the buffer & non-buffer together. This will
#'   be imposed after \code{areaMultiplier}.
#' @return a data.table with year, fire ID, buffer status, and pixelID
#' @export
#' @importFrom data.table rbindlist data.table set
#' @importFrom terra buffer patches
rasterFireBufferDT <- function(fireRaster, year, flammableRTM, bufferForFireRaster, areaMultiplier, minSize, verb = 1) {

  #identify this year's fires
  origIndex <- which(fireRaster[] == year)
  if (length(origIndex) == 0) {
    return(NULL) #TODO: what should happen here? empty data.table/
  } else{
    fireRaster[] <- NA
    fireRaster[origIndex] <- 1
    #classify into discrete fires via buffer
    fireRaster <- buffer(fireRaster, width = bufferForFireRaster, background = NA)
    #doublecheck that no patches are lost... why are patch values greater than # of patches
    fireRaster <- terra::patches(fireRaster)
    fireIDs <- fireRaster[origIndex]

    #mask to the original pixels
    fireRaster[] <- NA
    fireRaster[origIndex] <- fireIDs
    fireBufferDT <- bufferToAreaRast(fireIDraster = fireRaster,
                                     flammableRTM = flammableRTM,
                                     areaMultiplier = areaMultiplier,
                                     minSize = minSize,
                                     verb = verb)
    return(fireBufferDT)
  }
}

#' create a variable sized buffer around a set of pixels belonging to the same fire ID
#' @param fireIDraster a \code{SpatRaster} with values representing distinct fires in a year
#' @param flammableRTM @template flammableRTM
#' @param areaMultiplier A scalar that will buffer \code{areaMultiplier * fireSize}
#' @param minSize The absolute minimum size of the buffer & non-buffer together. This will
#'   be imposed after \code{areaMultiplier}.
#' @return a data.table with fire ID, buffer status, and pixelID
#' @export
#' @importFrom data.table rbindlist data.table set
#' @importFrom terra buffer res values
#' @importFrom SpaDES.tools spread
bufferToAreaRast <- function(fireIDraster, areaMultiplier, minSize, flammableRTM, verb = 1) {

  am <- if (is(areaMultiplier, "function")) {
    areaMultiplier
  } else {
    function(x) areaMultiplier * x
  }


  ##########
  loci <- which(!is.na(fireIDraster[]))
  ids <- as.integer(values(fireIDraster, data.frame = FALSE)[loci])

  initialDf <- data.table(loci, ids, id = seq(ids))

  fireSize <- initialDf[, list(
    actualSize = .N,
    # simSize = .N,# needed for numIters
    goalSize = asInteger(pmax(minSize, am(.N)))
  ), by = "ids"]

  out <- list()
  simSizes <- initialDf[, list(simSize = .N), by = "ids"]
  simSizes <- fireSize[simSizes, on = "ids"]
  spreadProb <- 1
  it <- 1L
  maxIts <- 100L ## TODO: what's reasonable here???

  temp <- fireIDraster
  fireIDraster <- raster(fireIDraster)
  while ((length(loci) > 0) & (it <= maxIts)) {
    dups <- duplicated(loci)
    df <- data.table(loci = loci[!dups], ids = ids[!dups], id = seq_along(ids[!dups]))
    r1 <- SpaDES.tools::spread(fireIDraster,
                               loci = df$loci, iterations = 1,
                               spreadProb = flammableRTM, quick = TRUE,
                               returnIndices = TRUE)
    df <- df[r1, on = "id"]
    simSizes <- df[, list(simSize = .N), by = "ids"]
    simSizes <- fireSize[simSizes, on = "ids"]
    bigger <- simSizes$simSize > simSizes$goalSize


    if (any(bigger)) {
      idsBigger <- simSizes$ids[bigger]
      names(idsBigger) <- idsBigger
      out1 <- lapply(idsBigger, function(idBig) {
        wh <- which(df$ids %in% idBig)
        if (as.integer(verb) >= 2) {
          print(paste(
            "  Fire id:,", idBig, "finished. Num pixels in buffer:",
            simSizes[ids == idBig]$goalSize - simSizes[ids == idBig]$actualSize,
            ", in fire:", simSizes[ids == idBig]$actualSize
          ))
        }
        lastIters <- !df[wh]$active
        needMore <- simSizes[ids == idBig]$goalSize - sum(lastIters)
        if (needMore > 0) {
          dt <- try(rbindlist(list(
            df[wh][lastIters],
            df[wh][sample(which(df[wh]$active), needMore)]
          )))
        } else {
          dt <- df[wh][lastIters][sample(sum(lastIters), simSizes[ids == idBig]$goalSize)]
        }
        if (is(dt, "try-error")) browser() # stop("try error here")
        dtOut <- dt[, list(buffer = 0L, pixelID = indices, ids)]

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
      loci <- df$indices
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

  return(out3)
}

#' create a list of annual ignition points based on fire raster
#' @param fireBufferListDT a list of \code{data.table} with columns \code{buffer} (1 = burned),
#' \code{id} (unique fire ID), and \code{pixelID}
#' @param flammableRTM @template flammableRTM
#' @return a list of \code{SpatialPoints}
#' @export
#' @importFrom data.table rbindlist data.table set
#' @importFrom terra xyFromCell
#' @importFrom sf st_nearest_points st_centroid st_multipoint st_crs st_as_sf
rasterFireSpreadPoints <- function(fireBufferListDT) {
  #WIP
}
