utils::globalVariables(c(
  "pixelID", "fireID", "ids", "id", ".N", "x", "y"
))

#' this is a wrapper to simplify caching of lapply with bufferForFireRaster.
#' Years are iteratively processed by \code{makeFireID}.
#' @param years numeric fire years
#' @param fireRaster a \code{SpatRaster} with values representing fire years
#' @template flammableRTM
#' @param bufferForFireRaster buffer size used to group discrete patches of burned pixels as
#' belonging to the same fire
#' @param areaMultiplier A scalar that will buffer \code{areaMultiplier * fireSize}
#' @param verb Logical or numeric related to how much verbosity is printed. \code{FALSE} or
#'  \code{0} is none. \code{TRUE} or \code{1} is some. \code{2} is much more.
#' @param minSize The absolute minimum size of the buffer & non-buffer together. This will
#'   be imposed after \code{areaMultiplier}.
#' @param cores number of processor cores to use
#' @return a list of data.table named by year, with cols \code{ids}, \code{buffer},
#' and \code{pixelID}
#' @export
#' @importFrom parallelly availableCores
rasterFireBufferDT <- function(years, fireRaster, flammableRTM, bufferForFireRaster, areaMultiplier,
                               minSize = 5000, verb = 1, cores = 1) {

  maxCores <- parallelly::availableCores(constraints = "connections", omit = 1)
  cores <- min(min(length(years), cores), maxCores)
  if (cores > 1) {
    out <- parallel::mclapply(years, FUN = makeFireIDs, fireRaster = fireRaster,
                              flammableRTM = flammableRTM, bufferForFireRaster = bufferForFireRaster,
                              areaMultiplier = areaMultiplier, minSize = minSize, verb = verb)
  } else {
    fireBufferListDT <- lapply(years, makeFireIDs, fireRaster = fireRaster, flammableRTM = flammableRTM,
                               bufferForFireRaster = bufferForFireRaster,
                               areaMultiplier = areaMultiplier, minSize = minSize, verb = verb)
  }
  names(fireBufferListDT) <- paste0("year", years)
  return(fireBufferListDT)
}

#' identify each year's individual fires and buffer them accordingly
#' @param year numeric fire year
#' @param fireRaster a \code{SpatRaster} with values representing fire years
#' @template flammableRTM
#' @param bufferForFireRaster buffer size used to group discrete patches of burned pixels as
#' belonging to the same fire
#' @param areaMultiplier A scalar that will buffer \code{areaMultiplier * fireSize}
#' @param verb Logical or numeric related to how much verbosity is printed. \code{FALSE} or
#'   \code{0} is none. \code{TRUE} or \code{1} is some. \code{2} is much more.
#' @param minSize The absolute minimum size of the buffer & non-buffer together. This will
#'   be imposed after \code{areaMultiplier}.
#' @return a data.table with fire ID, buffer status, and pixelID
#' @export
#' @importFrom data.table rbindlist data.table set setcolorder
#' @importFrom terra buffer patches
makeFireIDs <- function(year, fireRaster, flammableRTM, bufferForFireRaster, areaMultiplier, minSize = 5000, verb = 1) {

  #identify this year's fires
  origIndex <- which(fireRaster[] == year)
  if (length(origIndex) == 0) {
    return(NULL) #TODO: what should happen here? empty data.table
  } else{
    fireRaster[] <- NA
    fireRaster[origIndex] <- 1
    #classify into discrete fires via buffer
    fireRaster <- buffer(fireRaster, width = bufferForFireRaster, background = NA)
    #TODO: doublecheck that no patches are lost... why are patch values greater than # of patches
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
    #match the polygonal order
    setcolorder(fireBufferDT, c("pixelID", "buffer", "ids"))
    return(fireBufferDT)
  }
}

#' create a variable sized buffer around a set of pixels belonging to the same fire ID
#' @param fireIDraster a \code{SpatRaster} with values representing distinct fires in a year
#' @param flammableRTM @template flammableRTM
#' @param areaMultiplier A scalar that will buffer \code{areaMultiplier * fireSize}
#' @param minSize The absolute minimum size of the buffer & non-buffer together. This will
#'   be imposed after \code{areaMultiplier}.
#' @param verb Logical or numeric related to how much verbosity is printed. \code{FALSE} or
#'  \code{0} is none. \code{TRUE} or \code{1} is some. \code{2} is much more.
#' @return a data.table with fire ID, buffer status, and pixelID
#' @export
#' @importFrom data.table rbindlist data.table set setnames
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
  emptyDT <- data.table(pixelID = integer(0), buffer = integer(0), ids = integer(0))

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
#' @param fireBufferDT a \code{data.table} with columns \code{buffer} (1 = burned),
#' \code{id} (unique fire ID), and \code{pixelID}
#' @param flammableRTM @template flammableRTM
#' @return a list of \code{sf} point objects
#' @export
#' @importFrom data.table rbindlist data.table set
#' @importFrom terra xyFromCell rast res
#' @importFrom sf st_nearest_points st_centroid st_multipoint st_crs  "st_crs<-" st_as_sf
#' @importFrom data.table as.data.table
rasterFireSpreadPoints <- function(fireBufferDT, flammableRTM) {

  if (inherits(flammableRTM, "RasterLayer")) {
    flammableRTM <- rast(flammableRTM)
  }

  burnedPoints <- fireBufferDT[buffer == 1,]
  burnLocs <- xyFromCell(flammableRTM, cell = burnedPoints$pixelID)
  burnLocs <- as.data.table(burnLocs)
  burnLocs$ids <- burnedPoints$ids

  #TODO - we can use st_centroid without the need to join if we are able to
  burnPoints <- st_as_sf(burnLocs, coords = c("x", "y"))
  st_crs(burnPoints) <- st_crs(flammableRTM)

  burnCentroidDT <- burnLocs[, .(x = mean(x), y = mean(y)), .(ids)]

  burnCentroids <- st_as_sf(burnCentroidDT, coords = c("x", "y"))
  st_crs(burnCentroids) <- st_crs(flammableRTM)

  #there is likely a faster way?
  ids <- unique(burnCentroidDT$id)
  firePoints <- lapply(ids, FUN = function(x, p = burnPoints, c = burnCentroids){
    p <- p[p$ids == x,]
    c <- c[c$ids == x,]
    #gives nearest y to each x - so nearest point to each centroid
    nearest <- sf::st_nearest_feature(x = c, y = p)
    p <- p[nearest,]
    return(p)
  })

  firePoints <- rbindlist(firePoints)
  #calculate size in hectares
  multiplier <- prod(res(flammableRTM))/10000
  burnSizes <-fireBufferDT[buffer == 1, .(POLY_HA = .N * multiplier, size = .N), .(ids)]
  firePoints <- firePoints[burnSizes, on = c("ids")]
  setnames(firePoints, old = "ids", new = "FIRE_ID") #following fireSenseUtils::makeLociList
  firePoints <- st_as_sf(firePoints)
  return(firePoints)
}
