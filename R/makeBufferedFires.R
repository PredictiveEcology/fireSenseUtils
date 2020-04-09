utils::globalVariables(c(".SD", "N"))

#' Create buffers around fires
#'
#' @param fireLocationsPolys DESCRIPTION NEEDED
#' @param rasterToMatch DESCRIPTION NEEDED
#' @param lowerTolerance DESCRIPTION NEEDED
#' @param upperTolerance DESCRIPTION NEEDED
#' @param verbose DESCRIPTION NEEDED
#' @param useParallel DESCRIPTION NEEDED
#'
#' @return DESCRIPTION NEEDED
#'
#' @export
#' @importFrom crayon green red yellow
#' @importFrom data.table data.table
#' @importFrom fasterize fasterize
#' @importFrom SpaDES.tools adj
#' @importFrom future.apply future_lapply
#' @importFrom raster crs raster setValues stack
#' @importFrom reproducible projectInputs
#' @importFrom sf st_as_sf
makeBufferedFires <- function(fireLocationsPolys,
                              rasterToMatch,
                              lowerTolerance = 3.8,
                              upperTolerance = 4.2,
                              verbose = getOption("verbose", TRUE),
                              useParallel = FALSE) {
  t1 <- Sys.time()
  # fireLocationsPolys: list of each year of SpatialPolygonsDataFrame with fire polygons
  # lowerTolerance: lower tolerance for buffer to be different from fire points (i.e. 0.8, 20% lower)
  # upperTolerance: higher tolerance for buffer to be different from fire points (i.e. 1.2, 20% higher)
  data.table::setDTthreads(1)
  if (!quickPlot::isRstudioServer() && isTRUE(useParallel) && is(plan(), "sequential"))
    future:::plan(multiprocess(workers = min(detectCores(), length(fireLocationsPolys))))


  fun <- ifelse(useParallel, future.apply::future_lapply, lapply)
  historicalFire <- do.call(what = fun, args = list(
    X = names(fireLocationsPolys),
    FUN = function(yr){
      # Projection is not the same, so I need to convert the polygon

      data.table::setDTthreads(1)
      fireLocationsPoly <- reproducible::projectInputs(
        x = fireLocationsPolys[[yr]],
        targetCRS = crs(rasterToMatch)
      )
      sf_fPY <- sf::st_as_sf(fireLocationsPoly)
      firePolyRas <- fasterize::fasterize(sf = sf_fPY, raster = raster(rasterToMatch), field = "NFIREID")
      names(firePolyRas) <- yr
      # Do the calculation for each fire
      fireIDS <- unique(firePolyRas[!is.na(firePolyRas)])
      allFires <- lapply(X = fireIDS, FUN = function(fireID){
        valsFireRas  <- which(firePolyRas[] == fireID)
        adj <- SpaDES.tools::adj(firePolyRas,
                                 valsFireRas,
                                 directions = 8,
                                 pairs = FALSE, match.adjacent = TRUE, cutoff.for.data.table = Inf)
        tb <- data.table(V1 = c(0, 1), N = c(1, 2))
        stateToCheck <- paste0("isFALSE(isTRUE(tb[1, N] < tb[2, N]*upperTolerance) & ",
                               "isTRUE(tb[1, N] > tb[2, N]*lowerTolerance))")
        while (eval(parse(text = stateToCheck))) {
          adj <- SpaDES.tools::adj(firePolyRas, adj, directions = 8, pairs = FALSE, match.adjacent = TRUE, cutoff.for.data.table = Inf)
          rasBuffer <- raster(firePolyRas)
          suppressWarnings(rasBuffer[adj] <- 0)
          suppressWarnings(rasBuffer[valsFireRas] <- 1)
          tb <- data.table(table(rasBuffer[]))
          perc <- round((tb[1, N] / tb[2, N]) * 100, 0)
          direcion <- ifelse(perc > 100, "larger", "smaller")
          perc <- ifelse(direcion == "larger", perc - 100, 100 - perc)
          if (verbose) {
            if (all(direcion == "larger", perc / 100 > upperTolerance - 1)) {
              message(crayon::red(paste0(
                "Buffered area is ", perc, "% ", direcion, " than fires.",
                " Outside of bounds, returning for ", yr
              )))
            } else {
              if (eval(parse(text = stateToCheck))) {
                message(crayon::yellow(paste0(
                  "Buffered area is ", perc, "% ", direcion, " than fires.",
                  " Trying again for ", yr
                )))
              } else {
                message(crayon::green(paste0(
                  "Buffered area is ", perc,
                  "% ", direcion, " than fires. Within bounds for ",
                  yr
                )))
              }
            }
          }
          if (tb[1, N] > tb[2, N] * upperTolerance) break
        }
        if (is.null(rasBuffer)){
          print("NULL raster? Debug")
          browser()
        }
        attr(rasBuffer, "buffer") <- adj
        return(rasBuffer)
      })
      # Convert to data table to speed up putting the rasters back together
      adjAll <- unlist(lapply(allFires, attr, which = "buffer"))
      allFires <- lapply(allFires, function(ras){
        attr(ras, "buffer") <- NULL
        return(ras)
      })
      stk <- stack(allFires)
      stkDT <- as.data.table(stk[])
      stkDT[, fires := rowSums(.SD, na.rm = TRUE)]
      rasBuffer <- setValues(raster(firePolyRas), stkDT$fires)
      # Put NA's back
      suppressWarnings(rasBuffer[rasBuffer > 1] <- 1) #TODO Just making sure all fires are 1 (the same pixel might have burned                                     more then one time). This whould have already happened so here it could be a test
      fires  <- which(rasBuffer[] == 1)
      suppressWarnings(rasBuffer[] <- NA)
      suppressWarnings(rasBuffer[adjAll] <- 0) # Buffer
      suppressWarnings(rasBuffer[fires] <- 1) # Fires
      return(rasBuffer)
    }))
  message(crayon::green(paste0("Finished fire buffering")))
  print(Sys.time() - t1)
  names(historicalFire) <- names(fireLocationsPolys)
  if (verbose) {
    print(Sys.time() - t1)
  }
  return(historicalFire)
}

#' Simplify buffered fires
#'
#' DESCRIPTION NEEDED
#'
#' @param fireBuffered DESCRIPTION NEEDED
#'
#' @return DESCRIPTION NEEDED
#'
#' @export
#' @importFrom data.table data.table
simplifyFireBuffered <- function(fireBuffered) {
  lapply(fireBuffered, function(r) {
    ras <- raster(r)
    nonNA <- which(!is.na(r[]))
    ras[r[] == 0] <- 0L
    ras[r[] == 1] <- 1L
    data.table(buffer = ras[][nonNA], pixelID = nonNA)
  })
}
