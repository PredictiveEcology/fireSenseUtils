utils::globalVariables(c(
  ".BY", ".I", "indices", "numNeighs"
))

#' Cleaning up the polygon points
#'
#' Mostly this is about 2 things:
#' 1. remove fires that were so small that they take less than 1 pixel so they are
#' not in the \code{buff} object but are in the \code{cent} object.
#' 2. the centroid cell is in a buffer or otherwise nonburnable cell (e.g., water).
#' For 1) remove these from the centroid data.
#' For 2) this function will search in the neighbourhood for the next closest pixel
#' that has at least 7 available neighbours that can burn. If not, remove these.
#'
#' @param cent List of points as \code{SpatialPointsDataFrame}
#' @param idCol The column name as a character string with the fire ids.
#'   Defaults to \code{"FIRE_ID"}.
#' @param buff List of \code{data.table} objects with 3 columns, "buffer" which is 1 (in the fire)
#'   or 0 (in a buffer), \code{ids} which are the fire ids which MUST match the ids
#'   in the \code{cent}.
#' @param ras The raster that created the \code{pixelIDs} in the \code{buff}.
#'
#' @export
#' @importFrom data.table as.data.table data.table setkeyv set
#' @importFrom pemisc rasterToMatch
#' @importFrom purrr pmap
#' @importFrom sf st_as_sf st_crs "st_crs<-" st_transform st_coordinates
#' @importFrom LandR .compareCRS
#' @importFrom SpaDES.tools distanceFromEachPoint spread2
#' @importFrom terra cellFromXY xyFromCell
#' @importFrom utils head tail
harmonizeBufferAndPoints <- function(cent, buff, ras, idCol = "FIRE_ID") {
  purrr::pmap(list(
    cent = cent,
    buff = buff
  ),
  .f = function(cent, buff, fireIDcol = idCol) {
    if (!nrow(buff) > 0) { # cent can be >1 row while buff = 0, if poly is small
      return(NULL)
    }

    if (!.compareCRS(ras, cent)) {
      cent <- st_project(cent, st_crs(ras))
    }

    whToUse <- cent[[fireIDcol]] %in% buff$ids
    idsNotInBuffer <- cent[[fireIDcol]][!whToUse]
    if (NROW(idsNotInBuffer) > 0) {
      polyCentroids <- cent[whToUse, ]
    } else {
      polyCentroids <- cent
    }
    inOrigFire <- buff[buffer == 1, ]
    centDT <- data.table(
      pixelID = cellFromXY(st_coordinates(polyCentroids), object = ras),
      ids = polyCentroids[[fireIDcol]]
    )

    #for rbindlist, need to ensure col order
    colOrders <- names(cent)

    notInAFire <- centDT[!inOrigFire, on = c("pixelID")]
    if (NROW(notInAFire)) {

      inAFire <- buff[buffer == 1]
      fr <- cbind(xyFromCell(ras, inAFire$pixelID),
                  id = inAFire$ids, pixelID = inAFire$pixelID
      )
      from <- cbind(id = notInAFire$ids, xyFromCell(ras, notInAFire$pixelID))
      dfep <- distanceFromEachPoint(from, fr)
      dfep <- as.data.table(dfep)
      if (nrow(dfep) == 0) {
        return(NULL)
      }

      ## TODO: make sure it is not surrounded by NAs
      setkeyv(dfep, c("id", "dists"))
      i <- 1
      replacementCentroids <- dfep[,
                                   list(centroidIndex = {
                                     if (.N > 1) {
                                       notFound <- TRUE
                                       iter <- 1
                                       out1 <- maxSoFar <- integer()
                                       while (notFound) {
                                         spr <- SpaDES.tools::spread2(
                                           start = tail(head(pixelID, iter * 20), 20), landscape = ras,
                                           spreadProb = 1, iterations = 1, asRaster = FALSE,
                                           allowOverlap = TRUE)

                                         out <- spr[, list(numNeighs = sum(values(ras, mat = FALSE)[pixels],
                                                                           na.rm = TRUE)),
                                                    by = "initialPixels"][, numNeighs := numNeighs - 1]
                                         notFound <- (!any(out$numNeighs > 6))
                                         if (!notFound) {
                                           ind <- min(which(out$numNeighs > 6))
                                           out1 <- .I[ind]
                                         } else {
                                           iter <- iter + 1
                                           print(paste(.BY, ":", iter))
                                           ind <- which.max(out$numNeighs)
                                           maxSoFar <- c(maxSoFar, out$numNeighs[ind])
                                           out1 <- c(out1, .I[ind])
                                           # if (i == 1) browser()
                                           # if (.BY[[1]] == "706") browser()
                                           if (iter * 20 > .N) {
                                             notFound <- FALSE
                                             ind1 <- which.max(maxSoFar)
                                             out1 <- out1[ind1]
                                           }
                                         }
                                       }
                                     } else {
                                       out1 <- .I[1L]
                                     }
                                     out1
                                   }),
                                   by = "id"
      ]
      replacementCentroids <- dfep[replacementCentroids$centroidIndex]

      spOrig <- as.data.table(polyCentroids[match(replacementCentroids$id, polyCentroids[[fireIDcol]]), ])
      spOrig <- spOrig[
        match(replacementCentroids$id, spOrig[[fireIDcol]]),
        .SD, .SDcol = grep("^x|y$|coords", names(spOrig), value = TRUE, invert = TRUE)]

      #the join does not work with named argument for column, so assign to temp
      #we could use merge but I prefer this approach
      set(spOrig, NULL, "tempFoo", spOrig[[fireIDcol]])
      sp <- spOrig[replacementCentroids, on = c("tempFoo" = "id")]
      sp[, tempFoo := NULL]
      sp <- st_as_sf(sp, coords = c("x", "y"))
      st_crs(sp) <- st_crs(ras)

      #subset original out
      whichToKeep <- !polyCentroids[[fireIDcol]] %in% sp[[fireIDcol]]
      polyCentroids <- polyCentroids[whichToKeep,]
      set(sp, NULL, c("dists", "pixelID"), NULL)

      polyCentroids <- rbind(polyCentroids, sp)

    }
    #col order must be preserved

    polyCentroids <- setcolorder(polyCentroids, colOrders)
  }
  )
}
