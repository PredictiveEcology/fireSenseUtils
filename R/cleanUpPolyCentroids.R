utils::globalVariables(c(
  ".BY", ".I", "indices", "numNeighs"
))

#' Cleaning up the polygon points
#'
#' Mostly this is about 2 things: 1) remove fires that were so small that they take less
#' than 1 pixel so they are not in the \code{buff} object but are in the \code{cent}
#' object. 2) the centroid cell is in a buffer or otherwise nonburnable cell (e.g., water).
#' For 1) remove these from the centroid data. For 2) this function will search
#' in the neighbourhood for the next closest pixel that
#' has at least 7 available neighbours that can burn
#'
#' @param cent List of points as \code{SpatialPointsDataFrame}
#' @param idCol The column name as a character string with the fire ids. Defaults to
#'   \code{"NFIREID"}
#' @param buff List of \code{data.table} objects with 3 columns, "buffer" which is 1 (in the fire)
#'   or 0 (in a buffer), \code{ids} which are the fire ids which MUST match the ids
#'   in the \code{cent}
#' @param ras The raster that created the \code{pixelIDs} in the \code{buff}
#'
#' @export
#' @importFrom data.table as.data.table data.table setkeyv
#' @importFrom pemisc rasterToMatch
#' @importFrom purrr pmap
#' @importFrom sp SpatialPointsDataFrame spTransform
#' @importFrom SpaDES.tools distanceFromEachPoint spread
#' @importFrom raster cellFromXY compareCRS crs xyFromCell
#' @importFrom utils head tail
harmonizeBufferAndPoints <- function(cent, buff, ras, idCol = "FIRE_ID") {
  purrr::pmap(list(
    cent = cent,
    buff = buff
  ),
  .f = function(cent, buff) {
    if (!nrow(buff) > 0) { #cent can be >1 row while buff = 0, if poly is small
      return(NULL)
    }
    if (!compareCRS(crs(ras), crs(cent))) {
      cent <- sp::spTransform(cent, crs(ras))
    }

    whToUse <- cent[[idCol]] %in% buff$ids
    idsNotInBuffer <- cent[[idCol]][!whToUse]
    if (NROW(idsNotInBuffer) > 0) {
      polyCentroids <- cent[whToUse, ]
    } else {
      polyCentroids <- cent
    }
    inOrigFire <- buff[buff$buffer == 1, ]
    centDT <- data.table(
      pixelID = cellFromXY(spTransform(polyCentroids, crs(ras)), object = ras),
      ids = polyCentroids[[idCol]]
    )
    notInAFire <- centDT[!inOrigFire, on = c("pixelID")]
    if (NROW(notInAFire)) {
      inAFire <- buff[buffer == 1]
      fr <- cbind(xyFromCell(ras, inAFire$pixelID),
        id = inAFire$ids, pixelID = inAFire$pixelID
      )
      from <- cbind(id = notInAFire$ids, xyFromCell(ras, notInAFire$pixelID))
      dfep <- distanceFromEachPoint(from, fr)
      dfep <- as.data.table(dfep)
      # Make sure it is not surrounded by NAs

      setkeyv(dfep, c("id", "dists"))
      i <- 1
      replacementCentroids <- dfep[,
        list(centroidIndex = {
          if (.N > 1) {
            # if (i == 1) browseri <- 1
            notFound <- TRUE
            iter <- 1
            out1 <- maxSoFar <- integer()
            while (notFound) {
              spr <- spread(
                loci = tail(head(pixelID, iter * 20), 20), ras, spreadProb = 1, iterations = 1,
                allowOverlap = TRUE, returnIndices = TRUE
              )
              out <- spr[, list(numNeighs = sum(ras[][indices], na.rm = TRUE)), by = "id"][, numNeighs := numNeighs - 1]
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

      spOrig <- as.data.frame(polyCentroids[match(replacementCentroids$id, polyCentroids[[idCol]]), ])
      spOrig <- spOrig[
        match(replacementCentroids$id, spOrig[[idCol]]),
        grep("^x|y$|coords", names(spOrig), value = TRUE, invert = TRUE)
      ]
      sp <- SpatialPointsDataFrame(replacementCentroids[, c("x", "y")],
        spOrig,
        proj4string = crs(ras)
      )

      suppressWarnings({
        browser()
        polyCentroids <- rbind(polyCentroids[-match(replacementCentroids$id, polyCentroids[[idCol]]), ], sp)
      })
    }
    centDT2 <- data.table(
      pixelID = cellFromXY(spTransform(polyCentroids, crs(ras)), object = ras),
      ids = polyCentroids[[idCol]]
    )
    notInAFire <- centDT2[!inOrigFire, on = c("pixelID")]

    # if (NROW(notInAFire) > 0) browser()
    #fires will be rejected if centroid is outside

    polyCentroids
  })
}
