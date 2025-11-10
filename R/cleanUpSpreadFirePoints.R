utils::globalVariables(c(
  "dists", "geometry", "isFlammable"
))

#' Ensure fire points are located on flammable pixels inside a fire polygon
#' Intended to be run using `Map`
#'
#' @param firePoints a `sf` points object representing annual ignitions
#'
#' @template flammableRTM
#'
#' @param bufferDT a data.table of burned cells, output from `bufferToArea`
#'
#' @return a list of harmonized points and polygons
#'
#' @export
#' @importFrom SpaDES.tools distanceFromEachPoint
#' @importFrom data.table as.data.table is.data.table data.table
#' @importFrom sf st_as_sf st_crs
#' @importFrom terra extract xyFromCell
cleanUpSpreadFirePoints <- function(firePoints, bufferDT, flammableRTM, idCol) {
  if (!is.data.table(bufferDT)) {
    bufferDT <- as.data.table(bufferDT)
  }
  FlamPoints <- as.data.table(extract(flammableRTM, firePoints, cells = TRUE))
  setnames(FlamPoints, c("ID", "isFlammable", "cells"))
  FlamPoints[, isFlammable := as.numeric(as.character(isFlammable))] # otherwise factor = 1 and 2
  FlamPoints[is.na(isFlammable), isFlammable := 0]
  if (any(FlamPoints$isFlammable == 0)) {
    badPoints <- FlamPoints[isFlammable == 0]
    badStarts <- firePoints[firePoints[[idCol]] %in% FlamPoints[isFlammable == 0]$ID, ]
    FlamPoints <- FlamPoints[isFlammable == 1, ] # keep the good ones
    polysWBadStarts <- bufferDT[ids %in% badPoints$ID, ]
    cells <- polysWBadStarts[polysWBadStarts$buffer == 1, ]
    flammableInPolys <- extract(flammableRTM, cells$pixelID)[, 1] == 1
    cells <- cells[flammableInPolys, ]

    badPolys <- setdiff(polysWBadStarts$ids, cells$ids) # truly bad; no flammable pixels anywhere in fire
    if (any(flammableInPolys, na.rm = TRUE)) {
      xyPolys <- cbind(
        id = cells$ids,
        pixelID = cells$pixelID,
        xyFromCell(flammableRTM, cells$pixelID)
      )
      xyPoints <- cbind(
        id = badPoints$ID,
        # pixelID = badStartsPixels[, "cells"],
        xyFromCell(flammableRTM, badPoints$cells)
      )
      dd <- as.data.table(distanceFromEachPoint(to = xyPolys, from = xyPoints))
      nearestPixels <- dd[, .SD[which.min(dists)], by = "id"]
      nearestPixels <- nearestPixels[, .(x, y, id)]
      badStarts <- as.data.table(badStarts)
      badStarts[, geometry := NULL]
      onVec <- c("id") |> setNames(idCol)
      badStarts <- badStarts[nearestPixels, on = onVec]
      newStarts <- st_as_sf(badStarts, coords = c("x", "y"), crs = st_crs(firePoints))
      firePoints <- firePoints[!firePoints[[idCol]] %in% newStarts[[idCol]], ]
      firePoints <- rbind(firePoints, newStarts)
    }
    # 2. rm bad points - those without ANY flammable pixels
    if (length(badPolys) > 0) {
      bufferDT <- bufferDT[!ids %in% badPolys]
      firePoints <- firePoints[!firePoints[[idCol]] %in% badPolys, ]
    }
  }
  list(SpatialPoints = firePoints, FireBuffered = bufferDT)
}
