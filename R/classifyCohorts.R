globalVariables(c(
  "age", "B", "BperClass", "FuelClass", "foo", "Leading", "LeaderValue",
  "maxAge", "NspeciesWithMaxB", "pixelIndex", "totalBiomass"
))

#' Classify `pixelGroups` by flammability
#'
#' @template cohortData
#'
#' @template pixelGroupMap
#'
#' @template sppEquiv
#'
#' @template sppEquivCol
#'
#' @param landcoverDT Optional table of non-forest land cover classes and pixel indices.
#'                    It will override pixel values in `cohortData`, if supplied.
#'
#' @template flammableRTM
#'
#' @param cutoffForYoungAge age at and below which pixels are considered 'young'
#'
#' @param fuelClassCol the column in `sppEquiv` that describes unique fuel classes
#' 
#' @return a `SpatRaster` of biomass by fuel class as determined by `fuelClassCol` and `cohortData`.
#'
#' @export
#' @inheritParams fireSenseCovariatesCreate
#' @importFrom data.table copy setkey
#' @importFrom LandR asInteger
#' @importFrom SpaDES.tools rasterizeReduced
#' @importFrom terra values rast
#'
cohortsToFuelClasses <- function(cohortData, pixelGroupMap, flammableRTM, landcoverDT = NULL,
                                 sppEquiv, sppEquivCol, cutoffForYoungAge, fuelClassCol = "FuelClass",
                                 requiredFuelClasses) {
  # cD <- copy(cohortData)
  joinCol <- c(fuelClassCol, eval(sppEquivCol))
  sppEquivSubset <- unique(sppEquiv[, .SD, .SDcols = joinCol])
  cD <- cohortData[sppEquivSubset, on = c("speciesCode" = sppEquivCol)]
  setnames(cD, old = fuelClassCol, new = "FuelClass") # so we don't have to use eval, which trips up some dt
  # data.table needs an argument for which column names are kept during join
  cD[, maxAge := max(age), .(pixelGroup)]
  cD[maxAge <= cutoffForYoungAge, FuelClass := youngAgeTxt]
  cD[, maxAge := NULL]
  cD <- cD[, .(BperClass = asInteger(sum(B))), by = c("FuelClass", "pixelGroup")]

  # youngAge is better treated as a binary cover variable than continuous measure of biomass
  cD[FuelClass == youngAgeTxt, BperClass := 1]

  # Fix zero age, zero biomass
  classes <- sort(unique(cD$FuelClass))
  
  # st <- profvis::profvis({
  
  pgmVals <- list(pixelGroup = values(pixelGroupMap, mat = FALSE), 
                  pixelId = seq(ncell(pixelGroupMap))) |> 
    setDT() |> na.omit()
  aa <- pgmVals[cD, on = "pixelGroup", allow.cartesian=TRUE]
  bb <- split(aa, by = "FuelClass")
  flamVals <- values(flammableRTM, mat = FALSE)
  flamValsGood <- !is.na(flamVals)
  cc <- Map(r = bb, function(r) {
    ras <- rastFromDF(r[, .(pixelId, BperClass)], rasTemplate = pixelGroupMap)
    rasVals <- values(ras, mat = FALSE)
    rasVals[flamValsGood & is.na(rasVals)] <- 0
    ras <- setValues(x = ras, values = rasVals)
    ras
    })
  
  noFuelForRequiredClass <- character()
  if (!is.null(requiredFuelClasses))
    noFuelForRequiredClass <- setdiff(requiredFuelClasses, names(cc))

  # This is where a species disappears from the map: create a map of zeros
  if (length(noFuelForRequiredClass)) {
    for (fuel in noFuelForRequiredClass) {
      ras <- Copy(pixelGroupMap)
      rasVals <- values(ras, mat = FALSE)
      rasVals[rasVals > 0] <- 0
      cc[[fuel]] <- setValues(x = ras, values = rasVals)
    }
    
  }
  dd <- rast(cc)
  classList <- dd[[order(names(dd))]]
  # })
  
  
  # st2 <- profvis::profvis({
  # classList <- lapply(classes, makeRastersFromCD,
  #   flammableRTM = flammableRTM,
  #   pixelGroupMap = pixelGroupMap,
  #   cohortData = cD
  # )
  # 
  # classList <- rast(classList)
  # })
  if (!is.null(landcoverDT)) {
    # find rows that aren't empty i.e. have non-forest landcover
    landcoverDT[, foo := rowSums(.SD, na.rm = TRUE), .SDcols = setdiff(names(landcoverDT), nonNFColNamesTxt)]
    # terra needs protection from zero-length index
    if (nrow(landcoverDT[foo > 0, ]) > 0) {
      classList[landcoverDT[foo > 0]$pixelID] <- 0 # must be 0
    }
    landcoverDT[, foo := NULL]
  }

  names(classList) <- classes
  
  # Need to confirm that there was at least 1 youngAge ... sometime there are none e.g., with 9.2.1 plains
  if (!youngAgeTxt %in% names(classList)) {
    ya <- as.int(is.na(cc[[1]]) )
    vals <- values(ya, mat = FALSE)
    ya[vals == 1L] <- NA
    names(ya) <- youngAgeTxt
    classList <- c(classList, ya)
  }
  return(classList)
}

#' Put `cohortData` back into a `SpatRaster` with some extra details
#'
#' @param class `fuelClass` from  `sppEquiv`
#'
#' @template flammableRTM
#'
#' @template pixelGroupMap
#'
#' @template cohortData
#'
#' @return a `SpatRaster` with values equal to `class` biomass (B)
#'
#' @importFrom data.table data.table
#' @importFrom SpaDES.tools rasterizeReduced
#' @importFrom terra values setValues rast
makeRastersFromCD <- function(class, cohortData, flammableRTM, pixelGroupMap) {
  cohortDataFB <- cohortData[FuelClass == class, ]
  ras <- rasterizeReduced(
    reduced = cohortDataFB,
    fullRaster = pixelGroupMap,
    newRasterCols = "BperClass",
    mapcode = "pixelGroup"
  )
  ## fuel class is 0 and not NA if absent entirely
  ## to prevent NAs following aggregation during ignitionFit
  flamVals <- values(flammableRTM, mat = FALSE)
  rasVals <- values(ras, mat = FALSE)
  rasVals[!is.na(flamVals) & is.na(rasVals)] <- 0
  ras <- setValues(x = ras, values = rasVals)
  return(ras)
}

#' Modify `cohortData` with burn column
#'
#' @param year length-two vector giving temporal period used to subset `firePolys`. Closed interval.
#'
#' @param pixelGroupMap either a `SpatRaster` with `pixelGroups` or list of `SpatRasters` named by year.
#'
#' @param cohortData either a `cohortData` object or list of `cohortData` objects named by year.
#'
#' @param firePolys the output of [getFirePolygons()] with `YEAR` column.
#'
#' @return `cohortData` modified with burn status
#'
#' @importFrom LandR addPixels2CohortData
#' @importFrom terra rasterize
buildCohortBurnHistory <- function(cohortData, pixelGroupMap, firePolys, year) {
  ## build fire raster
  firePolys <- do.call(rbind, firePolys)
  firePolys <- firePolys[firePolys$YEAR >= min(year) & firePolys$YEAR <= max(year), ]
  fireRas <- rasterize(firePolys, pixelGroupMap, field = "YEAR", fun = min)
  cdLong <- addPixels2CohortData(cohortData, pixelGroupMap)
  cdLong[, burned := fireRas[pixelIndex]]
  return(cdLong)
}

