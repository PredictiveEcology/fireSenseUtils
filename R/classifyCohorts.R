globalVariables(c(
  "age", "B", "BperClass", "FuelClass", "foo", "Leading", "LeaderValue",
  "totalBiomass"
))

#' Classify \code{pixelGroups} by flammability
#'
#' @template cohortData
#' @template pixelGroupMap
#' @template sppEquiv
#' @template sppEquivCol
#' @param yearCohort the year the \code{cohortData} represents
#' @template flammableRTM
#' @param cutoffForYoungAge age at and below which pixels are considered 'young'
#' @param fuelClassCol the column in sppEquiv that describes unique fuel classes
#'
#' @return a raster stack of fuel classes defined by leading species and \code{sppEquiv$fuelClass}
#'
#' @export
#' @importFrom data.table copy setkey
#' @importFrom SpaDES.tools rasterizeReduced
#' @importFrom raster getValues nlayers raster stack
#'
cohortsToFuelClasses <- function(cohortData, yearCohort, pixelGroupMap, flammableRTM,
                                 sppEquiv, sppEquivCol, cutoffForYoungAge, fuelClassCol = "FuelClass") {
  cD <- copy(cohortData)
  joinCol <- c(fuelClassCol, eval(sppEquivCol))
  sppEquivSubset <- unique(sppEquiv[, .SD, .SDcols = joinCol])
  cD <- cD[sppEquivSubset, on = c('speciesCode' = sppEquivCol)]
  setnames(cD, old = fuelClassCol, new = "FuelClass") #so we don't have to use eval, which trips up some dt
  #data.table needs an argument for which column names are kept during join

  cD[age <= cutoffForYoungAge, FuelClass := "youngAge"]
  if (!"totalBiomass" %in% names(cD))
    cD[, totalBiomass := asInteger(sum(B)), by = c("pixelGroup")]

  cD[, BperClass := sum(B), by = c("FuelClass", "pixelGroup")]

  # Fix zero age, zero biomass
  cD[, Leading := BperClass == max(BperClass), .(pixelGroup)]
  cD <- unique(cD[Leading == TRUE, .(FuelClass, pixelGroup)]) # unique due to species
  cD[, N := .N, .(pixelGroup)]

  #In the event of a tie, we randomly pick a fuel class
  ties <- cD[N > 1]
  noTies <- cD[N == 1]
  if (nrow(ties) > 1) {
    ties$foo <- sample(x = 1:nrow(ties), size = nrow(ties))
    setkey(ties, foo)
    ties <- ties[!duplicated(ties[, .(pixelGroup)])]
    ties[, foo := NULL]
    cD <- rbind(ties, noTies)
  } else {
    cD <- noTies
  }
  classes <- sort(unique(cD$FuelClass))
  classList <- lapply(classes, makeRastersFromCD,
                      flammableRTM =  flammableRTM,
                      pixelGroupMap = pixelGroupMap,
                      cohortData = cD)

  classList <- stack(classList)
  names(classList) <- classes
  return(classList)
}

#' put cohortData back into raster with some extra details
#' @param class fuelClass from  \code{sppEquiv}
#' @template flammableRTM
#' @template pixelGroupMap
#' @template cohortData
#' @return a raster
#' @importFrom SpaDES.tools rasterizeReduced
#' @importFrom data.table data.table
makeRastersFromCD <- function(class, cohortData, flammableRTM, pixelGroupMap) {

  cohortDataub <- cohortData[FuelClass == class, ]
  cohortDataub[,  LeaderValue := 1]
  ras <- rasterizeReduced(reduced = cohortDataub, fullRaster = pixelGroupMap,
                          newRasterCols = "LeaderValue",
                          mapcode = "pixelGroup")
  #fuel class is 0 and not NA if absent entirely
  #to prevent NAs following aggregation to 25 km pixels
  ras[!is.na(flammableRTM[]) & is.na(ras[])] <- 0
  return(ras)

}
