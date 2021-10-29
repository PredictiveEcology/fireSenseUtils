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
#'
#' @return a raster stack of fuel classes defined by leading species and \code{sppEquiv$fuelClass}
#'
#' @export
#' @importFrom data.table copy setkey
#' @importFrom SpaDES.tools rasterizeReduced
#' @importFrom raster getValues nlayers raster stack
#'
cohortsToFuelClasses <- function(cohortData, yearCohort, pixelGroupMap, flammableRTM,
                                 sppEquiv, sppEquivCol, cutoffForYoungAge) {
  cohortData <- copy(cohortData)
  joinCol <- c('FuelClass', eval(sppEquivCol))
  sppEquivSubset <- unique(sppEquiv[, .SD, .SDcols = joinCol])

  cohortData <- cohortData[sppEquivSubset, on = c('speciesCode' = sppEquivCol)]
  #data.table needs an argument for which column names are kept during join

  cohortData[age <= cutoffForYoungAge, FuelClass := "youngAge"]
  if (!"totalBiomass" %in% names(cohortData))
    cohortData[, totalBiomass := asInteger(sum(B)), by = c("pixelGroup")]

  cohortData[, BperClass := sum(B), by = c("FuelClass", "pixelGroup")]

  # Fix zero age, zero biomass
  # testthat::expect_true(NROW(cohortData) == NROW(na.omit(cohortData)))
  cohortData[, Leading := BperClass == max(BperClass), .(pixelGroup)]
  cohortData <- unique(cohortData[Leading == TRUE, .(FuelClass, pixelGroup)]) #unique due to species
  cohortData[, N := .N, .(pixelGroup)]

  #In the event of a tie, we randomly pick a fuel class
  ties <- cohortData[N > 1]
  noTies <- cohortData[N == 1]
  if (nrow(ties) > 1) {
    ties$foo <- sample(x = 1:nrow(ties), size = nrow(ties))
    setkey(ties, foo)
    ties <- ties[!duplicated(ties[, .(pixelGroup)])]
    ties[, foo := NULL]
    cohortData <- rbind(ties, noTies)
  } else {
    cohortData <- noTies
  }
  classes <- sort(unique(cohortData$FuelClass))
  classList <- lapply(classes, makeRastersFromCD,
                      flammableRTM =  flammableRTM,
                      pixelGroupMap = pixelGroupMap,
                      cohortData = cohortData)

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
