globalVariables(c(
  "age", "B", "BperClass", "burnClass", "foo", "Leading", "LeaderValue", "propBurnClassFire",
  "totalBiomass"
))

#' Classify \code{pixelGroups} by flammability
#'
#' @template cohortData
#' @template pixelGroupMap
#' @template sppEquiv
#' @template sppEquivCol
#' @param yearCohort the year the \code{cohortData} represents
#' @param flammableMap  binary map of flammable pixels - see \code{LandR::defineFlammable}
#'
#' @return a raster stack of fuel classes defined by leading species and \code{sppEquiv$fuelClass}
#'
#' @export
#' @importFrom data.table copy setkey
#' @importFrom SpaDES.tools rasterizeReduced
#' @importFrom raster getValues nlayers raster stack
#'
classifyCohortsFireSenseSpread <- function(cohortData, yearCohort, pixelGroupMap, flammableMap,
                                           sppEquiv, sppEquivCol) {

  cohortData <- copy(cohortData)
  joinCol <- c('FuelClass', eval(sppEquivCol))
  sppEquivSubset <- sppEquiv[, .SD, .SDcols = joinCol]

  cohortData <- cohortData[sppEquivSubset, on = c('speciesCode' = sppEquivCol)]
  #data.table needs an argument for which column names are kept during join
  setnames(cohortData, 'FuelClass', 'burnClass')
  cohortData[age < 15, burnClass := "class1"]
  if (!"totalBiomass" %in% names(cohortData))
    cohortData[, totalBiomass := asInteger(sum(B)), by = c("pixelGroup")]

  #Assertion - no longer accurate with no class 5
  # if (!((NROW(cohortData[is.na(totalBiomass) & burnClass != "class5", ]) == 0) &
  #       (NROW(cohortData[!is.na(totalBiomass) & burnClass == "class5", ]) == 0))) {
  #   stop('there is a problem setting the burn class. contact module developers')
  # }

  cohortData[, BperClass := sum(B), by = c("burnClass", "pixelGroup")]

  # Fix zero age, zero biomass
  # testthat::expect_true(NROW(cohortData) == NROW(na.omit(cohortData)))
  cohortData[, Leading := BperClass == max(BperClass), .(pixelGroup)]
  cohortData <- unique(cohortData[Leading == TRUE, .(burnClass, pixelGroup)]) #unique due to species
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
  classes <-sort(unique(cohortData$burnClass))
  classList <- lapply(classes, function(cls){
    cohortDataub <- cohortData[burnClass == cls, ]
    cohortDataub[,  LeaderValue := 1]
    ras <- rasterizeReduced(reduced = cohortDataub, fullRaster = pixelGroupMap,
                            newRasterCols = "LeaderValue",
                            mapcode = "pixelGroup")
    #need to make sure fuel class is 0 and not NA if absent entirely
    #to prevent NAs following aggregation to 25 km pixels
    ras[!is.na(flammableMap[]) & flammableMap[] == 1 & is.na(ras[])] <- 0
    return(ras)
  })

  classList <- stack(classList)
  names(classList) <- classes
  return(classList)
}
