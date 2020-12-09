globalVariables(c(
  "age", "B", "BperClass", "burnClass", "propBurnClassFire", "totalBiomass"
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
#' @return DESCRIPTION NEEDED when this finally runs
#'
#' @export
#' @importFrom data.table copy
#' @importFrom SpaDES.tools rasterizeReduced
#' @importFrom raster getValues nlayers raster
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

  cohortData[, propBurnClassFire := BperClass/totalBiomass]
  # Fix 0/0
  cohortData[is.na(propBurnClassFire), propBurnClassFire := 0]
  # testthat::expect_true(NROW(cohortData) == NROW(na.omit(cohortData)))

  # Remove speciesCode so I can remove duplicates (i.e. different species that make the same class)
  toRemove <- names(cohortData)[!names(cohortData) %in% c("pixelGroup", "burnClass", "propBurnClassFire")]
  cohortData[, c(toRemove) := NULL]
  cohortData <- unique(cohortData)

  classList <- lapply(paste0("class", 1:4), function(cls){
    cohortDataub <- cohortData[burnClass == cls, ]
    ras <- rasterizeReduced(reduced = cohortDataub, fullRaster = pixelGroupMap,
                            newRasterCols = "propBurnClassFire",
                            mapcode = "pixelGroup")
    names(ras) <- paste0(cls, "_", yearCohort)
    return(ras)
  })

  # Identify non-forested pixels (non-ice/water/rocks) as class5
  # Pixels that are *NOT* NA in the RTM when this has been NA'ed for water, ice, and rocks, and
  # ARE NA in the pixelGroupMap are the pixels that are class5
  names(classList) <- paste0("fuelclass", 1:4)
  return(classList)
}
