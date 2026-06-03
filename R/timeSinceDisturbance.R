utils::globalVariables(c(
  "sumRows"
))
#' Prepare a time since disturbance map from stand age and fire data
#'
#' Combines an initial stand age map with disturbance (fire) history to produce a
#' time-since-disturbance (TSD) raster. Stand age is trusted for most pixels, but
#' for pixels whose stand age is unreliable (e.g. non-forest land cover) the fire
#' history is used instead: a recorded recent burn sets the TSD, while the absence
#' of one marks the pixel as old (`cutoffForYoungAge + 1`).
#'
#' Which pixels to update from fire history, and which pixels are flammable, can be
#' supplied in one of two ways:
#' * directly, via `pixToUpdate` and `flammablePixels` (general purpose); or
#' * via `lcc`, the `landcoverDT` produced by `fireSense_dataPrepFit`, from which
#'   both are derived (non-forest pixels are updated, and `lcc$pixelID` defines the
#'   flammable mask). This is retained for backwards compatibility.
#'
#' @param standAgeMap initial stand age map
#'
#' @param firePolys list of `spatialPolygon` objects comprising annual fires.
#'   `fireRaster` will supersede `firePolys` if provided.
#'
#' @param fireRaster a `RasterLayer` with values representing fire years.
#'
#' @param year the year represented by `standAge`.
#'
#' @param lcc Optional `data.table` with landcover values, i.e., `landcoverDT`
#'   (see [makeLandcoverDT()]). When supplied, `pixToUpdate` and `flammablePixels`
#'   are derived from it unless they are passed explicitly. Specific to
#'   `fireSense_dataPrepFit`; for general use prefer `pixToUpdate` /
#'   `flammablePixels`.
#'
#' @param pixToUpdate Optional integer vector of `pixelID`s whose age should be
#'   taken from fire history rather than `standAgeMap` (e.g. pixels with no reliable
#'   stand age). If `NULL` and `lcc` is supplied, these are the non-forest pixels in
#'   `lcc` (those with a positive sum across the non-forest landcover columns).
#'
#' @param flammablePixels Optional integer vector of flammable `pixelID`s. Pixels
#'   not in this set are set to `NA` in the output (non-flammable). If `NULL` and
#'   `lcc` is supplied, this is `lcc$pixelID`. If `NULL` and `lcc` is not supplied,
#'   no pixels are masked to `NA`.
#'
#' @inheritParams castCohortData
#'
#' @return a `SpatRaster` with values representing time since disturbance
#'
#' @export
#' @importFrom data.table data.table as.data.table
#' @importFrom terra values rast setValues rasterize vect set.names
makeTSD <- function(year, firePolys = NULL, fireRaster = NULL,
                    standAgeMap, lcc = NULL, cutoffForYoungAge = 15,
                    pixToUpdate = NULL, flammablePixels = NULL) {
  if (!is.null(fireRaster)) {
    baseYear <- rast(fireRaster)
    baseYear <- setValues(baseYear, year)
    initialTSD <- baseYear - fireRaster
    initialTSD[initialTSD < 0] <- cutoffForYoungAge + 1
    ## these pixels burn in the future - can't infer prior disturbance
  } else if (!is.null(firePolys)) {
    ## get particular fire polys in format that can be fasterized
    polysNeeded <- firePolys[names(firePolys) %in% paste0("year", c(year - cutoffForYoungAge - 1):year - 1)]
    polysNeeded <- polysNeeded[sapply(polysNeeded, length) > 0]
    # polysNeeded <- vect(polysNeeded) #terrarize
    # polysNeeded <- do.call(rbind, polysNeeded)
    polysNeeded <- Reduce(rbind, polysNeeded)
    ## create background raster with TSD
    initialTSD <- rasterize(polysNeeded,
      y = standAgeMap,
      background = year - cutoffForYoungAge - 1,
      field = "YEAR", fun = "max"
    ) |> terra::mask(standAgeMap)
    initialTSD <- year - initialTSD
  } else {
    stop("Please provide either firePolys or fireRaster")
  }

  ## Derive `pixToUpdate` (pixels to age from fire history) and `flammablePixels`
  ## (the non-NA mask) from `lcc` when they are not supplied directly. These are
  ## the operations specific to `fireSense_dataPrepFit`'s `landcoverDT`; explicit
  ## arguments take precedence so the function can be used without an `lcc`.
  if (!is.null(lcc)) {
    nfLCC <- names(lcc)[!names(lcc) %in% nonNFColNamesTxt]
    lcc[, sumRows := rowSums(.SD, na.rm = TRUE), .SDcols = nfLCC]
    if (is.null(pixToUpdate)) {
      pixToUpdate <- lcc[sumRows > 0]$pixelID |> unique()
    }
    lcc[, sumRows := NULL]
    if (is.null(flammablePixels)) {
      flammablePixels <- lcc$pixelID
    }
  }

  standAgeVals <- data.table(pixelID = 1:ncell(standAgeMap), age = values(standAgeMap, mat = FALSE))
  # data.table::setnames(standAgeVals, c("pixelID", "age"))

  # standAgeVals <- values(standAgeMap, mat = FALSE)
  ## these have no disturbance history but are apparently young

  falseYoungs <- standAgeVals[c(pixelID %in% pixToUpdate & age <= cutoffForYoungAge) |
                                c(pixelID %in% pixToUpdate & is.na(age))]$pixelID
  ## disturbance history suggests young

  ## Read TSD as a plain numeric vector once. `which()` drops NA entries so a
  ## missing fireRaster value (no recorded disturbance) is not flagged as young
  ## and trueYoungs / trueAges stay length-locked for the := assignment.
  tsdAtPix <- terra::values(initialTSD, mat = FALSE)[pixToUpdate]
  youngPos <- which(tsdAtPix <= cutoffForYoungAge)
  trueYoungs <- pixToUpdate[youngPos]
  trueAges <- tsdAtPix[youngPos]
  standAgeVals[pixelID %in% falseYoungs, age := cutoffForYoungAge + 1]
  ## note that by doing this second, pixels in both groups are correctly set to trueYoung

  standAgeVals[pixelID %in% trueYoungs, age := trueAges]
  if (!is.null(flammablePixels)) {
    standAgeVals[!pixelID %in% flammablePixels, age := NA] #these should be NA - not flammable
  }

  standAgeMap <- setValues(standAgeMap, standAgeVals$age)
  set.names(standAgeMap, paste0("timeSinceDisturbance", year))

  return(standAgeMap)
}

#' Iteratively calculate `youngAge` column in FS covariates
#'
#' @param standAgeMap template `SpatRaster`
#' @param years the years over which to iterate
#' @param fireBufferedListDT data.table containing non-annual burn and buffer `pixelID`s
#' @param annualCovariates list of data.table objects with `pixelID`
#' @inheritParams castCohortData
#'
#' @return a raster layer with unified stand age and time-since-disturbance values
#'
#' @export
#' @importFrom data.table data.table
#' @importFrom terra rast setValues values
calcYoungAge <- function(years, annualCovariates, standAgeMap, fireBufferedListDT,
                         cutoffForYoungAge = 15) {
  # this is safest way to subset given the NULL year
  yearsIsCorrectNaming <- all(years %in% names(annualCovariates))
  if (yearsIsCorrectNaming %in% FALSE) {
    years <- sapply(years, function(yr) grep(pattern = yr, years, value = TRUE))
  }
  for (year in years) {
    ann <- annualCovariates[[year]] ## no copy made
    fires <- fireBufferedListDT[[year]]
    if (!is.null(fires)) {
      ageVals <- values(standAgeMap, mat = FALSE)
      set(ann, NULL, "youngAge", as.integer(ageVals[ann$pixelID] <= cutoffForYoungAge) |
        is.na(ageVals[ann$pixelID])) ## cannot have NAs
      burnedPix <- fires$pixelID[fires$buffer == 1]
      ageVals[burnedPix] <- 0
      standAgeMap <- setValues(standAgeMap, values = ageVals)
    }
    standAgeMap <- setValues(standAgeMap, values = values(standAgeMap, mat = FALSE) + 1)
  }
  return(annualCovariates)
}
