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
                    standAgeMap, lcc = NULL, cutoffForYoungAge = fireSenseYoungAgeCutoff,
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
    initialTSD <- if (is.null(polysNeeded)) {
      ## no fire in the young-age window anywhere in the study area (e.g. a low-fire ELF, or an
      ## early data year whose window predates the fire record): nothing burned recently
      setValues(rast(standAgeMap), year - cutoffForYoungAge - 1) |> terra::mask(standAgeMap)
    } else {
      rasterize(polysNeeded,
        y = standAgeMap,
        background = year - cutoffForYoungAge - 1,
        field = "YEAR", fun = "max"
      ) |> terra::mask(standAgeMap)
    }
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

  falseYoungs <- standAgeVals[c(pixelID %in% pixToUpdate & isYoungAge(age, cutoffForYoungAge)) |
                                c(pixelID %in% pixToUpdate & is.na(age))]$pixelID
  ## disturbance history suggests young

  ## Read TSD as a plain numeric vector once. `which()` drops NA entries so a
  ## missing fireRaster value (no recorded disturbance) is not flagged as young
  ## and trueYoungs / trueAges stay length-locked for the := assignment.
  tsdAtPix <- terra::values(initialTSD, mat = FALSE)[pixToUpdate]
  youngPos <- which(isYoungAge(tsdAtPix, cutoffForYoungAge))
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

#' `youngAge` of pixels at one fire year
#'
#' Time since disturbance at `year` is the smaller of the data year's time since disturbance
#' aged to `year` (`tsd + (year - dataYear)`) and the years since the last fire before `year`
#' (`year - fireYear`), taken over **all** fires in `firePixelsByYear`, not only the fires being
#' fitted. A pixel is young when that is at or below `cutoffForYoungAge`. A pixel with `NA`
#' time since disturbance (non-flammable or unknown) is never young, whether or not a fire covers
#' it. A fire in `year` itself does not count: fire years are compared with the state before
#' they burn.
#'
#' It works on pixel IDs, so the same call serves the buffered pixels of a spread fit and every
#' pixel of the landscape (ignition, validation). Only fires in
#' `(year - cutoffForYoungAge):(year - 1)` can make a pixel young, so only those are read.
#'
#' @param tsd `data.table` with columns `pixelID` and `tsd`, the time since disturbance at
#'   `dataYear` (for example from [makeTSD()]).
#' @param dataYear the year `tsd` describes.
#' @param year the fire year.
#' @param firePixelsByYear list named by fire year (`"2001"`, or `"year2001"`) of the `pixelID`s
#'   that burned that year (see [firePixelsByYear()]). Years with no fire may be absent.
#' @param pixelID integer, the pixels to return, in this order. Default: all rows of `tsd`.
#' @inheritParams castCohortData
#'
#' @return integer vector, 1 for young and 0 otherwise, one per `pixelID`.
#'
#' @export
youngAgeAtYear <- function(tsd, dataYear, year, firePixelsByYear = list(),
                           cutoffForYoungAge = fireSenseYoungAgeCutoff,
                           pixelID = tsd$pixelID) {
  age <- tsd$tsd[match(pixelID, tsd$pixelID)] + (year - dataYear)
  names(firePixelsByYear) <- gsub("[^0-9]", "", names(firePixelsByYear))
  for (fy in (year - cutoffForYoungAge):(year - 1)) {
    burned <- firePixelsByYear[[as.character(fy)]]
    if (length(burned)) {
      pos <- which(pixelID %in% burned & !is.na(age))
      age[pos] <- pmin(age[pos], year - fy)
    }
  }
  as.integer(isYoungAge(age, cutoffForYoungAge))
}

#' Pixel IDs burned in each year
#'
#' Turns fire polygons or a fire-year raster into the per-year pixel lists
#' [youngAgeAtYear()] reads. Polygons are rasterized at cell centres, as [makeTSD()] does.
#'
#' @param firePolys list named by year (`"year2001"`) of `SpatVector`s, or `NULL` entries for years
#'   without fires.
#' @param fireRaster `SpatRaster` whose values are fire years; an alternative to `firePolys`. It
#'   holds one fire per pixel, so earlier fires on the same pixel are not seen.
#' @param template `SpatRaster` giving the pixels; required with `firePolys`.
#'
#' @return list named by year (`"2001"`) of integer `pixelID`s; years with no burned pixel are
#'   absent.
#'
#' @export
#' @importFrom terra rasterize values
firePixelsByYear <- function(firePolys = NULL, fireRaster = NULL, template = NULL) {
  if (!is.null(fireRaster)) {
    v <- terra::values(fireRaster, mat = FALSE)
    return(split(which(!is.na(v)), v[!is.na(v)]))
  }
  if (is.null(firePolys) || is.null(template)) stop("Provide firePolys and template, or fireRaster")
  out <- lapply(firePolys, function(p) {
    if (is.null(p) || !length(p)) return(integer(0))
    r <- terra::rasterize(p, template, field = 1, background = NA)
    which(!is.na(terra::values(r, mat = FALSE)))
  })
  names(out) <- gsub("[^0-9]", "", names(firePolys))
  out[lengths(out) > 0]
}

#' Iteratively calculate `youngAge` column in FS covariates
#'
#' Deprecated: use [youngAgeAtYear()]. This one treats `NA` ages as young, only sees the fires
#' in `fireBufferedListDT` (those being fitted), and needs a raster copy per year.
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
                         cutoffForYoungAge = fireSenseYoungAgeCutoff) {
  .Deprecated("youngAgeAtYear")
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
      set(ann, NULL, "youngAge", as.integer(isYoungAge(ageVals[ann$pixelID], cutoffForYoungAge)) |
        is.na(ageVals[ann$pixelID])) ## cannot have NAs
      burnedPix <- fires$pixelID[fires$buffer == 1]
      ageVals[burnedPix] <- 0
      standAgeMap <- setValues(standAgeMap, values = ageVals)
    }
    standAgeMap <- setValues(standAgeMap, values = values(standAgeMap, mat = FALSE) + 1)
  }
  return(annualCovariates)
}

#' Is an age young?
#'
#' The single rule for `youngAge`: an age is young when it is at or below
#' `cutoffForYoungAge` (`age <= cutoffForYoungAge`). `NA` is never young. Every
#' place in fireSenseUtils that decides young or not young calls this function, so
#' fitting and prediction cannot disagree at the cutoff.
#'
#' @param age numeric vector of ages (stand age or time since disturbance), in years.
#' @param cutoffForYoungAge age at and below which a pixel is young.
#'
#' @return logical vector the length of `age`; `FALSE` where `age` is `NA`.
#'
#' @export
isYoungAge <- function(age, cutoffForYoungAge = fireSenseYoungAgeCutoff) {
  age <= cutoffForYoungAge & !is.na(age)
}
