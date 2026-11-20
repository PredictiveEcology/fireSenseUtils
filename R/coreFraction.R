#' The share of each coarse pixel's flammable cells that are core, by fire year
#'
#' `fireSense_dataPrepFit` adds it to the ignition covariates of the fit and `fireSense_dataPrepPredict` to the covariates
#' it predicts, with these functions, so both define the core alike.
#'
#' The ignition table has a row for every coarse pixel (`igAggFactor`), buffer included, but the ignitions it
#' counts start in the core of `rasterToMatchLargeELF`. `coreFraction` lets the ignition fit weight a pixel by
#' the part of it that can hold a counted ignition. It is the mean, over the pixel's flammable cells (the
#' cells of value 1 in the data year's binary `flammableRTM`; 0 and `NA` are not flammable), of
#' `inELFCore()`, aggregated with the same `terra::aggregate(fact =, fun = mean)` as the coarse
#' covariates. A multi-ELF categorical raster has no buffer (all `1`); without a raster everything is core.
#'
#' @param years named list, by data year, of the fire years (`year<year>`) that use that data year's vegetation.
#' @param flammableRTMs list of `SpatRaster`, one per data year, in the order of `years`.
#' @param rasterToMatchLargeELF `SpatRaster` read by `inELFCore()`, or `NULL`.
#' @param fact the aggregation factor, `igAggFactor`.
#' @return a `data.table` with `pixelID` (the coarse cell), `year` (numeric) and `coreFraction` ([coreFractionTxt]), a row for every
#'   coarse pixel with flammable cells, for every fire year.
#' @export
#' @rdname coreFraction
coreFractionByYear <- function(years, flammableRTMs, rasterToMatchLargeELF, fact) {
  byDataYear <- Map(function(fireYears, flammableRTM) {
    flammable <- terra::values(flammableRTM, mat = FALSE) %in% 1
    core <- inELFCore(seq_len(terra::ncell(flammableRTM)), rasterToMatchLargeELF,
                                      grid = flammableRTM)
    fineCore <- terra::setValues(terra::rast(flammableRTM, nlyrs = 1), ifelse(flammable, as.numeric(core), NA_real_))
    coarse <- terra::values(terra::aggregate(fineCore, fact = fact, fun = mean, na.rm = TRUE), mat = FALSE)
    pixelID <- which(!is.na(coarse))
    out <- data.table::data.table(pixelID = rep(pixelID, length(fireYears)),
                                  year = rep(as.numeric(gsub(yearTxt, "", fireYears)),
                                             each = length(pixelID)),
                                  fraction = rep(coarse[pixelID], length(fireYears)))
    data.table::setnames(out, "fraction", coreFractionTxt)
  }, years, flammableRTMs)
  data.table::rbindlist(byDataYear)
}

#' Add `coreFraction` to the ignition covariates
#'
#' @param covariates `data.table` with `pixelID` and `year`, as made by [mergePreparedCovs()].
#' @return `covariates` with the column `coreFraction` added (by reference); `ignitions` is not changed.
#' @export
#' @rdname coreFraction
addCoreFraction <- function(covariates, years, flammableRTMs, rasterToMatchLargeELF, fact) {
  fractions <- coreFractionByYear(years, flammableRTMs, rasterToMatchLargeELF, fact)
  key <- function(d) d$year * 1e10 + d$pixelID  ## exact in a double: pixelID < 1e10
  fraction <- fractions[[coreFractionTxt]][match(key(covariates), key(fractions))]
  if (anyNA(fraction))
    stop("coreFraction is missing for ", sum(is.na(fraction)),
         " ignition covariate rows: their coarse pixel has no flammable cell")
  data.table::set(covariates, NULL, coreFractionTxt, fraction)
  covariates
}
