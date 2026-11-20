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
#' @return `coreFractionByYear()`: a `data.table` with `pixelID` (the coarse cell), `year` (numeric) and `coreFraction` ([coreFractionTxt]), a row for every
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

#' @param covariates `data.table` with `pixelID` and `year`, as made by [mergePreparedCovs()].
#' @param fractions the result of `coreFractionByYear()` for these arguments, if the caller has it already
#'   (`fireSense_dataPrepPredict` keeps it from year to year, see [coreFractionKey()]); made if not given.
#' @param onMissing what to do with a covariate row whose coarse pixel has no flammable cell, so no `coreFraction`:
#'   `"stop"` (the default; the fit's data must have one for every row) or `"zero"`: the row gets 0 (no core) and one
#'   warning gives the count.
#' @return `addCoreFraction()`: `covariates` with the column `coreFraction` added (by reference); `ignitions` is not changed.
#' @export
#' @rdname coreFraction
addCoreFraction <- function(covariates, years, flammableRTMs, rasterToMatchLargeELF, fact,
                            onMissing = c("stop", "zero"),
                            fractions = coreFractionByYear(years, flammableRTMs, rasterToMatchLargeELF, fact)) {
  onMissing <- match.arg(onMissing)
  key <- function(d) d$year * 1e10 + d$pixelID  ## exact in a double: pixelID < 1e10
  fraction <- fractions[[coreFractionTxt]][match(key(covariates), key(fractions))]
  if (anyNA(fraction)) {
    msg <- paste0("coreFraction is missing for ", sum(is.na(fraction)),
                  " ignition covariate rows: their coarse pixel has no flammable cell")
    if (onMissing == "stop") stop(msg)
    warning(msg, "; they get 0 (no core)", call. = FALSE)
    fraction[is.na(fraction)] <- 0
  }
  data.table::set(covariates, NULL, coreFractionTxt, fraction)
  covariates
}

#' A key of the inputs of `coreFractionByYear()`
#'
#' Two calls with the same key give the same fractions (apart from the year labels), so a caller that
#' makes them every year (`fireSense_dataPrepPredict`) can keep the last result and recompute only when the key
#' changes: when the flammable cells, the core raster, the grid or `fact` change.
#'
#' @return A digest.
#' @export
#' @rdname coreFraction
coreFractionKey <- function(flammableRTMs, rasterToMatchLargeELF, fact) {
  geometry <- function(r) list(dim(r), as.vector(terra::ext(r)), terra::crs(r))
  .robustDigest(list(
    fact, lapply(flammableRTMs, function(r) list(geometry(r), terra::values(r, mat = FALSE) %in% 1)),
    if (!is.null(rasterToMatchLargeELF))
      list(geometry(rasterToMatchLargeELF), terra::values(rasterToMatchLargeELF, mat = FALSE),
           terra::levels(rasterToMatchLargeELF))))
}
