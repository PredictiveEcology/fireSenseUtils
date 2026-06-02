
#' Prepare and cache lightning rasters from Google Drive sources
#'
#' @description
#' Downloads, reads, aggregates, and caches multiple lightning-related rasters
#' (e.g., daily flashes, flash density, positive CG counts) using
#' `reproducible::prepInputs()`, `terra::aggregate()`, and `reproducible::Cache()`.
#' The function builds a cache key that combines local processing context
#' (`rasterToMatch`, `igAggFactor`, and the lightning reader) with the remote
#' hashes of the Google Drive files, so cached results are invalidated whenever
#' either local parameters or upstream files change.
#'
#' @param rtm `terra::SpatRaster` or raster-like object used as the target
#'   template (i.e., an alignment/extent/projection reference) for reading
#'   lightning rasters. 
#' @param igAggFactor `numeric(1)` aggregation factor forwarded to
#'   `terra::aggregate(fact = ...)`, typically used to coarsen lightning rasters
#'   for ignition modeling.
#' @param dPath Character. Destination path forwarded to
#'   `reproducible::prepInputs(destinationPath = ...)`; the directory where
#'   downloaded lightning rasters are stored / cached on disk.
#'
#' @details
#' The function maintains an internal mapping of Google Drive file IDs for
#' four lightning products:
#'
#' - `lightningDays`
#' - `lightningDensity`
#' - `positiveCG`
#' - `positiveCGdensity`
#'
#' For each product, it:
#'
#' 1. Resolves the Drive ID to a human URL, obtains remote metadata (hash),
#' 2. Calls `reproducible::prepInputs()` with a reading function that should
#'    return a raster aligned to `sim$rasterToMatch`,
#' 3. Aggregates the raster with `terra::aggregate(fact = P(sim)$igAggFactor)`,
#' 4. Caches the result with `reproducible::Cache()` under a product-specific
#'    function label, seeding the cache with both local and remote digests.
#'
#' A final `Cache()` call wraps the list of per-product results to further
#' stabilize caching across the set.
#'
#' @return
#' A named list of aggregated lightning rasters (class `list` of
#' `terra::SpatRaster`), with names matching:
#' `c("lightningDays", "lightningDensity", "positiveCG", "positiveCGdensity")`.
#'
#' @section Caching and reproducibility:
#' The cache key includes:
#' - The remote hash of each Google Drive resource (so cloud-side updates bust
#'   the cache),
#' - The local context digest of `rasterToMatch`, `igAggFactor`, and the reader,
#'   via `.robustDigest`.
#'
#' This design minimizes recomputation while remaining sensitive to upstream and
#' local changes.
#'
#' @note
#' If `prepInputs()` expects a function for `fun`, prefer:
#' ```
#' fun = function(...) rld(targetFile, to = sim$rasterToMatch)
#' ```
#' instead of calling `rld(...)` immediately.
#'
#' @examples
#' \dontrun{
#' out <- prepare_LightningData(
#'   rtm = sim$rasterToMatch,
#'   igAggFactor = P(sim)$igAggFactor,
#'   rld = readLightningData
#' )
#' out$lightningDensity
#' }
#'
#' @seealso
#' [reproducible::prepInputs()], [reproducible::Cache()], [terra::aggregate()]
#'
#' @references
#' - **User-provided function source** (this documentation derives from the code
#'   you supplied).
#' - Reproducible workflows with the `reproducible` package and raster handling
#'   with `terra` as described in their package manuals and vignettes.
#'
#' @keywords data-preparation caching reproducibility raster lightning
#' @importFrom terra aggregate
#' @export
prepare_LightningData <- function(rtm, igAggFactor, dPath) {
  lightningUrls <- list(lightningDays = "1jeKJquhVJsesoNk2EPP1QZkttX3Zwp5c",
                        lightningDensity = "12fnhfKtER-JXkl06M4_yZ3GvpZWtQlIr" ,
                        positiveCG = "1bn6cQ23tvPicFLHn1tz4Z3AqDzJI4r60",
                        positiveCGdensity = "1GNixhXj1Ex1jT0tWXfhmxef-dX3ze1a4")
  digCE <- .robustDigest(list(rtm = rtm, 
                              igAggFactor = igAggFactor,
                              rld = fireSenseUtils::readLightningData))
  .getRemoteMetadata <- utils::getFromNamespace("getRemoteMetadata", "reproducible")
  .googledriveIDtoHumanURL <- utils::getFromNamespace("googledriveIDtoHumanURL", "reproducible")
  digURLs <- Map(url = lightningUrls, function(url)
    .getRemoteMetadata(url = .googledriveIDtoHumanURL(url),
                       isGDurl = TRUE)$remoteHash)
  Map(url = lightningUrls, nam = names(lightningUrls), digURL = digURLs,
                           function(url, nam, digURL) {
                             {
                               prepInputs(url = url,
                                          fun = readLightningData(targetFile,
                                                                  to = rtm),
                                          destinationPath = dPath, useCache = FALSE)  |>
                                 terra::aggregate(fact = igAggFactor)} |>
                               Cache(.functionName = paste0("prepInputs_lightning_", nam),
                                     omitArgs = c("...", "x"), # x comes from terra::aggregate and is undefined at call; so returns different each time
                                     .cacheExtra = append(list(url = digURL), digCE))
                           }) |> Cache(.functionName = "prepInputs_lightning",
                                       .cacheExtra = append(digURLs, digCE))

}


#' Aggregate and cache ignition-related climate rasters
#'
#' @description
#' Aggregates a list of ignition-relevant climate rasters to a coarser spatial
#' resolution using `terra::aggregate()` (mean across source cells) and caches
#' the resulting list with `reproducible::Cache()`. The cache key is augmented
#' with a robust digest of the input rasters and a caller-provided upstream
#' digest.
#'
#' @param ignitionClimateList A named or unnamed `list` of raster-like objects
#'   (typically `terra::SpatRaster`) representing historical climate variables
#'   used in ignition modeling.
#' @param fact `numeric(1)` aggregation factor passed to
#'   `terra::aggregate(fact = ...)`.
#' @param useCache `logical` (default `TRUE`). If `TRUE`, wrap the per-element
#'   aggregation in `reproducible::Cache()` so that re-running with the same
#'   inputs returns the cached result. If `FALSE`, force recomputation.
#' @param digest A digest object (e.g., `list`, atomic vector) representing
#'   upstream processing context, which will be appended to the input-derived
#'   digest when constructing the cache key.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Computes a robust digest of `ignitionClimateList` using
#'         `.robustDigest()` to capture changes in input rasters.
#'   \item Appends that digest to the user-supplied `digest` to form a composite
#'         cache key (`dig3`).
#'   \item Aggregates each raster via `terra::aggregate(fact = fact, fun = mean)`.
#'   \item Wraps the aggregated list in `reproducible::Cache()` with
#'         `omitArgs = "x"` (to exclude the ephemeral list element passed to
#'         `lapply`) and `.cacheExtra = dig3` to drive cache invalidation from
#'         the composite digest.
#' }
#'
#' All rasters are processed independently via `lapply()`, preserving list order
#' and names.
#'
#' @return
#' A `list` of aggregated climate rasters (typically `terra::SpatRaster`), in
#' one-to-one correspondence with `ignitionClimateList`.
#'
#' @section Caching and reproducibility:
#' Cache invalidation is controlled by `.cacheExtra = dig3`, which combines:
#' \itemize{
#'   \item a robust digest of the input rasters (`.robustDigest(ignitionClimateList)`), and
#'   \item the caller-provided `digest` (e.g., pipeline/module context).
#' }
#'
#' @section Assumptions and required environment:
#' The symbols `.robustDigest()` and `reproducible::Cache()` must be available.
#' No validation is performed on input types. Aggregation is controlled by the `fact` parameter.
#'
#' @examples
#' \dontrun{
#' # Example (requires `terra` rasters and `reproducible` cache setup):
#' upstreamDigest <- .robustDigest(list(module = "Fire", version = "1.0.0"))
#' ignitionClimateCoarse <- prepare_ignitionClimate(
#'   ignitionClimateList = sim$historicalClimateRasters[
#'     sim$climateVariablesForFire$ignition
#'   ],
#'   fact = 4,
#'   digest = upstreamDigest
#' )
#' }
#'
#' @seealso
#' [terra::aggregate()], [reproducible::Cache()]
#'
#' @references
#' - `terra` package reference manual: `aggregate()` (spatial raster aggregation).
#' - `reproducible` package reference manual and vignettes: `Cache()` and `.cacheExtra`.
#'
#' @keywords climate ignition raster aggregation caching reproducibility
#' @importFrom terra aggregate
#' @export
prepare_ignitionClimate <- function(ignitionClimateList, fact, digest = NULL, useCache = TRUE) {
  # ignitionClimateList <- sim$historicalClimateRasters[sim$climateVariablesForFire$ignition]
  digClimate <- .robustDigest(ignitionClimateList)
  dig3 <- append(digest, digClimate)
  ignitionClimateCoarse <- lapply(X = ignitionClimateList, FUN = terra::aggregate,
                                  fact = fact, fun = mean) |>
    Cache(.functionName = "aggregate_historicalClimateRasters_to_coarse",
          omitArgs = "x", .cacheExtra = dig3, useCache = useCache)
}



#' Create and aggregate fuel covariate rasters to a coarse resolution
#'
#' @description
#' Creates fuel covariate layers from tabular inputs using
#' `fireSenseUtils:::fireSenseCovariatesCreate()`, converts the resulting data
#' frame into a raster using `SpaDES.tools::rastFromDF()`, and aggregates the
#' raster to a coarser spatial resolution with `terra::aggregate()` using the
#' mean.
#'
#' @param ... Arguments forwarded directly to
#'   `fireSenseUtils:::fireSenseCovariatesCreate()`, defining the inputs and
#'   options required to construct fuel covariates.
#' @param rasTemplate A raster-like object (typically `terra::SpatRaster`) used
#'   as a spatial template (extent, resolution, projection) when converting the
#'   fuel covariate data frame to a raster.
#' @param fact `numeric(1)` aggregation factor passed to
#'   `terra::aggregate(fact = ...)`, controlling the spatial coarsening scale.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Generates a fuel covariate data frame using
#'         `fireSenseUtils:::fireSenseCovariatesCreate(...)`.
#'   \item Converts the data frame to a raster aligned with `rasTemplate` via
#'         `SpaDES.tools::rastFromDF()`.
#'   \item Aggregates the raster to a coarser resolution using
#'         `terra::aggregate(fact = fact, fun = mean)`.
#' }
#'
#' The aggregation is applied across all raster layers (if present), computing
#' the mean of contributing cells.
#'
#' @return
#' A `terra::SpatRaster` containing aggregated fuel covariate layers at the
#' coarser spatial resolution defined by `fact`.
#'
#' @section Assumptions and required environment:
#' The function assumes that the following are available:
#' \itemize{
#'   \item `fireSenseUtils:::fireSenseCovariatesCreate()` (non-exported API),
#'   \item `SpaDES.tools::rastFromDF()`,
#'   \item the `terra` package.
#' }
#'
#' No input validation is performed on the structure of the data returned by
#' `fireSenseCovariatesCreate()`.
#'
#' @note
#' This function uses a non-exported function from `fireSenseUtils` via `:::`,
#' which may change without notice across package versions.
#'
#' @examples
#' \dontrun{
#' fuelCovsCoarse <- prepare_FuelCovsCoarse(
#'   fuels = fuelTable,
#'   rasTemplate = templateRaster,
#'   fact = 4
#' )
#' }
#'
#' @seealso
#' [terra::aggregate()], [SpaDES.tools::rastFromDF()]
#'
#' @references
#' - `terra` package reference manual: raster aggregation with `aggregate()`.
#' - `SpaDES.tools` package reference manual: `rastFromDF()`.
#' - `fireSenseUtils` package source code: `fireSenseCovariatesCreate()` (internal API).
#'
#' @keywords fuel covariates raster aggregation
#' @importFrom terra aggregate
#' @export
prepare_FuelCovsCoarse <- function(..., rasTemplate, fact) {
  fcn <- fireSenseCovariatesCreate(...)
  fuelCovs <- SpaDES.tools::rastFromDF(fcn, rasTemplate = rasTemplate)
  terra::aggregate(fuelCovs, fact = fact, fun = mean)
}




#' Build ignition covariates table by stacking/extracting rasters and joining lightning
#'
#' For each `year` (or year period) provided, stacks and extracts covariates from
#' coarse fuel rasters (split into forested vs. non-forested LCC groups) and climate,
#' at ignition point locations, then merges results across years, filters unusable
#' pixels, appends lightning summaries, and orders columns. Results are optionally
#' cached.
#'
#' @param years A vector or list of year identifiers passed to
#'   [fireSenseUtils::stackAndExtract()]. Can be a list matching the structure of
#'   `fuelCovsCoarse`/`ignitionClimateCoarse`.
#' @param fuelCovsCoarse A list (possibly nested by year) of named raster-layer
#'   collections (e.g., `terra::SpatRaster` or lists of SpatRaster) representing
#'   fuel-related covariates. Names must overlap with LCC group names used to
#'   split forested vs. non-forested classes.
#' @param ignitionFirePoints Ignition point locations (e.g., `sf`/`sp` points or
#'   a \code{data.frame} convertible by \code{stackAndExtract}) used for extraction.
#' @param nonForestedLCCGroups A named vector or list whose names identify the LCC
#'   (land-cover) classes considered non-forested. These names are used to split
#'   `fuelCovsCoarse` into `LCC` (non-forested) and `fuel` (forested/other) inputs.
#' @param ignitionClimateCoarse A list (possibly by year) of coarse climate rasters
#'   (e.g., `terra::SpatRaster`) to include in the stack-and-extract process.
#'   The names of this list are used as `userTags` in caching and are later treated
#'   as climate column names when computing row-wise cover sums.
#' @param lightningMap A named `terra::SpatRaster` (or a named list coercible to one)
#'   containing one or more lightning-derived layers (e.g., `"lightningDays"`).
#'   These values are looked up by `pixelID` and appended as columns.
#' @param digest A cache key supplement (character or list) forwarded to
#'   [reproducible::Cache()] via `.cacheExtra` to ensure cache correctness for
#'   large inputs omitted from the hash.
#' @param useCache `logical` (default `TRUE`). If `TRUE`, wrap the Map call in
#'   [reproducible::Cache()] with `omitArgs = c("fuel", "LCC")`, `.cacheExtra = digest`,
#'   and `userTags = names(ignitionClimateCoarse)`.
#'
#' @details
#' **Processing steps:**
#' \enumerate{
#'   \item Calls [fireSenseUtils::stackAndExtract()] via `Map()` for each entry in `years`,
#'         splitting `fuelCovsCoarse` per year into:
#'         \itemize{
#'           \item `fuel`: all layers not in `names(nonForestedLCCGroups)`
#'           \item `LCC`: layers whose names are in `names(nonForestedLCCGroups)`
#'         }
#'         Additional `MoreArgs` include `climate = ignitionClimateCoarse` and
#'         `fires = ignitionFirePoints`.
#'   \item Optionally caches the list of per-year extractions with
#'         [reproducible::Cache()] (function name tagged as `"stackAndExtract"`).
#'   \item Combines the list to a single `data.table` via [data.table::rbindlist()].
#'   \item Removes rows for which the sum across all cover-type columns is zero
#'         (computed as `rowSums(.SD)` over all columns except `ignitions`,
#'         `yearTxt` (from `fireSenseUtils`), `cell`, and the climate variable
#'         names).
#'   \item Renames the `cell` column to `pixelID`, coerces `year` to numeric,
#'         and appends lightning columns by indexing `terra::values(terra::rast(lightningMap))`
#'         with `pixelID`.
#'   \item Optionally reorders columns to place `pixelID`, `ignitions`,
#'         climate columns (i.e., `names(ignitionClimateCoarse)`), and `youngAgeTxt`
#'         (if it exists) first.
#' }
#'
#' **Assumptions / Requirements:**
#' - `fireSenseUtils::stackAndExtract()` must accept arguments `years`, `fuel`, `LCC`,
#'   `climate`, and `fires`, and return a data.table-like object with at least
#'   `cell`, `ignitions`, and `year` (or `yearTxt`) columns.
#' - `fireSenseUtils::yearTxt` is used as a column name to exclude from cover sums;
#'   ensure it exists/aligns with the data produced by `stackAndExtract()`.
#' - `youngAgeTxt` is referenced when establishing column order but is **not**
#'   defined in this scope; if unavailable, it is silently omitted.
#' - `pixelID` is assumed to be a valid row index into `lightningMap` raster values.
#'
#' **Lightning layer choice:** The code appends all layers provided in `lightningMap`.
#' If a single best-performing layer is known (e.g., `"lightningDays"`), subset
#' `lightningMap` before calling this function.
#'
#' **Caching:** Large raster inputs are excluded from the cache key via
#' `omitArgs = c("fuel", "LCC")`; the caller-supplied `digest` should capture
#' the salient identity of those omitted objects to avoid cache collisions.
#'
#' @return A `data.table` of ignition covariates with:
#' \itemize{
#'   \item `pixelID` (formerly `cell`)
#'   \item `ignitions`
#'   \item Climate covariate columns (names derived from `ignitionClimateCoarse`)
#'   \item Cover/class covariates retained from `stackAndExtract()` (filtered to rows with non-zero cover sum)
#'   \item Lightning-derived columns appended from `lightningMap`
#'   \item `year` coerced to numeric
#' }
#'
#' @section Potential pitfalls:
#' - If `pixelID` indexes are not aligned with `lightningMap` cell indices,
#'   appended lightning values will be incorrect; ensure consistent indexing/resolution.
#' - If `fireSenseUtils::yearTxt` is absent in the returned data, the exclusion in
#'   `.SDcols` is harmless, but confirm that the intended year column exists for analysis.
#' - If `youngAgeTxt` is not defined, it will be ignored when reordering columns.
#'
#' @examples
#' \dontrun{
#' # Pseudocode, assuming you have compatible inputs prepared:
#' yrs <- list(`2001` = 2001L, `2002` = 2002L)
#' # fuelCovsCoarse[[year]] is a named list or SpatRaster with LCC class names
#' # ignitionClimateCoarse[[year]] is a SpatRaster of climate layers
#' # ignitionFirePoints are point locations (sf/sp/data.frame accepted by stackAndExtract)
#' # lightningMap is a named SpatRaster (e.g., "lightningDays")
#'
#' out <- mergePreparedCovs(
#'   years = yrs,
#'   fuelCovsCoarse = fuelCovsCoarse,
#'   ignitionFirePoints = ignitionFirePoints,
#'   nonForestedLCCGroups = c("Shrub", "Grass", "Bare"),
#'   ignitionClimateCoarse = ignitionClimateCoarse,
#'   lightningMap = lightningMap["lightningDays"],
#'   digest = list(seed = 1L, cfg = "v1"),
#'   useCache = TRUE
#' )
#' data.table::str(out)
#' }
#'
#' @seealso fireSenseUtils::stackAndExtract, reproducible::Cache, data.table::rbindlist,
#'   terra::rast, terra::values
#'
#' @importFrom data.table rbindlist set setnames setcolorder
#' @importFrom reproducible Cache
#' @importFrom terra rast values
#' @export
mergePreparedCovs <- function(years, fuelCovsCoarse, ignitionFirePoints, nonForestedLCCGroups,
                              ignitionClimateCoarse, lightningMap, digest, useCache = TRUE) {
  
  fireSense_ignitionCovariates <- Map(
    f = fireSenseUtils::stackAndExtract,
    years = years, # list(pre2012, post2012),
    fuel = Map(lrc = fuelCovsCoarse, function(lrc) lrc[[setdiff(names(lrc), names(nonForestedLCCGroups))]]),
    LCC = Map(lrc = fuelCovsCoarse, function(lrc) lrc[[names(nonForestedLCCGroups)]]), # list(LCCras$year2010, LCCras$year2020),
    MoreArgs = list(climate = ignitionClimateCoarse,
                    fires = ignitionFirePoints)
  ) |>
    Cache(omitArgs = c("fuel", "LCC"), 
          .cacheExtra = digest, useCache = useCache,
          .functionName = "stackAndExtract",
          userTags = names(ignitionClimateCoarse)
    )
  
  
  fireSense_ignitionCovariates <- rbindlist(fireSense_ignitionCovariates)
  
  ## remove any pixels that are 0 for all classes
  fireSense_ignitionCovariates[, coverSums := rowSums(.SD),
                               .SDcols = setdiff(names(fireSense_ignitionCovariates),
                                                 c(names(ignitionClimateCoarse), "cell", "ignitions", fireSenseUtils::yearTxt))]
  fireSense_ignitionCovariates <- fireSense_ignitionCovariates[coverSums > 0]
  set(fireSense_ignitionCovariates, NULL, "coverSums", NULL)
  
  ## rename cells to pixelID - though aggregated raster is not saved
  setnames(fireSense_ignitionCovariates, old = "cell", new = "pixelID")
  fireSense_ignitionCovariates[, year := as.numeric(year)]
  
  # add lightning
  # https://www.tandfonline.com/doi/full/10.1080/07055900.2020.1845117
  
  # Eliot tested using all 4 and there was one clear winner basedon on variable importance
  #   lightningDays -- in ELF 4.3
  # lightningMap <- sim$lightningMaps["lightningDays"]
  
  set(fireSense_ignitionCovariates, NULL, names(lightningMap),
      as.data.frame(terra::values(terra::rast(lightningMap))[fireSense_ignitionCovariates[["pixelID"]],]))
  
  # ## for random effect
  # if (grepl("xgb", Par$modelAlgorithm) %in% FALSE) {
  #   # ranEffs <- "fireSenseUtils::yearTxt"
  #   set(fireSense_ignitionCovariates, NULL, ranEffsLabel, as.character(fireSense_ignitionCovariates$year))
  # }
  firstCols <- c("pixelID", "ignitions", names(ignitionClimateCoarse), youngAgeTxt)
  firstCols <- firstCols[firstCols %in% names(fireSense_ignitionCovariates)]
  setcolorder(fireSense_ignitionCovariates, neworder = firstCols)
  
  fireSense_ignitionCovariates

}




#' Prepare and (optionally) cache scaled covariates for modeling
#'
#' Wraps [rescaleCovariates()] to rescale features according to the modeling
#' algorithm, and optionally caches results based on a robust digest of the
#' unscaled input data. Returns a list containing the (possibly) rescaled
#' covariates, the modeling formula (if any), per-variable rescalers (if any),
#' the x-axis variable name(s), and the digest used for caching.
#'
#' @param unscaledData A `data.frame` or `data.table` of covariates. Must
#'   include predictor columns and may include response columns such as
#'   `ignitions` and/or `escapes` (see Details).
#' @param algorithm `character`. Name of the model algorithm. If it contains
#'   `"xgb"`, covariates are standardized via [base::scale()] (centering and
#'   scaling); otherwise, covariates are rescaled by order of magnitude to
#'   roughly the \[0, 10\] range (see Details).
#' @param rescaleVars `logical`. If `TRUE`, perform rescaling/standardization
#'   according to `algorithm`; if `FALSE`, covariates are returned unchanged.
#' @param useCache `logical` (default `TRUE`). If `TRUE`, results are cached
#'   using [reproducible::Cache()] with an extra key derived from a digest of
#'   `unscaledData`.
#'
#' @details
#' This function:
#' * Computes a digest of `unscaledData` (if `useCache = TRUE`) using an internal
#'   `.robustDigest()` to form a `.cacheExtra` key for [reproducible::Cache()].
#' * Calls [rescaleCovariates()] with `formula = NULL` and `family = NULL`
#'   (non-xgboost code-paths were removed).
#' * Caches the result while omitting the large arguments `"covariates"` and
#'   `"formula"` from cache key construction (via `omitArgs`), relying on the
#'   explicit digest instead.
#'
#' **Expected columns**:
#' - `ignitions`: if present, is preserved through standardization (xgboost path).
#' - `escapes`: if present, is preserved through standardization (xgboost path).
#' - `year` or `yr`: used to infer `xvar` (see `rescaleCovariates()`).
#'
#' @return A `list` with elements:
#' \describe{
#'   \item{covariates}{A `data.table` with rescaled/standardized covariates.}
#'   \item{formula}{`NULL` (the non-xgboost formula path is disabled).}
#'   \item{ignitionRescalers}{Named numeric vector of per-variable 10-based
#'     rescalers (non-xgboost path) or `NULL` for xgboost / no-rescale.}
#'   \item{xvar}{Character vector indicating which year-like column(s) were found.}
#'   \item{digestOfData}{The digest used for caching, or `NULL` if `useCache = FALSE`.}
#' }
#'
#' @seealso [rescaleCovariates()], [reproducible::Cache()]
#'
#' @examples
#' \dontrun{
#' library(data.table)
#' dt <- data.table(pixelID = 1:3,
#'                  year = c(2001L, 2002L, 2003L),
#'                  ignitions = c(0, 1, 0),
#'                  x1 = c(0.01, 2.3, 15))
#'
#' out <- prepareCovariatesOuter(unscaledData = dt,
#'                               algorithm = "xgb_classifier",
#'                               rescaleVars = TRUE,
#'                               useCache = FALSE)
#' str(out)
#' }
#'
#' @export
prepareCovariatesOuter <- function(unscaledData, algorithm, rescaleVars,
                                   useCache = TRUE) {
  # covariatesHere <- igOrEscNames(igOrEsc, post = "Covariates") # paste0("fireSense_", igOrEsc, "Covariates")
  # objsNeeded <- c(covariatesHere)
  # objsNeeded <- c()
  
  # Eliot removed the ability to not use xgboost -- not maintained and less good
  # if (grepl("xgb", algorithm) %in% FALSE) {
  #   formulaHere <- igOrEscNames(igOrEsc, post = "Formula") # paste0("fireSense_", igOrEsc, "Formula")
  #   objsNeeded <- c(objsNeeded, formulaHere)
  #   modelFamily <- igOrEscNames(igOrEsc, pre = "", post = "Family") # paste0(igOrEsc, "Family")
  #   family <- P(sim)[[modelFamily]]
  #   formulaHere <- sim[[formulaHere]] # pull from simList
  # } else {
  formulaHere <- NULL
  family <- NULL
  # }
  
  # This takes time for large datasets
  digestOfData <- if (isTRUE(useCache)) 
    .robustDigest(list(unscaledData = unscaledData)) else NULL
  
  data <- rescaleCovariates(formula = formulaHere,
                            covariates = unscaledData,
                            rescaleVars = rescaleVars,
                            modelAlgorithm = algorithm) |>
    Cache(omitArgs = c("covariates", "formula"), .cacheExtra = digestOfData, useCache = useCache)
  
  data$digestOfData <- digestOfData
  # if (identical(igOrEsc, "escape")) {
  #   # Escape should not have lightning
  #   # if don't explicitly copy, then Cache above returns the "lightning"-removed data.table
  #   data$covariates <- data.table::copy(data$covariates)
  #   set(data$covariates, NULL, "lightning", NULL)
  # }
  data
}



#' Rescale/standardize covariates for modeling
#'
#' Applies either order-of-magnitude rescaling to approximately confine features
#' to the \[0, 10\] range (non-xgboost) or standardization via [base::scale()]
#' (xgboost path). Preserves key response columns when present.
#'
#' @param formula A model formula or `NULL`. Non-xgboost paths would parse and
#'   validate it, but the surrounding code path currently sets `formula = NULL`
#'   and disables non-xgboost modeling in the caller.
#' @param covariates A `data.frame` or `data.table` of covariates and possibly
#'   responses (e.g., `ignitions`, `escapes`).
#' @param rescaleVars `logical`. If `TRUE`, features are rescaled/standardized
#'   according to `modelAlgorithm`. If `FALSE`, returned unchanged (with
#'   `ignitionRescalers = NULL`).
#' @param modelAlgorithm `character`. If the string contains `"xgb"`, apply
#'   [base::scale()] (centering and scaling) to all columns except
#'   `yearChar`, `ignitions`, and `escapes`. Otherwise, log-10 order-of-magnitude
#'   rescaling is applied to all non-excluded covariates whose maxima fall
#'   outside \[0, 10\].
#'
#' @details
#' **Year variable inference**: If any of `c("year", "yr")` match column names (case-insensitive),
#' `xvar` is set to the intersection. Otherwise, the code references `rows` (a TODO),
#' which is undefined here. This appears to be a bug or placeholder and will error
#' unless `rows` exists in scope; consider providing a default (e.g., `NULL`) or removing.
#'
#' **Non-xgboost path**:
#'   - Excludes `pixelID`, `ignitions`, `escapes`, `year`, and `yearChar` from rescaling.
#'   - Computes per-variable `rescalers <- 10^(floor(log10(abs(max))))` for variables
#'     whose maxima fall outside \[0, 10\].
#'   - Calls `rescaleVarsByMagnitude(covariates, ignitionRescalers)` (must exist in scope).
#'
#' **xgboost path**:
#'   - Standardizes numeric columns except `yearChar`, `ignitions`, and `escapes`.
#'   - Preserves original `ignitions` and (if present) `escapes` columns.
#'   - Attaches centering/scaling attributes under `attr(covariates, "scaleData")`.
#'
#' **External references**:
#'   - Uses `ignitionsTxt` and `escapesTxt` to
#'     refer to the expected response column names. If these are not available,
#'     ensure compatible strings are provided (see exported constants below).
#'
#' @return A `list` with elements:
#' \describe{
#'   \item{covariates}{A `data.table` after rescaling/standardization.}
#'   \item{formula}{The input `formula` (often `NULL` in current usage).}
#'   \item{ignitionRescalers}{Named numeric vector of per-variable rescalers (non-xgb) or `NULL`.}
#'   \item{xvar}{Character vector with detected year-like column(s) or (bug) `rows`.}
#' }
#'
#' @examples
#' \dontrun{
#' library(data.table)
#' dt <- data.table(ignitions = c(0,1,0),
#'                  escapes = c(0,0,1),
#'                  year = 2001:2003,
#'                  x1 = c(0.05, 12, 100),
#'                  x2 = c(5, 6, 7))
#'
#' # XGBoost-like standardization
#' out_xgb <- rescaleCovariates(NULL, dt, rescaleVars = TRUE, modelAlgorithm = "xgb_tree")
#' str(out_xgb)
#'
#' # Non-xgboost magnitude rescaling
#' out_other <- rescaleCovariates(NULL, dt, rescaleVars = TRUE, modelAlgorithm = "glm")
#' str(out_other)
#' }
#'
#' @seealso [base::scale()]
#'
#' @export
rescaleCovariates <- function(formula, covariates, rescaleVars, modelAlgorithm) {
  
  covariates <- copy(setDT(covariates))
  if (any(c("year", "yr") %in% tolower(names(covariates)))) {
    xvar <- intersect(c("year", "yr"), tolower(names(covariates)))
  } else {
    xvar <- rows #TODO what is this?
  }
  
  
  if (grepl("xgb", modelAlgorithm) %in% FALSE) {
    formula <- as.formula(formula, env = .GlobalEnv)
    terms <- terms.formula(formula)
    
    if (attr(terms, "response")) {
      y <- formula[[2L]]
    } else {
      stop("Incomplete formula, the LHS is missing.")
    }
  }
  
  # if (!is.data.table(covariates))
  #   covariates <- as.data.table(covariates)
  
  if (rescaleVars) {
    if (grepl("xgb", modelAlgorithm) %in% FALSE) {
      # rescalers <- abs(sapply(covariates[, .SD, .SDcols = toRescale], FUN = max))
      message("Variables outside of [0,10] range will be rescaled to [0,10]")
      toRescale <- setdiff(names(covariates),
                           c("pixelID", ignitionsTxt, escapesTxt, "year", "yearChar"))
      rescalers <- sapply(covariates[, .SD, .SDcols = toRescale], max)
      needRescale <- sapply(rescalers, FUN = function(x) !inRange(x, 0, 10))
      cols <- names(rescalers)[which(needRescale)]
      message("rescaling the following variables: ", paste(cols, collapse = ", "))
      ignitionRescalers <- 10^(floor(log10(abs(rescalers[cols])))) # if range is 0,1, need + 1
      covariates <- rescaleVarsByMagnitude(covariates, ignitionRescalers)
    } else {
      
      SDcols <- setdiff(colnames(covariates), c("yearChar", ignitionsTxt, "escapes"))
      scaledData <- scale(covariates[, ..SDcols])
      origIgnitions <- covariates[[ignitionsTxt]]
      origEscapes <- covariates[["escapes"]]
      centeringData <- attributes(scaledData)
      covariates <- as.data.table(scaledData)
      set(covariates, NULL, ignitionsTxt, origIgnitions)
      if (!is.null(origEscapes))
        set(covariates, NULL, "escapes", origEscapes)
      
      # covariates <- covariates[,
      #                          append(
      #                            list(# yearChar = yearChar,
      #                              ignitions = ignitions
      #                              #, nfLCC_100 = nfLCC_100,
      #                              #, nfLCC_50_80 = nfLCC_50_80,
      #                              #, youngAge = youngAge
      #                            ),
      #                            lapply(.SD, scale, scale = TRUE)),
      #                          .SDcols = setdiff(colnames(covariates), c("yearChar", ignitionsTxt))
      #                          #.SDcols = c("CMDsm", "Betu_pap", "Pc_gl.Lr_la", "Pice_mar",
      #                          #"Pn_co.Pn_ba", "Pp_ba.Pp_tr", "lightning")
      # ]
      # cols <- grep("V1", value = TRUE, colnames(covariates))
      # setnames(covariates, old = cols, new = gsub(".V1", "", cols))
      # if (escapesTxt %in% colnames(covariates)) {
      #   set(covariates, NULL, escapesTxt, covariates[[escapesTxt]])
      # }
      
      setattr(covariates, name = "scaleData", value = centeringData)
      ignitionRescalers <- NULL #so that fire fireSense_IgnitionFit can add it
      
    }
  } else {
    ignitionRescalers <- NULL #so that fire fireSense_IgnitionFit can add it
  }
  
  return(list(covariates = covariates,
              formula = formula,
              ignitionRescalers = ignitionRescalers,
              xvar = xvar))
}

#' Construct standardized FireSense object names
#'
#' Creates a standardized name by combining a `pre` prefix, a normalized
#' `igOrEsc` token (case-transformed), and a `post` suffix.
#'
#' @param igOrEsc `character`. A token such as `"ignition"` or `"escape"`.
#' @param pre `character` scalar. Prefix string (default `"fireSense_"`).
#' @param post `character` scalar. Suffix string to append.
#' @param case One of `c("lower", "camel", "sentence", "title")`. Controls how
#'   `igOrEsc` is transformed before concatenation. Partial matching on the
#'   first letters is used: e.g., `"cam"`, `"sen"`, `"tit"`.
#'
#' @return `character(1)` constructed as `paste0(pre, transformed(igOrEsc), post)`.
#'
#' @examples
#' igOrEscNames("ignition", post = "Covariates")
#' igOrEscNames("escape", pre = "fs_", post = "_Formula", case = "camel")
#'
#' @export
igOrEscNames <- function(igOrEsc, pre = "fireSense_", post, case = c("lower", "camel", "sentence", "title")) {
  if (startsWith(tolower(case[1]), prefix = "cam"))
    igOrEsc <- .camelCase(igOrEsc)
  if (startsWith(tolower(case[1]), prefix = "sen") || startsWith(tolower(case[1]), prefix = "tit"))
    igOrEsc <- tools::toTitleCase(igOrEsc) # only has one word, so OK
  paste0(pre, igOrEsc, post)
}

# Internal: capitalize the first letter of each element. Used by
# igOrEscNames(case = "camel") to produce names like "Ignition" from
# "ignition" without depending on a function from another package.
.camelCase <- function(x) sub("^([a-z])", "\\U\\1", x, perl = TRUE)

