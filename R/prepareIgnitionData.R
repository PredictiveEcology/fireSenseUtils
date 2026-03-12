
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
  digURLs <- Map(url = lightningUrls, function(url) 
    reproducible:::getRemoteMetadata(url = reproducible:::googledriveIDtoHumanURL(url), 
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
  fcn <- fireSenseUtils:::fireSenseCovariatesCreate(...)
  fuelCovs <- SpaDES.tools::rastFromDF(fcn, rasTemplate = rasTemplate)
  terra::aggregate(fuelCovs, fact = fact, fun = mean)
}





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
                                                 c(names(ignitionClimateCoarse), "cell", "ignitions", fireSenseUtils::yearChar))]
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
  #   # ranEffs <- "fireSenseUtils::yearChar"
  #   set(fireSense_ignitionCovariates, NULL, ranEffsLabel, as.character(fireSense_ignitionCovariates$year))
  # }
  firstCols <- c("pixelID", "ignitions", names(ignitionClimateCoarse), youngAgeName)
  firstCols <- firstCols[firstCols %in% names(fireSense_ignitionCovariates)]
  setcolorder(fireSense_ignitionCovariates, neworder = firstCols)
  
  fireSense_ignitionCovariates

}




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
                           c("pixelID", fireSenseUtils:::ignitionsTxt, fireSenseUtils:::escapesTxt, "year", "yearChar"))
      rescalers <- sapply(covariates[, .SD, .SDcols = toRescale], max)
      needRescale <- sapply(rescalers, FUN = function(x) !inRange(x, 0, 10))
      cols <- names(rescalers)[which(needRescale)]
      message("rescaling the following variables: ", paste(cols, collapse = ", "))
      ignitionRescalers <- 10^(floor(log10(abs(rescalers[cols])))) # if range is 0,1, need + 1
      covariates <- rescaleVarsByMagnitude(covariates, ignitionRescalers)
    } else {
      
      SDcols <- setdiff(colnames(covariates), c("yearChar", fireSenseUtils:::ignitionsTxt, "escapes"))
      scaledData <- scale(covariates[, ..SDcols])
      origIgnitions <- covariates[[fireSenseUtils:::ignitionsTxt]]
      origEscapes <- covariates[["escapes"]]
      centeringData <- attributes(scaledData)
      covariates <- as.data.table(scaledData)
      set(covariates, NULL, fireSenseUtils:::ignitionsTxt, origIgnitions)
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
      #                          .SDcols = setdiff(colnames(covariates), c("yearChar", fireSenseUtils:::ignitionsTxt))
      #                          #.SDcols = c("CMDsm", "Betu_pap", "Pc_gl.Lr_la", "Pice_mar",
      #                          #"Pn_co.Pn_ba", "Pp_ba.Pp_tr", "lightning")
      # ]
      # cols <- grep("V1", value = TRUE, colnames(covariates))
      # setnames(covariates, old = cols, new = gsub(".V1", "", cols))
      # if (fireSenseUtils:::escapesTxt %in% colnames(covariates)) {
      #   set(covariates, NULL, fireSenseUtils:::escapesTxt, covariates[[fireSenseUtils:::escapesTxt]])
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

igOrEscNames <- function(igOrEsc, pre = "fireSense_", post, case = c("lower", "camel", "sentence", "title")) {
  if (startsWith(tolower(case[1]), prefix = "cam"))
    igOrEsc <- camelCase(igOrEsc)
  if (startsWith(tolower(case[1]), prefix = "sen") || startsWith(tolower(case[1]), prefix = "tit"))
    igOrEsc <- tools::toTitleCase(igOrEsc) # only has one word, so OK
  paste0(pre, igOrEsc, post)
}


ignitionsTxt <- "ignitions"

escapesTxt <- "escapes"

