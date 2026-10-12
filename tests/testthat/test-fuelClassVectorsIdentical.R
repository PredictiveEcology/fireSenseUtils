## .fuelClassVectors() was rewritten for speed (integer class index, one grouped max and one grouped
## sum, instead of a join on sppEquiv, split() and a match() per class). It must return exactly what
## the old code returned. `.fuelClassVectorsReference` below is a verbatim copy of the old function
## (e571ebd), kept in this file only as the reference for that comparison.

.fuelClassVectorsReference <- function(cohortData, pixelGroupMap, flammableRTM, landcoverDT,
                              sppEquiv, sppEquivCol, cutoffForYoungAge, fuelClassCol,
                              requiredFuelClasses) {
  joinCol <- c(fuelClassCol, eval(sppEquivCol))
  sppEquivSubset <- unique(sppEquiv[, .SD, .SDcols = joinCol])

  ## `unique()` above is over BOTH columns, but the join below keys on the species column
  ## ALONE. A species carrying two different `fuelClassCol` values therefore survives as
  ## two rows, and every `cohortData` row for it is multiplied -- data.table then stops
  ## with an opaque "Join results in N rows; more than nrow(x)+nrow(i)" that names neither
  ## the species nor this function. In LandR::sppEquivalencies_CA, `Pseu_men` (Douglas-fir)
  ## carries both "DgFrPoPine" and "CedrMplOther", so every study area containing
  ## Douglas-fir failed here while areas without it were fine.
  ##
  ## Naming the species rather than picking one: which fuel class a species belongs to
  ## decides what fuel the model sees, so silently choosing would quietly change the
  ## science. This is a fixable problem in the `sppEquiv` table.
  dupKeys <- sppEquivSubset[, .N, by = c(sppEquivCol)][get("N") > 1L]
  if (NROW(dupKeys)) {
    offenders <- dupKeys[[sppEquivCol]]
    detail <- sppEquivSubset[get(sppEquivCol) %in% offenders]
    setorderv(detail, c(sppEquivCol, fuelClassCol))
    stop("`sppEquiv` maps ", length(offenders), " species to more than one ", fuelClassCol,
         ", so joining it to `cohortData` on `", sppEquivCol, "` would multiply every ",
         "cohort of those species:\n",
         paste0("  ", detail[[sppEquivCol]], " -> ", detail[[fuelClassCol]], collapse = "\n"),
         "\nGive each of those species a single ", fuelClassCol, " in `sppEquiv`.",
         call. = FALSE)
  }

  cD <- cohortData[sppEquivSubset, on = c("speciesCode" = sppEquivCol)]
  setnames(cD, old = fuelClassCol, new = "FuelClass") # so we don't have to use eval, which trips up some dt
  # data.table needs an argument for which column names are kept during join
  cD[, maxAge := max(age), .(pixelGroup)]
  cD[isYoungAge(maxAge, cutoffForYoungAge), FuelClass := youngAgeTxt]
  cD[, maxAge := NULL]
  cD <- cD[, .(BperClass = asInteger(sum(B))), by = c("FuelClass", "pixelGroup")]

  # youngAge is better treated as a binary cover variable than continuous measure of biomass
  cD[FuelClass == youngAgeTxt, BperClass := 1]

  ## Everything below is vector arithmetic over cells. It used to build one SpatRaster per
  ## fuel class with rastFromDF() and stack them, then fireSenseCovariatesCreate() converted
  ## the stack straight back to a table: ~6 s of the ~10 s call on 6.7M cells.
  pgv <- values(pixelGroupMap, mat = FALSE)
  flamValsGood <- !is.na(values(flammableRTM, mat = FALSE))
  cD <- cD[!is.na(cD$pixelGroup)]
  cc <- lapply(split(cD, by = "FuelClass"), function(r) {
    ## a cell takes the class's biomass of its pixelGroup; NA where the class is absent
    rasVals <- as.numeric(r$BperClass[match(pgv, r$pixelGroup)])
    rasVals[flamValsGood & is.na(rasVals)] <- 0
    rasVals
  })

  noFuelForRequiredClass <- character()
  if (!is.null(requiredFuelClasses))
    noFuelForRequiredClass <- setdiff(requiredFuelClasses, names(cc))

  # This is where a species disappears from the map: create a map of zeros
  for (fuel in noFuelForRequiredClass) {
    ## 0 wherever a present class would be: every cell with a pixelGroup, and every flammable cell
    rasVals <- as.numeric(pgv)
    rasVals[!is.na(rasVals)] <- 0
    rasVals[flamValsGood & is.na(rasVals)] <- 0
    cc[[fuel]] <- rasVals
  }
  ## No tree fuel classes at all (no tree species and none required): there is nothing to
  ## stack. Carry an empty list until the youngAge layer below gives the stack its first layer.
  if (length(cc)) cc <- cc[order(names(cc))]

  if (!is.null(landcoverDT) && length(cc)) {
    # find rows that aren't empty i.e. have non-forest land cover
    hasLC <- rowSums(landcoverDT[, .SD, .SDcols = setdiff(names(landcoverDT), nonNFColNamesTxt)],
                     na.rm = TRUE) > 0
    ## must be 0
    if (any(hasLC)) cc <- lapply(cc, function(v) { v[landcoverDT$pixelID[hasLC]] <- 0; v })
  }

  # Need to confirm that there was at least 1 youngAge ... sometime there are none e.g., with 9.2.1 plains
  if (!youngAgeTxt %in% names(cc)) {
    ## every fuel-class layer is NA exactly where pixelGroupMap is, so with no layers the
    ## map itself is the template
    template <- if (length(cc)) cc[[1]] else pgv
    cc[[youngAgeTxt]] <- ifelse(is.na(template), NA_integer_, 0L)
  }
  cc
}

fcvInputs <- function(young = FALSE, fuelClassNA = FALSE) {
  ## 4 x 4 cells; pixelGroup 6 is in the map but has no cohorts; cell 12 is NA (outside the map);
  ## cell 16 has a pixelGroup but is not flammable
  pgm <- terra::rast(nrows = 4, ncols = 4, vals = c(1, 1, 2, 2, 3, 3, NA, 4, 4, 5, 5, NA, 1, 2, 6, 5))
  flam <- terra::rast(pgm, vals = ifelse(is.na(terra::values(pgm, mat = FALSE)), NA, 1))
  flam[16] <- NA
  sppEquiv <- data.table(LandR = c("A", "B", "C", "D", "E"),
                         FuelClass = c("x", "y", "x", "z", "w")) # E: no cohorts, so class w is empty
  if (fuelClassNA) sppEquiv$FuelClass[2] <- NA
  cd <- data.table(pixelGroup = c(1L, 1L, 2L, 3L, 3L, 4L, 5L, 5L, 7L),
                   speciesCode = factor(c("A", "B", "A", "C", "B", "A", "D", "D", "Q"),
                                        levels = c("A", "B", "C", "D", "Q")), # Q is not in sppEquiv
                   age = c(50L, 60L, 10L, 70L, 80L, 5L, 3L, NA, 90L),
                   B = c(100L, 200L, 50L, 5L, 6L, 7L, 8L, 9L, 10L))
  if (young) cd$age <- pmin(cd$age, 12L) # every pixelGroup is young
  list(cohortData = cd, pixelGroupMap = pgm, flammableRTM = flam, sppEquiv = sppEquiv,
       sppEquivCol = "LandR", fuelClassCol = "FuelClass")
}

test_that(".fuelClassVectors returns exactly what the old implementation returned", {
  withr::local_package("terra")
  withr::local_package("data.table")
  lc <- data.table(pixelID = c(1L, 5L, 8L), nfLCC_a = c(1, 0, 1), nfLCC_b = c(0, 0, 0))
  for (young in c(FALSE, TRUE)) for (naClass in c(FALSE, TRUE)) {
    a <- fcvInputs(young, naClass)
    for (landcoverDT in list(NULL, lc)) for (cutoff in c(15, 40)) for (req in list(NULL, c("x", "y", "z", "w", "q"))) {
      args <- c(a, list(landcoverDT = landcoverDT, cutoffForYoungAge = cutoff, requiredFuelClasses = req))
      expect_identical(do.call(.fuelClassVectors, args), do.call(.fuelClassVectorsReference, args),
                       info = paste("allYoung", young, "NAclass", naClass, "landcoverDT", !is.null(landcoverDT),
                                    "cutoff", cutoff, "required", length(req)))
    }
  }
})

test_that(".fuelClassVectors matches the old implementation with character speciesCode and no cohorts", {
  withr::local_package("terra")
  withr::local_package("data.table")
  a <- fcvInputs()
  a$cohortData$speciesCode <- as.character(a$cohortData$speciesCode)
  args <- c(a, list(landcoverDT = NULL, cutoffForYoungAge = 15, requiredFuelClasses = c("x", "y")))
  expect_identical(do.call(.fuelClassVectors, args), do.call(.fuelClassVectorsReference, args))

  args$cohortData <- args$cohortData[integer(0), , drop = FALSE]
  expect_identical(suppressWarnings(do.call(.fuelClassVectors, args)),
                   suppressWarnings(do.call(.fuelClassVectorsReference, args)))
})

test_that("the .robustDigest of the cache key is computed only when caching is on", {
  withr::local_package("terra")
  withr::local_package("data.table")
  nCalls <- 0L
  testthat::local_mocked_bindings(
    .robustDigest = function(...) { nCalls <<- nCalls + 1L; "dig" }, .package = "reproducible")

  clim <- list(MDC = terra::rast(nrows = 4, ncols = 4, vals = 1:16))
  prepare_ignitionClimate(clim, fact = 2, useCache = FALSE)
  expect_identical(nCalls, 0L)

  ## the fireSenseCovariatesCreate() branch that digests landcoverDT and flammableRTM
  pgm <- terra::rast(nrows = 2, ncols = 2, vals = 0L)
  covsArgs <- list(
    cohortData = data.table(pixelGroup = integer(), speciesCode = character(), age = integer(), B = integer()),
    pixelGroupMap = pgm, flammableRTM = terra::rast(pgm, vals = 1),
    sppEquiv = data.table(LandR = character(), FuelClass = character()),
    landcoverDT = data.table(pixelID = 1:4, nfLCC_40 = c(1, 0, 0, 1)),
    fuelClassCol = "FuelClass", sppEquivCol = "LandR", missingLCCgroup = "nfLCC_40",
    nonForestedLCCGroups = c(nfLCC_40 = 40), nonForest_timeSinceDisturbance = c(5, 30, NA, 200),
    cutoffForYoungAge = 15, nonForestCanBeYoungAge = TRUE, studyAreaName = "test")
  suppressWarnings(do.call(fireSenseCovariatesCreate, c(covsArgs, list(useCache = FALSE))))
  expect_identical(nCalls, 0L)
})

test_that("with caching on the cache key still holds the robust digest", {
  withr::local_package("terra")
  clim <- list(MDC = terra::rast(nrows = 4, ncols = 4, vals = 1:16))
  expect_identical(.digestIfCaching(clim, TRUE), reproducible::.robustDigest(clim))
  expect_identical(.digestIfCaching(clim, "overwrite"), reproducible::.robustDigest(clim))
  expect_null(.digestIfCaching(clim, FALSE))

  ## .cacheExtra of the Cache() call is upstream digest + the climate digest, as before
  keys <- list()
  testthat::local_mocked_bindings(
    Cache = function(..., .cacheExtra = NULL, useCache = TRUE) { keys[[length(keys) + 1L]] <<- .cacheExtra; list(...)[[1]] },
    .package = "fireSenseUtils")
  prepare_ignitionClimate(clim, fact = 2, digest = "up", useCache = TRUE)
  expect_identical(keys[[1]], c("up", reproducible::.robustDigest(clim)))
})
