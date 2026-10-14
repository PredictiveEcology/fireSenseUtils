globalVariables(c(
  "age", "B", "BperClass", "FuelClass", "foo", "Leading", "LeaderValue",
  "maxAge", "NspeciesWithMaxB", "pixelIndex", "totalBiomass"
))

#' Classify `pixelGroups` by flammability
#'
#' @template cohortData
#'
#' @template pixelGroupMap
#'
#' @template sppEquiv
#'
#' @template sppEquivCol
#'
#' @param landcoverDT Optional table of non-forest land cover classes and pixel indices.
#'                    It will override pixel values in `cohortData`, if supplied.
#'
#' @template flammableRTM
#'
#' @param cutoffForYoungAge age at and below which pixels are considered 'young'
#'
#' @param fuelClassCol the column in `sppEquiv` that describes unique fuel classes
#'
#' @param asTable if `TRUE`, return a `data.table` with `pixelID` and one column per layer,
#'   restricted to cells where at least one layer is not `NA` -- exactly
#'   `as.data.table(as.data.frame(<the SpatRaster>, cells = TRUE))` with `cell` renamed to
#'   `pixelID`, without building the raster. Used by [fireSenseCovariatesCreate()].
#'
#' @return a `SpatRaster` of biomass by fuel class as determined by `fuelClassCol` and `cohortData`
#'   (or a `data.table`, see `asTable`).
#'
#' @export
#' @inheritParams fireSenseCovariatesCreate
#' @importFrom data.table copy setkey
#' @importFrom LandR asInteger
#' @importFrom SpaDES.tools rasterizeReduced
#' @importFrom terra as.int values rast
#'
cohortsToFuelClasses <- function(cohortData, pixelGroupMap, flammableRTM, landcoverDT = NULL,
                                 sppEquiv, sppEquivCol, cutoffForYoungAge, fuelClassCol = fireSenseFuelClassCol,
                                 requiredFuelClasses, asTable = FALSE) {
  cc <- .fuelClassVectors(cohortData, pixelGroupMap, flammableRTM, landcoverDT, sppEquiv,
                          sppEquivCol, cutoffForYoungAge, fuelClassCol, requiredFuelClasses)
  if (asTable) {
    keep <- Reduce(`|`, lapply(cc, function(v) !is.na(v)))
    ids <- which(keep)
    out <- setDT(c(list(pixelID = ids), lapply(cc, function(v) v[ids])))
    return(out)
  }
  classList <- rast(pixelGroupMap, nlyrs = length(cc))
  values(classList) <- do.call(cbind, cc)
  names(classList) <- names(cc)
  classList
}

## One numeric vector (length ncell(pixelGroupMap)) per fuel-class layer, in layer order
.fuelClassVectors <- function(cohortData, pixelGroupMap, flammableRTM, landcoverDT,
                              sppEquiv, sppEquivCol, cutoffForYoungAge, fuelClassCol,
                              requiredFuelClasses) {
  joinCol <- c(fuelClassCol, eval(sppEquivCol))
  ## An NA fuel class means "not a fuel" (a species with no cohorts when the model was fitted has
  ## FuelClass = NA in the spread ledger, see fireSense_dataPrepFit R/sppEquivFuelClass.R). Drop
  ## those species like species absent from `sppEquiv`, and an NA from `requiredFuelClasses`;
  ## otherwise a layer named NA reaches the covariate table.
  sppEquivSubset <- unique(sppEquiv[, .SD, .SDcols = joinCol])
  sppEquivSubset <- sppEquivSubset[!is.na(sppEquivSubset[[fuelClassCol]])]
  requiredFuelClasses <- requiredFuelClasses[!is.na(requiredFuelClasses)]

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

  ## Everything below works on integer indices instead of the 10M-row right join, the grouped
  ## max/sum by (FuelClass, pixelGroup), `split()` and a `match()` per class that this function
  ## used to run: same values, ~3x faster. `speciesCode` is mapped to a fuel-class index through
  ## its factor levels, and rows whose species is not in `sppEquiv` are dropped (as the join did).
  spp <- sppEquivSubset[[sppEquivCol]]
  fcOfSpp <- sppEquivSubset[[fuelClassCol]]
  classes <- unique(fcOfSpp)
  youngK <- length(classes) + 1L # index of the youngAge class
  clsOfSpp <- match(fcOfSpp, classes)
  sc <- cohortData$speciesCode
  cls <- if (is.factor(sc)) clsOfSpp[match(levels(sc), spp)][as.integer(sc)] else clsOfSpp[match(sc, spp)]
  pg <- cohortData$pixelGroup
  age <- cohortData$age
  B <- cohortData$B
  keep <- which(!is.na(cls) & !is.na(pg))
  if (length(keep) < length(pg)) {
    pg <- pg[keep]; cls <- cls[keep]; age <- age[keep]; B <- B[keep]
  }

  ## a pixelGroup is young when its oldest cohort is (an NA age is never young)
  maxPG <- if (length(pg)) max(pg) else 0L
  isYoungPG <- logical(maxPG)
  if (length(pg)) {
    maxAge <- data.table(pg = pg, age = age)[, .(maxAge = max(age)), by = pg]
    isYoungPG[maxAge$pg[isYoungAge(maxAge$maxAge, cutoffForYoungAge)]] <- TRUE
  }
  cls[isYoungPG[pg]] <- youngK

  ## biomass per (pixelGroup, class): one grouped sum (GForce) on a single integer key
  perClass <- data.table(k = (pg - 1L) * youngK + cls, B = B)[, .(BperClass = sum(B)), by = k]
  pgOfRow <- (perClass$k - 1L) %/% youngK + 1L
  clsOfRow <- (perClass$k - 1L) %% youngK + 1L
  BperClass <- asInteger(perClass$BperClass)
  # youngAge is better treated as a binary cover variable than continuous measure of biomass
  BperClass[clsOfRow == youngK] <- 1L

  ## Everything below is vector arithmetic over cells. It used to build one SpatRaster per
  ## fuel class with rastFromDF() and stack them, then fireSenseCovariatesCreate() converted
  ## the stack straight back to a table: ~6 s of the ~10 s call on 6.7M cells.
  pgv <- values(pixelGroupMap, mat = FALSE)
  flamValsGood <- !is.na(values(flammableRTM, mat = FALSE))
  pgIdx <- pgv
  pgIdx[!is.na(pgIdx) & (pgIdx < 1 | pgIdx > maxPG)] <- NA
  cc <- list()
  for (k in sort(unique(clsOfRow))) {
    ## a cell takes the class's biomass of its pixelGroup; NA where the class is absent
    w <- clsOfRow == k
    vals <- rep(NA_real_, maxPG)
    vals[pgOfRow[w]] <- as.numeric(BperClass[w])
    rasVals <- vals[pgIdx]
    rasVals[flamValsGood & is.na(rasVals)] <- 0
    cc[[if (k == youngK) youngAgeTxt else classes[k]]] <- rasVals
  }

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

#' Put `cohortData` back into a `SpatRaster` with some extra details
#'
#' @param class `fuelClass` from  `sppEquiv`
#'
#' @template flammableRTM
#'
#' @template pixelGroupMap
#'
#' @template cohortData
#'
#' @return a `SpatRaster` with values equal to `class` biomass (B)
#'
#' @importFrom data.table data.table
#' @importFrom SpaDES.tools rasterizeReduced
#' @importFrom terra values setValues rast
makeRastersFromCD <- function(class, cohortData, flammableRTM, pixelGroupMap) {
  cohortDataFB <- cohortData[FuelClass == class, ]
  ras <- rasterizeReduced(
    reduced = cohortDataFB,
    fullRaster = pixelGroupMap,
    newRasterCols = "BperClass",
    mapcode = "pixelGroup"
  )
  ## fuel class is 0 and not NA if absent entirely
  ## to prevent NAs following aggregation during ignitionFit
  flamVals <- values(flammableRTM, mat = FALSE)
  rasVals <- values(ras, mat = FALSE)
  rasVals[!is.na(flamVals) & is.na(rasVals)] <- 0
  ras <- setValues(x = ras, values = rasVals)
  return(ras)
}

#' Modify `cohortData` with burn column
#'
#' @param year length-two vector giving temporal period used to subset `firePolys`. Closed interval.
#'
#' @param pixelGroupMap either a `SpatRaster` with `pixelGroups` or list of `SpatRasters` named by year.
#'
#' @param cohortData either a `cohortData` object or list of `cohortData` objects named by year.
#'
#' @param firePolys the output of [getFirePolygons()] with `YEAR` column.
#'
#' @return `cohortData` modified with burn status
#'
#' @importFrom LandR addPixels2CohortData
#' @importFrom terra rasterize
buildCohortBurnHistory <- function(cohortData, pixelGroupMap, firePolys, year) {
  ## build fire raster
  firePolys <- do.call(rbind, firePolys)
  firePolys <- firePolys[firePolys$YEAR >= min(year) & firePolys$YEAR <= max(year), ]
  fireRas <- rasterize(firePolys, pixelGroupMap, field = "YEAR", fun = min)
  cdLong <- addPixels2CohortData(cohortData, pixelGroupMap)
  cdLong[, burned := fireRas[pixelIndex]]
  return(cdLong)
}

