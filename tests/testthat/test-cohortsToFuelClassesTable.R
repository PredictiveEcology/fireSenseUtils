## cohortsToFuelClasses() used to build one SpatRaster per fuel class (rastFromDF) and
## fireSenseCovariatesCreate() converted the stack back to a table. It now works on cell
## vectors, and `asTable = TRUE` returns the table directly. This must change nothing: the
## pre-change implementation is kept below as the reference.

skip_if_not_installed("terra")

## --- reference: cohortsToFuelClasses() as it was before the table rewrite ---
oldCohortsToFuelClasses <- function(cohortData, pixelGroupMap, flammableRTM, landcoverDT = NULL,
                                 sppEquiv, sppEquivCol, cutoffForYoungAge, fuelClassCol = fireSenseFuelClassCol,
                                 requiredFuelClasses) {
  joinCol <- c(fuelClassCol, eval(sppEquivCol))
  sppEquivSubset <- unique(sppEquiv[, .SD, .SDcols = joinCol])

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
  cD[, maxAge := max(age), .(pixelGroup)]
  cD[isYoungAge(maxAge, cutoffForYoungAge), FuelClass := youngAgeTxt]
  cD[, maxAge := NULL]
  cD <- cD[, .(BperClass = asInteger(sum(B))), by = c("FuelClass", "pixelGroup")]

  cD[FuelClass == youngAgeTxt, BperClass := 1]

  classes <- sort(unique(cD$FuelClass))
  
  
  pgmVals <- list(pixelGroup = values(pixelGroupMap, mat = FALSE), 
                  pixelId = seq(ncell(pixelGroupMap))) |> 
    setDT() |> na.omit()
  aa <- pgmVals[cD, on = "pixelGroup", allow.cartesian=TRUE]
  bb <- split(aa, by = "FuelClass")
  flamVals <- values(flammableRTM, mat = FALSE)
  flamValsGood <- !is.na(flamVals)
  cc <- Map(r = bb, function(r) {
    ras <- rastFromDF(r[, .(pixelId, BperClass)], rasTemplate = pixelGroupMap)
    rasVals <- values(ras, mat = FALSE)
    rasVals[flamValsGood & is.na(rasVals)] <- 0
    ras <- setValues(x = ras, values = rasVals)
    ras
    })
  
  noFuelForRequiredClass <- character()
  if (!is.null(requiredFuelClasses))
    noFuelForRequiredClass <- setdiff(requiredFuelClasses, names(cc))

  if (length(noFuelForRequiredClass)) {
    for (fuel in noFuelForRequiredClass) {
      ras <- Copy(pixelGroupMap)
      rasVals <- values(ras, mat = FALSE)
      rasVals[rasVals > 0] <- 0
      cc[[fuel]] <- setValues(x = ras, values = rasVals)
    }
    
  }
  if (length(cc)) {
    dd <- rast(cc)
    classList <- dd[[order(names(dd))]]
  } else {
    classList <- NULL
  }
  
  
  
  if (!is.null(landcoverDT) && !is.null(classList)) {
    landcoverDT[, foo := rowSums(.SD, na.rm = TRUE), .SDcols = setdiff(names(landcoverDT), nonNFColNamesTxt)]
    if (nrow(landcoverDT[foo > 0, ]) > 0) {
      classList[landcoverDT[foo > 0]$pixelID] <- 0 # must be 0
    }
    landcoverDT[, foo := NULL]
  }

  
  if (!youngAgeTxt %in% names(classList)) {
    template <- if (is.null(classList)) pixelGroupMap else classList[[1]]
    ya <- as.int(is.na(template))
    vals <- values(ya, mat = FALSE)
    ya[vals == 1L] <- NA
    names(ya) <- youngAgeTxt
    classList <- if (is.null(classList)) ya else c(classList, ya)
  }
  return(classList)
}

environment(oldCohortsToFuelClasses) <- asNamespace("fireSenseUtils")

fuelFixture <- function(n = 30, seed = 1, landcover = TRUE, young = TRUE, required = TRUE, noTrees = FALSE) {
  set.seed(seed)
  pg <- terra::rast(nrows = n, ncols = n, vals = sample(c(0:30, NA), n * n, TRUE, prob = c(rep(1, 31), 3)))
  fl <- terra::rast(pg, vals = sample(c(1, NA), n * n, TRUE, prob = c(9, 1)))
  sp <- c("A", "B", "C", "D")
  sE <- data.table::data.table(LandR = sp, FuelClass = c("x", "y", "y", "z"))
  cd <- data.table::CJ(pixelGroup = 1:30, speciesCode = sp)
  cd <- cd[sample(nrow(cd), 70)]
  cd$age <- sample(c(0:10, 40:150), nrow(cd), TRUE)
  cd$B <- sample(c(1:900, NA), nrow(cd), TRUE)
  if (noTrees) { cd <- cd[integer(0), , drop = FALSE]; sE <- sE[integer(0), , drop = FALSE] }
  ldt <- NULL
  if (landcover) {
    ids <- sample(n * n, 60)
    ldt <- data.table::data.table(pixelID = ids, nfLCC_1 = rbinom(60, 1, 0.5), nfLCC_2 = rbinom(60, 1, 0.5))
  }
  list(cohortData = cd, pixelGroupMap = pg, flammableRTM = fl, landcoverDT = ldt, sppEquiv = sE,
       sppEquivCol = "LandR", cutoffForYoungAge = if (young) 15 else -1, fuelClassCol = "FuelClass",
       requiredFuelClasses = if (required && !noTrees) c("x", "y", "z", "w") else NULL)
}

test_that("cohortsToFuelClasses gives the same raster and table as the pre-change implementation", {
  cases <- list(
    list(), list(landcover = FALSE), list(required = FALSE), list(young = FALSE),
    list(seed = 2), list(seed = 3, landcover = FALSE, required = FALSE),
    list(noTrees = TRUE), list(noTrees = TRUE, landcover = FALSE)
  )
  for (cs in cases) {
    args <- do.call(fuelFixture, cs)
    lcBefore <- data.table::copy(args$landcoverDT)

    ref <- suppressWarnings(do.call(oldCohortsToFuelClasses, args))
    refDT <- data.table::as.data.table(as.data.frame(ref, cells = TRUE))
    data.table::setnames(refDT, "cell", "pixelID")

    newRas <- suppressWarnings(do.call(cohortsToFuelClasses, args))
    newDT <- suppressWarnings(do.call(cohortsToFuelClasses, c(args, list(asTable = TRUE))))

    expect_identical(names(newRas), names(ref))
    expect_identical(terra::values(newRas), terra::values(ref))
    expect_true(terra::compareGeom(newRas, ref, stopOnError = FALSE))
    expect_identical(newDT, refDT)
    expect_identical(args$landcoverDT, lcBefore) # not left with a scratch column
  }
})

test_that("fireSenseCovariatesCreate with a fuelClassTable gives the same covariates, and leaves the table alone", {
  args <- fuelFixture(seed = 4, landcover = FALSE)
  ## every pixel flammable and in a pixelGroup, so no pixel is NA in the covariates
  set.seed(4)
  terra::values(args$pixelGroupMap) <- sample(1:30, terra::ncell(args$pixelGroupMap), TRUE)
  terra::values(args$flammableRTM) <- 1
  nPix <- terra::ncell(args$pixelGroupMap)
  nf <- sample(nPix, 50)
  ldt <- data.table::data.table(pixelID = 1:nPix, nfLCC_1 = 0L, nfLCC_2 = 0L)
  data.table::set(ldt, nf[1:25], "nfLCC_1", 1L)
  data.table::set(ldt, nf[26:50], "nfLCC_2", 1L)
  tsd <- terra::rast(args$pixelGroupMap, vals = sample(0:60, nPix, TRUE))
  mk <- function(mode, ...) fireSenseCovariatesCreate(
    cohortData = args$cohortData, pixelGroupMap = args$pixelGroupMap, flammableRTM = args$flammableRTM,
    sppEquiv = args$sppEquiv, landcoverDT = data.table::copy(ldt), fuelClassCol = "FuelClass",
    requiredFuelClasses = args$requiredFuelClasses, sppEquivCol = "LandR", missingLCCgroup = "nfLCC_1",
    nonForestedLCCGroups = list(nfLCC_1 = 1, nfLCC_2 = 2), nonForest_timeSinceDisturbance = tsd,
    cutoffForYoungAge = 15, nonForestCanBeYoungAge = TRUE, studyAreaName = "t", useCache = FALSE,
    fuelCovariates = mode, ...)
  tbl <- cohortsToFuelClasses(cohortData = args$cohortData, pixelGroupMap = args$pixelGroupMap,
                              flammableRTM = args$flammableRTM, landcoverDT = data.table::copy(ldt),
                              sppEquiv = args$sppEquiv, sppEquivCol = "LandR", fuelClassCol = "FuelClass",
                              requiredFuelClasses = args$requiredFuelClasses, cutoffForYoungAge = 15,
                              asTable = TRUE)
  tblBefore <- data.table::copy(tbl)
  for (mode in c("species", "domSecWetland")) {
    expect_identical(suppressWarnings(mk(mode, fuelClassTable = tbl)), suppressWarnings(mk(mode)))
    expect_identical(tbl, tblBefore)
  }
})
