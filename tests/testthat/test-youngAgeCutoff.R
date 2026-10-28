## youngAge is "age <= cutoffForYoungAge" everywhere (makeTSD, cohortsToFuelClasses,
## youngAgeAtYear, fireSenseYoungAgeCutoff); calcNonForestYoungAge() used "<".

test_that("calcNonForestYoungAge: a non-forest pixel whose age equals the cutoff is young", {
  withr::local_package("terra")
  withr::local_package("data.table")

  LCC <- rast(nrows = 1, ncols = 4, vals = 1)
  names(LCC) <- "LCC1"
  dt <- data.table(pixelID = 1:4, NF1 = 1)
  ages <- c(14, 15, 16, 200)

  out <- calcNonForestYoungAge(dt, NFTSD = ages, LCCras = LCC, cutoffForYoungAge = 15)
  expect_equal(as.vector(values(out$youngAge)), c(1, 1, 0, 0))
  expect_equal(as.vector(values(out$LCC1)), c(0, 0, 1, 1))
})

test_that("fireSenseCovariatesCreate(youngAge = TRUE): non-forest pixel at the cutoff is young", {
  withr::local_package("terra")
  withr::local_package("data.table")

  pixelGroupMap <- rast(nrows = 2, ncols = 2, vals = 0L)
  flammableRTM <- rast(pixelGroupMap, vals = 1)
  noCohorts <- data.table(pixelGroup = integer(), speciesCode = character(),
                          age = integer(), B = integer())
  covs <- suppressWarnings(
    fireSenseCovariatesCreate(
      cohortData = noCohorts, pixelGroupMap = pixelGroupMap, flammableRTM = flammableRTM,
      sppEquiv = data.table(LandR = character(), FuelClass = character()),
      landcoverDT = data.table(pixelID = 1:4, nfLCC_40 = c(1, 1, 1, 1)),
      fuelClassCol = "FuelClass", sppEquivCol = "LandR",
      missingLCCgroup = "nfLCC_40", nonForestedLCCGroups = c(nfLCC_40 = 40),
      nonForest_timeSinceDisturbance = c(14, 15, 16, 200),
      cutoffForYoungAge = 15, nonForestCanBeYoungAge = TRUE,
      studyAreaName = "test", useCache = FALSE
    )
  )
  setkey(covs, pixelID)
  expect_equal(covs$youngAge, c(1, 1, 0, 0))
})

test_that("cohortsToFuelClasses works when terra is not attached", {
  skip_if_not_installed("callr")
  pkgPath <- normalizePath(file.path(testthat::test_path(), "..", ".."))
  res <- callr::r(function(pkgPath) {
    ## a source tree has R/*.R files; an installed package (R CMD check) has only R/<pkg>.rdb
    if (length(list.files(file.path(pkgPath, "R"), pattern = "[.]R$"))) {
      pkgload::load_all(pkgPath, quiet = TRUE)
    } else {
      loadNamespace("fireSenseUtils")
    }
    stopifnot(!"package:terra" %in% search())
    pgm <- terra::rast(nrows = 2, ncols = 2, vals = 0L)
    cD <- data.table::data.table(pixelGroup = integer(), speciesCode = character(),
                                 age = integer(), B = integer())
    out <- suppressWarnings(fireSenseUtils::cohortsToFuelClasses(
      cohortData = cD, pixelGroupMap = pgm, flammableRTM = terra::rast(pgm, vals = 1),
      sppEquiv = data.table::data.table(LandR = character(), FuelClass = character()),
      sppEquivCol = "LandR", cutoffForYoungAge = 15, fuelClassCol = "FuelClass",
      requiredFuelClasses = character()
    ))
    names(out)
  }, list(pkgPath), libpath = .libPaths())
  expect_true("youngAge" %in% res)
})

test_that("isYoungAge: age equal to the cutoff is young, NA is not", {
  expect_identical(isYoungAge(c(14, 15, 16, NA), 15), c(TRUE, TRUE, FALSE, FALSE))
  expect_identical(isYoungAge(15), TRUE)  ## default cutoff is fireSenseYoungAgeCutoff
  expect_identical(isYoungAge(fireSenseYoungAgeCutoff + 1), FALSE)
})

test_that("isYoungAge is the only place cutoffForYoungAge is compared", {
  rDir <- file.path(testthat::test_path(), "..", "..", "R")
  skip_if_not(dir.exists(rDir), "package sources not available")
  nm <- "(cutoffForYoungAge|fireSenseYoungAgeCutoff)"
  hits <- unlist(lapply(list.files(rDir, pattern = "\\.R$", full.names = TRUE), function(f) {
    ln <- readLines(f, warn = FALSE)
    ln <- ln[!grepl("^\\s*#", ln)]
    ln <- gsub("<-|->", " ASSIGN ", ln)  ## assignments are not comparisons
    ln[grepl(paste0("[<>]=?\\s*", nm, "|", nm, "\\s*[<>]=?"), ln)]
  }))
  ## the one comparison lives in isYoungAge()
  hits <- hits[!grepl("age <= cutoffForYoungAge & !is.na(age)", hits, fixed = TRUE)]
  expect_length(hits, 0L)
})

## The young non-forest layer is built as a vector. The raster-per-class path it replaced
## (putBackIntoRaster() then the youngAge layer of calcNonForestYoungAge()) gave the same values.
youngNFLandscape <- function(ages) {
  tmpl <- terra::rast(nrows = 4, ncols = 4, vals = 1)
  tmpl[c(1, 6, 16)] <- NA                                    # cells outside the study area
  list(
    flammableRTM = tmpl,
    ## two non-forest classes in the cells that are not NA in the template
    landcoverDT = data.table::data.table(
      pixelID = c(2:5, 7:15),
      nfLCC_40 = c(1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0),
      nfLCC_50 = c(0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1)
    ),
    NFTSD = ages
  )
}
oldYoungAgeValues <- function(x, cutoff) {
  lcc <- putBackIntoRaster(setdiff(names(x$landcoverDT), "pixelID"), data.table::copy(x$landcoverDT),
                           x$flammableRTM)
  terra::values(calcNonForestYoungAge(data.table::copy(x$landcoverDT), x$NFTSD, lcc, cutoff)$youngAge,
                mat = FALSE)
}

test_that("the vector youngAge equals the youngAge layer of the raster path, with and without young pixels", {
  for (ages in list(c(3, 20, 5, 15, 16, 2, 99, 1, 40, 14, 15, 80, 7, 8, 9, 10),   # some young
                    rep(200, 16),                                                   # none young
                    rep(1, 16),                                                     # all young
                    c(NA, 3, NA, 5, 15, NA, 2, 99, 4, 100, 3, 3, 3, 3, 3, NA))) {  # NA ages
    x <- youngNFLandscape(ages)
    pix <- youngNonForestPixelIDs(x$landcoverDT, x$NFTSD, cutoffForYoungAge = 15)
    expect_identical(nonForestYoungAgeValues(x$flammableRTM, pix), oldYoungAgeValues(x, 15))
    expect_identical(names(x$landcoverDT), c("pixelID", "nfLCC_40", "nfLCC_50"))  # not modified
  }
})

test_that("fireSenseCovariatesCreate flags young non-forest without making a raster layer per class", {
  withr::local_package("terra")
  withr::local_package("data.table")
  local_mocked_bindings(putBackIntoRaster = function(...) stop("a raster per class is not needed"))
  create <- function(ages) {
    x <- youngNFLandscape(ages)
    pgm <- rast(x$flammableRTM, vals = 0L)
    cD <- data.table(pixelGroup = integer(), speciesCode = character(), age = integer(), B = integer())
    covs <- suppressWarnings(fireSenseCovariatesCreate(
      cohortData = cD, pixelGroupMap = pgm, flammableRTM = x$flammableRTM,
      sppEquiv = data.table(LandR = character(), FuelClass = character()),
      landcoverDT = x$landcoverDT, fuelClassCol = "FuelClass", sppEquivCol = "LandR",
      missingLCCgroup = "nfLCC_40", nonForestedLCCGroups = c(nfLCC_40 = 40, nfLCC_50 = 50),
      nonForest_timeSinceDisturbance = ages, cutoffForYoungAge = 15, nonForestCanBeYoungAge = TRUE,
      studyAreaName = "test", useCache = FALSE
    ))
    setkey(covs, pixelID)
    list(covs = covs, x = x)
  }
  ages <- c(3, 20, 5, 15, 16, 2, 99, 1, 40, 14, 15, 80, 7, 8, 9, 10)
  r <- create(ages)
  expect_equal(r$covs$youngAge, as.integer(ages[r$covs$pixelID] <= 15))
  expect_gt(sum(r$covs$youngAge), 0)
  expect_equal(sum(create(rep(200, 16))$covs$youngAge), 0)
})

test_that("fireSenseCovariatesCreate follows land cover and burnable pixels that change between calls", {
  withr::local_package("terra")
  withr::local_package("data.table")
  ages <- c(3, 20, 5, 15, 16, 2, 99, 1, 40, 14, 15, 80, 7, 8, 9, 10)
  x <- youngNFLandscape(ages)
  pgm <- rast(x$flammableRTM, vals = 0L)
  cD <- data.table(pixelGroup = integer(), speciesCode = character(), age = integer(), B = integer())
  create <- function(landcoverDT, flammableRTM) {
    covs <- suppressWarnings(fireSenseCovariatesCreate(
      cohortData = cD, pixelGroupMap = pgm, flammableRTM = flammableRTM,
      sppEquiv = data.table(LandR = character(), FuelClass = character()),
      landcoverDT = landcoverDT, fuelClassCol = "FuelClass", sppEquivCol = "LandR",
      missingLCCgroup = "nfLCC_40", nonForestedLCCGroups = c(nfLCC_40 = 40, nfLCC_50 = 50),
      nonForest_timeSinceDisturbance = ages, cutoffForYoungAge = 15, nonForestCanBeYoungAge = TRUE,
      studyAreaName = "test", useCache = FALSE
    ))
    setkey(covs, pixelID)
    covs
  }
  first <- create(x$landcoverDT, x$flammableRTM)
  ## other modules change land cover: pixels 2-5 become burnable but no longer in landcoverDT, and
  ## pixels 10-15 are no longer burnable
  lc2 <- x$landcoverDT[pixelID > 5]
  rtm2 <- x$flammableRTM
  rtm2[10:15] <- NA
  lc2 <- lc2[pixelID < 10]
  second <- create(lc2, rtm2)
  expect_identical(first$pixelID, x$landcoverDT$pixelID)
  expect_identical(second$pixelID, lc2$pixelID)
  expect_equal(second$youngAge, as.integer(ages[lc2$pixelID] <= 15))
  expect_false(identical(first$youngAge, second$youngAge))
})

