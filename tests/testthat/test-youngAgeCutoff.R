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
