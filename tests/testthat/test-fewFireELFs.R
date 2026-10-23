## fewFireELFs() on the toy ELF maps of helper-toyELFs.R. The two fireregimetools loaders read
## national shapefiles, so they are replaced by toy fire records; everything else is the real
## code (fewFireELFs() and the counting, status, plan and merge). prepInputs() is replaced by a function
## that evaluates the loader call it is given, and Cache() by a pass-through.
##
## Fire records, 2001-2003, thresholds minNaturalIgnitions = 3 and minFirePolygons = 1 (minEscapes 0
## unless a test sets it):
##
##   ELF     natural ignitions counted      escapes   polygons counted     status
##   3.1.1   1  (col 2, 2001)               1         1 (cols 1-2, 2001)   few (1 < 3)
##   3.1.2   3  (col 5: 2001, 2002, 2003)   2         1 (col 5, 2002)      ok (the 2003 fire is 100 ha,
##                                                                          below the 2500 ha pixel)
##   5.1     1  (col 8, 2002)               1         0                    zero (no polygon)
##   1.1     not counted: ecozone 1
##
## Not counted: a human-caused point (col 2), a point from 1999 (col 5), a point in 1.1 (col 11),
## and a polygon of half a cell in 5.1 (1250 ha, not above one 2500 ha pixel).
## So 3.1.1 merges with 3.1.2, its only neighbour with the same base "3.1" (1 + 3 = 4 >= 3
## ignitions, 1 + 1 = 2 >= 1 polygons), and 5.1 has no neighbour sharing base "5": skipped.

## The fire records are the toy ones (`points` replaces the NFDB points), for the calling test.
local_toyFireRecords <- function(loaderArgs = NULL, cache = FALSE, points = toyFirePoints(), env = parent.frame()) {
  local_mocked_bindings(
    fireRecordShapefile = function(url, destinationPath, ...) file.path(destinationPath, sub("\\.zip$", ".shp", basename(url))),
    .env = env
  )
  local_mocked_bindings(
    prepInputs = function(url, targetFile, destinationPath, fun, ...) eval(fun, list(targetFile = targetFile)),
    .package = "reproducible", .env = env
  )
  if (!cache) local_mocked_bindings(Cache = function(x, ...) x, .package = "reproducible", .env = env)
  local_mocked_bindings(
    load_nfdb_points = function(nfdb_shp, study_area, fire_years = NULL, min_size_ha = 1) {
      if (!is.null(loaderArgs)) loaderArgs$points <- list(nfdb_shp = nfdb_shp, study_area = study_area,
                                                          fire_years = fire_years, min_size_ha = min_size_ha)
      points
    },
    load_nbac_polys = function(nbac_shp, study_area, fire_years = NULL, min_size_ha = 1) {
      if (!is.null(loaderArgs)) loaderArgs$polys <- list(nbac_shp = nbac_shp, study_area = study_area,
                                                         fire_years = fire_years, min_size_ha = min_size_ha)
      toyFirePolys()
    },
    .package = "fireregimetools", .env = env
  )
}

runFewFire <- function(..., minEscapes = 0, loaderArgs = NULL, cache = FALSE) {
  local_toyFireRecords(loaderArgs = loaderArgs, cache = cache)
  fewFireELFs(toyELFs(), fireYears = 2001:2003, pixelAreaHa = 2500, destinationPath = tempdir(),
              nfdbUrl = "https://example.org/points.zip", nbacUrl = "https://example.org/polys.zip",
              minEscapes = minEscapes, ...)
}

test_that("fire counts and status are those of the toy fire records", {
  out <- suppressMessages(runFewFire(minNaturalIgnitions = 3, minFirePolygons = 1))
  st <- as.data.frame(out$status)
  st <- st[order(st$ELF), ]
  expect_identical(st$ELF, c("3.1.1", "3.1.2", "5.1"))        # 1.1 is arctic: not counted
  expect_identical(st$naturalIgnitions, c(1L, 3L, 1L))
  expect_identical(st$escapes, c(1L, 2L, 1L))
  expect_identical(st$firePolygons, c(1L, 1L, 0L))
  expect_identical(st$status, c("few", "ok", "zero"))
})

test_that("the thin ELF is merged with the neighbour sharing its base; the other is skipped", {
  out <- suppressMessages(runFewFire(minNaturalIgnitions = 3, minFirePolygons = 1))
  plan <- out$plan
  expect_identical(plan$action, c("merge", "skip"))
  expect_identical(plan$ELF, c("3.1.1_2", NA))
  expect_identical(plan$members, list(c("3.1.1", "3.1.2"), "5.1"))
  expect_identical(plan$naturalIgnitions, c(4L, 1L))           # 1 + 3; 1
  expect_identical(plan$firePolygons, c(2L, 0L))               # 1 + 1; 0
  expect_identical(out$excluded, "5.1")
})

test_that("the returned maps hold the merged ELF in place of its members", {
  out <- suppressMessages(runFewFire(minNaturalIgnitions = 3, minFirePolygons = 1))
  expect_setequal(names(out$ELFs$rasWhole), c("5.1", "1.1", "3.1.1_2"))
  expect_setequal(names(out$ELFs$rasCentered), c("5.1", "1.1", "3.1.1_2"))
  merged <- out$ELFs$rasWhole[["3.1.1_2"]]
  expect_identical(names(merged), "3.1.1_2")
  ## row 1 of the merged map: core over cols 1-6, the 3.1.2 buffer at col 7, 0 beyond
  expect_identical(as.vector(merged[1, ])[[1]], c(2, 2, 2, 2, 2, 2, 1, 0, 0, 0, 0, 0))
  ## the untouched ELFs are returned as they came in
  expect_identical(terra::values(out$ELFs$rasWhole[["5.1"]]), terra::values(toyELFs()$rasWhole[["5.1"]]))
})

test_that("each decision is announced once", {
  msgs <- capture_messages(runFewFire(minNaturalIgnitions = 3, minFirePolygons = 1))
  msgs <- grep("^fireSense_ELFs: ", msgs, value = TRUE)
  expect_identical(
    trimws(msgs),
    c("fireSense_ELFs: merge 3.1.1 + 3.1.2 -> 3.1.1_2 (too few fires; merged with 3.1.2, the longest shared border)",
      "fireSense_ELFs: skip 5.1 (too few fires; no neighbour shares its base)"))
})

test_that("the thresholds given are the ones used, for the status and for the plan", {
  ## with 5 ignitions needed, 3.1.2 (3) is thin too, and 3.1.1 + 3.1.2 = 4 is still too few
  out <- suppressMessages(runFewFire(minNaturalIgnitions = 5, minFirePolygons = 1))
  st <- as.data.frame(out$status)
  expect_identical(st$status[order(st$ELF)], c("few", "few", "zero"))
  expect_identical(out$plan$action, c("skip", "skip"))
  expect_identical(sort(out$excluded), c("3.1.1", "3.1.2", "5.1"))
  expect_setequal(names(out$ELFs$rasWhole), c("3.1.1", "3.1.2", "5.1", "1.1")) # nothing merged

  ## with 1 ignition and 0 polygons needed only 5.1 is left out, because a "zero" ELF is thin
  ## whatever the thresholds
  lax <- suppressMessages(runFewFire(minNaturalIgnitions = 1, minFirePolygons = 0))
  expect_identical(lax$excluded, "5.1")
  expect_identical(lax$plan$action, "skip")
})

test_that("an ELF with enough ignitions and polygons but too few escapes is few, and merged to reach them", {
  ## 3.1.2 has 3 ignitions (enough) but 2 escapes
  out <- suppressMessages(runFewFire(minNaturalIgnitions = 3, minFirePolygons = 1, minEscapes = 3))
  st <- as.data.frame(out$status)
  expect_identical(st$status[order(st$ELF)], c("few", "few", "zero"))
  expect_identical(out$plan$action, c("merge", "skip"))
  expect_identical(out$plan$escapes, c(3L, 1L))
  ## 3.1.1 + 3.1.2 = 3 escapes is still short of 4: neither is fitted
  out <- suppressMessages(runFewFire(minNaturalIgnitions = 3, minFirePolygons = 1, minEscapes = 4))
  expect_identical(out$plan$action, c("skip", "skip"))
  expect_setequal(out$excluded, c("3.1.1", "3.1.2", "5.1"))
})

test_that("escapeSizeHa decides which ignitions escape", {
  ## the toy points are 5000 ha (one is 100 ha); none reaches 6000 ha
  big <- suppressMessages(runFewFire(minNaturalIgnitions = 3, minFirePolygons = 1,
                                     escapeSizeHa = 6000))
  expect_identical(sum(big$status$escapes), 0L)
})

test_that("the defaults are 50 ignitions, 50 polygons, 5 escapes of 50 ha", {
  out <- suppressMessages(runFewFire(minEscapes = 5))
  ## nothing in the toy records reaches 50, so every counted ELF is thin and none can be merged
  expect_identical(sort(out$excluded), c("3.1.1", "3.1.2", "5.1"))
  expect_identical(formals(fewFireELFs)$minNaturalIgnitions, 50)
  expect_identical(formals(fewFireELFs)$minFirePolygons, 50)
  expect_identical(formals(fewFireELFs)$minEscapes, 5)
  expect_identical(formals(fewFireELFs)$escapeSizeHa, 50)
})

test_that("the clipped fire polygons are cached: the second call reads no polygons", {
  withr::local_options(reproducible.cachePath = tempfile("cache"), reproducible.useCache = TRUE)
  loads <- new.env()
  loads$n <- 0L
  local_mocked_bindings(
    fireRecordShapefile = function(url, destinationPath, ...) file.path(destinationPath, sub("\\.zip$", ".shp", basename(url)))
  )
  local_mocked_bindings(
    prepInputs = function(url, targetFile, destinationPath, fun, ...) eval(fun, list(targetFile = targetFile)),
    .package = "reproducible"
  )
  local_mocked_bindings(
    load_nfdb_points = function(...) toyFirePoints(),
    load_nbac_polys = function(...) { loads$n <- loads$n + 1L; toyFirePolys() },
    .package = "fireregimetools"
  )
  run <- function(years = 2001:2003, minNat = 3) suppressMessages(fewFireELFs(
    toyELFs(), fireYears = years, pixelAreaHa = 2500, destinationPath = tempdir(),
    nfdbUrl = "https://example.org/points.zip", nbacUrl = "https://example.org/polys.zip",
    minNaturalIgnitions = minNat, minFirePolygons = 1, minEscapes = 0))
  first <- run()
  expect_identical(loads$n, 1L)
  second <- run(minNat = 4)                   # a change of plan rule, not of the data
  expect_identical(loads$n, 1L)
  run(years = 2001:2002)                                   # other years: the table is rebuilt
  expect_identical(loads$n, 2L)
})

test_that("the fire records are asked for over the whole map, every size, the years given", {
  seen <- new.env()
  suppressMessages(runFewFire(minNaturalIgnitions = 3, minFirePolygons = 1, loaderArgs = seen))
  expect_identical(seen$points$nfdb_shp, "points.shp")
  expect_identical(seen$polys$nbac_shp, "polys.shp")
  expect_identical(seen$points$fire_years, 2001:2003)
  expect_identical(seen$polys$fire_years, 2001:2003)
  expect_identical(seen$points$min_size_ha, 0)   # fires of every size, as the fit reads them
  expect_identical(seen$polys$min_size_ha, 0)
  expect_equal(unname(as.vector(terra::ext(seen$points$study_area))), c(0, 60000, 0, 20000))
  expect_equal(unname(as.vector(terra::ext(seen$polys$study_area))), c(0, 60000, 0, 20000))
})

## The size rule (fireSenseUtils::ELFsizePlan()) runs after the fire-count plan, on the regions that remain.
## The toy cores are 300 km2 each, so a huge minimum makes them all small. Land cover is one class everywhere.
toyLandCoverFile <- function(code = 210) {
  f <- tempfile(fileext = ".tif")
  terra::writeRaster(terra::setValues(toyGrid(), code), f)
  f
}
## thresholds that leave 3.1.1 and 3.1.2 as they are in the fire-count plan (5.1 has no polygon: zero, skipped)
sizeArgs <- list(minNaturalIgnitions = 1, minFirePolygons = 0, minEscapes = 0, maxBurnRatio = Inf)

test_that("small regions are merged with a similar neighbour after the fire-count plan", {
  out <- suppressMessages(do.call(runFewFire, c(sizeArgs, list(landCoverFile = toyLandCoverFile(),
                                                             minRegionAreaKm2 = 1e6))))
  expect_identical(out$plan$action, c("skip", "merge"))
  merge <- out$plan[out$plan$action == "merge", ]
  expect_identical(merge$ELF, "3.1.1_2")
  expect_identical(merge$members[[1]], c("3.1.1", "3.1.2"))
  expect_match(merge$reason, "^smaller than 1000000 km2; merged with 3\\.1\\.[12] \\(land cover distance 0\\.00")
  expect_identical(merge$naturalIgnitions, 4L)           # counts of the original members, from the status
  expect_identical(out$excluded, "5.1")                  # the skipped ELF is not a partner and not merged
  expect_setequal(names(out$ELFs$rasWhole), c("5.1", "1.1", "3.1.1_2"))
})

test_that("the size rule is off without a land-cover map or when minRegionAreaKm2 is NA", {
  none <- suppressMessages(do.call(runFewFire, c(sizeArgs, list(minRegionAreaKm2 = 1e6))))
  expect_identical(none$plan$action, "skip")
  off <- suppressMessages(do.call(runFewFire, c(sizeArgs, list(landCoverFile = toyLandCoverFile(),
                                                             minRegionAreaKm2 = NA))))
  expect_identical(off$plan$action, "skip")
  expect_identical(formals(fewFireELFs)$minRegionAreaKm2, 35000)
  expect_identical(formals(fewFireELFs)$maxLandCoverDist, 0.35)
  expect_identical(formals(fewFireELFs)$maxBurnRatio, 6)
})

test_that("a size merge replaces the fire-count merge it contains, and a skipped ELF stays skipped", {
  fire <- data.table::data.table(
    action = c("merge", "skip"), ELF = c("3.1.1_2", NA), members = list(c("3.1.1", "3.1.2"), "5.1"),
    naturalIgnitions = c(4L, 1L), escapes = c(4L, 1L), firePolygons = c(2L, 0L), reason = c("few", "few"))
  status <- data.table::data.table(ELF = c("3.1.1", "3.1.2", "3.1.3", "5.1"), naturalIgnitions = c(1L, 3L, 7L, 1L),
                                   escapes = c(1L, 3L, 7L, 1L), firePolygons = c(1L, 1L, 2L, 0L))
  size <- data.table::data.table(action = "merge", ELF = "3.1.1_2_3", members = list(c("3.1.1", "3.1.2", "3.1.3")),
                                 reason = "smaller than 35000 km2; merged with 3.1.3", areaKm2 = 1)
  plan <- combinePlans(fire, size, status)
  expect_identical(plan$action, c("skip", "merge"))
  expect_identical(plan$ELF, c(NA, "3.1.1_2_3"))
  expect_identical(plan$naturalIgnitions, c(1L, 11L))
  expect_identical(combinePlans(fire, size[0, ], status), fire)
})

test_that("sizeRuleLandCover() downloads through LandR's SCANFI mirror, not Google Drive", {
  skip_if_not(exists(".withScanfiRemap", asNamespace("LandR")))
  withr::local_options(reproducible.urlRemap = NULL, reproducible.inputPaths = NULL, LandR.scanfiMirror = TRUE)
  remapped <- NULL
  local_mocked_bindings(
    preProcess = function(url, targetFile, ...) {
      remapped <<- getOption("reproducible.urlRemap")(url, filename = targetFile)
      list(targetFilePath = targetFile)
    },
    .package = "reproducible"
  )
  sizeRuleLandCover(withr::local_tempdir())
  expect_match(remapped, "^https://object-arbutus\\.")
})

test_that("sizeRuleLandCover() returns a local copy's path as a download returns it", {
  ## Cache keys on this path: through a symlink it must not differ from the normalized download path
  withr::local_options(reproducible.inputPaths = NULL)
  real <- withr::local_tempdir()
  link <- file.path(withr::local_tempdir(), "link")
  skip_if_not(file.symlink(real, link))
  file.create(file.path(real, LandR:::.scanfiLCCFAOSource(2020, "V2")$targetFile))
  expect_identical(sizeRuleLandCover(link),
                   normalizePath(file.path(real, LandR:::.scanfiLCCFAOSource(2020, "V2")$targetFile)))
})

test_that("ELFfireRates gives the yearly natural fires and area burned in each core, after the merges", {
  ## an extra 2002 point in column 4: core of 3.1.2, buffer of 3.1.1; it counts for 3.1.2 only
  local_toyFireRecords(points = rbind(toyFirePoints(), toyPoints(1, 4, 2002L, "L")))
  plan <- data.table::data.table(action = "merge", ELF = "3.1.1_2", members = list(c("3.1.1", "3.1.2")),
                                 reason = "too few fires")
  out <- ELFfireRates(toyELFs(), plan, fireYears = 2001:2003, pixelAreaHa = 2500, destinationPath = tempdir(),
                      nfdbUrl = "https://example.org/points.zip", nbacUrl = "https://example.org/polys.zip")
  out <- as.data.frame(out)[order(out$ELF), ]
  expect_identical(out$ELF, c("3.1.1_2", "5.1"))                 # merged; 1.1 is arctic: not counted
  ## 3.1.1: 1 point; 3.1.2: 3 + the column-4 point; 5.1: 1 -- over 3 years
  expect_equal(out$naturalFiresPerYear, c(5, 1) / 3)
  ## 5000 ha in each of 3.1.1 and 3.1.2, 1250 ha in 5.1 (geodesic areas, so not exact)
  expect_equal(out$haBurnedPerYear, c(10000, 1250) / 3, tolerance = 0.01)
})
