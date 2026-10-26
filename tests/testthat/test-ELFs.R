## makeELFs() and the functions around it, on a small made-up map: three ecoprovinces side by
## side in EPSG:3978 (4.1 has a 20 x 20 km island just off 4.2's coast; 4.3 is twice the size of
## the others). The ecostratification downloads are stubbed with those polygons.

sq <- function(x0, x1, y0 = 0, y1 = 3e5)
  sprintf("POLYGON ((%1$s %3$s, %2$s %3$s, %2$s %4$s, %1$s %4$s, %1$s %3$s))", x0, x1, y0, y1)

makeTestEcos <- function() {
  provs <- c("4.1", "4.2", "4.3")
  withIsland <- paste0("MULTIPOLYGON (((0 0, 300000 0, 300000 300000, 0 300000, 0 0)), ",
                       "((400000 310000, 420000 310000, 420000 330000, 400000 330000, 400000 310000)))")
  ecos <- terra::vect(c(withIsland, sq(3e5, 6e5), sq(6e5, 12e5)), crs = "EPSG:3978")
  ecos$ECOZONE <- "4"
  ecos$ECOREGION <- ecos$ECOPROVINC <- ecos$ECODISTRIC <- provs
  ecos
}

makeTestFRU <- function(ecos) {
  x <- terra::rasterize(ecos, terra::rast(ecos, resolution = 5000), field = "ECOPROVINC")
  x <- as.numeric(x) + 1
  names(x) <- "FRU"
  x
}

## bufferOut() drops output directories under tempdir(), so the test's must be outside it.
localELFsDir <- function(env = parent.frame()) {
  dp <- withr::local_tempdir(tmpdir = dirname(tempdir()), .local_envir = env)
  withr::local_options(reproducible.useCache = FALSE,
                       reproducible.cachePath = file.path(dp, "cache"), .local_envir = env)
  dp
}

cellAt <- function(r, x, y) terra::extract(r, cbind(x, y))[[1]]

test_that("makeELFs builds one buffered ELF per ecoprovince, splitting the oversized one", {
  skip_if_not_installed("dismo")
  skip_if_not_installed("deldir")
  ecos <- makeTestEcos()
  dp <- localELFsDir()
  local_mocked_bindings(prepInputs = function(url, ...) ecos, .package = "reproducible")
  set.seed(1)

  ## split_poly()'s sf::st_intersection() warns about attributes; unrelated to the result
  expect_warning(
    out <- makeELFs(makeTestFRU(ecos), destinationPath = dp, useCache = FALSE, maxArea = 1e11),
    "spatially constant")

  expect_named(out, c("rasCentered", "rasWhole", "poly"))
  expect_identical(names(out$rasWhole), c("4.1", "4.2", "4.3.1", "4.3.2"))
  expect_identical(names(out$rasCentered), names(out$rasWhole))
  ## each ELF: core (2) and buffer (1) polygon
  expect_identical(sort(unique(out$poly$ID)), names(out$rasWhole))
  expect_setequal(out$poly$buffer, 1:2)

  w <- out$rasWhole
  expect_identical(cellAt(w$`4.1`, 150000, 150000), 2)      # its own core
  expect_identical(cellAt(w$`4.1`, 310000, 150000), 1)      # 10 km into 4.2: buffer
  expect_identical(cellAt(w$`4.1`, 450000, 150000), 0)      # 150 km into 4.2: outside
  expect_identical(cellAt(w$`4.2`, 290000, 150000), 1)      # buffer reaches back into 4.1

  ## the island is too small to be its own piece of 4.1: it moves to 4.2, whose buffer covers it
  expect_identical(cellAt(w$`4.1`, 410000, 320000), 0)
  expect_identical(cellAt(w$`4.2`, 410000, 320000), 2)

  ## the cores partition the land: every land cell is core in exactly one ELF. 4.1's convex hull
  ## (with its island) reaches into 4.2, and must not take 4.2's cells; projecting each core back
  ## from its centred map must not leave gaps or overlaps along shared borders.
  cores <- sapply(w, function(r) terra::values(r, mat = FALSE) %in% 2)
  ## land: the ecoprovinces' cells (makeELFs()'s focal() hole filling also grows the coast a little)
  land <- !is.na(terra::values(terra::rasterize(ecos, w[[1]]), mat = FALSE))
  expect_identical(sum(rowSums(cores) > 1), 0L)                     # no overlaps
  expect_identical(sum(land & rowSums(cores) == 0), 0L)             # no gaps
  expect_false(any(rowSums(cores)[!land] > 0))

  written <- basename(list.files(dp, recursive = TRUE, pattern = "^(r|ca)_.*tif$"))
  expect_setequal(written, paste0(c("r_", "ca_"), rep(names(w), each = 2), ".tif"))
})

test_that("makeELFs works when its cached rasters come back from disk", {
  skip_if_not_installed("dismo")
  skip_if_not_installed("deldir")
  ## a cache restore reads a raster back from a .tif, which names the first column of its categories
  ## "value" instead of "ID"; the ecoprovince names must still be found
  ecos <- makeTestEcos()
  dp <- localELFsDir()
  local_mocked_bindings(prepInputs = function(url, ...) ecos, .package = "reproducible")
  fromDisk <- function(x) {
    if (inherits(x, "SpatRaster")) {
      f <- tempfile(tmpdir = dp, fileext = ".tif")
      terra::writeRaster(x, f)
      terra::rast(f)
    } else if (is.list(x) && !is.data.frame(x)) {
      structure(lapply(x, fromDisk), names = names(x))
    } else x
  }
  local_mocked_bindings(Cache = function(x, ...) fromDisk(x), .package = "reproducible")
  local_mocked_bindings(Cache = function(x, ...) fromDisk(x))
  set.seed(1)
  out <- suppressWarnings(makeELFs(makeTestFRU(ecos), destinationPath = dp, useCache = FALSE, maxArea = 1e11))
  expect_identical(names(out$rasWhole), c("4.1", "4.2", "4.3.1", "4.3.2"))
})

test_that("makeELFs accepts fire regime polygons (sf) as well as a raster", {
  ecos <- makeTestEcos()[2:3, ]
  ecos$ECOPROVINC <- ecos$ECOREGION <- ecos$ECODISTRIC <- c("4.1", "4.2")
  fru <- sf::st_as_sf(ecos)
  fru$FRU <- c(1, 2)
  dp <- localELFsDir()
  local_mocked_bindings(prepInputs = function(url, ...) ecos, .package = "reproducible")

  out <- makeELFs(fru, destinationPath = dp, useCache = FALSE)

  expect_identical(names(out$rasWhole), c("4.1", "4.2"))
  expect_identical(cellAt(out$rasWhole$`4.2`, 900000, 150000), 2)
})

test_that("makeELFs stops when the ecoprovince raster has no categories", {
  ecos <- makeTestEcos()
  ecos$ECOPROVINC <- c(4.1, 4.2, 4.3)
  dp <- localELFsDir()
  local_mocked_bindings(prepInputs = function(url, ...) ecos, .package = "reproducible")
  expect_error(makeELFs(makeTestFRU(ecos), destinationPath = dp, useCache = FALSE),
               "ecoprovince raster has no categories")
})

## Two ELFs sharing a border: A's core is columns 1-4, B's is 7-10, each with a 2-column buffer.
makeTestELFRasters <- function() {
  r <- terra::rast(nrows = 4, ncols = 10, xmin = 0, xmax = 10000, ymin = 0, ymax = 4000,
                   crs = "EPSG:3978")
  a <- terra::setValues(r, rep(c(2, 2, 2, 2, 1, 1, 0, 0, 0, 0), 4))
  b <- terra::setValues(r, rep(c(0, 0, 0, 0, 1, 1, 2, 2, 2, 2), 4))
  list(rasWhole = list(A = a, B = b))
}

test_that("ELFsInStudyArea labels cores, shared buffer and background", {
  ras <- makeTestELFRasters()
  poly <- do.call(rbind, unname(Map(r = ras$rasWhole, nam = names(ras$rasWhole), function(r, nam) {
    v <- terra::as.polygons(terra::classify(r, cbind(0, NA)))
    v$ID <- nam
    v
  })))
  studyArea <- terra::vect(sq(0, 10000, 0, 4000), crs = "EPSG:3978")

  out <- ELFsInStudyArea(studyArea, inputPath = tempdir(), ELFsRaster = ras, ELFsPolygon = poly)

  expect_true(terra::is.factor(out$rast))
  labs <- terra::levels(out$rast)[[1]]
  expect_setequal(labs$ELFind, c("none", "A", "B"))
  row1 <- as.character(terra::as.data.frame(out$rast[1, ], na.rm = FALSE)[[1]])
  expect_identical(row1[1:4], rep("A", 4))
  expect_identical(row1[7:10], rep("B", 4))
  expect_true(all(row1[5:6] %in% c("A", "B")))   # buffers of both: labelled, not "none"
  expect_setequal(out$poly$ID, c("A", "B"))
})

## One row of 12 cells: A's core is cells 1-3 and its buffer 4-8; B's core is 10-12, buffer 5-9.
overlapELFs <- function() {
  r <- terra::rast(nrows = 1, ncols = 12, xmin = 0, xmax = 12000, ymin = 0, ymax = 1000,
                   crs = "EPSG:3978")
  a <- terra::setValues(r, c(2, 2, 2, 1, 1, 1, 1, 1, 0, 0, 0, 0))
  b <- terra::setValues(r, c(0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2))
  poly <- terra::vect(c(sq(0, 8000, 0, 1000), sq(4000, 12000, 0, 1000)), crs = "EPSG:3978")
  poly$ID <- c("A", "B")
  list(ras = list(rasWhole = list(A = a, B = b)), poly = poly)
}
labelsOf <- function(r) as.character(terra::as.data.frame(r, na.rm = FALSE)[[1]])

test_that("ELFsInStudyArea gives a shared buffer cell to the ELF with the nearest core", {
  x <- overlapELFs()
  studyArea <- terra::vect(sq(0, 12000, 0, 1000), crs = "EPSG:3978")
  out <- ELFsInStudyArea(studyArea, tempdir(), ELFsRaster = x$ras, ELFsPolygon = x$poly)
  ## cells 5-6 are nearer A's core, 7-8 nearer B's; 4 and 9 are buffer of one ELF only
  expect_identical(labelsOf(out$rast), c(rep("A", 6), rep("B", 6)))
})

test_that("ELFsInStudyArea works when the study area overlaps one ELF", {
  x <- overlapELFs()
  studyArea <- terra::vect(sq(0, 3000, 0, 1000), crs = "EPSG:3978")   # A's core only
  poly <- x$poly[x$poly$ID == "A"]
  out <- ELFsInStudyArea(studyArea, tempdir(), ELFsRaster = x$ras, ELFsPolygon = poly)
  expect_identical(labelsOf(out$rast), c(rep("A", 8), rep("none", 4)))
})

test_that("ELFsInStudyArea builds the ELF polygons when none are given", {
  ras <- makeTestELFRasters()
  poly <- terra::vect(c(sq(0, 6000, 0, 4000), sq(4000, 10000, 0, 4000)), crs = "EPSG:3978")
  poly$ID <- c("A", "B")
  called <- NULL
  local_mocked_bindings(
    ELFtemplateRaster = function(inputPath) {
      called <<- inputPath
      ras$rasWhole$A
    },
    makeELFs = function(x, ...) list(poly = poly),
    Cache = function(x, ...) x
  )
  studyArea <- terra::vect(sq(0, 10000, 0, 4000), crs = "EPSG:3978")

  out <- ELFsInStudyArea(studyArea, inputPath = "here", ELFsRaster = ras)

  expect_identical(called, "here")
  expect_setequal(out$poly$ID, c("A", "B"))
})

test_that("ELFtemplateRaster aggregates the downloaded raster by 8 and writes it", {
  dp <- localELFsDir()
  src <- terra::rast(nrows = 16, ncols = 16, xmin = 0, xmax = 16, ymin = 0, ymax = 16,
                     vals = 1, crs = "EPSG:3978")
  local_mocked_bindings(getRemoteMetadata = function(...) list(remoteHash = "abc"),
                        .package = "reproducible")
  local_mocked_bindings(prepInputs = function(url, destinationPath, ...) src)

  out <- ELFtemplateRaster(dp)

  expect_equal(dim(out)[1:2], c(2, 2))
  expect_identical(normalizePath(terra::sources(out)),   # macOS: /var is /private/var
                   normalizePath(file.path(dp, "rastTemplate_Canada.tif")))
})

test_that("ELFtemplateRaster does not rewrite the template another run may be reading", {
  ## runs share inputPath; a rewrite while another run read the file gave
  ## "[trim] only cells with NA found" in fireSense_ELFs (6.6.1 hindcast, 2026-10-25)
  ## Under covr the second call misses the cache (the instrumented code is digested), so the file is rewritten;
  ## R CMD check runs it uninstrumented on every platform.
  testthat::skip_on_covr()
  dp <- localELFsDir()
  src <- terra::rast(nrows = 16, ncols = 16, xmin = 0, xmax = 16, ymin = 0, ymax = 16,
                     vals = 1, crs = "EPSG:3978")
  local_mocked_bindings(getRemoteMetadata = function(...) list(remoteHash = "abc"),
                        .package = "reproducible")
  local_mocked_bindings(prepInputs = function(url, destinationPath, ...) src)
  withr::local_options(reproducible.cachePath = withr::local_tempdir(), reproducible.useCache = TRUE)
  f <- file.path(dp, "rastTemplate_Canada.tif")

  first <- ELFtemplateRaster(dp)
  old <- as.POSIXct("2020-01-01", tz = "UTC")
  Sys.setFileTime(f, old)
  second <- ELFtemplateRaster(dp)

  expect_equal(as.numeric(file.mtime(f)), as.numeric(old))
  ## a cache hit may return the cache's own copy of the file, not inputPath's: what matters is the file is not rewritten
  expect_equal(terra::values(second, mat = FALSE), terra::values(first, mat = FALSE))
})

test_that("ELFtemplateRaster defaults to the SCANFI v2 Douglas-fir file on Arbutus", {
  expect_identical(formals(ELFtemplateRaster)$url,
                   "https://object-arbutus.alliancecan.ca/529feda02d404836b0e91b1d360b4df0:predictiveecology/SCANFI_v2/2020/SCANFI_spsCC_PSEU_MEN_2020_v2_20260119.tif")
  dp <- localELFsDir()
  src <- terra::rast(nrows = 16, ncols = 16, xmin = 0, xmax = 16, ymin = 0, ymax = 16,
                     vals = 1, crs = "EPSG:3978")
  got <- NULL
  local_mocked_bindings(getRemoteMetadata = function(...) list(remoteHash = "abc"),
                        .package = "reproducible")
  local_mocked_bindings(prepInputs = function(url, destinationPath, ...) {
    got <<- url
    src
  })
  ELFtemplateRaster(dp)
  expect_identical(got, formals(ELFtemplateRaster)$url)
})

test_that("ELFtemplateRaster returns a 1/NA mask of where the aggregated source has data", {
  dp <- localELFsDir()
  v <- rep(37, 16 * 16)                       # cover values, not 1
  src <- terra::rast(nrows = 16, ncols = 16, xmin = 0, xmax = 16, ymin = 0, ymax = 16,
                     vals = v, crs = "EPSG:3978")
  src[1:8, 1:8] <- NA                         # top-left aggregation block has no data
  src[9:16, 9:16] <- 0                        # zero cover is still data
  local_mocked_bindings(getRemoteMetadata = function(...) list(remoteHash = "abc"),
                        .package = "reproducible")
  local_mocked_bindings(prepInputs = function(url, destinationPath, ...) src)

  out <- ELFtemplateRaster(dp)

  expect_equal(terra::values(out, mat = FALSE), c(NA, 1, 1, 1))
})

test_that("ELFtemplateRaster downloads the url it is given, and caches per url", {
  dp <- localELFsDir()
  src <- terra::rast(nrows = 16, ncols = 16, xmin = 0, xmax = 16, ymin = 0, ymax = 16,
                     vals = 1, crs = "EPSG:3978")
  got <- character()
  local_mocked_bindings(getRemoteMetadata = function(...) list(remoteHash = "abc"),
                        .package = "reproducible")
  local_mocked_bindings(prepInputs = function(url, destinationPath, ...) {
    got <<- c(got, url)
    src
  })
  withr::local_options(reproducible.cachePath = withr::local_tempdir(), reproducible.useCache = TRUE)
  ELFtemplateRaster(dp, url = "https://example.org/a.tif")
  ELFtemplateRaster(dp, url = "https://example.org/b.tif")
  expect_identical(unique(got), c("https://example.org/a.tif", "https://example.org/b.tif"))
  ## covr rewrites function bodies to count lines, which changes what Cache() digests, so the
  ## "served from cache" count only holds outside coverage runs
  testthat::skip_on_covr()
  ELFtemplateRaster(dp, url = "https://example.org/a.tif")   # cached
  expect_identical(got, c("https://example.org/a.tif", "https://example.org/b.tif"))
})

test_that("moveSliversToOtherELFs reports pixels no other ELF covers", {
  r <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2,
                   vals = c(2, 2, 0, 0), crs = "EPSG:3978")
  lost <- list(data.table::data.table(pixelID = 3:4, value = c(2, 2)))
  expect_message(out <- moveSliversToOtherELFs(lost, ca = list(r), i = 1, r = list(r)),
                 "From ELF 1, lost 4 isolated pixels")
  expect_identical(terra::values(out$ca[[1]], mat = FALSE), c(2, 2, 0, 0))
})

## runELFs(): the module run is stubbed; the simList carries the ELF names.
localRunELFs <- function(excluded = NULL, env = parent.frame()) {
  sim <- SpaDES.core::simInit()
  rasList <- setNames(as.list(1:4), c("1.1", "4.1", "4.2", "6.1"))
  sim$ELFs <- list(rasCentered = rasList, rasWhole = rasList)
  sim$spreadFitPreRun <- setNames(data.frame(c("2.1", "4.1", "4.2")), polygonIDTxt)
  sim$ELFsExcluded <- excluded
  SpaDES.core::outputs(sim) <- data.frame(objectName = c("ELFs", "ELFsPostMerge", "studyAreas"),
                                          file = c("ELF_polygons.rds", "ELFs_postMerge.rds", "studyAreas4.rds"))
  sim@outputs$file <- sub("\\.rds$", ".png", sim@outputs$file)  # the pngs Plots() saves
  local_mocked_bindings(simInitAndSpades2 = function(l) sim, .package = "SpaDES.core", .env = env)
  invisible(sim)
}

test_that("runELFs returns fitted, all, or map ELFs without the excluded ones", {
  prj <- list(modules = c("fireSense_ELFs", "other"),
              paths = list(modulePath = withr::local_tempdir()),
              params = list(fireSense_ELFs = list()))
  localRunELFs(excluded = "4.2")

  expect_identical(runELFs(prj), "4.1")                         # 2.1 is Arctic, 4.2 has no fire
  expect_identical(runELFs(prj, whatOut = "all"), c("4.1", "6.1"))
  maps <- runELFs(prj, whatOut = "maps")
  expect_named(maps$rasWhole, c("4.1", "6.1"))
})

test_that("runELFs accepts partial whatOut matching (\"fit\" / \"allNames\")", {
  prj <- list(modules = "fireSense_ELFs", paths = list(modulePath = withr::local_tempdir()),
              params = list(fireSense_ELFs = list()))
  localRunELFs(excluded = "4.2")

  expect_identical(runELFs(prj, whatOut = "fit"), "4.1")
  expect_identical(runELFs(prj, whatOut = "allNames"), c("4.1", "6.1"))
})

test_that("runELFs maps path keeps all ELFs when none are excluded (Arctic excluded as always)", {
  prj <- list(modules = "fireSense_ELFs", paths = list(modulePath = withr::local_tempdir()),
              params = list(fireSense_ELFs = list()))
  localRunELFs(excluded = NULL)

  maps <- runELFs(prj, whatOut = "maps")
  expect_named(maps$rasWhole, c("4.1", "4.2", "6.1"))   # 1.1 is Arctic, excluded
})

test_that("runELFs runs the module on every call, with no Cache, Drive hash or png upload", {
  prj <- list(modules = "fireSense_ELFs", paths = list(modulePath = withr::local_tempdir()),
              params = list(fireSense_ELFs = list()))
  sim <- localRunELFs()
  runs <- 0L
  seen <- NULL
  local_mocked_bindings(simInitAndSpades2 = function(l) { runs <<- runs + 1L; seen <<- l; sim },
                        .package = "SpaDES.core")
  local_mocked_bindings(getRemoteMetadata = function(...) stop("runELFs must not read Drive"),
                        .package = "reproducible")
  local_mocked_bindings(drive_update = function(...) stop("runELFs must not upload"),
                        .package = "googledrive")

  runELFs(prj)
  runELFs(prj)

  expect_identical(runs, 2L)
  expect_null(seen$params$fireSense_ELFs$hashSpreadFitRemoteFile)
})

test_that("mergeAndSplitRas gives the same split whatever RNG state or kind the caller has", {
  skip_if_not_installed("dismo")
  skip_if_not_installed("deldir")
  ecos <- makeTestEcos()
  dp <- localELFsDir()
  tmpl <- terra::rast(ecos, resolution = 5000)
  ecopRseg <- terra::rast(lapply(seq_len(nrow(ecos)), function(i) terra::rasterize(ecos[i], tmpl, field = 1)))
  names(ecopRseg) <- ecos$ECOPROVINC
  split <- function() {
    out <- suppressWarnings(mergeAndSplitRas(ecopRseg, ecos, maxArea = 1e11, useCache = FALSE,
                                             destinationPath = dp))
    terra::values(out[["4.3"]], mat = FALSE)
  }

  withr::with_seed(1, ref <- split())
  expect_gt(length(unique(stats::na.omit(ref))), 1)       # 4.3 really was split
  withr::with_seed(99, expect_identical(split(), ref))
  withr::with_seed(5, expect_identical(split(), ref), .rng_kind = "Wichmann-Hill")

  ## the caller's stream is restored, so the draw after the call does not depend on the call
  set.seed(7); expected <- stats::runif(1)
  set.seed(7); split()
  expect_identical(stats::runif(1), expected)
})

test_that("makeELFs reads the ecostratification as a SpatVector whatever reproducible.shapefileRead says", {
  skip_if_not_installed("dismo")
  skip_if_not_installed("deldir")
  ecos <- makeTestEcos()
  dp <- localELFsDir()
  funs <- list()
  local_mocked_bindings(prepInputs = function(url, ..., fun = NULL) {
    funs[[length(funs) + 1]] <<- fun
    ecos
  }, .package = "reproducible")
  suppressWarnings(makeELFs(makeTestFRU(ecos), destinationPath = dp, useCache = FALSE, maxArea = 1e11))
  expect_length(funs, 4)
  expect_true(all(vapply(funs, identical, logical(1), terra::vect)))
})

test_that("makeELFs on the real ecostratification gives the ELF maps of the SCANFI v2 template", {
  skip_on_cran()
  skip_on_ci()   # downloads the SCANFI template and four ecostratification shapefiles; takes many minutes
  skip_if_offline("sis.agr.gc.ca")
  skip_if_not_installed("dismo")
  skip_if_not_installed("deldir")
  skip_if_not_installed("digest")
  dp <- withr::local_tempdir(tmpdir = dirname(tempdir()))
  ## cold cache and a caller RNG state unlike the campaign's; the result must not depend on either
  withr::local_options(reproducible.cachePath = file.path(dp, "cache"),
                       reproducible.destinationPathShared = file.path(dp, "shared"),
                       reproducible.useCloud = FALSE)
  withr::local_seed(456)
  runif(1000)
  ## one call only: a second makeELFs() in a fresh cache can fail with "no categories" (reproducible)
  out <- makeELFs(ELFtemplateRaster(dp), desiredBuffer = 20000, destinationPath = dp)
  ## the ca_<id>.tif files of ELFs_final are out$rasWhole
  f <- vapply(out$rasWhole, terra::sources, character(1))
  ids <- sub("^ca_(.*)\\.tif$", "\\1", basename(f))
  f <- f[order(ids)]
  sums <- vapply(f, function(p) {
    r <- terra::rast(p)
    digest::digest(list(round(as.vector(terra::ext(r)), 2), dim(r), terra::values(r, mat = FALSE)),
                   algo = "xxhash64")
  }, character(1))
  names(sums) <- sort(ids)
  expect_length(sums, 79)
  ## "ec1a815c01d736e6" was the 2026-10-05 campaign's maps (SCANFI v1.2 template); "3dd5c361ed1acdda" the
  ## SCANFI v2 ones (82 ELFs) before the cores tiled the map: the convex hull of a multi-polygon province took
  ## its neighbours' cells, and 1.2.1, 1.2.2 and 2.3.3 were made of hull cells outside the land template only
  expect_identical(digest::digest(sums, algo = "xxhash64"), "f6bfdbfd1c6bf726")
})
