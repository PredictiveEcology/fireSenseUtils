## The gate must reproduce the counts that the 2026-09-11 investigation validated against the fit
## itself (3.2.2 = 0 natural ignitions, 10.1 = few), and it must not remove an ELF that has fire.

makeTestELFs <- function() {
  r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10000,
                   ymin = 0, ymax = 10000, crs = "EPSG:3978")
  busy <- terra::setValues(r, 2)            # fire everywhere
  empty <- terra::setValues(r, 2)           # same footprint, no fire put in it
  thin <- terra::setValues(r, 2)
  out <- c(busy, empty, thin)
  names(out) <- c("busy", "empty", "thin")
  out
}

## n points in the ELF footprint, spread over `years`.
makePoints <- function(n, years, cause = "N", sizeHa = 100) {
  if (n == 0) {
    v <- terra::vect(cbind(-1e6, -1e6), crs = "EPSG:3978")
    v$YEAR <- 1985L
    v$CAUSE <- cause
    v$SIZE_HA <- sizeHa[1]
    return(v)
  }
  xy <- cbind(seq(500, 9500, length.out = n), rep(5000, n))
  v <- terra::vect(xy, crs = "EPSG:3978")
  v$YEAR <- as.integer(rep_len(years, n))
  v$CAUSE <- cause
  v$SIZE_HA <- rep_len(sizeHa, n)
  v
}

## n square polygons well above one pixel, inside the footprint.
makePolys <- function(n, years) {
  if (n == 0) {
    e <- terra::vect("POLYGON ((-1e6 -1e6, -1e6 -9e5, -9e5 -9e5, -1e6 -1e6))", crs = "EPSG:3978")
    e$YEAR <- 1985L
    return(e)
  }
  ps <- lapply(seq_len(n), function(i) {
    x0 <- 500 + (i - 1) * 100
    terra::vect(sprintf("POLYGON ((%f 4000, %f 4900, %f 4900, %f 4000, %f 4000))",
                        x0, x0, x0 + 800, x0 + 800, x0),
                crs = "EPSG:3978")
  })
  v <- do.call(rbind, ps)
  v$YEAR <- as.integer(rep_len(years, n))
  v
}

test_that("ELFfireCounts counts points and polygons per ELF and year", {
  elfs <- makeTestELFs()
  years <- 1985:1989
  counts <- ELFfireCounts(elfs, makePoints(10, years), makePolys(5, years),
                          fireYears = years, pixelAreaHa = 5.76)

  expect_s3_class(counts, "data.table")
  expect_setequal(counts$ELF, names(elfs))
  expect_identical(nrow(counts), 3L * length(years))
  ## Every layer has the same footprint, so every layer sees the same fires.
  expect_identical(sum(counts$naturalIgnitions), 30L)
  expect_true(all(counts$firePolygons >= 0))
})

test_that("ELFfireCounts counts escapes: natural ignitions that escapedFires() keeps", {
  elfs <- makeTestELFs()
  years <- 1985:1989
  ## sizes 1, 5, 49.9, 50, 200 ha; pixel 5.76 ha. Only the 50 and 200 ha fires escape.
  counts <- ELFfireCounts(elfs, makePoints(10, years, sizeHa = c(1, 5, 49.9, 50, 200)),
                          makePolys(5, years), fireYears = years, pixelAreaHa = 5.76,
                          escapeSizeHa = 50)
  expect_identical(sum(counts$naturalIgnitions), 30L)
  expect_identical(sum(counts$escapes), 12L)
  ## a lower threshold keeps more; 5 ha is below the pixel area's 5.76 so 1 and 5 never escape
  low <- ELFfireCounts(elfs, makePoints(10, years, sizeHa = c(1, 5, 49.9, 50, 200)),
                       makePolys(5, years), fireYears = years, pixelAreaHa = 5.76, escapeSizeHa = 0)
  expect_identical(sum(low$escapes), 18L)
  ## human-caused fires never count as escapes
  human <- ELFfireCounts(elfs, makePoints(10, years, cause = "H"), makePolys(5, years),
                         fireYears = years, pixelAreaHa = 5.76)
  expect_identical(sum(human$escapes), 0L)
})

test_that("escapedFires keeps fires above one pixel and at least escapeSizeHa", {
  fires <- data.frame(id = 1:5, SIZE_HA = c(1, 5.76, 6, 50, 500))
  expect_identical(escapedFires(fires, escapeSizeHa = 50, pixSizeHa = 5.76)$id, 4:5)
  expect_identical(escapedFires(fires, escapeSizeHa = 0, pixSizeHa = 5.76)$id, 3:5)
  expect_identical(escapedFires(fires, 0, 0, sizeCol = "SIZE_HA")$id, 1:5)
  v <- terra::vect(cbind(1:5, 1:5), crs = "EPSG:3978")
  v$SIZE_HA <- fires$SIZE_HA
  expect_equal(nrow(escapedFires(v, 50, 5.76)), 2)
})

test_that("ELFfireOverlay clips fire polygons to each ELF's fitting area and core, once", {
  r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10000, ymin = 0, ymax = 10000, crs = "EPSG:3978")
  ## "half": core on the left 5 columns, buffer on the right 5; "far": a footprint no polygon reaches
  half <- terra::setValues(r, rep(rep(c(2, 1), each = 5), 10))
  far <- terra::setValues(r, 0)
  far[1, 1] <- 2
  elfs <- c(half, far)
  names(elfs) <- c("half", "far")
  ## 800 x 900 m = 72 ha: one polygon in the core, one across the core/buffer edge (x 4600-5400), one in the buffer
  polys <- terra::vect(c(
    "POLYGON ((500 4000, 500 4900, 1300 4900, 1300 4000, 500 4000))",
    "POLYGON ((4600 4000, 4600 4900, 5400 4900, 5400 4000, 4600 4000))",
    "POLYGON ((7000 4000, 7000 4900, 7800 4900, 7800 4000, 7000 4000))",
    "POLYGON ((500 4000, 500 4900, 1300 4900, 1300 4000, 500 4000))"), crs = "EPSG:3978")
  polys$YEAR <- c(1985L, 1985L, 1986L, 1999L)               # the 1999 polygon is outside fireYears
  ov <- ELFfireOverlay(elfs, polys, fireYears = 1985:1989)
  expect_s3_class(ov, "data.table")
  expect_named(ov, c("ELF", "zone", "YEAR", "ha"))
  expect_setequal(ov$ELF, "half")                             # nothing reaches "far"
  fit <- ov[ov$zone == "fit", ]
  core <- ov[ov$zone == "core", ]
  expect_equal(sort(fit$ha), c(72, 72, 72), tolerance = 1e-2)
  expect_equal(sort(core$ha), c(36, 72), tolerance = 1e-2)    # the edge polygon is half in the core
  expect_identical(sort(fit$YEAR), c(1985L, 1985L, 1986L))
  ## no polygon at all: an empty table with the same columns
  none <- ELFfireOverlay(elfs, polys[0, ], fireYears = 1985:1989)
  expect_identical(nrow(none), 0L)
  expect_named(none, names(ov))
})

test_that("ELFfireCounts and ELFregionStats read a given overlay and give the same answer as without", {
  elfs <- makeTestELFs()
  years <- 1985:1989
  pts <- makePoints(10, years)
  polys <- makePolys(5, years)
  ov <- ELFfireOverlay(elfs, polys, years)
  expect_identical(ELFfireCounts(elfs, pts, polys, years, pixelAreaHa = 5.76),
                   ELFfireCounts(elfs, pts, polys, years, pixelAreaHa = 5.76, overlay = ov))
  ## a polygon that does not clip to more than one pixel is not counted, whatever the overlay holds
  expect_identical(sum(ELFfireCounts(elfs, pts, polys, years, pixelAreaHa = 5000, overlay = ov)$firePolygons), 0L)

  lcc <- terra::setValues(elfs[[1]], 210)
  direct <- ELFregionStats(elfs, lcc, firePolys = polys, fireYears = years)
  viaOverlay <- ELFregionStats(elfs, lcc, fireYears = years, overlay = ov)
  expect_identical(viaOverlay, direct)
  ## a merged region's burn is the sum over its members; its burn rate is that over its own area
  merged <- max(elfs[[c("busy", "empty")]])
  names(merged) <- "busy_empty"
  st <- ELFregionStats(merged, lcc, fireYears = years, overlay = ov,
                       members = list(busy_empty = c("busy", "empty")))
  one <- direct$burnRate[1]
  expect_equal(st$burnRate, one * 2 / 1, tolerance = 1e-6)   # twice the burned area over the same core area
})

test_that("ELFfireCounts counts only natural causes", {
  elfs <- makeTestELFs()
  years <- 1985:1989
  human <- ELFfireCounts(elfs, makePoints(10, years, cause = "H"), makePolys(5, years),
                         fireYears = years, pixelAreaHa = 5.76)
  expect_identical(sum(human$naturalIgnitions), 0L)
  ## Polygons are all-cause, as in the spread fit, so they still count.
  expect_gt(sum(human$firePolygons), 0L)
})

test_that("ELFfireCounts ignores fires outside the fitted years", {
  elfs <- makeTestELFs()
  counts <- ELFfireCounts(elfs, makePoints(10, 1985:1989), makePolys(5, 1985:1989),
                          fireYears = 2000:2004, pixelAreaHa = 5.76)
  expect_identical(sum(counts$naturalIgnitions), 0L)
  expect_identical(sum(counts$firePolygons), 0L)
})

test_that("a polygon smaller than one pixel is not counted", {
  elfs <- makeTestELFs()
  years <- 1985L
  ## 800 x 900 m = 72 ha; a 5000 ha threshold excludes it.
  counts <- ELFfireCounts(elfs, makePoints(1, years), makePolys(1, years),
                          fireYears = years, pixelAreaHa = 5000)
  expect_identical(sum(counts$firePolygons), 0L)
})

test_that("ELFfitStatus separates zero, few and ok", {
  counts <- data.table::data.table(
    ELF = rep(c("noIgnitions", "noPolys", "thin", "busy"), each = 2),
    year = rep(1985:1986, 4),
    naturalIgnitions = c(0L, 0L, 5L, 5L, 3L, 4L, 60L, 70L),
    escapes = c(0L, 0L, 5L, 5L, 3L, 4L, 60L, 70L),
    firePolygons = c(4L, 4L, 0L, 0L, 2L, 2L, 80L, 90L)
  )
  status <- ELFfitStatus(counts, minNaturalIgnitions = 50, minFirePolygons = 50)

  expect_identical(status$status[status$ELF == "noIgnitions"], "zero")
  expect_identical(status$status[status$ELF == "noPolys"], "zero")
  expect_identical(status$status[status$ELF == "thin"], "few")
  expect_identical(status$status[status$ELF == "busy"], "ok")
  expect_identical(status$escapes[status$ELF == "busy"], 130L)
  ## Counts are summed over the whole window, not judged per year.
  expect_identical(status$naturalIgnitions[status$ELF == "busy"], 130L)
  expect_identical(status$yearsWithFire[status$ELF == "busy"], 2L)
})

test_that("an ELF with many small fires but one escape is few", {
  ## ELF 7.3: plenty of ignitions and polygons, 1 escape, so the escape fit has 1 positive.
  counts <- data.table::data.table(
    ELF = rep(c("smallFires", "bigFires"), each = 2), year = rep(1985:1986, 2),
    naturalIgnitions = c(40L, 40L, 40L, 40L), escapes = c(1L, 0L, 10L, 15L),
    firePolygons = c(40L, 40L, 40L, 40L)
  )
  status <- ELFfitStatus(counts)
  expect_identical(status$status[status$ELF == "smallFires"], "few")
  expect_identical(status$status[status$ELF == "bigFires"], "ok")
  expect_identical(ELFfitStatus(counts, minEscapes = 1)$status[1:2], c("ok", "ok"))
})

test_that("ELFsExcluded returns only zero ELFs", {
  status <- data.table::data.table(
    ELF = c("a", "b", "c"),
    status = c("zero", "few", "ok")
  )
  expect_identical(ELFsExcluded(status), "a")
})

test_that("a few-fire ELF is not excluded", {
  ## 10.1 has fire but not much: it must be flagged, never removed.
  counts <- data.table::data.table(
    ELF = "10.1", year = 1985:1986,
    naturalIgnitions = c(16L, 16L), escapes = c(16L, 16L), firePolygons = c(57L, 58L)
  )
  status <- ELFfitStatus(counts)
  expect_identical(status$status, "few")
  expect_identical(ELFsExcluded(status), character(0))
})

## Flammable area is a REPORT, not a gate (Eliot, 2026-09-11: no floor; SpreadFit may try).
## The cases that matter are the ones VLCE2 leaves unclassified, which measure 0.

makeLCC <- function(codes) {
  r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10000,
                   ymin = 0, ymax = 10000, crs = "EPSG:3978")
  terra::setValues(r, rep_len(codes, terra::ncell(r)))
}

test_that("ELFflammableArea counts only flammable codes", {
  elfs <- makeTestELFs()
  ## half coniferous (210, flammable), half water (20, not)
  fa <- ELFflammableArea(elfs, makeLCC(c(210, 20)))

  expect_s3_class(fa, "data.table")
  expect_setequal(fa$ELF, names(elfs))
  expect_true(all(abs(fa$flammableFraction - 0.5) < 0.01))
  expect_true(all(fa$flammablePixels > 0))
})

test_that("an ELF outside the land-cover extent measures zero, and is not an error", {
  ## VLCE2 code 0 over everything: this is 8.2, 10.3.2 and 3.2.4.
  fa <- ELFflammableArea(makeTestELFs(), makeLCC(0))
  expect_true(all(fa$flammablePixels == 0))
  expect_true(all(fa$flammableFraction == 0))
  ## No stop(), no exclusion: SpreadFit is still allowed to try.
  expect_identical(nrow(fa), 3L)
})

test_that("nonflammableLCC is honoured", {
  elfs <- makeTestELFs()
  ## Treat coniferous as nonflammable and nothing is left.
  fa <- ELFflammableArea(elfs, makeLCC(210), nonflammableLCC = c(0, 210))
  expect_true(all(fa$flammablePixels == 0))
  ## With the default, all of it burns.
  fa2 <- ELFflammableArea(elfs, makeLCC(210))
  expect_true(all(fa2$flammableFraction == 1))
})

test_that("flammable area is reported in hectares of the land-cover grid", {
  elfs <- makeTestELFs()
  fa <- ELFflammableArea(elfs, makeLCC(210))
  ## 1000 x 1000 m pixels = 100 ha each, 100 cells fully inside the ELF footprint.
  expect_equal(fa$flammableAreaHa[1], fa$flammablePixels[1] * 100)
})

## Track A step 6: an ELF with no SCANFI tree species must NOT be excluded. fireSense fits
## nonForest fuel classes, and the dataPrep chain now carries a zero-layer speciesLayers
## through, so the old noSCANFIELFs regex would silently withhold fittable ELFs.

test_that("the exclusion regex keeps no-tree ELFs and still drops Arctic and zero-fire ones", {
  ## Rebuilt exactly as runELFs() builds it (R/ELFs.R, excludeELFs).
  buildExclude <- function(zeroFireELFs = character(0)) {
    paste(c(fireSenseUtils:::.ELFsNotFittedPattern(),
            if (length(zeroFireELFs))
              paste0("^", gsub(".", "\\.", zeroFireELFs, fixed = TRUE), "$")),
          collapse = "|")
  }

  ## 3.2.1 and 3.2.4 have no SCANFI species; they must survive now.
  keep <- c("3.2.1", "3.2.4", "3.2.5", "3.3.2", "10.1", "6.2.2")
  expect_false(any(grepl(buildExclude(), keep)))

  ## Arctic is still excluded, and zero-fire ELFs still are when named.
  expect_true(all(grepl(buildExclude(), c("1.2.3", "2.5.4", "10.2", "10.2.1", "10.2.2"))))
  expect_false(any(grepl(buildExclude(), c("10.20", "10.21.1", "10.1.1", "12.1", "11.2"))))
  expect_true(all(grepl(buildExclude(c("3.2.2", "3.3.1")), c("3.2.2", "3.3.1"))))

  ## Anchoring: 3.2.2 must not also catch 3.2.25.
  expect_false(grepl(buildExclude(c("3.2.2")), "3.2.25"))
})
