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
makePoints <- function(n, years, cause = "N") {
  if (n == 0) {
    v <- terra::vect(cbind(-1e6, -1e6), crs = "EPSG:3978")
    v$YEAR <- 1985L
    v$CAUSE <- cause
    return(v)
  }
  xy <- cbind(seq(500, 9500, length.out = n), rep(5000, n))
  v <- terra::vect(xy, crs = "EPSG:3978")
  v$YEAR <- as.integer(rep_len(years, n))
  v$CAUSE <- cause
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
    firePolygons = c(4L, 4L, 0L, 0L, 2L, 2L, 80L, 90L)
  )
  status <- ELFfitStatus(counts, minNaturalIgnitions = 50, minFirePolygons = 50)

  expect_identical(status$status[status$ELF == "noIgnitions"], "zero")
  expect_identical(status$status[status$ELF == "noPolys"], "zero")
  expect_identical(status$status[status$ELF == "thin"], "few")
  expect_identical(status$status[status$ELF == "busy"], "ok")
  ## Counts are summed over the whole window, not judged per year.
  expect_identical(status$naturalIgnitions[status$ELF == "busy"], 130L)
  expect_identical(status$yearsWithFire[status$ELF == "busy"], 2L)
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
    naturalIgnitions = c(16L, 16L), firePolygons = c(57L, 58L)
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
    paste(c("^1\\.|^2\\.",
            if (length(zeroFireELFs))
              paste0("^", gsub(".", "\\.", zeroFireELFs, fixed = TRUE), "$")),
          collapse = "|")
  }

  ## 3.2.1 and 3.2.4 have no SCANFI species; they must survive now.
  keep <- c("3.2.1", "3.2.4", "3.2.5", "3.3.2", "10.1", "6.2.2")
  expect_false(any(grepl(buildExclude(), keep)))

  ## Arctic is still excluded, and zero-fire ELFs still are when named.
  expect_true(all(grepl(buildExclude(), c("1.2.3", "2.5.4"))))
  expect_true(all(grepl(buildExclude(c("3.2.2", "3.3.1")), c("3.2.2", "3.3.1"))))

  ## Anchoring: 3.2.2 must not also catch 3.2.25.
  expect_false(grepl(buildExclude(c("3.2.2")), "3.2.25"))
})
