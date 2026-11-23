## modelledFires(): the one place that picks the fires of the modelled causes and drops prescribed burns.
## Values are those of the real files: NBAC PRESCRIBED is "true" or NA; NFDB points have CAUSE H-PB (to 20250519) and
## PRESCRIBED "PB" / FIRE_TYPE "PB" on N, U and H rows (20260811).

nbacPolys <- function() {
  sq <- function(x0) sprintf("POLYGON ((%f 0, %f 0, %f 1000, %f 1000, %f 0))", x0, x0 + 1000, x0 + 1000, x0, x0)
  p <- terra::vect(vapply(c(0, 2000, 4000, 6000, 8000), sq, ""), crs = "EPSG:3978")
  p$YEAR <- 2000L
  p$FIRECAUS <- c("Natural", "Natural", "Human", "Undetermined", "Natural")
  p$PRESCRIBED <- c(NA, "true", "true", NA, NA)
  p$id <- 1:5
  p
}

test_that("modelledFires drops a prescribed natural-cause NBAC polygon, and one of an unmodelled cause", {
  p <- nbacPolys()
  expect_identical(modelledFires(p, "natural", "nbac")$id, c(1L, 5L))             # id 2 prescribed, id 3 human, 4 unknown
  expect_identical(modelledFires(p, c("natural", "unknown"), "nbac")$id, c(1L, 4L, 5L))
  expect_identical(modelledFires(p, c("natural", "unknown", "human"), "nbac")$id, c(1L, 4L, 5L)) # prescribed human too
  expect_identical(modelledFires(p, source = "nbac")$id, c(1L, 4L, 5L))           # no cause cut: prescribed only
})

test_that("modelledFires drops NFDB prescribed burns: CAUSE H-PB, PRESCRIBED PB, FIRE_TYPE PB", {
  pts <- terra::vect(cbind(1:6, 1), crs = "EPSG:3978")
  pts$CAUSE <- c("N", "H-PB", "N", "U", "U", "H")
  pts$PRESCRIBED <- c(NA, "PB", "PB", NA, "No", "0")
  pts$FIRE_TYPE <- c("Fire", "PB", "Fire", "Prescribed Burn", "Fire", "Fire")
  pts$id <- 1:6
  expect_identical(modelledFires(pts, c("natural", "unknown", "human"), "nfdb")$id, c(1L, 5L, 6L))
  expect_identical(modelledFires(pts, c("natural", "unknown"), "nfdb")$id, c(1L, 5L))
  expect_identical(modelledFires(pts, "human", "nfdb")$id, 6L)
})

test_that("modelledFires: a missing cause column is an error, NULL stays NULL, explicit codes are used", {
  p <- nbacPolys()
  p$FIRECAUS <- NULL
  expect_error(modelledFires(p, "natural", "nbac"), "FIRECAUS")
  expect_identical(modelledFires(p, source = "nbac")$id, c(1L, 4L, 5L)) # no cause asked: no cause needed (id 3 prescribed)
  expect_null(modelledFires(NULL, "natural"))
  expect_identical(modelledFires(nbacPolys(), codes = "Undetermined", source = "nbac")$id, 4L)
})

test_that("a prescribed fire is in no observed record: hindcastObserved, escapes and ELF counts", {
  core <- terra::setValues(terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10000, ymin = 0, ymax = 10000,
                                       crs = "EPSG:3978"), 2)
  flam <- terra::setValues(terra::rast(core), 1)
  sq <- function(x0, y0, w) sprintf("POLYGON ((%f %f, %f %f, %f %f, %f %f, %f %f))",
                                    x0, y0, x0 + w, y0, x0 + w, y0 + w, x0, y0 + w, x0, y0)
  polys <- terra::vect(c(sq(500, 500, 1000), sq(3500, 500, 1000), sq(6500, 500, 1000)), crs = "EPSG:3978")
  polys$YEAR <- 2000L
  polys$FIRECAUS <- c("Natural", "Natural", "Human")
  polys$PRESCRIBED <- c(NA, "true", "true")
  pts <- terra::vect(cbind(c(1000, 4000), 1000), crs = "EPSG:3978")
  pts$YEAR <- 2000L
  pts$CAUSE <- c("N", "N")
  pts$SIZE_HA <- 100
  pts$PRESCRIBED <- c(NA, "PB")
  obs <- hindcastObserved(core, pts, polys, 2000, pixelAreaHa = 5.76, flammableRTM = flam)$counts
  expect_identical(obs$ignitions, 1L)
  expect_equal(obs$areaBurnedHa, 100, tolerance = 0.03)          # the natural, not prescribed one
  expect_equal(obs$areaBurnedHaAllCause, 100, tolerance = 0.03)  # all-cause has no prescribed burn either
  obs2 <- hindcastObserved(core, pts, polys, 2000, pixelAreaHa = 5.76, flammableRTM = flam,
                           fireCauses = "natural", escapeSizeHa = 50, coreStart = TRUE)$counts
  expect_identical(obs2$ignitions, 1L)
  expect_identical(obs2$escapes, 1L)
  expect_equal(NROW(escapePolys(polys, 50, 5.76, causes = fireCauseCodes("natural", "nbac"))), 1)
  cnt <- ELFfireCounts(stats::setNames(core, "e"), pts, polys, 2000, pixelAreaHa = 5.76)
  expect_equal(cnt$naturalIgnitions, 1)
  expect_equal(sum(cnt$firePolygons), 1)                         # the polygon count has no prescribed fire either
})
