## One table says which NFDB CAUSE and NBAC FIRECAUS values are which cause. The codes below are those of the
## real releases: NFDB 20260811 (H, N, U) and NBAC 20260513 (Human, Natural, Undetermined).

test_that("fireCauseCodes: natural is today's L/N and Natural; unknown is U and Undetermined", {
  expect_identical(fireCauseCodes("natural", "nfdb"), c("L", "N"))
  expect_identical(fireCauseCodes("natural", "nbac"), "Natural")
  nf <- fireCauseCodes(c("natural", "unknown"), "nfdb")
  nb <- fireCauseCodes(c("natural", "unknown"), "nbac")
  expect_true(all(c("N", "U", NA) %in% nf) && !"H" %in% nf)
  expect_true(all(c("Natural", "Undetermined", NA) %in% nb) && !"Human" %in% nb)
  real <- c("H", "N", "U")
  expect_identical(real[real %in% fireCauseCodes("unknown", "nfdb")], "U")
  real <- c("Human", "Natural", "Undetermined")
  expect_identical(real[real %in% fireCauseCodes("unknown", "nbac")], "Undetermined")
  expect_error(fireCauseCodes("lightning"), "unknown fire cause group")
})

test_that("hindcastObserved and ELFfireCounts: fireCauses = NULL is as before, unknown adds the unknown-cause fires", {
  core <- terra::setValues(terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10000, ymin = 0, ymax = 10000,
                                       crs = "EPSG:3978"), 2)
  pts <- terra::vect(cbind(c(1000, 2000, 3000, 4000), 5500), crs = "EPSG:3978")
  pts$YEAR <- 2000L
  pts$CAUSE <- c("N", "U", "H", "U")
  pts$SIZE_HA <- 10
  sq <- function(x0, y0, w) sprintf("POLYGON ((%f %f, %f %f, %f %f, %f %f, %f %f))",
                                    x0, y0, x0 + w, y0, x0 + w, y0 + w, x0, y0 + w, x0, y0)
  polys <- terra::vect(c(sq(500, 500, 1000), sq(3500, 500, 1000), sq(6500, 500, 1000)), crs = "EPSG:3978")
  polys$YEAR <- 2000L
  polys$FIRECAUS <- c("Natural", "Undetermined", "Human")
  flam <- terra::setValues(terra::rast(core), 1)
  obs <- function(...) hindcastObserved(core, pts, polys, 2000, pixelAreaHa = 5.76, flammableRTM = flam, ...)$counts

  expect_identical(obs()$ignitions, 1L)                                     # N only, as before
  expect_identical(obs(fireCauses = "natural")$ignitions, 1L)
  expect_identical(obs(fireCauses = c("natural", "unknown"))$ignitions, 3L) # N, U, U
  expect_equal(obs()$areaBurnedHa, 100, tolerance = 0.03)
  expect_equal(obs(fireCauses = c("natural", "unknown"))$areaBurnedHa, 200, tolerance = 0.03)

  cnt <- function(...) ELFfireCounts(stats::setNames(core, "e"), pts, polys, 2000, pixelAreaHa = 5.76, ...)$naturalIgnitions
  expect_equal(cnt(), 1)
  expect_equal(cnt(fireCauses = c("natural", "unknown")), 3)
})
