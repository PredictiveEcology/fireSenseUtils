## The hindcast compares simulated rep-years with what was observed in the ELF core. The toy world is a
## 10 x 10 grid of 1 km cells (100 ha each) whose left five columns are the core.

toyCore <- function() {
  r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10000, ymin = 0, ymax = 10000, crs = "EPSG:3978")
  terra::setValues(r, ifelse(terra::colFromCell(r, seq_len(terra::ncell(r))) <= 5, 2, NA))
}

## Year 2000: natural points of 10, 60 and 5 ha in the core, one natural point outside the core, one human
## point in the core; a natural point in 2001.
toyPoints <- function() {
  v <- terra::vect(cbind(c(1000, 2000, 3000, 8000, 1500, 1500), rep(5500, 6)), crs = "EPSG:3978")
  v$YEAR <- c(2000L, 2000L, 2000L, 2000L, 2000L, 2001L)
  v$CAUSE <- c("N", "L", "N", "N", "H", "N")
  v$SIZE_HA <- c(10, 60, 5, 80, 70, 90)
  v
}

## Year 2000: a natural polygon wholly in the core (about 100 ha), a natural one half in the core (about 50 ha
## of 100 ha), a human one in the core (about 200 ha); 2001: a natural one in the core.
toyPolys <- function() {
  sq <- function(x0, y0, w, h) sprintf("POLYGON ((%f %f, %f %f, %f %f, %f %f, %f %f))",
                                       x0, y0, x0 + w, y0, x0 + w, y0 + h, x0, y0 + h, x0, y0)
  v <- terra::vect(c(sq(500, 500, 1000, 1000), sq(4500, 3000, 1000, 1000), sq(500, 7000, 1000, 2000),
                     sq(500, 3000, 1000, 1000)), crs = "EPSG:3978")
  v$YEAR <- c(2000L, 2000L, 2000L, 2001L)
  v$FIRECAUS <- c("Natural", "Natural", "Human", "Natural")
  v
}

test_that("hindcastObserved counts ignitions, escapes and core-clipped area with the fit's definitions", {
  o <- hindcastObserved(toyCore(), toyPoints(), toyPolys(), year = 2000, escapeSizeHa = 50, pixelAreaHa = 5.76)
  cn <- o$counts
  expect_identical(cn$year, 2000L)
  expect_identical(cn$ignitions, 3L)            # N, L, N in the core: not the outside one, the human one or 2001
  expect_identical(cn$escapes, 1L)              # only the 60 ha point reaches 50 ha
  expect_equal(cn$areaBurnedHa, 100 + 50, tolerance = 0.03)
  expect_equal(cn$areaBurnedHaAllCause, 100 + 50 + 200, tolerance = 0.03)
  expect_identical(cn$nFires50, 2L)             # 100 ha and the 50 ha half
  expect_identical(cn$nFires50AllCause, 3L)
  expect_equal(sort(o$sizes$ha[o$sizes$cause == "natural"]), c(50, 100), tolerance = 0.03)
  expect_identical(sum(o$sizes$cause == "allCause"), 3L)
  expect_equal(NROW(o$polys), 2)           # the year's natural polygons, unclipped, for the map
  expect_true(all(o$polys$YEAR == 2000L & o$polys$FIRECAUS == "Natural"))
})

test_that("hindcastObserved takes the list of annual polygons the simulation holds", {
  p <- toyPolys()
  asList <- list(year2000 = p[p$YEAR == 2000L, ], year2001 = NULL, year2002 = p[p$YEAR == 2001L, ])
  a <- hindcastObserved(toyCore(), toyPoints(), asList, year = 2000, pixelAreaHa = 5.76)
  b <- hindcastObserved(toyCore(), toyPoints(), p, year = 2000, pixelAreaHa = 5.76)
  expect_equal(a$counts, b$counts)
})

test_that("hindcastObserved reports a year without polygons as zero area, and a missing cause column as NA", {
  o <- hindcastObserved(toyCore(), toyPoints(), toyPolys(), year = 1999, pixelAreaHa = 5.76)
  expect_identical(o$counts$ignitions, 0L)
  expect_identical(o$counts$areaBurnedHa, 0)
  expect_identical(NROW(o$sizes), 0L)

  p <- toyPolys()
  p$FIRECAUS <- NULL
  expect_warning(o2 <- hindcastObserved(toyCore(), toyPoints(), p, year = 2000, pixelAreaHa = 5.76), "FIRECAUS")
  expect_true(is.na(o2$counts$areaBurnedHa))
  expect_equal(o2$counts$areaBurnedHaAllCause, 350, tolerance = 0.03)
})

test_that("hindcastSimulated clips each fire to the core", {
  r <- toyCore()
  id <- terra::setValues(terra::rast(r), NA_integer_)
  id[c(1, 2, 3)] <- 1L                          # row 1, columns 1-3: in the core, 300 ha
  id[c(15, 16)] <- 2L                           # row 2, columns 5-6: one pixel in the core, one out
  id[c(20, 30)] <- 3L                           # column 10: outside the core
  ## ignitions in cells 1 (core, escaped), 2 (core), 6 (outside the core), 15 (core, escaped), 20 (outside)
  ig <- data.table::data.table(pixelID = c(1L, 2L, 6L, 15L, 20L), escaped = c(TRUE, FALSE, TRUE, TRUE, FALSE))
  s <- hindcastSimulated(id, r, ig, escapeSizeHa = 150)
  expect_identical(s$sizes$fire_id, c(1L, 2L))
  expect_equal(s$sizes$ha, c(300, 100))
  expect_identical(s$counts$ignitions, 3L)
  expect_identical(s$counts$escapes, 2L)
  expect_equal(s$counts$areaBurnedHa, 400)
  expect_identical(s$counts$nFires50, 1L)           # only the 300 ha fire reaches 150 ha in the core

  none <- hindcastSimulated(terra::setValues(terra::rast(r), NA_integer_), r)
  expect_identical(NROW(none$sizes), 0L)
  expect_identical(none$counts$areaBurnedHa, 0)
  expect_true(is.na(none$counts$ignitions))
})

test_that("hindcastSimulated places pixels by their centres when the core is on another grid", {
  r <- toyCore()                                  # 1 km cells
  fine <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 10000, ymin = 0, ymax = 10000, crs = "EPSG:3978")
  id <- terra::setValues(fine, NA_integer_)
  id[c(1, 2, 20)] <- 1L                           # row 1: columns 1-2 in the core, column 20 outside
  s <- hindcastSimulated(id, r)
  expect_equal(s$sizes$ha, 2 * 25)                # two 500 m pixels
})

test_that("the hindcast figures are built from the simulated reps and the observed year", {
  skip_if_not_installed("ggplot2")
  o <- hindcastObserved(toyCore(), toyPoints(), toyPolys(), year = 2000, pixelAreaHa = 5.76)
  set.seed(1)
  reps <- data.table::data.table(ignitions = rpois(50, 3), escapes = rpois(50, 1),
                                 areaBurnedHa = rexp(50, 1 / 150), nFires50 = rpois(50, 1))
  g <- plotHindcastDistributions(reps, o, title = "toy 2000")
  expect_s3_class(g, "ggplot")
  b <- ggplot2::ggplot_build(g)
  expect_identical(length(unique(b$data[[1]]$PANEL)), 4L)
  vl <- b$data[[2]]                               # observed lines: 4 natural + 2 all-cause
  expect_identical(NROW(vl), 6L)
  expect_true(3 %in% vl$xintercept)               # observed ignitions

  sizes <- c(60, 120, 80, 400)
  gs <- plotHindcastSizes(sizes, o)
  expect_s3_class(gs, "ggplot")
  expect_identical(NROW(ggplot2::ggplot_build(gs)$plot$data), 4L + 2L + 3L)   # >= 50 ha: 4 sim, 2 natural, 3 all

  burn <- terra::setValues(terra::rast(toyCore()), sample(0:50, 100, replace = TRUE))
  gm <- plotHindcastMap(burn, nReps = 50, firePolys = toyPolys()[1:2, ], coreRaster = toyCore(), title = "toy")
  expect_s3_class(gm, "ggplot")
  expect_identical(NROW(ggplot2::ggplot_build(gm)$data[[1]]), 50L)            # the core's 50 cells
  expect_identical(length(gm$layers), 2L)
})
