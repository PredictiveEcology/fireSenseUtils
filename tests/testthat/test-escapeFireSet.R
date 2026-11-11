## Eliot, 2026-11-10: the buffer of an ELF only removes edge effects. A fire belongs to the ELF if and only if its
## start (the polygon centroid moved to the nearest burnable pixel of the fire) is in the core; it then counts at its
## whole extent, core and buffer, in the fit, the escape model and the hindcast. One polygon is one escape.

## A 20 x 20 grid of 1 km cells (100 ha): columns 1-10 are the core (2), 11-20 the buffer (1); column 3 cannot burn.
escLarge <- function() {
  r <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20000, ymin = 0, ymax = 20000, crs = "EPSG:3978")
  terra::setValues(r, ifelse(terra::colFromCell(r, seq_len(terra::ncell(r))) <= 10, 2, 1))
}
escFlammable <- function() {
  r <- escLarge()
  terra::setValues(r, ifelse(terra::colFromCell(r, seq_len(terra::ncell(r))) == 3, 0, 1))
}
## the core alone, as rasterToMatchELF is in a hindcast: NA outside it
escCoreOnly <- function() terra::ifel(escLarge() == 2, 2, NA)

## in units of 100 m
escSq <- function(x0, y0, w, h) {
  x0 <- 100 * x0; y0 <- 100 * y0; w <- 100 * w; h <- 100 * h
  sprintf("POLYGON ((%f %f, %f %f, %f %f, %f %f, %f %f))", x0, y0, x0 + w, y0, x0 + w, y0 + h, x0, y0 + h, x0, y0)
}
## 2000: A starts in the core and burns into the buffer (24 pixels, 6 of them in the buffer); B starts in the buffer
## and burns into the core (6 core pixels); C is in the core with its centroid on the unburnable column (9 pixels,
## 6 burnable); D is human, in the core (4 pixels); E is natural and in the core but below the escape size (4
## pixels). 2001: F, as A.
escPolys <- function() {
  v <- terra::vect(c(escSq(42, 22, 80, 30), escSq(82, 122, 80, 30), escSq(12, 162, 26, 30),
                     escSq(42, 82, 20, 20), escSq(62, 172, 20, 20), escSq(42, 22, 80, 30)),
                   crs = "EPSG:3978")
  v$YEAR <- c(rep(2000L, 5), 2001L)
  v$FIRECAUS <- c("Natural", "Natural", "Natural", "Human", "Natural", "Natural")
  v$FIRE_ID <- c(1:5, 1L)
  v$POLY_HA <- terra::expanse(v, unit = "ha")
  v
}
escDef <- function() escapeDefinition(500, "natural")

test_that("escapeFireSet: one start per escape polygon, a buffer start is not a core start", {
  s <- escapeFireSet(escPolys(), escFlammable(), escLarge(), escDef(), idCol = "FIRE_ID")
  expect_identical(s$YEAR, c(2000L, 2000L, 2000L, 2001L))
  expect_identical(s$FIRE_ID, c(1L, 2L, 3L, 1L))         # not the human D nor the small E
  expect_identical(s$core, c(TRUE, FALSE, TRUE, TRUE))    # B starts in the buffer though it burns into the core
  ## A, which burns into the buffer, at its whole extent: 24 pixels, not the 18 in the core
  expect_equal(s$burnableHa[1], 2400)
  expect_equal(s$polyHa[1], 2400, tolerance = 0.001)
  ## the fit set: the core starts only
  expect_identical(s$FIRE_ID[s$core & s$YEAR == 2000L], c(1L, 3L))
})

test_that("a centroid on an unburnable pixel is moved to a burnable pixel inside its polygon", {
  p <- escPolys()
  expect_identical(terra::extract(escFlammable(), terra::centroids(p[3, ]))[[2]], 0) # on the unburnable column
  s <- escapeFireSet(p, escFlammable(), escLarge(), escDef(), idCol = "FIRE_ID")
  c3 <- s[s$FIRE_ID == 3L, ]
  expect_identical(terra::values(escFlammable(), mat = FALSE)[c3$pixelID], 1)
  expect_true(terra::relate(c3, p[3, ], "intersects")[1, 1])
  expect_equal(c3$burnableHa, 600)
})

test_that("starts given (the spread fit's) are used as they are", {
  p <- escPolys()
  starts <- terra::vect(cbind(c(5500, 15500, 3500), c(2500, 12500, 16500)), crs = "EPSG:3978")
  starts$YEAR <- 2000L
  starts$FIRE_ID <- 1:3
  s <- escapeFireSet(p, escFlammable(), escLarge(), escDef(), starts = starts, idCol = "FIRE_ID")
  expect_identical(s$FIRE_ID, 1:3)                         # F (2001) has no start here
  expect_equal(terra::crds(s), terra::crds(starts), ignore_attr = TRUE)
  expect_identical(s$core, c(TRUE, FALSE, TRUE))
  expect_error(escapeFireSet(p, escFlammable(), escLarge(), escDef(), starts = starts), "idCol")
})

## NFDB points of 2000: 1 in A (matched), 2 in B's core part (matched to a buffer start), 3 in the core in no
## polygon, 4 in the buffer in no polygon, 5 human in the core
escPoints <- function() {
  v <- terra::vect(cbind(c(7000, 9000, 3000, 18500, 5000), c(3000, 13500, 10000, 18500, 9000)), crs = "EPSG:3978")
  v$YEAR <- 2000L
  v$CAUSE <- c("N", "N", "N", "L", "H")
  v$SIZE_HA <- c(2000, 2000, 10, 10, 400)
  v
}

test_that("escape-model ignitions: core points, a matched point replaced by its fire's start, no double count", {
  p <- escPolys()
  p <- p[p$YEAR == 2000L, ]
  s <- escapeFireSet(p, escFlammable(), escLarge(), escDef(), idCol = "FIRE_ID")
  ig <- escapeIgnitionPoints(escPoints(), p, s, escLarge(), escDef(), pixSizeHa = 100)
  ## point 3 (not escaped), then the starts of A and C (escaped); not point 1, whose fire is A's start, nor point 2,
  ## whose fire starts in the buffer, nor the buffer point or the human one
  expect_identical(ig$escaped, c(FALSE, TRUE, TRUE))
  expect_equal(terra::crds(ig)[1, ], c(x = 3000, y = 10000))
  expect_equal(terra::crds(ig)[2:3, ], terra::crds(s[s$core, ]), ignore_attr = TRUE)
  expect_identical(ig$CAUSE, c("N", "N", "N"))
  expect_equal(ig$SIZE_HA[2:3], s$polyHa[s$core])
  ## one escape per polygon: A has its start and nothing else
  inA <- terra::relate(ig, p[1, ], "intersects")[, 1]
  expect_identical(sum(inA), 1L)
})

test_that("hindcastObserved(coreStart = TRUE): buffer-start fires are out, core-start fires count whole", {
  o <- hindcastObserved(escCoreOnly(), escPoints(), escPolys(), year = 2000, escapeSizeHa = 500,
                        flammableRTM = escFlammable(), fireCauses = "natural", coreStart = TRUE)
  cn <- o$counts
  expect_identical(cn$ignitions, 3L)                       # point 3, and the starts of A and C
  expect_identical(cn$escapes, 2L)                         # A and C, one each
  expect_equal(cn$areaBurnedHa, 2400 + 600 + 400)          # A whole, C, E; not B
  expect_equal(cn$areaBurnedHaAllCause, 2400 + 600 + 400 + 400) # and the human D
  expect_identical(cn$nFires50, 2L)                        # A and C reach 500 ha
  expect_equal(sort(o$sizes$ha[o$sizes$cause == "natural"]), c(400, 600, 2400))
  expect_equal(cn$areaBurnedHaPolygon, 2400 + 780 + 400, tolerance = 0.001)
  expect_identical(sum(o$points$escaped), 2L)
  ## the same data clipped to the core, as before: B's core part in, A cut to its core part (areas of the polygon
  ## pieces on burnable pixels)
  old <- hindcastObserved(escCoreOnly(), escPoints(), escPolys(), year = 2000, escapeSizeHa = 500,
                          flammableRTM = escFlammable(), fireCauses = "natural")
  expect_equal(sort(old$sizes$ha[old$sizes$cause == "natural"]), c(400, 480, 540, 1740), tolerance = 0.01)
  expect_error(hindcastObserved(escCoreOnly(), escPoints(), escPolys(), year = 2000, flammableRTM = escFlammable(),
                                coreStart = TRUE), "needs fireCauses")
})

## two simulated fires: 1 starts at column 7 (core), 6 of its 10 pixels in the buffer; 2 starts at column 13 (buffer)
## and has 2 core pixels. `ig` holds their ignitions, both escaped.
escBurn <- function() {
  burnID <- terra::setValues(terra::rast(escLarge()), NA_integer_)
  cell <- function(row, col) terra::cellFromRowCol(burnID, row, col)
  fire1 <- cell(5, 7:16)
  fire2 <- cell(10, 13:9)
  burnID[fire1] <- 1L
  burnID[fire2] <- 2L
  list(burnID = burnID, ig = data.table::data.table(pixelID = c(fire1[1], fire2[1]), escaped = TRUE))
}

test_that("hindcastSimulated(coreStart = TRUE): a buffer ignition's fire is out, a core one counts whole", {
  b <- escBurn()
  burnID <- b$burnID
  ig <- b$ig
  s <- hindcastSimulated(burnID, escCoreOnly(), ig, escapeSizeHa = 500, coreStart = TRUE)
  expect_identical(s$sizes$fire_id, 1L)
  expect_equal(s$sizes$ha, 1000)
  expect_identical(c(s$counts$ignitions, s$counts$escapes), c(1L, 1L))
  expect_equal(s$counts$areaBurnedHa, 1000)
  old <- hindcastSimulated(burnID, escCoreOnly(), ig, escapeSizeHa = 500)
  expect_equal(old$sizes$ha, c(400, 200))
  expect_error(hindcastSimulated(burnID, escCoreOnly(), coreStart = TRUE), "needs ignitions")
})

test_that("escapeDefinition is the ledger record, and checkEscapeDefinition stops on a mismatch", {
  d <- escapeDefinition()
  expect_identical(d$escapeSizeHa, 50)
  expect_identical(d$fireCauses, "natural,unknown")
  expect_identical(escapeDefinition(50, c("unknown", "natural")), d)
  expect_error(escapeDefinition(50, "lightning"), "unknown fire cause group")
  expect_error(escapeDefinition(-1), "0 or more")
  ledgerRow <- cbind(data.table::data.table(objFunVal = 1), d)
  expect_true(checkEscapeDefinition(ledgerRow, d))
  expect_error(checkEscapeDefinition(ledgerRow, escapeDefinition(100)), "differs from this run's")
  expect_error(checkEscapeDefinition(ledgerRow, escapeDefinition(50, "natural"), what = "ledger row 3"),
               "escape definition of ledger row 3")
  expect_identical(checkEscapeDefinition(data.table::data.table(objFunVal = 1), d), NA) # an older ledger row
})

## the multi-ELF rasterToMatchLargeELF of fireSense_ELFs: every cell labelled with its ELF (codes from 0), no buffer
escCategorical <- function() {
  v <- terra::vect(c(escSq(0, 0, 100, 200), escSq(100, 0, 100, 200)), crs = "EPSG:3978")
  v$ELFind <- c("4.3", "6.6.2")
  terra::rasterize(v, terra::rast(escLarge()), field = "ELFind")
}

test_that("inELFCore: the core of a one-ELF raster, of a core mask and of a multi-ELF (categorical) raster", {
  cells <- terra::cellFromRowCol(escLarge(), 1, c(1, 10, 11, 20))
  expect_identical(inELFCore(cells, escLarge()), c(TRUE, TRUE, FALSE, FALSE))
  expect_identical(inELFCore(cells, escLarge(), "mask"), rep(TRUE, 4))        # a 1 is core in a mask
  expect_identical(inELFCore(cells, escCoreOnly(), "mask"), c(TRUE, TRUE, FALSE, FALSE))
  expect_identical(inELFCore(cells, escCategorical()), rep(TRUE, 4))         # code 0 ("4.3") too
  pts <- terra::vect(terra::xyFromCell(escLarge(), cells), crs = "EPSG:3978")
  expect_identical(inELFCore(pts, escLarge()), c(TRUE, TRUE, FALSE, FALSE))
  expect_identical(inELFCore(sf::st_as_sf(pts), escLarge()), c(TRUE, TRUE, FALSE, FALSE))
  ## cells of another grid are placed by their centres: one 10 km cell per ELF zone
  coarse <- terra::aggregate(terra::rast(escLarge()), 10)
  expect_identical(inELFCore(1:4, escLarge(), grid = coarse), c(TRUE, FALSE, TRUE, FALSE))
  expect_identical(inELFCore(1:3, NULL), rep(TRUE, 3))
  expect_identical(inELFCore(integer(), escLarge()), logical())
})

test_that("a multi-ELF (categorical) raster is core everywhere, in the starts and in both hindcasts", {
  s <- escapeFireSet(escPolys(), escFlammable(), escCategorical(), escDef(), idCol = "FIRE_ID")
  expect_true(all(s$core))
  b <- escBurn()
  sim <- hindcastSimulated(b$burnID, escCategorical(), b$ig, escapeSizeHa = 500, coreStart = TRUE)
  expect_identical(sim$sizes$fire_id, 1:2)
  expect_identical(c(sim$counts$ignitions, sim$counts$escapes), c(2L, 2L))
  o <- hindcastObserved(escCategorical(), escPoints(), escPolys(), year = 2000, escapeSizeHa = 500,
                        flammableRTM = escFlammable(), fireCauses = "natural", coreStart = TRUE)
  expect_equal(o$counts$areaBurnedHa, 2400 + 2400 + 600 + 400)  # A, B (in the core now), C and E
  expect_identical(o$counts$escapes, 3L)                         # A, B and C
})

test_that("the hindcasts read a one-ELF raster (2 core, 1 buffer) as such when told, and its buffer is not core", {
  b <- escBurn()
  expect_identical(hindcastSimulated(b$burnID, escLarge(), b$ig, 500, coreStart = TRUE, coreIs = "ELF"),
                   hindcastSimulated(b$burnID, escCoreOnly(), b$ig, 500, coreStart = TRUE))
  obs <- function(core, ...) hindcastObserved(core, escPoints(), escPolys(), year = 2000, escapeSizeHa = 500,
                                              flammableRTM = escFlammable(), fireCauses = "natural", ...)
  expect_identical(obs(escLarge(), coreStart = TRUE, coreIs = "ELF")$counts, obs(escCoreOnly(), coreStart = TRUE)$counts)
  expect_identical(obs(escLarge(), coreIs = "ELF")$counts, obs(escCoreOnly())$counts)
})

test_that("checkEscapeDefinition checks one row per ELF, skips rows without a definition, names the ELF", {
  d <- escapeDefinition()
  rows <- data.table::data.table(polygonID = c("4.3", "6.6.2", "9.1"), escapeSizeHa = c(50, NA, 50),
                                 fireCauses = c(d$fireCauses, NA, d$fireCauses))
  expect_true(checkEscapeDefinition(rows, d))
  rows$escapeSizeHa[3] <- 100
  expect_error(checkEscapeDefinition(rows, d), "escape definition of the fit of ELF 9.1 \\(escapeSizeHa 100")
  expect_identical(checkEscapeDefinition(rows[2], d), NA)           # its only row has none
  ## an old fit's definition (escapeSizeHa NA) is recorded as such, and checks nothing
  old <- escapeDefinition(NA)
  expect_identical(old$escapeSizeHa, NA_real_)
  expect_identical(checkEscapeDefinition(rows, old), NA)
  expect_error(escapeDefinition("50"), "one number")
})

test_that("escapeDefinitionOf: the definition of a run's modules, which must agree", {
  mp <- withr::local_tempdir()
  for (m in c("escA", "escB")) writeLines(sprintf(c(
    'defineModule(sim, list(name = "%1$s", description = "", keywords = "", authors = person("a"),',
    '  childModules = character(0), version = list(%1$s = "0.0.1"), timeframe = as.POSIXlt(c(NA, NA)),',
    '  timeunit = "year", citation = list(), documentation = list(), reqdPkgs = list(),',
    '  parameters = rbind(defineParameter("escapeSizeHa", "numeric", 50, NA, NA, ""),',
    '                     defineParameter("fireCauses", "character", "natural", NA, NA, "")),',
    '  inputObjects = bindrows(), outputObjects = bindrows(createsOutput("def_%1$s", "data.table", ""))))',
    'doEvent.%1$s <- function(sim, eventTime, eventType) {',
    '  sim$def_%1$s <- fireSenseUtils::escapeDefinitionOf(sim)',
    '  invisible(sim)',
    '}'), m), file.path(mp, m, paste0(m, ".R")) |> (\(f) { dir.create(dirname(f)); f })())
  run <- function(a, b) suppressMessages(SpaDES.core::simInitAndSpades(
    modules = c("escA", "escB"), paths = list(modulePath = mp), times = list(start = 0, end = 0),
    params = list(escA = a, escB = b)))
  sim <- run(list(escapeSizeHa = 100), list(escapeSizeHa = 100))
  expect_identical(sim$def_escA, escapeDefinition(100, "natural"))
  expect_identical(sim$def_escB, sim$def_escA)
  expect_identical(run(list(escapeSizeHa = NA_real_), list(escapeSizeHa = NA_real_))$def_escA$escapeSizeHa, NA_real_)
  expect_error(run(list(escapeSizeHa = 100), list(escapeSizeHa = 50)), "multiple values for escapeSizeHa")
})

test_that("snappedCentroids() moves a centroid that is off its fire and keeps the `.poly` id", {
  ras <- hbpRaster()
  ## a U-shaped fire: its centroid is in the gap, which is not burnable
  uCells <- fireCells(ras, 5:10, 5:10)
  gap <- fireCells(ras, 7:8, 5:8)
  ras[gap] <- NA
  uCells <- setdiff(uCells, gap)
  r <- terra::rast(ras); r[] <- NA; r[uCells] <- 1L
  fire <- terra::as.polygons(r, dissolve = TRUE) |> terra::aggregate()
  fire$.poly <- 1L
  fire$YEAR <- 2001L
  cl <- data.table::data.table(.poly = 1L, pixelID = uCells)

  out <- snappedCentroids(fire, cl, ras)

  expect_identical(out$.poly, 1L)
  expect_true(terra::cellFromXY(ras, as.matrix(out[, c("x", "y")])) %in% uCells)
})
