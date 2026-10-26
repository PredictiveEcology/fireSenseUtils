## 2026-10-25 (Eliot): "The escape model is about the transition from point to polygon." An escape is an
## ignition point that became a burned-area polygon of at least escapeSizeHa, the fires the spread model is
## fitted to, not a point whose reported size reached escapeSizeHa. One polygon is one escape.

sq <- function(x0, y0, w) sprintf("POLYGON ((%f %f, %f %f, %f %f, %f %f, %f %f))",
                                  x0, y0, x0 + w, y0, x0 + w, y0 + w, x0, y0 + w, x0, y0)

## 2000: a 100 ha polygon holding two points, a 20 ha polygon holding one, a 100 ha polygon with no point;
## 2001: a 100 ha polygon.
toyEscPolys <- function() {
  v <- terra::vect(c(sq(0, 0, 1000), sq(2000, 2000, 447.2), sq(4000, 4000, 1000), sq(6000, 0, 1000)),
                   crs = "EPSG:3978")
  v$YEAR <- c(2000L, 2000L, 2000L, 2001L)
  v
}

## 2000: two points in the big polygon (10 and 80 ha reported), a 60 ha point in the 20 ha polygon, a point
## inside the 2001 polygon, a 200 ha point in no polygon.
toyEscPoints <- function() {
  v <- terra::vect(cbind(c(500, 600, 2200, 6500, 9000), c(500, 600, 2200, 500, 9000)), crs = "EPSG:3978")
  v$YEAR <- 2000L
  v$SIZE_HA <- c(10, 80, 60, 30, 200)
  v
}

test_that("an ignition escaped if it is in a polygon of the same year of at least escapeSizeHa, one per polygon", {
  e <- escapedIgnitions(toyEscPoints(), toyEscPolys(), escapeSizeHa = 50, pixSizeHa = 6.25)
  ## the big polygon's larger point; not the point in the 20 ha polygon, nor the one in another year's
  ## polygon, nor the 200 ha point outside every polygon; the big polygon without a point adds nothing
  expect_identical(e, c(FALSE, TRUE, FALSE, FALSE, FALSE))
  ## (the big polygon without a point becomes an ignition through polygonIgnitions(), tested below)
  ## the point-size rule it replaces called points 2, 3 and 5 escapes
  expect_identical(which(toyEscPoints()$SIZE_HA >= 50), c(2L, 3L, 5L))
})

test_that("a polygon must also be larger than one pixel, and a lower threshold admits smaller polygons", {
  expect_false(any(escapedIgnitions(toyEscPoints(), toyEscPolys(), escapeSizeHa = 50, pixSizeHa = 150)))
  expect_identical(escapedIgnitions(toyEscPoints(), toyEscPolys(), escapeSizeHa = 10, pixSizeHa = 6.25),
                   c(FALSE, TRUE, TRUE, FALSE, FALSE))
})

test_that("a point is one escape even where same-year polygons overlap", {
  p <- toyEscPolys()
  both <- rbind(p, p[1, ])                     # the big polygon twice
  e <- escapedIgnitions(toyEscPoints(), both, escapeSizeHa = 50, pixSizeHa = 6.25)
  expect_identical(e, c(TRUE, TRUE, FALSE, FALSE, FALSE))  # the second copy takes the other point
})

test_that("polygons may be a list of annual SpatVectors with NULLs, sf, or in another CRS", {
  p <- toyEscPolys()
  ref <- escapedIgnitions(toyEscPoints(), p, 50, 6.25)
  asList <- list(year2000 = p[p$YEAR == 2000L, ], year1999 = NULL, year2001 = p[p$YEAR == 2001L, ])
  expect_identical(escapedIgnitions(toyEscPoints(), asList, 50, 6.25), ref)
  expect_identical(escapedIgnitions(sf::st_as_sf(toyEscPoints()), sf::st_as_sf(p), 50, 6.25), ref)
  expect_identical(escapedIgnitions(toyEscPoints(), terra::project(p, "EPSG:3347"), 50, 6.25), ref)
})

test_that("no polygons or no points give no escapes", {
  expect_identical(escapedIgnitions(toyEscPoints(), NULL, 50, 6.25), logical(5))
  expect_identical(escapedIgnitions(toyEscPoints(), list(NULL), 50, 6.25), logical(5))
  expect_identical(escapedIgnitions(toyEscPoints()[0, ], toyEscPolys(), 50, 6.25), logical(0))
})

test_that("a point just outside a polygon matches it (NFDB rounding), but one inside comes first", {
  poly <- terra::vect(sq(0, 0, 1000), crs = "EPSG:3978")       # 100 ha
  poly$YEAR <- 2000L
  near <- terra::vect(cbind(1300, 500), crs = "EPSG:3978")     # 300 m outside, 80 ha reported
  near$YEAR <- 2000L
  near$SIZE_HA <- 80
  expect_true(escapedIgnitions(near, poly, 50, 6.25))
  expect_false(escapedIgnitions(near, poly, 50, 6.25, maxDistM = 0))
  expect_false(escapedIgnitions(near, poly, 50, 6.25, maxDistM = 200))
  inside <- terra::vect(cbind(500, 500), crs = "EPSG:3978")    # inside, smaller
  inside$YEAR <- 2000L
  inside$SIZE_HA <- 5
  expect_identical(escapedIgnitions(rbind(near, inside), poly, 50, 6.25), c(FALSE, TRUE))
})

## Eliot, 2026-10-25: a natural-cause polygon of at least escapeSizeHa with no natural ignition point is added as
## an ignition that escaped, at the fire's start point (its centroid moved into the fire, as the spread fit has it).
test_that("a natural polygon with no point adds one ignition at its start point, and that ignition escapes", {
  polys <- toyEscPolys()
  polys$FIRECAUS <- c("Natural", "Natural", "Natural", "Human")
  polys$FIRE_ID <- 1:4
  ## start points as the spread fit gives them: one per polygon, inside it
  starts <- terra::vect(cbind(c(400, 2100, 4100, 6500), c(400, 2100, 4900, 500)), crs = "EPSG:3978")
  starts$YEAR <- polys$YEAR
  starts$FIRE_ID <- 1:4
  pts <- toyEscPoints()
  pts$CAUSE <- "N"
  added <- polygonIgnitions(pts, polys, starts, idCol = "FIRE_ID", escapeSizeHa = 50, pixSizeHa = 6.25)
  ## only polygon 3: polygon 1 has points, polygon 2 is 20 ha, polygon 4 (2001) is human-caused
  expect_equal(NROW(added), 1)
  expect_equal(terra::crds(added)[1, ], c(x = 4100, y = 4900))
  expect_identical(added$YEAR, 2000L)
  expect_identical(added$CAUSE, "N")
  expect_equal(added$SIZE_HA, 100, tolerance = 0.01)
  ## bound to the ignitions: one more ignition and one more escape
  all <- rbind(pts, added)
  expect_identical(escapedIgnitions(all, polys, 50, 6.25), c(FALSE, TRUE, FALSE, FALSE, FALSE, TRUE))
  ## and adding again adds nothing: the polygon is now matched
  expect_equal(NROW(polygonIgnitions(all, polys, list(year2000 = starts[1:3, ], year2001 = starts[4, ]),
                                         "FIRE_ID", 50, 6.25)), 0)
  ## a polygon without a start point (dropped by the spread fit) adds nothing
  expect_equal(NROW(polygonIgnitions(pts, polys, starts[-3, ], "FIRE_ID", 50, 6.25)), 0)
  ## no polygons: an empty result of the same shape
  none <- polygonIgnitions(pts, NULL, starts, "FIRE_ID", 50, 6.25)
  expect_s4_class(none, "SpatVector")
  expect_equal(NROW(none), 0)
  expect_true(all(c("YEAR", "CAUSE", "SIZE_HA") %in% names(none)))
})
