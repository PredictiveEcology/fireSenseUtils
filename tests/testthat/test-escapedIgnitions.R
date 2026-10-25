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
