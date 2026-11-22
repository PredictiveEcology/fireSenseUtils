## mosaicELFs(): per-ELF rasters (core plus buffer, each in its own projection) blended into one national grid, each ELF
## weighted by its distance inward from the outer edge of its footprint.

## Two overlapping ELFs on the toy 5 km grid (4 rows x 12 columns, EPSG:3978): A covers columns 1-8, B columns 5-12, so
## they overlap in columns 5-8 (20 km). A's raster has the value `a` in EPSG:3978, B's has `b` in UTM 15N.
overlapELFs <- function() {
  foot <- function(cols) {
    r <- toyGrid()
    col <- terra::colFromCell(r, seq_len(terra::ncell(r)))
    r[col %in% cols] <- 1
    r
  }
  rasWhole <- list(A = foot(1:8), B = foot(5:12))
  names(rasWhole$A) <- "A"; names(rasWhole$B) <- "B"
  list(rasWhole = rasWhole)
}

## a raster of `value` over an ELF's footprint (extent of its columns), in `crs`, at 240 m, with NA `holes` (extents)
elfRaster <- function(cols, value, crs, holes = list()) {
  e <- terra::ext((min(cols) - 1) * 5000, max(cols) * 5000, 0, 20000)
  r <- terra::rast(terra::ext(terra::project(terra::as.polygons(e, crs = "EPSG:3978"), crs)), resolution = 240, crs = crs)
  r[] <- value
  for (h in holes) {
    hole <- terra::project(terra::as.polygons(h, crs = "EPSG:3978"), crs)
    r <- terra::mask(r, hole, inverse = TRUE)
  }
  f <- tempfile(fileext = ".tif")
  terra::writeRaster(r, f, overwrite = TRUE)
  f
}

mosaicToy <- function(holes = list(), file = tempfile(fileext = ".tif")) {
  files <- list(sim = c(A = elfRaster(1:8, 10, "EPSG:3978", holes), B = elfRaster(5:12, 20, "EPSG:32615")))
  out <- mosaicELFs(files, overlapELFs(), filename = c(sim = file), res = 250)
  list(file = out[["sim"]], r = terra::rast(out[["sim"]]))
}

## the values along the middle of the map, west to east, at cell centres
rowValues <- function(r, y = 10000) {
  x <- seq(2000, 58000, by = 250)
  stats::setNames(terra::extract(r, cbind(x, y))[[1]], x)
}

test_that("mosaicELFs() keeps each ELF's value where only it reaches, and ramps between them across the overlap", {
  r <- mosaicToy()$r
  v <- rowValues(r)
  expect_equal(unname(v[as.character(c(2000, 5000, 10000))]), c(10, 10, 10))        # only A (west of column 5)
  expect_equal(unname(v[as.character(c(45000, 50000, 58000))]), c(20, 20, 20))      # only B (east of column 8)
  inOverlap <- v[as.numeric(names(v)) > 20000 & as.numeric(names(v)) < 40000]
  expect_false(anyNA(inOverlap))
  expect_true(all(diff(inOverlap) >= -1e-6))                                        # monotone: A's value to B's
  expect_gt(max(inOverlap), 19); expect_lt(min(inOverlap), 11)                      # reaches both ends
  expect_equal(v[["30000"]], 15, tolerance = 0.05)                                  # halfway by symmetry
  expect_true(any(inOverlap > 11 & inOverlap < 19))                                 # a ramp, not a step
})

test_that("mosaicELFs() writes a Float32 cloud-optimized GeoTIFF in EPSG:3978 with the cell size asked for", {
  m <- mosaicToy()
  expect_true(any(grepl("LAYOUT=COG", terra::describe(m$file))))
  expect_true(any(grepl("Float32", terra::describe(m$file))))
  expect_equal(terra::res(m$r), c(250, 250))
  expect_equal(terra::crs(m$r, describe = TRUE)$code, "3978")
})

test_that("an NA hole inside one ELF changes no weight: values elsewhere are as without the hole, the hole takes the other ELF's value", {
  plain <- mosaicToy()$r
  hole <- terra::ext(27500, 32500, 5000, 15000)                 # in the overlap, inside A's footprint
  holed <- mosaicToy(holes = list(hole))$r
  expect_equal(rowValues(holed)[["30000"]], 20)                  # only B has a value in the hole
  far <- terra::ext(0, 24000, 0, 20000)                          # A only, no hole: unchanged
  expect_equal(terra::values(terra::crop(holed, far)), terra::values(terra::crop(plain, far)), tolerance = 1e-5)
  ## two cells away from the hole, the blend is the one without the hole
  x <- c(24000, 36000)
  expect_equal(terra::extract(holed, cbind(x, 10000))[[1]], terra::extract(plain, cbind(x, 10000))[[1]], tolerance = 1e-3)
  ## a hole where no other ELF reaches stays NA
  hole2 <- terra::ext(2000, 6000, 5000, 15000)
  holed2 <- mosaicToy(holes = list(hole2))$r
  expect_true(is.na(terra::extract(holed2, cbind(4000, 10000))[[1]]))
})

test_that("mosaicELFs() blends every product with the same weights", {
  files <- list(sim = c(A = elfRaster(1:8, 10, "EPSG:3978"), B = elfRaster(5:12, 20, "EPSG:32615")),
                hist = c(A = elfRaster(1:8, 1, "EPSG:3978"), B = elfRaster(5:12, 3, "EPSG:32615")))
  out <- mosaicELFs(files, overlapELFs(), filename = c(sim = tempfile(fileext = ".tif"), hist = tempfile(fileext = ".tif")))
  s <- rowValues(terra::rast(out[["sim"]])); h <- rowValues(terra::rast(out[["hist"]]))
  inOverlap <- as.numeric(names(s)) > 20000 & as.numeric(names(s)) < 40000
  ## sim = 10 + 10 t and hist = 1 + 2 t for one blending share t
  expect_equal(unname(h[inOverlap]), unname(1 + 2 * (s[inOverlap] - 10) / 10), tolerance = 1e-4)
})
