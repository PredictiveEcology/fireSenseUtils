## bufferToArea() grows each buffer from its newest ring only. It must give what the whole-buffer
## spread2() version (helper-bufferToAreaSpread2.R) gave: the same size for every fire and the same
## cells in every complete ring. Only the random draws differ: for the cells of the last, partial
## ring, and for a cell that two fires reach in the same iteration.

squareFire <- function(col, row, half, id, n) { # half-width `half` cells; 1 km cells, n x n raster
  x <- (col - 0.5) * 1000; y <- (n - row + 0.5) * 1000; h <- half * 1000 # row 1 is the top row
  sf::st_sf(FIRE_ID = id, geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(x - h, y - h), c(x + h, y - h), c(x + h, y + h), c(x - h, y + h), c(x - h, y - h)))), crs = 3978))
}
gridRTM <- function(n) terra::rast(nrows = n, ncols = n, extent = c(0, n * 1000, 0, n * 1000),
                                   crs = "EPSG:3978", vals = 1L)
cellsOf <- function(b, id) sort(b$pixelID[b$ids == id])
## cells within `k` rings (Chebyshev distance) of a one-cell fire at (col, row)
squareCells <- function(rtm, col, row, k) {
  rc <- expand.grid(r = (row - k):(row + k), c = (col - k):(col + k))
  rc <- rc[rc$r >= 1 & rc$r <= terra::nrow(rtm) & rc$c >= 1 & rc$c <= terra::ncol(rtm), ]
  sort(terra::cellFromRowCol(rtm, rc$r, rc$c))
}
bothBuffers <- function(polys, rtm, am, minSize = 1, seed = 1) {
  set.seed(seed); old <- bufferToAreaSpread2(polys, rtm, areaMultiplier = am, field = "FIRE_ID", minSize = minSize)
  set.seed(seed); new <- bufferToArea(polys, rtm, areaMultiplier = am, field = "FIRE_ID", minSize = minSize)
  list(old = old, new = new)
}

test_that("a fire whose target is a complete ring gets the same cells as before", {
  rtm <- gridRTM(40)
  ## one-cell fires at (10, 10) and (30, 30); target 25 = 1 + 8 + 16 cells: rings 0-2, a 5 x 5 square
  polys <- rbind(squareFire(10, 10, 0.4, 1L, 40), squareFire(30, 30, 0.4, 2L, 40))
  b <- bothBuffers(polys, rtm, function(x) 25 * x)
  for (id in 1:2) {
    expect_equal(cellsOf(b$new, id), cellsOf(b$old, id))
  }
  expect_equal(cellsOf(b$new, 1L), squareCells(rtm, 10, 10, 2))
  expect_equal(cellsOf(b$new, 2L), squareCells(rtm, 30, 30, 2))
  expect_equal(b$new[ids == 1L & buffer == 1L]$pixelID, terra::cellFromRowCol(rtm, 10, 10))
})

test_that("a fire whose target ends in a partial ring has the same size and the same inner rings", {
  rtm <- gridRTM(40)
  polys <- rbind(squareFire(10, 10, 1.4, 1L, 40), squareFire(30, 28, 0.4, 2L, 40)) # 3 x 3 cells and 1 cell
  b <- bothBuffers(polys, rtm, function(x) 10 * x + 5) # targets 95 and 15: partial rings
  for (id in 1:2) {
    expect_identical(nrow(b$new[ids == id]), nrow(b$old[ids == id]))
  }
  expect_equal(nrow(b$new[ids == 1L]), 95L)
  expect_equal(nrow(b$new[ids == 2L]), 15L)
  ## fire 1: rings up to 9 x 9 = 81 cells are whole and the last 14 cells come from the 11 x 11 ring;
  ## fire 2: 3 x 3 = 9 cells are whole and 6 come from the 5 x 5 ring
  for (v in list(b$new, b$old)) {
    expect_true(all(squareCells(rtm, 10, 10, 4) %in% cellsOf(v, 1L)))
    expect_true(all(squareCells(rtm, 30, 28, 1) %in% cellsOf(v, 2L)))
    expect_true(all(cellsOf(v, 1L) %in% squareCells(rtm, 10, 10, 5)))
    expect_true(all(cellsOf(v, 2L) %in% squareCells(rtm, 30, 28, 2)))
  }
  expect_equal(anyDuplicated(b$new$pixelID), 0)
})

test_that("a fire that fills the landscape before its target gets the same cells as before", {
  rtm <- gridRTM(12)
  polys <- squareFire(3, 4, 0.4, 1L, 12)
  b <- bothBuffers(polys, rtm, function(x) 1000 * x) ## target 1000 > 144 cells
  expect_equal(cellsOf(b$new, 1L), seq_len(terra::ncell(rtm)))
  expect_equal(cellsOf(b$new, 1L), cellsOf(b$old, 1L))
})

test_that("a fire that finishes frees its cells for the fire beside it", {
  rtm <- gridRTM(12)
  ## fire 2 (one cell) reaches its target of 9 within two iterations; fire 1 (3 x 3 cells, target 1000)
  ## has by then grown against it, then grows over fire 2's cells to fill the whole raster
  polys <- rbind(squareFire(3, 6, 1.4, 1L, 12), squareFire(6, 6, 0.4, 2L, 12))
  b <- bothBuffers(polys, rtm, function(x) if (x == 1) 9 else 1000)
  expect_equal(nrow(b$new[ids == 2L]), 9L)
  ## a pixel is kept once, for the fire that finished first. Without the freed cells fire 1 never
  ## reaches the 16 cells fire 2 did not keep, and the raster is not full
  expect_equal(nrow(b$new), terra::ncell(rtm))
  expect_equal(nrow(b$old), terra::ncell(rtm))
  expect_equal(nrow(b$new[ids == 1L]), nrow(b$old[ids == 1L]))
})

test_that("many fires: fires out of each other's reach have the cells they had, all fires keep their burned cells", {
  n <- 200
  rtm <- gridRTM(n)
  set.seed(3)
  fires <- data.frame(id = 1:14, col = sample(10:190, 14), row = sample(10:190, 14),
                      half = sample(c(0.4, 1.4, 2.4), 14, replace = TRUE))
  polys <- do.call(rbind, lapply(seq_len(nrow(fires)), function(i)
    squareFire(fires$col[i], fires$row[i], fires$half[i], fires$id[i], n)))
  b <- bothBuffers(polys, rtm, function(x) 10 * x)
  expect_equal(anyDuplicated(b$new$pixelID), 0)
  expect_setequal(b$new$ids, fires$id)
  burnedIn <- terra::rasterize(polys, rtm, field = "FIRE_ID")
  expect_setequal(b$new[buffer == 1L]$pixelID, which(!is.na(terra::values(burnedIn, mat = FALSE))))
  ## a fire with no neighbour within 25 cells (its buffer reaches at most 8) cannot be contested:
  ## same size, and the same cells in every complete ring (ring = Chebyshev distance from its block)
  alone <- fires$id[vapply(seq_len(nrow(fires)), function(i)
    all(pmax(abs(fires$col[-i] - fires$col[i]), abs(fires$row[-i] - fires$row[i])) > 25), logical(1))]
  expect_gt(length(alone), 3)
  for (id in alone) {
    f <- fires[fires$id == id, ]
    ringOf <- function(cells) {
      rc <- terra::rowColFromCell(rtm, cells)
      pmax(0, abs(rc[, 1] - f$row) - floor(f$half), abs(rc[, 2] - f$col) - floor(f$half))
    }
    kNew <- cellsOf(b$new, id); kOld <- cellsOf(b$old, id)
    expect_length(kNew, length(kOld))
    outer <- max(ringOf(kNew))
    expect_equal(max(ringOf(kOld)), outer)
    expect_equal(kNew[ringOf(kNew) < outer], kOld[ringOf(kOld) < outer])
  }
})

test_that("a large fire's buffer is built from its frontier, not the whole buffer, each iteration", {
  skip_on_cran()
  rtm <- gridRTM(1000)
  fire <- sf::st_buffer(sf::st_sf(FIRE_ID = 1L, geometry = sf::st_sfc(sf::st_point(c(500000, 500000)), crs = 3978)), 30000)
  am <- function(x) 10 * x
  tNew <- system.time(new <- bufferToArea(fire, rtm, areaMultiplier = am, field = "FIRE_ID", minSize = 1))[["elapsed"]]
  tOld <- system.time(old <- bufferToAreaSpread2(fire, rtm, areaMultiplier = am, field = "FIRE_ID", minSize = 1))[["elapsed"]]
  expect_equal(nrow(new), nrow(old))
  expect_lt(tNew, tOld / 5)
})
