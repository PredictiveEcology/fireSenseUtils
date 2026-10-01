## The edge ring of a fire's buffer: buffer pixels with a queen neighbour outside THAT fire's own
## buffer. A simulated fire that burns one has reached the outer edge of where it may spread.

r10 <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10)
cellOf <- function(row, col) (row - 1L) * 10L + col

block <- function(rows, cols) as.vector(outer(rows, cols, cellOf))

test_that("interior pixels are not edge; pixels next to a hole, another fire's buffer or a removed pixel are", {
  ## fire 1: rows 2-8 x cols 2-8 with a hole at (5,5), which is fire 2's one-pixel buffer, and a
  ## hole at (7,7), a pixel removed as non-flammable (no row in the table at all)
  f1 <- setdiff(block(2:8, 2:8), c(cellOf(5, 5), cellOf(7, 7)))
  dt <- data.table::data.table(ids = c(rep(1L, length(f1)), 2L), pixelID = c(f1, cellOf(5, 5)),
                               buffer = 0L)
  dt[pixelID == cellOf(4, 4), buffer := 1L]
  out <- bufferEdge(dt, r10)
  expect_type(out, "logical")
  expect_length(out, nrow(dt))

  interior <- setdiff(block(3:7, 3:7), c(block(4:6, 4:6), block(6:7, 6:7)))
  expect_length(interior, 13L)
  expect_setequal(dt$pixelID[dt$ids == 1L & !out], interior)
  ## the outer ring of the block, and the pixels around each hole, are edge
  expect_true(all(out[dt$ids == 1L & dt$pixelID %in% block(2:8, c(2, 8))]))
  expect_true(out[dt$ids == 1L & dt$pixelID == cellOf(4, 4)])  # next to fire 2's pixel (and in the union's interior)
  expect_true(out[dt$ids == 1L & dt$pixelID == cellOf(6, 7)])  # next to the removed pixel
  ## a fire of one pixel is all edge
  expect_true(out[dt$ids == 2L])
})

test_that("the raster boundary counts as outside", {
  dt <- data.table::data.table(ids = 1L, pixelID = 1:100, buffer = 0L)
  out <- bufferEdge(dt, r10)
  expect_equal(sum(out), 36L)                                  # the 10 x 10 frame
  expect_false(any(out[dt$pixelID %in% block(2:9, 2:9)]))
  expect_true(all(out[dt$pixelID %in% c(block(1, 1:10), block(10, 1:10), block(1:10, 1), block(1:10, 10))]))
})

test_that("an edge is relative to the pixel's own fire: overlapping buffers, character ids", {
  ## fire "a" is rows 2-6 x cols 2-6; fire "b" is rows 3-5 x cols 3-5, wholly inside it, so b's
  ## pixels are in the union of a's buffer but at b's own edge; a's pixels around b are interior to a
  dt <- data.table::data.table(ids = c(rep("a", 25), rep("b", 9)),
                               pixelID = c(block(2:6, 2:6), block(3:5, 3:5)), buffer = 0L)
  out <- bufferEdge(dt, r10)
  aInterior <- dt$pixelID[dt$ids == "a" & !out]
  expect_setequal(aInterior, block(3:5, 3:5))
  bInterior <- dt$pixelID[dt$ids == "b" & !out]
  expect_equal(bInterior, cellOf(4, 4))
})

test_that("addBufferEdge() adds the column to each year's table, copies, and leaves a made one", {
  dt <- data.table::data.table(ids = 1L, pixelID = 1:100, buffer = 0L)
  lst <- list(year2001 = dt, year2002 = data.table::copy(dt))
  made <- addBufferEdge(lst, r10)
  expect_true(all(vapply(made, function(x) "edge" %in% names(x), logical(1))))
  expect_false("edge" %in% names(lst$year2001))
  expect_identical(made$year2001$edge, bufferEdge(dt, r10))
  again <- addBufferEdge(made, r10)
  expect_identical(again, made)
  expect_null(addBufferEdge(NULL, r10))
})
