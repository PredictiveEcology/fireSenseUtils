## Tests for removeBufferedFiresOutsideRTM

library(data.table)

test_that("removeBufferedFiresOutsideRTM: drops only pixels outside flammable RTM", {
  ## flammableRTM with NA at pixelID 4; pixelIDs 1:3 are flammable
  flammableRTM <- c(1, 1, 1, NA, 1)

  ## Two fires: pixels 1-3 entirely flammable; pixel 4 of fire 'B' is NA
  dt <- data.table(
    pixelID = c(1L, 2L, 3L, 4L, 5L),
    ids     = c("A", "A", "B", "B", "B"),
    buffer  = c(1L, 0L, 1L, 1L, 0L)
  )
  out <- removeBufferedFiresOutsideRTM(dt, flammableRTM)

  ## Only the NA-pixel row is dropped; the rest of fire B remains
  expect_false(4L %in% out$pixelID)
  expect_true(all(c(1L, 2L, 3L, 5L) %in% out$pixelID))
})

test_that("removeBufferedFiresOutsideRTM: internal 'flammable' column is stripped from output", {
  flammableRTM <- c(1, 1, 1)
  dt <- data.table(pixelID = 1:3, ids = "X", buffer = 1L)
  out <- removeBufferedFiresOutsideRTM(dt, flammableRTM)
  expect_false("flammable" %in% colnames(out))
  expect_setequal(colnames(out), c("pixelID", "ids", "buffer"))
})

test_that("removeBufferedFiresOutsideRTM: returns a data.table", {
  flammableRTM <- c(1, NA)
  dt <- data.table(pixelID = 1:2, ids = "X", buffer = 1L)
  out <- removeBufferedFiresOutsideRTM(dt, flammableRTM)
  expect_true(is.data.table(out))
})

test_that("removeBufferedFiresOutsideRTM: empty result when all pixels are NA", {
  flammableRTM <- c(NA, NA, NA)
  dt <- data.table(pixelID = 1:3, ids = "X", buffer = 1L)
  out <- removeBufferedFiresOutsideRTM(dt, flammableRTM)
  expect_equal(nrow(out), 0)
})
