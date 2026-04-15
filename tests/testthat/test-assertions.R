## Tests for assertion helpers in assertions.R

# ---------------------------------------------------------------------------
# chk_duplicatedStartPixels
# ---------------------------------------------------------------------------
test_that("chk_duplicatedStartPixels: no duplicates – returns input unchanged", {
  cells <- c(1L, 2L, 3L)
  sizes <- c(10L, 20L, 30L)
  out <- chk_duplicatedStartPixels(cells, sizes)
  expect_equal(out$loci,  cells)
  expect_equal(out$sizes, sizes)
})

test_that("chk_duplicatedStartPixels: single duplicate – keeps largest fire", {
  cells <- c(1L, 2L, 2L)
  sizes <- c(10L, 5L, 50L)   # pixel 2 appears twice; keep size=50
  out <- suppressWarnings(chk_duplicatedStartPixels(cells, sizes))
  expect_false(anyDuplicated(out$loci) > 0)
  idx2 <- which(out$loci == 2L)
  expect_equal(out$sizes[idx2], 50L)
})

test_that("chk_duplicatedStartPixels: issues a warning when duplicates exist", {
  expect_warning(
    chk_duplicatedStartPixels(c(1L, 1L), c(5L, 10L)),
    regexp = "No more than one fire"
  )
})

test_that("chk_duplicatedStartPixels: multiple duplicate pixels handled", {
  cells <- c(1L, 1L, 2L, 2L, 3L)
  sizes <- c(5L, 15L, 30L, 10L, 7L)
  out <- suppressWarnings(chk_duplicatedStartPixels(cells, sizes))
  expect_false(anyDuplicated(out$loci) > 0)
  expect_equal(out$sizes[out$loci == 1L], 15L)
  expect_equal(out$sizes[out$loci == 2L], 30L)
  expect_equal(out$sizes[out$loci == 3L],  7L)
})

test_that("chk_duplicatedStartPixels: returns a list with loci and sizes", {
  out <- chk_duplicatedStartPixels(c(1L, 2L), c(10L, 20L))
  expect_type(out, "list")
  expect_named(out, c("loci", "sizes"))
})

test_that("chk_duplicatedStartPixels: result lengths are equal", {
  cells <- c(5L, 5L, 5L, 9L)
  sizes <- c(1L, 100L, 50L, 20L)
  out <- suppressWarnings(chk_duplicatedStartPixels(cells, sizes))
  expect_equal(length(out$loci), length(out$sizes))
})

test_that("chk_duplicatedStartPixels: three-way tie keeps single largest", {
  cells <- c(4L, 4L, 4L)
  sizes <- c(3L, 9L, 6L)
  out <- suppressWarnings(chk_duplicatedStartPixels(cells, sizes))
  expect_length(out$loci, 1L)
  expect_equal(out$sizes, 9L)
})

test_that("chk_duplicatedStartPixels: size-1 vector, no duplicate", {
  out <- chk_duplicatedStartPixels(7L, 42L)
  expect_equal(out$loci,  7L)
  expect_equal(out$sizes, 42L)
})

test_that("chk_duplicatedStartPixels: order of non-duplicate entries is preserved", {
  cells <- c(10L, 20L, 30L)
  sizes <- c(1L, 2L, 3L)
  out <- chk_duplicatedStartPixels(cells, sizes)
  expect_equal(out$loci, cells)
})

test_that("chk_duplicatedStartPixels: all same pixel keeps exactly one entry", {
  cells <- rep(1L, 5)
  sizes <- c(1L, 3L, 2L, 5L, 4L)
  out <- suppressWarnings(chk_duplicatedStartPixels(cells, sizes))
  expect_length(out$loci, 1L)
  expect_equal(out$sizes, 5L)
})

test_that("chk_duplicatedStartPixels: equal sizes for duplicate – keeps one", {
  # which.max keeps first max, so first occurrence among tied winners stays
  cells <- c(1L, 1L)
  sizes <- c(10L, 10L)
  out <- suppressWarnings(chk_duplicatedStartPixels(cells, sizes))
  expect_length(out$loci, 1L)
  expect_equal(out$sizes, 10L)
})

test_that("chk_duplicatedStartPixels: non-duplicate pixels are not dropped", {
  cells <- c(1L, 2L, 2L, 3L)
  sizes <- c(5L, 10L, 20L, 15L)
  out <- suppressWarnings(chk_duplicatedStartPixels(cells, sizes))
  expect_true(1L %in% out$loci)
  expect_true(3L %in% out$loci)
})
