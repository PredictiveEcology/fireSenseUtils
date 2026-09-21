## cropToCells() gives SpaDES.tools::spread() a smaller landscape: spread() allocates
## landscape-length state on every call, and objFunInner() calls it Nreps times per fire year on a
## landscape where that year's pixels fill a small part (10-33% by bounding box on ELFs 5.3.1,
## 5.3.2 and 13.1). Measured with identical seeds: identical objective values, and an evaluation
## 1.6-2.7x faster on ELF 5.3.1.
## What could go wrong is a cell index landing on the wrong cell, or the crop changing what
## spread() draws. Both are tested against terra and spread() themselves, not against the mapping.

r <- terra::rast(nrows = 30, ncols = 40, xmin = 1000, xmax = 1400, ymin = 500, ymax = 800)
## a block of rows 10-15, columns 12-20
block <- as.integer(outer(12:20, (10:15 - 1L) * 40L, `+`))

test_that("toCrop() lands on the same place on the ground, by terra's own coordinates", {
  cr <- cropToCells(r, block)
  expect_equal(terra::xyFromCell(cr$r, cr$toCrop(block)), terra::xyFromCell(r, block))
  expect_equal(terra::res(cr$r), terra::res(r))
})

test_that("toFull() undoes toCrop(), for every cell of the crop", {
  cr <- cropToCells(r, block)
  all <- seq_len(cr$ncell)
  expect_identical(cr$toCrop(cr$toFull(all)), all)
  expect_equal(terra::xyFromCell(r, cr$toFull(all)), terra::xyFromCell(cr$r, all))
})

test_that("the crop is the bounding box plus one cell each side", {
  cr <- cropToCells(r, block)
  expect_identical(c(terra::nrow(cr$r), terra::ncol(cr$r)), c(6 + 2, 9 + 2))
  expect_identical(cr$ncell, as.integer(terra::ncell(cr$r)))
})

test_that("the margin stops at the landscape's edge", {
  cr <- cropToCells(r, c(1L, 42L)) # top-left corner
  expect_identical(c(terra::nrow(cr$r), terra::ncol(cr$r)), c(3, 3))
  expect_identical(cr$toCrop(1L), 1L)
  cr <- cropToCells(r, as.integer(terra::ncell(r))) # bottom-right corner
  expect_identical(cr$toFull(cr$ncell), as.integer(terra::ncell(r)))
})

test_that("both mappings keep cells in order", {
  cr <- cropToCells(r, block)
  expect_false(is.unsorted(cr$toCrop(sort(block)), strictly = TRUE))
  expect_false(is.unsorted(cr$toFull(seq_len(cr$ncell)), strictly = TRUE))
})

test_that("spread() on the crop burns the same cells as on the landscape, draw for draw", {
  skip_if_not_installed("SpaDES.tools")
  big <- terra::rast(nrows = 200, ncols = 300, xmin = 0, xmax = 300, ymin = 0, ymax = 200)
  ## burnable cells: two separate blocks, so fires run into the edge of what can burn
  pix <- as.integer(c(outer(50:120, (40:90 - 1L) * 300L, `+`), outer(200:260, (100:160 - 1L) * 300L, `+`)))
  set.seed(11); sp <- runif(length(pix), 0.18, 0.30)
  loci <- as.integer(c(85 + (65 - 1L) * 300L, 230 + (130 - 1L) * 300L))
  full <- numeric(terra::ncell(big)); full[pix] <- sp; full[loci] <- 1
  cr <- cropToCells(big, c(pix, loci))
  small <- numeric(cr$ncell); small[cr$toCrop(pix)] <- sp; small[cr$toCrop(loci)] <- 1
  expect_lt(cr$ncell, terra::ncell(big) / 2)

  for (seed in 1:3) {
    set.seed(seed)
    onFull <- SpaDES.tools::spread(big, loci = loci, spreadProb = full, maxSize = c(2000, 2000),
                                   returnIndices = TRUE, allowOverlap = FALSE, quick = TRUE)
    set.seed(seed)
    onCrop <- SpaDES.tools::spread(cr$r, loci = cr$toCrop(loci), spreadProb = small, maxSize = c(2000, 2000),
                                   returnIndices = TRUE, allowOverlap = FALSE, quick = TRUE)
    expect_gt(nrow(onFull), 50) # the fires did spread
    expect_identical(cr$toFull(onCrop$indices), onFull$indices)
    expect_identical(cr$toFull(onCrop$initialLocus), onFull$initialLocus)
  }
})
