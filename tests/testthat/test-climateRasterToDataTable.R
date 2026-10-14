test_that("climateRasterToDataTable merges any number of climate variables", {
  mkStack <- function(offset) {
    r <- terra::rast(nrows = 2, ncols = 2, nlyr = 2, vals = offset + 1:8)
    names(r) <- c("1985", "1986")
    r
  }
  allVars <- c("CMD", "CMDsm", "cumMDC", "MDC")
  for (n in 1:4) {
    vars <- allVars[seq_len(n)]
    rasters <- setNames(lapply(seq_len(n), function(i) mkStack(100 * i)), vars)
    out <- climateRasterToDataTable(rasters)
    expect_setequal(names(out), c("pixelID", "year", vars))
    expect_equal(nrow(out), 8L) # 4 cells x 2 years
    for (i in seq_len(n)) {
      expect_equal(sort(out[[vars[i]]]), 100 * i + 1:8)
    }
    # the variables of one pixel-year row come from the same pixel-year
    expect_true(all(out[[vars[n]]] - out[[vars[1]]] == 100 * (n - 1)))
  }
})
