## Tests for pooledFireSizes() in objFunSpread.R
##
## Regression: the adTest in objFunSpread() used to run inside the year-batch loop
## under `ii == 2`, pooling only the final batch's *simulated* sizes while comparing
## them against *every* year's observed sizes. The omitted years are the two with the
## largest area burned, so the observed sample kept an upper tail the simulated sample
## structurally could not have.

historicalFiresAboveMin <- list(
  year2001 = data.frame(cells = 1:2, size = c(500, 900)), ## a "large" year
  year2002 = data.frame(cells = 3:4, size = c(10, 20)),
  year2003 = data.frame(cells = 5:6, size = c(30, 40))
)

test_that("pooledFireSizes: both samples cover every year that ran", {
  out <- pooledFireSizes(
    fireSizesList = list(c(480, 870), c(12, 18, 33, 44)),
    yrsDoneList = list("year2001", c("year2002", "year2003")),
    historicalFiresAboveMin = historicalFiresAboveMin
  )
  ## the large-fire batch must appear on BOTH sides, not just the observed one
  expect_true(all(c(480, 870) %in% out$simulated))
  expect_true(all(c(500, 900) %in% out$observed))
  expect_equal(length(out$simulated), length(out$observed))
})

test_that("pooledFireSizes: observed is restricted to the years that ran", {
  ## only the small-fire batch ran; year2001 must be absent from BOTH samples
  out <- pooledFireSizes(
    fireSizesList = list(NULL, c(12, 18, 33, 44)),
    yrsDoneList = list(NULL, c("year2002", "year2003")),
    historicalFiresAboveMin = historicalFiresAboveMin
  )
  expect_false(any(c(500, 900) %in% out$observed))
  expect_equal(length(out$simulated), length(out$observed))
})

test_that("pooledFireSizes: no year contributes to only one side", {
  ## the invariant that broke -- an unbalanced batch set must not silently pass
  for (yrsDone in list("year2001", c("year2002", "year2003"),
                       c("year2001", "year2002", "year2003"))) {
    sizes <- lapply(historicalFiresAboveMin[yrsDone], function(x) x$size)
    out <- pooledFireSizes(
      fireSizesList = list(unlist(sizes)),
      yrsDoneList = list(yrsDone),
      historicalFiresAboveMin = historicalFiresAboveMin
    )
    expect_equal(length(out$simulated), length(out$observed))
  }
})

test_that("pooledFireSizes: empty accumulators give empty samples", {
  out <- pooledFireSizes(list(), list(), historicalFiresAboveMin)
  expect_length(out$simulated, 0)
  expect_length(out$observed, 0)
})
