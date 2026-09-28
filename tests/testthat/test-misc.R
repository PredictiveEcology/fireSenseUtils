## Tests for makeMutuallyExclusive and paramsSeparate edge cases

library(data.table)

# ---------------------------------------------------------------------------
# makeMutuallyExclusive
# ---------------------------------------------------------------------------
test_that("makeMutuallyExclusive: zeroes matched columns where cov1 is non-zero", {
  dt <- data.table(youngAge = c(0, 5, 0, 3),
                   vegPC1   = c(10, 20, 30, 40),
                   vegPC2   = c(1,  2,  3,  4))
  out <- makeMutuallyExclusive(dt)
  # rows 2 and 4 have youngAge != 0
  expect_equal(out$vegPC1[c(2, 4)], c(0, 0))
  expect_equal(out$vegPC2[c(2, 4)], c(0, 0))
  # rows 1 and 3 are untouched
  expect_equal(out$vegPC1[c(1, 3)], c(10, 30))
})

test_that("makeMutuallyExclusive: leaves rows alone where cov1 is zero", {
  dt <- data.table(youngAge = c(0, 0), vegPC1 = c(5, 10))
  out <- makeMutuallyExclusive(dt)
  expect_equal(out$vegPC1, c(5, 10))
})

test_that("makeMutuallyExclusive: custom mutuallyExclusiveCols argument", {
  dt <- data.table(fire = c(0, 1, 1, 0),
                   veg_a = c(3, 3, 3, 3),
                   veg_b = c(1, 1, 1, 1),
                   other = c(9, 9, 9, 9))
  out <- makeMutuallyExclusive(dt, mutuallyExclusiveCols = list("fire" = c("veg_")))
  expect_equal(out$veg_a[c(2, 3)], c(0, 0))
  expect_equal(out$veg_b[c(2, 3)], c(0, 0))
  expect_equal(out$other, c(9, 9, 9, 9))  # not matched
})

test_that("makeMutuallyExclusive: no columns match grep – no change", {
  dt <- data.table(youngAge = c(1, 2), colA = c(5, 6))
  out <- makeMutuallyExclusive(dt, mutuallyExclusiveCols = list("youngAge" = c("ZZZNOMATCH")))
  expect_equal(out$colA, c(5, 6))
})

test_that("makeMutuallyExclusive: returns a data.table", {
  dt <- data.table(youngAge = c(0, 1), vegPC1 = c(3, 4))
  out <- makeMutuallyExclusive(dt)
  expect_true(is.data.table(out))
})

test_that("makeMutuallyExclusive: modifies in place (same object)", {
  dt <- data.table(youngAge = c(1), vegPC1 = c(7))
  out <- makeMutuallyExclusive(dt)
  expect_true(identical(dt, out))
})

test_that("makeMutuallyExclusive: multiple grep patterns for one cov", {
  dt <- data.table(youngAge = c(0, 2),
                   vegPC1   = c(10, 10),
                   bio1     = c(5, 5))
  out <- makeMutuallyExclusive(dt,
    mutuallyExclusiveCols = list("youngAge" = c("vegPC", "bio")))
  expect_equal(out$vegPC1[2], 0)
  expect_equal(out$bio1[2],   0)
})

test_that("makeMutuallyExclusive: all cov1 zero – nothing changed", {
  dt <- data.table(youngAge = c(0, 0, 0), vegPC1 = c(1, 2, 3))
  out <- makeMutuallyExclusive(dt)
  expect_equal(out$vegPC1, c(1, 2, 3))
})

# ---------------------------------------------------------------------------
# makeMutuallyExclusive: youngAge is exclusive with everything, not zeroed itself
#
# fireSense_SpreadFit::spreadFitPrep() appends every non-annual column name to youngAge's own
# pattern list, so when youngAge itself is a non-annual column, one of those patterns is
# "youngAge" (see fireSense_SpreadFit's own tests). Column order (youngAge before or after the
# other covariates in the pattern list) must not matter.
# ---------------------------------------------------------------------------
test_that("makeMutuallyExclusive: youngAge stays 1 and is not zeroed by its own pattern (youngAge first)", {
  dt <- data.table(youngAge = c(1, 0, 1), nfLCC_40 = c(1, 1, 0), nfLCC_50 = c(0, 1, 1))
  out <- makeMutuallyExclusive(dt,
    mutuallyExclusiveCols = list(youngAge = c("youngAge", "nfLCC_40", "nfLCC_50")))
  expect_equal(out$youngAge, c(1, 0, 1))     # never zeroed by its own pattern
  expect_equal(out$nfLCC_40, c(0, 1, 0))     # zeroed on young rows only
  expect_equal(out$nfLCC_50, c(0, 1, 0))
})

test_that("makeMutuallyExclusive: youngAge stays 1 regardless of pattern order (youngAge last)", {
  dt <- data.table(nfLCC_40 = c(1, 1, 0), nfLCC_50 = c(0, 1, 1), youngAge = c(1, 0, 1))
  out <- makeMutuallyExclusive(dt,
    mutuallyExclusiveCols = list(youngAge = c("nfLCC_40", "nfLCC_50", "youngAge")))
  expect_equal(out$youngAge, c(1, 0, 1))
  expect_equal(out$nfLCC_40, c(0, 1, 0))
  expect_equal(out$nfLCC_50, c(0, 1, 0))
})

test_that("makeMutuallyExclusive: an earlier pattern zeroing a column does not blind a later pattern", {
  ## if whToZero were recomputed from dt[[cov1]] after an earlier pattern zeroed cov1 (the old
  ## bug), a pattern's own column would go quiet and later patterns would see no rows to zero
  dt <- data.table(youngAge = c(1, 0), youngAgeAlias = c(1, 0), nfLCC_40 = c(1, 1))
  out <- makeMutuallyExclusive(dt,
    mutuallyExclusiveCols = list(youngAge = c("youngAgeAlias", "nfLCC_40")))
  expect_equal(out$youngAge, c(1, 0))
  expect_equal(out$youngAgeAlias, c(0, 0))
  expect_equal(out$nfLCC_40, c(0, 1))  # still zeroed on the young row
})

# ---------------------------------------------------------------------------
# youngAgeExclusiveCols
# ---------------------------------------------------------------------------
test_that("youngAgeExclusiveCols: matches nfLCC_*, treedWetland and supplied fuel columns, not youngAge", {
  covNames <- c("youngAge", "BlkSprc", "nfLCC_40", "nfLCC_50_80", "treedWetland", "CMDsm")
  out <- youngAgeExclusiveCols(covNames, fuelCols = "BlkSprc")
  expect_named(out, "youngAge")
  expect_setequal(out$youngAge, c("BlkSprc", "nfLCC_40", "nfLCC_50_80", "treedWetland"))
  expect_false("CMDsm" %in% out$youngAge)     # climate is left alone
  expect_false("youngAge" %in% out$youngAge)
})

# ---------------------------------------------------------------------------
# paramsSeparate edge cases
# ---------------------------------------------------------------------------
test_that("paramsSeparate: parsModel equal to length gives all covPars", {
  par <- c(1.1, 2.2, 3.3)
  res <- paramsSeparate(par, parsModel = 3)
  expect_equal(res$covPars, par)
  expect_length(res$logisticPars, 0)
})

test_that("paramsSeparate: parsModel = 1 splits cleanly", {
  par <- c(0.5, 0.3, 0.1)
  res <- paramsSeparate(par, parsModel = 1)
  expect_equal(res$covPars,      0.1)
  expect_equal(res$logisticPars, c(0.5, 0.3))
})
