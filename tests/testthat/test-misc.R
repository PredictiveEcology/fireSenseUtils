## Tests for makeMutuallyExclusive, extractSpecial, and pw/oom edge cases

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
# extractSpecial
# ---------------------------------------------------------------------------
test_that("extractSpecial: returns list with variable and knot", {
  out <- extractSpecial(x, 5)
  expect_type(out, "list")
  expect_named(out, c("variable", "knot"))
})

test_that("extractSpecial: knot is a character string of the supplied value", {
  out <- extractSpecial(myVar, 10)
  expect_equal(out$knot, "10")
})

test_that("extractSpecial: errors when k is missing", {
  expect_error(extractSpecial(myVar), regexp = "knotName")
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

# ---------------------------------------------------------------------------
# oom additional edge cases
# ---------------------------------------------------------------------------
test_that("oom: returns 10 for exactly 1 (since ceiling(log10(1)) = 0, 10^0 = 1)", {
  # log10(1) = 0, ceiling(0) = 0, 10^0 = 1
  expect_equal(oom(1), 1)
})

test_that("oom: returns correct value for 0.01", {
  # log10(0.01) = -2, ceiling(-2) = -2, 10^-2 = 0.01
  expect_equal(oom(0.01), 0.01)
})

test_that("oom: consistent for large values", {
  expect_equal(oom(1e6), 1e6)
  expect_equal(oom(1.5e6), 1e7)
})

# ---------------------------------------------------------------------------
# pw edge cases
# ---------------------------------------------------------------------------
test_that("pw: works with negative knot", {
  expect_equal(pw(0, -3), 3)
  expect_equal(pw(-5, -3), 0)
})

test_that("pw: fractional values", {
  expect_equal(pw(1.5, 1.0), 0.5)
  expect_equal(pw(0.9, 1.0), 0)
})
