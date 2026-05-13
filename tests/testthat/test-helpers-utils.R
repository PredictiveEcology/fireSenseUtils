## Tests for utility helper functions: oom, pw, logMinB, paramsSeparate,
## dtReplaceNAwith0, rbetaBetween, toX1000, yearChar, youngAgeName

library(data.table)

# ---------------------------------------------------------------------------
# oom – order of magnitude
# ---------------------------------------------------------------------------
test_that("oom: returns 10 for values in (1, 10]", {
  expect_equal(oom(5),   10)
  expect_equal(oom(10),  10)
})

test_that("oom: returns 100 for values in (10, 100]", {
  expect_equal(oom(11),  100)
  expect_equal(oom(100), 100)
})

test_that("oom: returns 1 for values in (0, 1]", {
  expect_equal(oom(0.5), 1)
  expect_equal(oom(1),   1)
})

test_that("oom: handles negative numbers via abs()", {
  expect_equal(oom(-5),  oom(5))
  expect_equal(oom(-50), oom(50))
})

test_that("oom: returns 1000 for 500", {
  expect_equal(oom(500), 1000)
})

# ---------------------------------------------------------------------------
# pw – piecewise (hinge) function
# ---------------------------------------------------------------------------
test_that("pw: returns 0 when variable <= knot", {
  expect_equal(pw(3, 5),   0)
  expect_equal(pw(5, 5),   0)
  expect_equal(pw(-1, 0),  0)
})

test_that("pw: returns positive excess when variable > knot", {
  expect_equal(pw(7, 5),  2)
  expect_equal(pw(10, 3), 7)
})

test_that("pw: vectorised over variable", {
  res <- pw(c(1, 3, 5, 7), 4)
  expect_equal(res, c(0, 0, 1, 3))
})

test_that("pw: output is always non-negative", {
  expect_true(all(pw(c(-10, 0, 5, 10), 3) >= 0))
})

# ---------------------------------------------------------------------------
# logMinB
# ---------------------------------------------------------------------------
test_that("logMinB: values above minimum are log-transformed unchanged", {
  minimumB <- exp(log(100) - 1)  # same as in source
  big <- minimumB * 10
  expect_equal(logMinB(big), log(big))
})

test_that("logMinB: values below minimum are clamped then log-transformed", {
  minimumB <- exp(log(100) - 1)
  expect_equal(logMinB(0), log(minimumB))
  expect_equal(logMinB(1), log(minimumB))
})

test_that("logMinB: result is never less than log(minimumB)", {
  minimumB <- exp(log(100) - 1)
  vals <- c(0, 0.001, 1, 10, 100, 1000)
  expect_true(all(logMinB(vals) >= log(minimumB)))
})

test_that("logMinB: vectorised", {
  res <- logMinB(c(50, 100, 200))
  expect_length(res, 3)
  expect_true(all(is.finite(res)))
})

# ---------------------------------------------------------------------------
# paramsSeparate
# ---------------------------------------------------------------------------
test_that("paramsSeparate: splits correctly with 3 logistic + 2 cov pars", {
  par <- 1:5
  res <- paramsSeparate(par, parsModel = 2)
  expect_equal(res$covPars,      c(4L, 5L))
  expect_equal(res$logisticPars, c(1L, 2L, 3L))
})

test_that("paramsSeparate: all covariates case", {
  par <- 1:4
  res <- paramsSeparate(par, parsModel = 4)
  expect_equal(res$covPars,      1:4)
  expect_length(res$logisticPars, 0)
})

test_that("paramsSeparate: all logistic case", {
  par <- 1:4
  res <- paramsSeparate(par, parsModel = 0)
  expect_length(res$covPars, 0)
  expect_equal(res$logisticPars, 1:4)
})

test_that("paramsSeparate: single element", {
  res <- paramsSeparate(42, parsModel = 1)
  expect_equal(res$covPars,      42)
  expect_length(res$logisticPars, 0)
})

test_that("paramsSeparate: returns a list with two elements", {
  res <- paramsSeparate(1:6, parsModel = 3)
  expect_type(res, "list")
  expect_named(res, c("covPars", "logisticPars"))
})

# ---------------------------------------------------------------------------
# dtReplaceNAwith0
# ---------------------------------------------------------------------------
test_that("dtReplaceNAwith0: replaces NAs with 0 in all columns by default", {
  DT <- data.table(a = c(1, NA, 3), b = c(NA, 2, NA))
  out <- dtReplaceNAwith0(DT)
  expect_equal(out$a, c(1, 0, 3))
  expect_equal(out$b, c(0, 2, 0))
})

test_that("dtReplaceNAwith0: only touches specified colsToUse", {
  DT <- data.table(a = c(1, NA, 3), b = c(NA, 2, NA))
  out <- dtReplaceNAwith0(DT, colsToUse = "a")
  expect_equal(out$a, c(1, 0, 3))
  expect_true(is.na(out$b[1]))  # b untouched
})

test_that("dtReplaceNAwith0: no-op when no NAs present", {
  DT <- data.table(x = 1:3, y = 4:6)
  before <- copy(DT)
  out <- dtReplaceNAwith0(DT)
  expect_equal(out, before)
})

test_that("dtReplaceNAwith0: modifies in place (same object)", {
  DT <- data.table(a = c(NA, 2))
  out <- dtReplaceNAwith0(DT)
  expect_true(identical(DT, out))
})

test_that("dtReplaceNAwith0: handles all-NA column", {
  DT <- data.table(a = c(NA_real_, NA_real_))
  out <- dtReplaceNAwith0(DT)
  expect_equal(out$a, c(0, 0))
})

test_that("dtReplaceNAwith0: handles integer column with NA", {
  DT <- data.table(a = c(1L, NA_integer_, 3L))
  out <- dtReplaceNAwith0(DT)
  expect_equal(out$a[2], 0L)
})

# ---------------------------------------------------------------------------
# rbetaBetween
# ---------------------------------------------------------------------------
test_that("rbetaBetween: output is in [l, u]", {
  set.seed(42)
  out <- rbetaBetween(100, l = 2, u = 5, m = 3, shape1 = 3)
  expect_true(all(out >= 2 & out <= 5))
})

test_that("rbetaBetween: correct length", {
  set.seed(1)
  out <- rbetaBetween(50, l = 0, u = 1, m = 0.5, shape1 = 2)
  expect_length(out, 50)
})

test_that("rbetaBetween: explicit shape2 bypasses mean calculation", {
  set.seed(7)
  out <- rbetaBetween(20, l = 0, u = 1, m = 0.5, shape1 = 2, shape2 = 2)
  expect_true(all(out >= 0 & out <= 1))
})

test_that("rbetaBetween: mean of output near m (stochastic, loose check)", {
  set.seed(99)
  out <- rbetaBetween(1e4, l = 0, u = 10, m = 4, shape1 = 5)
  expect_equal(mean(out), 4, tolerance = 0.2)
})

# ---------------------------------------------------------------------------
# toX1000
# ---------------------------------------------------------------------------
test_that("toX1000: multiplies numeric columns by 1000 and coerces to integer", {
  df <- list(data.frame(pixelID = 1:3, cov1 = c(0.1, 0.2, 0.3)))
  out <- toX1000(df)
  expect_equal(out[[1]]$cov1, c(100L, 200L, 300L))
})

test_that("toX1000: pixelID column is omitted by default", {
  df <- list(data.frame(pixelID = 1:3, cov1 = c(0.5, 0.5, 0.5)))
  out <- toX1000(df)
  expect_equal(out[[1]]$pixelID, 1:3)  # unchanged
})

test_that("toX1000: custom omitCols respected", {
  df <- list(data.frame(id = 1:2, x = c(0.1, 0.2), y = c(0.3, 0.4)))
  out <- toX1000(df, omitCols = c("id", "x"))
  expect_equal(out[[1]]$id, 1:2)       # untouched
  expect_equal(out[[1]]$x, c(0.1, 0.2)) # untouched
  expect_equal(out[[1]]$y, c(300L, 400L))
})

test_that("toX1000: works on a list of multiple data.frames", {
  dfs <- list(
    yr1 = data.frame(pixelID = 1:2, cov = c(0.1, 0.2)),
    yr2 = data.frame(pixelID = 1:2, cov = c(0.3, 0.4))
  )
  out <- toX1000(dfs)
  expect_length(out, 2)
  expect_equal(out[["yr1"]]$cov, c(100L, 200L))
  expect_equal(out[["yr2"]]$cov, c(300L, 400L))
})

# ---------------------------------------------------------------------------
# yearChar and youngAgeName (exported constants)
# ---------------------------------------------------------------------------
test_that("yearChar is the string 'year'", {
  expect_equal(yearChar, "year")
})

test_that("youngAgeName is the string 'youngAge'", {
  expect_equal(youngAgeName, "youngAge")
})
