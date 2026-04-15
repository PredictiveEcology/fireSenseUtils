## Tests for logistic helper functions in helpers.R

# ---------------------------------------------------------------------------
# logistic4p
# ---------------------------------------------------------------------------
test_that("logistic4p: midpoint at x=0 with symmetric pars", {
  # par = c(min, max, hillSlope, asymmetry)
  # at x=0: min + (max-min)/(1+1)^asymmetry = 0 + 1/(2^1) = 0.5
  expect_equal(logistic4p(0, c(0, 1, 1, 1)), 0.5)
})

test_that("logistic4p: approaches upper asymptote as x -> Inf", {
  expect_equal(logistic4p(1e6, c(0, 1, 1, 1)), 1, tolerance = 1e-6)
})

test_that("logistic4p: approaches lower asymptote as x -> -Inf", {
  expect_equal(logistic4p(-1e6, c(0, 1, 1, 1)), 0, tolerance = 1e-6)
})

test_that("logistic4p: output shifts with non-zero lower asymptote", {
  result <- logistic4p(0, c(0.2, 0.8, 1, 1))
  expect_equal(result, 0.2 + (0.8 - 0.2) / 2, tolerance = 1e-10)
})

test_that("logistic4p: vectorised over x", {
  x <- c(-1e6, 0, 1e6)
  res <- logistic4p(x, c(0, 1, 1, 1))
  expect_length(res, 3)
  expect_true(res[1] < res[2] && res[2] < res[3])
})

test_that("logistic4p: hill slope controls steepness", {
  steep  <- logistic4p(0.5, c(0, 1, 10, 1))
  gentle <- logistic4p(0.5, c(0, 1,  1, 1))
  expect_true(steep > gentle)
})

# ---------------------------------------------------------------------------
# logistic5p
# ---------------------------------------------------------------------------
test_that("logistic5p: output is numeric", {
  expect_type(logistic5p(0, c(0, 1, 1, 1, 1)), "double")
})

test_that("logistic5p: approaches bounds at extremes", {
  expect_equal(logistic5p( 1e6, c(0, 1, 1, 1, 1)), 1, tolerance = 1e-5)
  expect_equal(logistic5p(-1e6, c(0, 1, 1, 1, 1)), 0, tolerance = 1e-5)
})

test_that("logistic5p: vectorised over x", {
  res <- logistic5p(c(-1, 0, 1), c(0, 1, 1, 2, 1))
  expect_length(res, 3)
  expect_true(all(is.finite(res)))
})

test_that("logistic5p: asymmetry parameter shifts midpoint", {
  sym   <- logistic5p(0, c(0, 1, 1, 1, 1))
  asym  <- logistic5p(0, c(0, 1, 1, 1, 2))
  expect_false(isTRUE(all.equal(sym, asym)))
})

# ---------------------------------------------------------------------------
# logistic3p
# ---------------------------------------------------------------------------
test_that("logistic3p: uses par1 as lower baseline", {
  # at x = -1e6, result -> par1
  expect_equal(logistic3p(-1e6, c(1, 1, 1), par1 = 0.05), 0.05, tolerance = 1e-5)
})

test_that("logistic3p: default par1 = 0.1", {
  expect_equal(logistic3p(-1e6, c(1, 1, 1)), 0.1, tolerance = 1e-5)
})

test_that("logistic3p: approaches par[1] as x -> Inf", {
  expect_equal(logistic3p(1e6, c(0.9, 1, 1)), 0.9, tolerance = 1e-5)
})

test_that("logistic3p: vectorised", {
  res <- logistic3p(c(-1, 0, 1), c(1, 1, 1))
  expect_length(res, 3)
  expect_true(res[1] < res[2] && res[2] < res[3])
})

# ---------------------------------------------------------------------------
# logistic2p
# ---------------------------------------------------------------------------
test_that("logistic2p: default par1 = 0.1, par4 = 0.5", {
  # at x = -Inf, -> 0.1
  expect_equal(logistic2p(-1e6, c(1, 1)), 0.1, tolerance = 1e-5)
})

test_that("logistic2p: approaches par[1] at x -> Inf", {
  expect_equal(logistic2p(1e6, c(0.95, 1)), 0.95, tolerance = 1e-5)
})

test_that("logistic2p: custom par1 and par4", {
  res <- logistic2p(0, c(1, 1), par1 = 0.2, par4 = 1)
  # = 0.2 + (1-0.2)/(1+exp(0)^(-1))^1 = 0.2 + 0.8/2 = 0.6
  expect_equal(res, 0.6, tolerance = 1e-10)
})

test_that("logistic2p: vectorised", {
  res <- logistic2p(c(-5, 0, 5), c(1, 1))
  expect_length(res, 3)
  expect_true(all(res >= 0.1 & res <= 1))
})

# ---------------------------------------------------------------------------
# logisticAll
# ---------------------------------------------------------------------------
test_that("logisticAll: dispatches 2p correctly", {
  mat     <- matrix(1, nrow = 1)
  covPars <- 0
  pars2p  <- c(0.9, 1)
  res <- logisticAll(pars2p, mat, covPars, lowerSpreadProb = 0.1)
  expect_equal(res, logistic2p(mat %*% covPars, pars2p, par1 = 0.1))
})

test_that("logisticAll: dispatches 3p correctly", {
  mat     <- matrix(1, nrow = 1)
  covPars <- 0
  pars3p  <- c(0.9, 1, 1)
  res <- logisticAll(pars3p, mat, covPars, lowerSpreadProb = 0.1)
  expect_equal(res, logistic3p(mat %*% covPars, pars3p, par1 = 0.1))
})

test_that("logisticAll: 4p raises error", {
  mat     <- matrix(1, nrow = 1)
  covPars <- 0
  pars4p  <- c(0, 0.9, 1, 1)
  expect_error(logisticAll(pars4p, mat, covPars, lowerSpreadProb = 0.1))
})

test_that("logisticAll: vectorised over mat rows", {
  mat     <- matrix(c(-2, -1, 0, 1, 2), ncol = 1)
  covPars <- 1
  pars2p  <- c(0.9, 1)
  res <- logisticAll(pars2p, mat, covPars, lowerSpreadProb = 0.1)
  expect_length(res, 5)
  expect_true(all(diff(res) > 0)) # monotone increasing
})

# ---------------------------------------------------------------------------
# logisticParamNames (constants)
# ---------------------------------------------------------------------------
test_that("logisticParamNames has correct names", {
  expect_named(logisticParamNames, c("2p", "3p", "4p", "5p"))
})

test_that("logisticParamNames element lengths match parameter counts", {
  expect_length(logisticParamNames[["2p"]], 2)
  expect_length(logisticParamNames[["3p"]], 3)
  expect_length(logisticParamNames[["4p"]], 4)
  expect_length(logisticParamNames[["5p"]], 5)
})

test_that("logisticParamNames entries are character vectors", {
  expect_true(all(vapply(logisticParamNames, is.character, logical(1))))
})
