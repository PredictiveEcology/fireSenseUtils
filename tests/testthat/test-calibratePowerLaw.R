## calibratePowerLaw(): p_c50 and gamma of the "powerLaw" link from runaway-curve simulations on a flammable map.

flammableMap <- function(n = 201L, breaks = FALSE) {
  r <- terra::rast(nrows = n, ncols = n, xmin = 0, xmax = n * 240, ymin = 0, ymax = n * 240, vals = 1)
  if (breaks) { # non-flammable lines every 10 cells, each with a 4-cell gap every 10 cells: fire must find the gaps
    m <- matrix(1, n, n)
    m[seq(5L, n, by = 10L), ] <- 0
    m[, seq(5L, n, by = 10L)] <- 0
    gaps <- sort(outer(1:4, seq(0L, n, by = 10L), "+"))
    gaps <- gaps[gaps <= n]
    m[seq(5L, n, by = 10L), gaps] <- 1
    m[gaps, seq(5L, n, by = 10L)] <- 1
    terra::values(r) <- as.vector(t(m))
  }
  r
}
pGridSmall <- seq(0.15, 0.40, by = 0.025)

test_that("a uniform map percolates near 0.25 and gives a positive gamma", {
  cal <- calibratePowerLaw(flammableMap(), pGrid = pGridSmall, nFires = 40, cores = 1, seed = 1)
  expect_named(cal, c("p50", "gamma", "R2", "table"))
  expect_gt(cal$p50, 0.23)
  expect_lt(cal$p50, 0.30)
  expect_gt(cal$gamma, 0.5)
  expect_lt(cal$gamma, 6)
  expect_gt(cal$R2, 0.8)
  expect_true(all(c("p", "meanHa", "shareRunaway") %in% names(cal$table)))
  expect_identical(cal$table$p, pGridSmall)
})

test_that("non-flammable breaks raise p_c50", {
  open <- calibratePowerLaw(flammableMap(), pGrid = pGridSmall, nFires = 40, seed = 1)
  broken <- calibratePowerLaw(flammableMap(breaks = TRUE), pGrid = pGridSmall, nFires = 40, seed = 1)
  expect_gt(broken$p50, open$p50)
})

test_that("the result depends on the seed only, not on the number of cores", {
  a <- calibratePowerLaw(flammableMap(101L), pGrid = pGridSmall, nFires = 20, cores = 1, seed = 7)
  b <- calibratePowerLaw(flammableMap(101L), pGrid = pGridSmall, nFires = 20, cores = 2, seed = 7)
  d <- calibratePowerLaw(flammableMap(101L), pGrid = pGridSmall, nFires = 20, cores = 1, seed = 8)
  expect_identical(a, b)
  expect_false(identical(a$table, d$table))
})

test_that("the caller's RNG state is left alone", {
  set.seed(3); before <- .Random.seed
  calibratePowerLaw(flammableMap(51L), pGrid = pGridSmall, nFires = 5, seed = 1)
  expect_identical(.Random.seed, before)
})

test_that("a grid that never reaches a runaway stops with a message", {
  expect_error(calibratePowerLaw(flammableMap(51L), pGrid = c(0.15, 0.16, 0.17, 0.18, 0.19), nFires = 10),
               "never reach")
})

test_that("powerLawPcBounds(): floor + 0.01 to p50 + margin, and p50 must be inside the outer limits", {
  expect_equal(powerLawPcBounds(0.26), c(spreadProbFloor + 0.01, 0.265))
  expect_equal(powerLawPcBounds(0.26, margin = 0.02)[2], 0.28)
  expect_equal(powerLawPcBounds(0.349)[2], 0.35)
  expect_error(powerLawPcBounds(0.4), "outside")
  expect_error(powerLawPcBounds(0.1), "outside")
})
