## The ledger's "best" parameter sets, chosen by replicated mean.
##
## FireSense, 2026-09-17: the ledger's "5 best" were five copies of one member -- fireSense_SpreadFit took
## the 5 GENERATIONS with the lowest best value, which all hold the same frozen best member -- and that
## member was itself a lucky draw: re-scored 10 times, DEoptim's best ranked 1st to 6th of 60 in eight fits.
## runDEoptim() now re-scores the final population, and bestByReplicatedMean() picks from the means.

## member i's true value is i; one member is recorded (and would be picked) as far better than it is
noisy <- function(par, shift = 0) par[["a"]] + shift + stats::rnorm(1, sd = 0.3)
pop <- cbind(a = c(5, 1, 4, 2, 3, 6), b = 0)

test_that("rescorePopulation() evaluates every member reps times, reproducibly", {
  s1 <- rescorePopulation(pop, noisy, reps = 4L, seed = 10L)
  expect_identical(nrow(s1), 24L)
  expect_identical(sort(unique(s1$member)), 1:6)
  expect_identical(s1, rescorePopulation(pop, noisy, reps = 4L, seed = 10L))
  expect_false(identical(s1$value, rescorePopulation(pop, noisy, reps = 4L, seed = 11L)$value))
})

test_that("rescorePopulation() passes fnArgs to the objective", {
  s <- rescorePopulation(pop, noisy, reps = 2L, fnArgs = list(shift = 100))
  expect_true(all(s$value > 90))
})

test_that("bestByReplicatedMean() returns the n lowest means, best first, with their parameters", {
  s <- rescorePopulation(pop, noisy, reps = 20L)
  b <- bestByReplicatedMean(pop, s, n = 3L)
  expect_identical(b$member, c(2L, 4L, 5L))                # true values 1, 2, 3
  expect_equal(b$params$a, c(1, 2, 3))
  expect_identical(names(b$params), c("a", "b"))
  expect_true(all(diff(b$objFunVal) > 0))
  expect_length(b$objFunValSD, 3L)
})

test_that("identical members count once, so the n returned are distinct", {
  dup <- rbind(pop, pop[2, , drop = FALSE], pop[2, , drop = FALSE])   # member 2 three times
  s <- rescorePopulation(dup, noisy, reps = 20L)
  b <- bestByReplicatedMean(dup, s, n = 5L)
  expect_identical(nrow(unique(b$params)), 5L)
  expect_identical(b$member[1], 2L)
  expect_length(bestByReplicatedMean(pop, rescorePopulation(pop, noisy, reps = 2L), n = 10L)$member, 6L)
})

test_that("rescorePopulation() runs on a cluster's workers", {
  skip_on_cran()
  cl <- parallel::makeCluster(2L)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  expect_identical(rescorePopulation(pop, noisy, reps = 3L, cl = cl, seed = 5L),
                   rescorePopulation(pop, noisy, reps = 3L, seed = 5L))
})

test_that("runDEoptim() re-scores the final population with the early stop off, after the fit", {
  src <- paste(deparse(runDEoptim), collapse = "\n")
  expect_match(src, "rescorePopulation", fixed = TRUE)
  expect_match(src, "thresh = Inf", fixed = TRUE)
  expect_lt(regexpr("DEoptimIterative2", src), regexpr("rescorePopulation", src))
  expect_identical(formals(runDEoptim)$rescoreReps, 10L)
})

test_that("the ledger carries covMinMax_spread, which prediction needs to rescale covariates", {
  expect_true("covMinMax_spread" %in% spreadFitAdditionalColNamesTxt)
})
