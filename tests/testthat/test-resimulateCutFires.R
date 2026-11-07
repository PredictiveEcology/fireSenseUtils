## A replicate fire cut at an open buffer edge is scored at the size it reached inside its buffer. With
## `annualDTx1000Whole` and `nonAnnualDTx1000Whole` given, every cut fire is re-simulated on the whole landscape
## from its ignition, with its replicate's year effect, until it burns an open edge pixel (rejection: given that
## the buffered fire is cut, the whole-landscape fire is a fire that burns an open edge pixel), so every
## size-based term scores its unbounded size. See helper-resimulateCut.R.

test_that("buffers over the whole landscape have no cut fires: the tables change nothing", {
  skip_if_not_installed("SpaDES.tools")
  off <- resimObjective(side = 60L)
  on <- resimObjective(side = 60L, whole = TRUE)
  expect_equal(off$runaways, 0)
  expect_equal(c(on$resimulated, on$capHit), c(0, 0))
  expect_identical(on[names(off)], off)
})

test_that("with the tables, cut fires are re-simulated and no longer scored as cut", {
  skip_if_not_installed("SpaDES.tools")
  off <- resimObjective()
  on <- resimObjective(whole = TRUE)
  expect_gt(off$runaways, 0)
  expect_equal(on$runaways, off$runaways)                  # the buffered cut count, as before
  expect_equal(on$resimulated + on$capHit, off$runaways)
  expect_gt(max(on$allFireSizes), max(off$allFireSizes))  # past the 11 x 11 buffer
  expect_lt(sum(on$crpsFire), Inf)
})

test_that("returnSims gives the buffered size, the re-simulated size and whether the fire is still cut", {
  skip_if_not_installed("SpaDES.tools")
  off <- resimObjective(returnSims = TRUE)
  on <- resimObjective(whole = TRUE, returnSims = TRUE)$sims
  expect_false(any(c("simBuffered", "cut") %in% names(off$sims)))
  expect_true(all(c("simBuffered", "cut") %in% names(on)))
  expect_identical(on$simBuffered, off$sims$sim)           # the same replicates, the same stream
  expect_gt(mean(on$sim), mean(on$simBuffered))          # a cut fire is redrawn on the whole landscape
  expect_true(any(on$sim > on$simBuffered))
})

test_that("a re-simulated cut fire is distributed as a whole-landscape fire that burns an open edge pixel", {
  skip_if_not_installed("SpaDES.tools")
  fx <- resimFixture()
  open <- addBufferEdge(list(fx$buf), resimLandscape)[[1]]
  openPix <- open$pixelID[open$open]
  N <- 4000L
  p <- 0.2
  ## the reference: the same fire on the whole landscape, kept when it burns an open edge pixel
  set.seed(11)
  ref <- vapply(seq_len(N), function(i) {
    o <- SpaDES.tools::spreadCpp(landscape = resimLandscape, loci = fx$fires$cells, spreadProb = rep(p, 3600))
    if (any(o$indices %in% openPix)) nrow(o) else NA_integer_
  }, integer(1))
  ref <- ref[!is.na(ref)]
  on <- resimObjective(p = p, whole = TRUE, Nreps = N, seed = 12, returnSims = TRUE, returnBurned = TRUE,
                       cutMaxTries = 200L)
  cutRep <- unique(on$burned[pixelID %in% openPix]$rep)           # the replicates cut in the buffer
  expect_gt(length(ref), 300)
  expect_gt(length(cutRep), 300)
  expect_identical(on$sims$sim[-cutRep], on$sims$simBuffered[-cutRep])   # an uncut fire keeps its size
  ## sizes are whole numbers; spread them over their unit so that the test has no ties
  untied <- function(x) x + stats::runif(length(x), -0.5, 0.5)
  set.seed(13)
  expect_gt(stats::ks.test(untied(on$sims$sim[cutRep]), untied(ref))$p.value, 0.01)
  ## and the buffered sizes of the same replicates are not that distribution: they are cut at the buffer
  expect_lt(stats::ks.test(untied(on$sims$simBuffered[cutRep]), untied(ref))$p.value, 0.01)
})

test_that("with the tables NULL the objective is what it was before the re-simulation existed", {
  skip_if_not_installed("SpaDES.tools")
  ## equal to the last bit (tolerance = 0); identical() would also compare data.table internals
  ## fixtures/objFunInnerBeforeResimulation.rds was made by the same calls on fireSenseUtils 0.2.3.9155 (c86d0a4)
  before <- readRDS(test_path("fixtures", "objFunInnerBeforeResimulation.rds"))
  expect_equal(resimObjective(yearSpreadSD = 0.4, Nreps = 30L, seed = 5), before$crps, tolerance = 0)
  expect_equal(resimObjective(Nreps = 30L, seed = 6), before$noEffect, tolerance = 0)
  sims <- resimObjective(yearSpreadSD = 0.4, Nreps = 30L, seed = 5, returnSims = TRUE, returnBurned = TRUE)
  expect_equal(sims, before$sims, tolerance = 0)
})

test_that("a replicate's year effect is the one its cut fire is re-simulated with", {
  skip_if_not_installed("SpaDES.tools")
  fx <- resimFixture()
  openPix <- addBufferEdge(list(fx$buf), resimLandscape)[[1]]
  openPix <- openPix$pixelID[openPix$open]
  spreadState <- data.table::data.table(rep = 1:2, initialLocus = fx$fires$cells, indices = fx$fires$cells)
  simSizes <- data.table::data.table(rep = 1:2, initialLocus = fx$fires$cells, N = 1L, cut = TRUE)
  run <- function(eps) resimulateCutFires(
    spreadState, simSizes, openPix, fx$fires$cells, eps,
    wholeLp = function() list(pixelID = 1:3600, lp = rep(0, 3600)), toP = function(z) stats::plogis(z - 1.4),
    r = resimLandscape, escapeMinPx = NULL, jumpTries = 0, jumpMeanDist = 0, maxTries = 500L)
  set.seed(2)
  big <- vapply(1:30, function(i) run(c(0.4, -0.4))$simSizes$N, numeric(2))   # plogis(-1) = 0.27, plogis(-1.8) = 0.14
  expect_gt(median(big[1, ]), median(big[2, ]))                              # replicate 1 had the larger effect
  out <- run(c(0.4, -0.4))
  expect_identical(out$simSizes$Nbuffered, c(1L, 1L))
  expect_identical(out$resimulated + out$capHit, 2L)
})

test_that("the whole-landscape tables stop when the cut fires are penalised on purpose, or when only one is given", {
  expect_error(checkCutFireTables(list(a = 1), list(b = 1), penaliseRunaways = TRUE), "penaliseRunaways = FALSE")
  expect_error(checkCutFireTables(list(a = 1), NULL, penaliseRunaways = FALSE), "go together")
  expect_false(checkCutFireTables(NULL, NULL, penaliseRunaways = TRUE))
  expect_true(checkCutFireTables(list(a = 1), list(b = 1), penaliseRunaways = FALSE))
  expect_error(checkCutFireTables(list(a = 1), list(b = 1), FALSE, annualDTx1000 = list(a = 1, c = 2)), "lacks the fire year")
})
