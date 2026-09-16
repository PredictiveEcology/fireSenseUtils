## A DEoptim generation costs the SLOWEST of its NP evaluations, not the median one: the
## generation is synchronous, so every worker waits for the last. Measured on ELF 4.1 during
## FireSense phase 2 (2026-09-16), p90 was 40.9 s against a max of 85.9 s -- a 2x spread inside
## the top decile alone. The tail sets the wall clock.
##
## .objfunSpreadFit() already has exactly one place where it can stop an evaluation early. Its
## work is two blocks -- `lrgSmallFireYears <- list(large = ..., small = ...)` -- and after the
## first it compares the accumulated SNLL against `thresh * numYrsDone`. That bound is STATIC for
## the whole fit, so it cannot tighten as the population improves: by late generations every
## trial that DEoptim will reject anyway still passes it, and nothing is cut.
##
## `pruneAbove` makes that same bound adaptive. The caller passes the worst value the current
## population would accept, and the check becomes min(static, pruneAbove).
##
## This is EXACT, not a heuristic. Block 2 contributes a non-negative SNLL, so block 1's
## accumulated value is a lower bound on the evaluation's final value. A trial pruned here would
## have finished at or above `pruneAbove`, which every parent already beats, so DEoptim's
## selection would have discarded it. The search trajectory is unchanged; only the time spent
## proving the trial is bad is saved.
##
## Two things below are therefore load-bearing for that argument and are asserted, not assumed:
## the default must be Inf (min(static, Inf) is today's behaviour exactly), and the bail must stay
## gated on `ii == 1`, because block 2 is the last block -- a lower bound there prunes nothing and
## the "non-negative remainder" reasoning no longer applies.

test_that(".objfunSpreadFit takes pruneAbove, defaulting to Inf", {
  ## `:::`: .objfunSpreadFit() is internal; DEoptim calls it through runDEoptim().
  fm <- formals(fireSenseUtils:::.objfunSpreadFit)
  expect_true("pruneAbove" %in% names(fm))
  expect_identical(eval(fm$pruneAbove), Inf)
})

test_that("the early-bail bound respects pruneAbove as well as the static thresh", {
  src <- paste(deparse(fireSenseUtils:::.objfunSpreadFit), collapse = "\n")
  expect_match(src, "min\\(thresh \\* numYrsDone, pruneAbove\\)")
})

test_that("the early bail stays restricted to the first block", {
  ## Pruning is only exact where a non-negative remainder is still to come, i.e. not in the
  ## final block. If this assertion ever fails, the exactness argument above has been broken.
  src <- paste(deparse(fireSenseUtils:::.objfunSpreadFit), collapse = "\n")
  expect_match(src, "ii == 1")
})
