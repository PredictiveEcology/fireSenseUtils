## A toy ledger in the national format: ELFs "1.1" and "1.2" have yearEffects, "1.3" has none (an older fit).
toyYears <- 2001:2008
toyEffects <- function(mult, priorSD = 0.5) data.frame(year = toyYears, observed = 5, expected = 5, priorSD = priorSD,
                                                      multiplier = mult, shape = if (priorSD > 0) 10 else NA_real_)
toyLedger <- function(dir, name = "fireSenseIgnitionParams_toy") {
  shared <- c(0.6, 1.5, 0.9, 2, 1.1, 0.7, 1.3, 0.8)
  x <- data.frame(polygonID = c("1.1", "1.2", "1.3"))
  x$fireSense_IgnitionFitted <- list(list(modelList = list(yearEffects = toyEffects(shared))),
                                     list(modelList = list(yearEffects = toyEffects(shared^2))),
                                     list(modelList = list()))
  x$fireSense_EscapeFitted <- list(list(modelList = list(yearEffects = toyEffects(rep(1, 8), priorSD = 0))),
                                   list(modelList = list(yearEffects = toyEffects(rev(shared)))),
                                   list(modelList = list()))
  f <- file.path(dir, paste0(name, ".rds")); saveRDS(x, f); f
}

test_that("readYearEffects() gives one long table and skips ELFs without yearEffects", {
  f <- toyLedger(withr::local_tempdir())
  expect_message(tab <- readYearEffects(f), "1.3")
  expect_setequal(unique(tab$ELF), c("1.1", "1.2"))
  expect_equal(nrow(tab), 2 * 2 * length(toyYears))
  expect_setequal(tab$process, c("ignition", "escape"))
  expect_identical(unique(tab$ledger), "fireSenseIgnitionParams_toy")
  expect_true(all(c("year", "observed", "expected", "priorSD", "multiplier", "shape") %in% names(tab)))
  expect_equal(unique(tab$priorSD[tab$ELF == "1.1" & tab$process == "escape"]), 0)
})

test_that("readYearEffects() stops on missing files, no files, and two rows for one ELF", {
  expect_error(readYearEffects(file.path(tempdir(), "nope.rds")), "not found")
  expect_error(readYearEffects(NA_character_), "at least one")
  f <- toyLedger(withr::local_tempdir()); x <- readRDS(f); saveRDS(rbind(x, x[1, ]), f)
  expect_error(readYearEffects(f), "more than one row for ELF")
})
