## readELFmetrics() collects the one-row files fireSense_summary writes per run, and drops runs made with other fits.

writeRun <- function(dir, runName, ELF, spreadFitFile, ignitionFitFile, ...) {
  dir.create(file.path(dir, runName), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(data.table::data.table(runName = runName, ELF = ELF, spreadFitFile = spreadFitFile,
                                            ignitionFitFile = ignitionFitFile, ...),
                     file.path(dir, runName, ELFmetricsFilename(runName)))
}

test_that("ELFmetricsFilename() names the file by the run", {
  expect_identical(ELFmetricsFilename("13.1_NRV_370_1"), "ELFmetrics_13.1_NRV_370_1.csv")
})

test_that("readELFmetrics() keeps each ELF's runs made with its newest fits, with ELF codes as text", {
  d <- withr::local_tempdir()
  yrs <- 1986:2024
  sp <- spreadFitFilenameFor(yrs)
  ig <- ignitionFitFilenameFor(yrs)
  writeRun(d, "4.10_NRV_370_1", "4.10", sp, ig, areaBurnedHaPerYear = 100)
  writeRun(d, "4.10_NRV_370_2", "4.10", sp, ig, areaBurnedHaPerYear = 300)
  writeRun(d, "4.10_NRV_370_3", "4.10", sp, ignitionFitFilenameFor(yrs, "_old20260101"), areaBurnedHaPerYear = 9)
  writeRun(d, "5.1_NRV_370_1", "5.1", spreadFitFilenameFor(yrs, "_test"), ignitionFitFilenameFor(yrs, "_test"))

  ## undated files are oldest, so 4.10's run with a dated (2026-01-01) ignition file is its newest; 5.1 has one run
  expect_message(m <- readELFmetrics(d), "dropped 2 of 4 runs made with other fits")
  expect_setequal(m$runName, c("4.10_NRV_370_3", "5.1_NRV_370_1"))
  expect_setequal(m$ELF, c("4.10", "5.1"))
  expect_equal(NROW(readELFmetrics(d, current = FALSE)), 4L)
  ## `files` pins the runs of those ledger files
  expect_setequal(suppressMessages(readELFmetrics(d, files = c(sp, ig)))$runName, c("4.10_NRV_370_1", "4.10_NRV_370_2"))
})

test_that("readELFmetrics() orders by when the fits were written, else by the date in the file names", {
  ## file names are labels now, so the fits' own writtenAt decides; older metrics files have only the names
  d <- withr::local_tempdir()
  yrs <- 1986:2024
  writeRun(d, "14.1_a", "14.1", spreadFitFilenameFor(yrs, "_gpscCause_20261029"),
           ignitionFitFilenameFor(yrs, "_gpscCause_20261029"))
  writeRun(d, "14.1_b", "14.1", spreadFitFilenameFor(yrs), ignitionFitFilenameFor(yrs))   # older: no date
  ## undated labels, but written after 2026-10-29
  writeRun(d, "6.2_a", "6.2", "fireSenseParams_myLabel.rds", "fireSenseIgnitionParams_myLabel.rds",
           spreadFitWrittenAt = "2026-11-05 10:00:00", ignitionFitWrittenAt = "2026-11-05 10:00:00")
  writeRun(d, "6.2_b", "6.2", spreadFitFilenameFor(yrs, "_x20261101"), ignitionFitFilenameFor(yrs, "_x20261101"),
           spreadFitWrittenAt = "2026-11-01 10:00:00", ignitionFitWrittenAt = "2026-11-01 10:00:00")
  m <- suppressMessages(readELFmetrics(d))
  expect_setequal(m$runName, c("14.1_a", "6.2_a"))
  ## a run whose files are not known (NA: fireSense_summary before it could read them) is never current
  writeRun(d, "9.3_a", "9.3", spreadFitFilenameFor(yrs, "_gpscCause_20261029"), NA_character_)
  expect_false("9.3_a" %in% suppressMessages(readELFmetrics(d))$runName)
})

test_that("readELFmetrics() returns NULL when there is no file", {
  expect_message(expect_null(readELFmetrics(withr::local_tempdir())), "no ELFmetrics")
})
