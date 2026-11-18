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

test_that("readELFmetrics() keeps the runs made with the current fits, with ELF codes as text", {
  d <- withr::local_tempdir()
  yrs <- 1986:2024
  sp <- spreadFitFilenameFor(yrs)
  ig <- ignitionFitFilenameFor(yrs)
  writeRun(d, "4.10_NRV_370_1", "4.10", sp, ig, areaBurnedHaPerYear = 100)
  writeRun(d, "4.10_NRV_370_2", "4.10", sp, ig, areaBurnedHaPerYear = 300)
  writeRun(d, "6.2_NRV_370_1", "6.2", sp, "fireSenseIgnitionParams_1986-2024_xgboost.rds", areaBurnedHaPerYear = 9)
  writeRun(d, "5.1_NRV_370_1", "5.1", spreadFitFilenameFor(yrs, "_test"), ignitionFitFilenameFor(yrs, "_test"),
           areaBurnedHaPerYear = 9)

  ## 6.2's ignition file is of another model; 5.1's test files are its newest of the current model
  expect_message(m <- readELFmetrics(d, fireYears = yrs), "dropped 1 of 4 runs made with other fits")
  expect_setequal(m$runName, c("4.10_NRV_370_1", "4.10_NRV_370_2", "5.1_NRV_370_1"))
  expect_setequal(m$ELF, c("4.10", "5.1"))

  expect_equal(NROW(readELFmetrics(d, fireYears = yrs, current = FALSE)), 4L)
  ## a test suffix's runs are current for that suffix
  expect_identical(suppressMessages(readELFmetrics(d, fireYears = yrs, suffix = "_test"))$runName, "5.1_NRV_370_1")
})

test_that("readELFmetrics() keeps each ELF's runs made with its newest fits, whatever their suffix", {
  ## fireSense_summary records the files "latest" read (any suffix); the newest by the date in the suffix is current
  d <- withr::local_tempdir()
  yrs <- 1986:2024
  writeRun(d, "14.1_a", "14.1", spreadFitFilenameFor(yrs, "_gpscCause_20261029"),
           ignitionFitFilenameFor(yrs, "_gpscCause_20261029"))
  writeRun(d, "14.1_b", "14.1", spreadFitFilenameFor(yrs), ignitionFitFilenameFor(yrs))   # older: no date
  writeRun(d, "6.2_a", "6.2", spreadFitFilenameFor(yrs, "_prod20261021"), ignitionFitFilenameFor(yrs, "_prod20261021"))
  ## a run over two ELFs, whose spread fits came from two files, is as old as its oldest file: its spread fits
  ## are as new as 6.2_a's, its ignition fit newer, so 6.2_a is dropped
  writeRun(d, "6.2_b", "6.2",
           paste(spreadFitFilenameFor(yrs, "_prod20261021"), spreadFitFilenameFor(yrs, "_x20261101"), sep = ", "),
           ignitionFitFilenameFor(yrs, "_x20261101"))
  m <- suppressMessages(readELFmetrics(d, fireYears = yrs))
  expect_setequal(m$runName, c("14.1_a", "6.2_b"))
  ## a suffix pins one experiment's files
  expect_identical(suppressMessages(readELFmetrics(d, fireYears = yrs, suffix = "_gpscCause_20261029"))$runName, "14.1_a")
})

test_that("readELFmetrics() returns NULL when there is no file", {
  expect_message(expect_null(readELFmetrics(withr::local_tempdir())), "no ELFmetrics")
})
