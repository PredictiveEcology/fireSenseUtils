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

  expect_message(m <- readELFmetrics(d, fireYears = yrs), "dropped 2 of 4 runs made with other fits")
  expect_setequal(m$runName, c("4.10_NRV_370_1", "4.10_NRV_370_2"))
  expect_identical(unique(m$ELF), "4.10")

  expect_equal(NROW(readELFmetrics(d, fireYears = yrs, current = FALSE)), 4L)
  ## a test suffix's runs are current for that suffix
  expect_identical(suppressMessages(readELFmetrics(d, fireYears = yrs, suffix = "_test"))$runName, "5.1_NRV_370_1")
})

test_that("readELFmetrics() returns NULL when there is no file", {
  expect_message(expect_null(readELFmetrics(withr::local_tempdir())), "no ELFmetrics")
})
