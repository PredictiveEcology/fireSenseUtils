## elfFitDates(): per ELF, the date of its latest fit in the spread and ignition ledgers (every ledger file of
## each folder, as latestSpreadFits() and latestIgnitionFits() read them). Drive is faked: no network.

spreadFile <- function(suffix) paste0("fireSenseParams_1985-2024", spreadFitFileTag, suffix, ".rds")
ignitionFile <- function(suffix) paste0(ignitionFitFilePrefix, "1985-2024", ignitionFitFileTag, suffix, ".rds")

## fakeDrive() needs names built from the model tags, which a list(name = ) call cannot take
drive <- function(...) {
  x <- list(...)
  fakeDrive(stats::setNames(lapply(x, `[[`, 2), vapply(x, `[[`, "", 1)), env = parent.frame())
}
f <- function(name, modified, ledger) list(name, list(modified = modified, ledger = ledger))

test_that("elfFitDates() takes each ELF's latest fit over several ledger files, and names the ledger", {
  d <- withr::local_tempdir()
  ## same fake folder for both: the prefixes tell the ledgers apart
  drive(
    f(spreadFile("_a"), "2026-10-01T10:00:00Z", ledgerRows(c("4.1", "4.3"), 1, writtenAt = "2026-10-01 08:00:00")),
    f(spreadFile("_b"), "2026-11-10T10:00:00Z", ledgerRows(c("4.1", "5.1"), 2, writtenAt = "2026-11-10 09:30:00")),
    f(ignitionFile("_c"), "2026-11-05T10:00:00Z", ledgerRows(c("4.1", "6.1"), 3, writtenAt = "2026-11-05 12:00:00")))
  out <- suppressMessages(elfFitDates("spreadFolder", "ignitionFolder", destinationPath = d))
  expect_identical(names(out), c("ELF", "spreadFitDate", "spreadLedger", "ignitionFitDate", "ignitionLedger"))
  expect_setequal(out$ELF, c("4.1", "4.3", "5.1", "6.1"))
  r <- function(id) out[out$ELF == id, ]
  expect_equal(r("4.1")$spreadFitDate, as.POSIXct("2026-11-10 09:30:00", tz = "UTC"))   # the newer file's row
  expect_identical(r("4.1")$spreadLedger, spreadFile("_b"))
  expect_equal(r("4.3")$spreadFitDate, as.POSIXct("2026-10-01 08:00:00", tz = "UTC"))
  expect_identical(r("4.1")$ignitionLedger, ignitionFile("_c"))
  expect_true(is.na(r("5.1")$ignitionFitDate) && is.na(r("5.1")$ignitionLedger))       # spread only
  expect_true(is.na(r("6.1")$spreadFitDate) && is.na(r("6.1")$spreadLedger))           # ignition only
})

test_that("a row with no writtenAt is dated by its file's name, and `ELFs` adds ELFs with no fit", {
  d <- withr::local_tempdir()
  drive(
    f(spreadFile("_old_20260925"), "2026-11-01T10:00:00Z", ledgerRows("4.1", 1)),
    f(spreadFile("_undated"), "2026-11-02T10:00:00Z", ledgerRows("4.3", 1)))
  out <- suppressMessages(elfFitDates("s", "i", destinationPath = d, ELFs = c("4.1", "4.3", "9.9")))
  expect_equal(out$spreadFitDate[out$ELF == "4.1"], as.POSIXct("2026-09-25", tz = "UTC"))
  expect_true(is.na(out$spreadFitDate[out$ELF == "4.3"]))       # no date anywhere: no date, though it has a fit
  expect_identical(out$spreadLedger[out$ELF == "4.3"], spreadFile("_undated"))
  expect_identical(out$ELF, c("4.1", "4.3", "9.9"))
  expect_true(all(is.na(out[out$ELF == "9.9", -1])))
})

test_that("elfFitDates() with no ledger files gives a table of NA", {
  d <- withr::local_tempdir()
  fakeDrive(list("somethingElse.rds" = list(modified = "2026-11-01T10:00:00Z", ledger = ledgerRows("4.1", 1))))
  out <- suppressMessages(elfFitDates("s", "i", destinationPath = d, ELFs = c("4.1", "4.3")))
  expect_identical(out$ELF, c("4.1", "4.3"))
  expect_true(all(is.na(out$spreadFitDate)) && all(is.na(out$ignitionFitDate)))
  expect_identical(NROW(suppressMessages(elfFitDates("s", "i", destinationPath = d))), 0L)
})

test_that(".ledgerRowDate() is writtenAt, else the date in the file name, else NA", {
  w <- as.numeric(as.POSIXct("2026-11-10 09:30:00", tz = "UTC"))
  f <- spreadFile("_x_20260925")
  expect_equal(.ledgerRowDate(c(w, NA, NA), c(f, f, spreadFile("")), spreadFitFileTag),
               c(w, as.numeric(as.POSIXct("2026-09-25", tz = "UTC")), NA))
  expect_length(.ledgerRowDate(numeric(0), character(0), spreadFitFileTag), 0L)
})
