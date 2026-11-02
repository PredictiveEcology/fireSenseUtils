shared <- "https://drive.google.com/drive/folders/1X9-mRjyLMNpgkP_cfqhbr_AQEPOsVCHf"

test_that("the ledger folders default to the shared fireSense folder", {
  withr::local_options(fireSenseUtils.spreadFitLedgerFolder = NULL,
                       fireSenseUtils.ignitionFitLedgerFolder = NULL)
  expect_identical(spreadFitLedgerFolder(), shared)
  expect_identical(ignitionFitLedgerFolder(), shared)
})

test_that("the spread option moves both ledgers; the ignition option moves only ignition", {
  withr::local_options(fireSenseUtils.spreadFitLedgerFolder = "spreadId",
                       fireSenseUtils.ignitionFitLedgerFolder = NULL)
  expect_identical(spreadFitLedgerFolder(), "spreadId")
  expect_identical(ignitionFitLedgerFolder(), "spreadId")
  withr::local_options(fireSenseUtils.ignitionFitLedgerFolder = "ignId")
  expect_identical(spreadFitLedgerFolder(), "spreadId")
  expect_identical(ignitionFitLedgerFolder(), "ignId")
})

test_that("an explicit folder wins over the options", {
  withr::local_options(fireSenseUtils.spreadFitLedgerFolder = "spreadId",
                       fireSenseUtils.ignitionFitLedgerFolder = "ignId")
  expect_identical(spreadFitLedgerFolder("x"), "x")
  expect_identical(ignitionFitLedgerFolder("y"), "y")
})

test_that("latestSpreadFits and latestIgnitionFits read the ledger folder by default", {
  withr::local_options(fireSenseUtils.spreadFitLedgerFolder = "spreadId",
                       fireSenseUtils.ignitionFitLedgerFolder = "ignId")
  expect_identical(formals(latestSpreadFits)$cloudFolderID, quote(spreadFitLedgerFolder()))
  expect_identical(formals(latestIgnitionFits)$cloudFolderID, quote(ignitionFitLedgerFolder()))
})

test_that("testing = TRUE gives the ledger folder's \"testing\" subfolder, for experiment fits", {
  withr::local_options(fireSenseUtils.spreadFitLedgerFolderTesting = NULL,
                       fireSenseUtils.ignitionFitLedgerFolderTesting = NULL)
  testing <- "https://drive.google.com/drive/folders/1nqSVuGRYxIvcgZMn7TIbn5HfRiUHR2xq"
  expect_identical(spreadFitLedgerFolder(testing = TRUE), testing)
  expect_identical(ignitionFitLedgerFolder(testing = TRUE), testing)   # one folder holds both ledgers
  expect_false(identical(spreadFitLedgerFolder(testing = TRUE), spreadFitLedgerFolder()))
  withr::local_options(fireSenseUtils.spreadFitLedgerFolderTesting = "tId")
  expect_identical(ignitionFitLedgerFolder(testing = TRUE), "tId")
})

## googledrive::drive_ls() lists a folder's own files only unless `recursive = TRUE`; the "testing"
## subfolder of the ledger folder is a row of type folder in that list, not its contents.
test_that("\"latest\" over the main folder lists it without recursing, so never reads \"testing\"", {
  expect_false(formals(googledrive::drive_ls)$recursive)
  d <- withr::local_tempdir()
  asked <- new.env()
  drive <- fakeDrive(list("fireSenseParams_1986-2024_linearFuel_esc50_20260925.rds" =
                            list(modified = "2026-10-01T10:00:00Z", ledger = ledgerRows("4.1", 1))))
  realLs <- googledrive::drive_ls
  testthat::local_mocked_bindings(drive_ls = function(path, ...) {
    asked$args <- list(...)
    rbind(data.frame(name = "fireSenseParams_1986-2024_linearFuel_esc50_20260925.rds", id = "id_1",
                     drive_resource = I(list(list(modifiedTime = "2026-10-01T10:00:00Z")))),
          data.frame(name = "testing", id = "id_t",
                     drive_resource = I(list(list(modifiedTime = "2026-11-02T10:00:00Z")))))
  }, .package = "googledrive")
  out <- latestSpreadFits("folder", d)
  expect_false(isTRUE(asked$args$recursive))
  expect_identical(out$polygonID, "4.1")
  expect_identical(drive$downloads, "fireSenseParams_1986-2024_linearFuel_esc50_20260925.rds")
})
