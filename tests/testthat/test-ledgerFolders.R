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
