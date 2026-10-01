test_that("ignitionFitFilenameFor() names a fit by its fire years and the model tag", {
  expect_identical(ignitionFitFilenameFor(1985:2024), "fireSenseIgnitionParams_1985-2024_xgboost.rds")
  expect_identical(ignitionFitFilenameFor(c(NA, 2001, 1990)), "fireSenseIgnitionParams_1990-2001_xgboost.rds")
  expect_error(ignitionFitFilenameFor(NA), "no years")
})

test_that("each polygon comes from the newest file that has it; other models' files are ignored", {
  d <- withr::local_tempdir()
  drive <- fakeDrive(list(
    "fireSenseIgnitionParams_1985-2024_xgboost.rds" = list(modified = "2026-09-20T10:00:00Z",
                                                      ledger = ledgerRows(c("4.1", "4.3"), c(1, 1))),
    "fireSenseIgnitionParams_1985-2025_xgboost.rds" = list(modified = "2026-09-24T10:00:00Z",
                                                      ledger = ledgerRows("4.1", 2)),
    ## an older model's file: has 5.1, but must never be read
    "fireSenseIgnitionParams_1985-2024.rds" = list(modified = "2026-09-25T10:00:00Z",
                                           ledger = ledgerRows(c("4.1", "5.1"), c(9, 9)))))
  out <- latestIgnitionFits("folder", d)
  expect_setequal(out$polygonID, c("4.1", "4.3"))
  expect_identical(out$objFunVal[out$polygonID == "4.1"], 2)      # the newer file's fit
  expect_identical(out$objFunVal[out$polygonID == "4.3"], 1)
  expect_identical(attr(out, "ignitionFitFiles")[["4.3"]], "fireSenseIgnitionParams_1985-2024_xgboost.rds")
  expect_false("fireSenseIgnitionParams_1985-2024.rds" %in% drive$downloads)
  expect_s3_class(sf::st_as_sf(out), "sf")                        # still a ledger CacheGeo can read
})

test_that("reading stops at the first file that has every polygon asked for, and reuses local copies", {
  d <- withr::local_tempdir()
  drive <- fakeDrive(list(
    "fireSenseIgnitionParams_1985-2024_xgboost.rds" = list(modified = "2026-09-20T10:00:00Z",
                                                      ledger = ledgerRows("4.3", 1)),
    "fireSenseIgnitionParams_1985-2025_xgboost.rds" = list(modified = "2026-09-24T10:00:00Z",
                                                      ledger = ledgerRows("4.1", 2))))
  out <- latestIgnitionFits("folder", d, polygonIDs = "4.1")
  expect_identical(out$polygonID, "4.1")
  expect_identical(drive$downloads, "fireSenseIgnitionParams_1985-2025_xgboost.rds")
  ## second call: the same file is on disk with Drive's MD5, so nothing is downloaded
  latestIgnitionFits("folder", d, polygonIDs = "4.1")
  expect_identical(drive$downloads, "fireSenseIgnitionParams_1985-2025_xgboost.rds")
  ## an ELF in no file: every file is read, and the result has only what exists
  out <- latestIgnitionFits("folder", d, polygonIDs = "9.9")
  expect_setequal(out$polygonID, c("4.1", "4.3"))
})

test_that("no matching file gives NULL", {
  d <- withr::local_tempdir()
  fakeDrive(list("fireSenseIgnitionParams.rds" = list(modified = "2026-09-20T10:00:00Z", ledger = ledgerRows("4.1", 1))))
  expect_null(suppressMessages(latestIgnitionFits("folder", d)))
})
