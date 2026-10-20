test_that("spreadFitFilenameFor() names a fit by its fire years and the model tag", {
  expect_identical(spreadFitFilenameFor(1985:2024), "fireSenseParams_1985-2024_linearFuel_esc50_20260925.rds")
  expect_identical(spreadFitFilenameFor(c(NA, 2001, 1990)), "fireSenseParams_1990-2001_linearFuel_esc50_20260925.rds")
  expect_error(spreadFitFilenameFor(NA), "no years")
})

test_that("undated and suffixed experiment files are not ledger files; the dated tag is", {
  d <- withr::local_tempdir()
  drive <- fakeDrive(list(
    "fireSenseParams_1985-2024_linearFuel_esc50.rds" = list(modified = "2026-10-12T10:00:00Z",
                                                      ledger = ledgerRows(c("4.1", "5.1"), c(9, 9))),
    "fireSenseParams_1985-2024_linearFuel_esc50_20260925_repStab_rep3.rds" = list(modified = "2026-10-12T11:00:00Z",
                                                      ledger = ledgerRows(c("4.1", "6.1"), c(8, 8))),
    "fireSenseParams_1985-2024_linearFuel_esc50_20260925.rds" = list(modified = "2026-09-01T10:00:00Z",
                                                      ledger = ledgerRows("4.1", 1))))
  out <- latestSpreadFits("folder", d)
  expect_identical(out$polygonID, "4.1")
  expect_identical(out$objFunVal, 1)
  expect_identical(drive$downloads, "fireSenseParams_1985-2024_linearFuel_esc50_20260925.rds")
  oldPattern <- "^fireSenseParams_.*_linearFuel_esc50\\.rds$"
  expect_false(grepl(oldPattern, spreadFitFilenameFor(1985:2024)))
})

test_that("each polygon comes from the newest file that has it; other models' files are ignored", {
  d <- withr::local_tempdir()
  drive <- fakeDrive(list(
    "fireSenseParams_1985-2024_linearFuel_esc50_20260925.rds" = list(modified = "2026-09-20T10:00:00Z",
                                                      ledger = ledgerRows(c("4.1", "4.3"), c(1, 1))),
    "fireSenseParams_1985-2025_linearFuel_esc50_20260925.rds" = list(modified = "2026-09-24T10:00:00Z",
                                                      ledger = ledgerRows("4.1", 2)),
    ## log fuel: newest of all, and has 5.1, but must never be read
    "fireSenseParams_1985-2024.rds" = list(modified = "2026-09-25T10:00:00Z",
                                           ledger = ledgerRows(c("4.1", "5.1"), c(9, 9))),
    ## the previous model (linear fuel, 1-pixel escape): has 6.1, but must never be read either
    "fireSenseParams_1985-2024_linearFuel.rds" = list(modified = "2026-09-25T11:00:00Z",
                                                      ledger = ledgerRows(c("4.1", "6.1"), c(8, 8)))))
  out <- latestSpreadFits("folder", d)
  expect_setequal(out$polygonID, c("4.1", "4.3"))
  expect_identical(out$objFunVal[out$polygonID == "4.1"], 2)      # the newer file's fit
  expect_identical(out$objFunVal[out$polygonID == "4.3"], 1)
  expect_identical(attr(out, "spreadFitFiles")[["4.3"]], "fireSenseParams_1985-2024_linearFuel_esc50_20260925.rds")
  expect_false("fireSenseParams_1985-2024.rds" %in% drive$downloads)
  expect_false("fireSenseParams_1985-2024_linearFuel.rds" %in% drive$downloads)
  expect_s3_class(sf::st_as_sf(out), "sf")                        # still a ledger CacheGeo can read
})

test_that("reading stops at the first file that has every polygon asked for", {
  d <- withr::local_tempdir()
  drive <- fakeDrive(list(
    "fireSenseParams_1985-2024_linearFuel_esc50_20260925.rds" = list(modified = "2026-09-20T10:00:00Z",
                                                      ledger = ledgerRows("4.3", 1)),
    "fireSenseParams_1985-2025_linearFuel_esc50_20260925.rds" = list(modified = "2026-09-24T10:00:00Z",
                                                      ledger = ledgerRows("4.1", 2))))
  out <- latestSpreadFits("folder", d, polygonIDs = "4.1")
  expect_identical(out$polygonID, "4.1")
  expect_identical(drive$downloads, "fireSenseParams_1985-2025_linearFuel_esc50_20260925.rds")
  ## an ELF in no file: every file is read, and the result has only what exists
  out <- latestSpreadFits("folder", d, polygonIDs = "9.9")
  expect_setequal(out$polygonID, c("4.1", "4.3"))
})

test_that("no matching file gives NULL", {
  d <- withr::local_tempdir()
  fakeDrive(list("fireSenseParams.rds" = list(modified = "2026-09-20T10:00:00Z", ledger = ledgerRows("4.1", 1))))
  expect_null(suppressMessages(latestSpreadFits("folder", d)))
})

test_that("each ledger file is read with reproducible::CacheGeoRead() from the Drive folder into destinationPath", {
  d <- withr::local_tempdir()
  nm <- "fireSenseParams_1985-2025_linearFuel_esc50_20260925.rds"
  drive <- fakeDrive(list(a = list(modified = "2026-09-20T10:00:00Z", ledger = ledgerRows("4.3", 1)),
                          b = list(modified = "2026-09-24T10:00:00Z", ledger = ledgerRows("4.1", 2))) |>
                       setNames(c("fireSenseParams_1985-2024_linearFuel_esc50_20260925.rds", nm)))
  out <- latestSpreadFits("folder", d, polygonIDs = "4.1")
  expect_identical(out$polygonID, "4.1")
  expect_length(drive$calls, 1L)
  ledger <- drive$calls[[1]]
  expect_identical(ledger$file, nm)                               # the newest file
  expect_identical(ledger$destinationPath, d)
  expect_identical(ledger$remote, "folder")
})

test_that("a suffix goes after the model tag; \"\" is the production name", {
  expect_identical(spreadFitFilenameFor(1986:2024, "_boundfix_rep81"),
                   "fireSenseParams_1986-2024_linearFuel_esc50_20260925_boundfix_rep81.rds")
  expect_identical(spreadFitFilenameFor(1986:2024, ""), spreadFitFilenameFor(1986:2024))
})

test_that("\"latest\" ignores suffixed files, and a suffixed run reads only its own file", {
  d <- withr::local_tempdir()
  prod <- spreadFitFilenameFor(1986:2024)
  test81 <- spreadFitFilenameFor(1986:2024, "_boundfix_rep81")
  drive <- fakeDrive(setNames(list(
    list(modified = "2026-09-20T10:00:00Z", ledger = ledgerRows("4.1", 1)),
    list(modified = "2026-10-15T10:00:00Z", ledger = ledgerRows(c("4.1", "5.1"), c(81, 81)))),
    c(prod, test81)))
  out <- latestSpreadFits("folder", d)
  expect_identical(out$polygonID, "4.1")
  expect_identical(out$objFunVal, 1)
  expect_identical(drive$downloads, prod)
  out <- latestSpreadFits("folder", d, suffix = "_boundfix_rep81")
  expect_setequal(out$polygonID, c("4.1", "5.1"))
  expect_identical(unique(out$objFunVal), 81)
})

## Drive unreachable, as on a compute node with no internet: drive_ls() fails to resolve the host
goOffline <- function(env = parent.frame()) {
  testthat::local_mocked_bindings(
    drive_ls = function(path, ...) stop("Could not resolve host: www.googleapis.com", call. = FALSE),
    .package = "googledrive", .env = env)
}

test_that("offline, the saved file list gives the same files and rows", {
  d <- withr::local_tempdir()
  drive <- fakeDrive(list(
    "fireSenseParams_1985-2024_linearFuel_esc50_20260925.rds" = list(modified = "2026-09-20T10:00:00Z",
                                                      ledger = ledgerRows(c("4.1", "4.3"), c(1, 1))),
    "fireSenseParams_1985-2025_linearFuel_esc50_20260925.rds" = list(modified = "2026-09-24T10:00:00Z",
                                                      ledger = ledgerRows("4.1", 2))))
  online <- latestSpreadFits("folder", d)
  expect_true(file.exists(file.path(d, ".ledgerListing_folder.rds")))
  goOffline()
  expect_message(offline <- latestSpreadFits("folder", d), "cannot be reached; using the file list saved")
  expect_identical(offline, online)
  expect_identical(drive$downloads[3:4], drive$downloads[1:2])    # the same files, in the same order
})

test_that("offline with no saved file list, the connection error is raised", {
  d <- withr::local_tempdir()
  goOffline()
  expect_error(latestSpreadFits("folder", d), "resolve host")
})

test_that("a Drive error that is not a connection failure is raised", {
  d <- withr::local_tempdir()
  testthat::local_mocked_bindings(drive_ls = function(path, ...) stop("Client error: (404) Not Found"),
                                  .package = "googledrive")
  expect_error(latestSpreadFits("folder", d), "404")
})
