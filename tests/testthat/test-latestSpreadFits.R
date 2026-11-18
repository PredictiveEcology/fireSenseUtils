test_that("spreadFitFilenameFor() names a fit by its fire years and the model tag", {
  expect_identical(spreadFitFilenameFor(1985:2024), "fireSenseParams_1985-2024_linearFuel_esc50_20260925.rds")
  expect_identical(spreadFitFilenameFor(c(NA, 2001, 1990)), "fireSenseParams_1990-2001_linearFuel_esc50_20260925.rds")
  expect_error(spreadFitFilenameFor(NA), "no years")
})

test_that("rows that record no fuel scale count only in files of the dated model tag", {
  d <- withr::local_tempdir()
  drive <- fakeDrive(list(
    "fireSenseParams_1985-2024_linearFuel_esc50.rds" = list(modified = "2026-10-12T10:00:00Z",
                                                      ledger = ledgerRows(c("4.1", "5.1"), c(9, 9))),
    "fireSenseParams_1985-2024_linearFuel_esc50_20260925_repStab_rep3.rds" = list(modified = "2026-10-12T11:00:00Z",
                                                      ledger = ledgerRows(c("4.1", "6.1"), c(8, 8))),
    "fireSenseParams_1985-2024_linearFuel_esc50_20260925.rds" = list(modified = "2026-09-01T10:00:00Z",
                                                      ledger = ledgerRows("4.1", 1))))
  expect_message(out <- latestSpreadFits("folder", d),
                 "skipped the fit of 5.1 in fireSenseParams_1985-2024_linearFuel_esc50.rds: its row records no fuel scale")
  ## both have no date in their suffix: the later modified is read first
  expect_setequal(out$polygonID, c("4.1", "6.1"))
  expect_identical(out$objFunVal[out$polygonID == "4.1"], 8)
  expect_false("fireSenseParams_1985-2024_linearFuel_esc50.rds" %in% attr(out, "spreadFitFiles"))
})

test_that("each polygon comes from its newest usable row; other models' rows are skipped", {
  d <- withr::local_tempdir()
  drive <- fakeDrive(list(
    "fireSenseParams_1985-2024_linearFuel_esc50_20260925.rds" = list(modified = "2026-09-20T10:00:00Z",
                                                      ledger = ledgerRows(c("4.1", "4.3"), c(1, 1))),
    "fireSenseParams_1985-2025_linearFuel_esc50_20260925.rds" = list(modified = "2026-09-24T10:00:00Z",
                                                      ledger = ledgerRows("4.1", 2)),
    ## log fuel: newest of all, and has 5.1, but its rows are never used
    "fireSenseParams_1985-2024.rds" = list(modified = "2026-09-25T10:00:00Z",
                                           ledger = ledgerRows(c("4.1", "5.1"), c(9, 9))),
    ## the previous model (linear fuel, 1-pixel escape): has 6.1, but its rows are never used either
    "fireSenseParams_1985-2024_linearFuel.rds" = list(modified = "2026-09-25T11:00:00Z",
                                                      ledger = ledgerRows(c("4.1", "6.1"), c(8, 8)))))
  out <- suppressMessages(latestSpreadFits("folder", d))
  expect_setequal(out$polygonID, c("4.1", "4.3"))
  expect_identical(out$objFunVal[out$polygonID == "4.1"], 2)      # the newer usable file's fit
  expect_identical(out$objFunVal[out$polygonID == "4.3"], 1)
  expect_identical(attr(out, "spreadFitFiles")[["4.3"]], "fireSenseParams_1985-2024_linearFuel_esc50_20260925.rds")
  expect_s3_class(sf::st_as_sf(out), "sf")                        # still a ledger CacheGeo can read
})

test_that("a row's recorded model decides, not its file's name: fuel scale and escape definition", {
  ## the file-name tag was kept by hand; a row now says what model it is (spreadFuelScale, escapeDefinition())
  d <- withr::local_tempdir()
  rows <- function(ids, value, scale, size = 50, writtenAt) {
    r <- ledgerRows(ids, value, writtenAt)
    r[[spreadFitFuelScaleTxt]] <- scale
    cbind(r, as.data.frame(escapeDefinition(size)))
  }
  fakeDrive(list(
    "fireSenseParams_1985-2024_linearFuel_esc50_20260925.rds" = list(modified = "2026-10-01T10:00:00Z",
                                                      ledger = rows(c("4.1", "5.1", "6.1"), 1, "linear", writtenAt = "2026-10-01")),
    ## a new label, no tag: usable, because its rows record the model
    "fireSenseParams_myExperiment.rds" = list(modified = "2026-11-01T10:00:00Z",
                                                      ledger = rows("4.1", 2, "linear", writtenAt = "2026-11-01")),
    "fireSenseParams_logFuelTrial.rds" = list(modified = "2026-11-02T10:00:00Z",
                                                      ledger = rows("5.1", 3, "log", writtenAt = "2026-11-02")),
    "fireSenseParams_esc100.rds" = list(modified = "2026-11-03T10:00:00Z",
                                                      ledger = rows("6.1", 4, "linear", size = 100, writtenAt = "2026-11-03"))))
  msgs <- testthat::capture_messages(out <- latestSpreadFits("folder", d, escapeDefinition = escapeDefinition(50)))
  val <- function(id) out$objFunVal[out$polygonID == id]
  expect_identical(val("4.1"), 2)
  expect_identical(val("5.1"), 1)       # the newer log-fuel fit is skipped
  expect_identical(val("6.1"), 1)       # the newer 100 ha fit is skipped
  expect_true(any(grepl("5.1 in fireSenseParams_logFuelTrial.rds: its fuel covariates are on the log scale", msgs)))
  expect_true(any(grepl("6.1 in fireSenseParams_esc100.rds: the escape definition", msgs)))
  ## without a definition to compare, the escape size is not checked
  unchecked <- suppressMessages(latestSpreadFits("folder", d))
  expect_identical(unchecked$objFunVal[unchecked$polygonID == "6.1"], 4)
})

test_that("reading stops once each polygon asked for has a row newer than the next file", {
  d <- withr::local_tempdir()
  drive <- fakeDrive(list(
    "fireSenseParams_1985-2024_linearFuel_esc50_20260925.rds" = list(modified = "2026-09-20T10:00:00Z",
                                                      ledger = ledgerRows("4.3", 1, "2026-09-20 09:00")),
    "fireSenseParams_1985-2025_linearFuel_esc50_20260925.rds" = list(modified = "2026-09-24T10:00:00Z",
                                                      ledger = ledgerRows("4.1", 2, "2026-09-24 09:00"))))
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
  drive <- fakeDrive(list(a = list(modified = "2026-09-20T10:00:00Z", ledger = ledgerRows("4.3", 1, "2026-09-20")),
                          b = list(modified = "2026-09-24T10:00:00Z", ledger = ledgerRows("4.1", 2, "2026-09-24"))) |>
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

test_that("\"latest\" reads every suffix, newest date in the suffix first; a suffix reads only its own files", {
  d <- withr::local_tempdir()
  prod <- spreadFitFilenameFor(1986:2024)                                  # no suffix: no date, oldest
  test81 <- spreadFitFilenameFor(1986:2024, "_boundfix_rep81")             # no date in the suffix either
  gpsc <- spreadFitFilenameFor(1986:2024, "_gpscCause_20261029")
  bufI <- spreadFitFilenameFor(1986:2024, "_bufIndep20261024")
  drive <- fakeDrive(setNames(list(
    list(modified = "2026-11-01T10:00:00Z", ledger = ledgerRows(c("4.1", "7.1"), c(1, 1))),
    list(modified = "2026-10-15T10:00:00Z", ledger = ledgerRows(c("4.1", "5.1", "7.1"), c(81, 81, 81))),
    list(modified = "2026-10-29T10:00:00Z", ledger = ledgerRows("4.1", 29)),
    ## modified after the gpsc file (a fit added later), but its suffix date is older
    list(modified = "2026-11-02T10:00:00Z", ledger = ledgerRows(c("4.1", "5.1"), c(24, 24)))),
    c(prod, test81, gpsc, bufI)))
  out <- latestSpreadFits("folder", d)
  val <- function(id) out$objFunVal[out$polygonID == id]
  expect_identical(val("4.1"), 29)      # 20261029 beats 20261024, whatever the modified times
  expect_identical(val("5.1"), 24)
  expect_identical(val("7.1"), 1)       # undated files: the later modified first
  ## files are read by modified time (no row is newer than its file), and the date in the suffix
  ## when that is later; no polygonIDs: every file
  expect_identical(drive$downloads, c(bufI, prod, gpsc, test81))
  ## reading stops once every ELF asked for is found
  drive$downloads <- character()
  out <- latestSpreadFits("folder", d, polygonIDs = c("4.1", "5.1"))
  ## gpsc's ceiling (its 2026-10-29 modified time) is after 4.1's row in bufI (20261024), so it is
  ## read; test81's (2026-10-15) is before every row kept, so reading stops there
  expect_identical(drive$downloads, c(bufI, prod, gpsc))
  ## a suffix reads only its own files
  out <- latestSpreadFits("folder", d, suffix = "_boundfix_rep81")
  expect_setequal(out$polygonID, c("4.1", "5.1", "7.1"))
  expect_identical(unique(out$objFunVal), 81)
})

## A file's date is when its campaign was named, not when each of its fits was written: a fit added
## later to an older campaign's file lost to an older fit in a newer-named file. reproducible's
## CacheGeoWrite() now stamps each row with writtenAt.
test_that("each polygon's newest row by writtenAt wins, whatever its file's name; unstamped rows by the file's date", {
  d <- withr::local_tempdir()
  old <- spreadFitFilenameFor(1986:2024, "_crps_20261028")
  new <- spreadFitFilenameFor(1986:2024, "_gpscCause_20261101")
  legacy <- spreadFitFilenameFor(1986:2024, "_prod20261021")
  drive <- fakeDrive(setNames(list(
    ## an older campaign's file, with a refit of 4.1 added on 2026-11-05
    list(modified = "2026-11-05T10:00:00Z",
         ledger = rbind(ledgerRows("4.1", 5, "2026-11-05 09:00"), ledgerRows("5.1", 28, "2026-10-28 09:00"))),
    list(modified = "2026-11-01T10:00:00Z", ledger = ledgerRows(c("4.1", "5.1"), c(1, 1), "2026-11-01 09:00")),
    ## written before the stamp: its rows date from 20261021
    list(modified = "2026-10-21T10:00:00Z", ledger = ledgerRows(c("4.1", "5.1", "7.1"), c(21, 21, 21)))),
    c(old, new, legacy)))
  out <- latestSpreadFits("folder", d)
  val <- function(id) out$objFunVal[out$polygonID == id]
  expect_identical(val("4.1"), 5)       # the refit, though its file's name is older
  expect_identical(val("5.1"), 1)       # 2026-11-01 beats 2026-10-28, in the same two files
  expect_identical(val("7.1"), 21)      # only in the unstamped file
  expect_identical(unname(attr(out, "spreadFitFiles")[c("4.1", "5.1", "7.1")]), c(old, new, legacy))

  ## asked for 4.1 only: its row (2026-11-05) is newer than the next file, so one file is read
  drive$downloads <- character()
  out <- latestSpreadFits("folder", d, polygonIDs = "4.1")
  expect_identical(drive$downloads, old)
  expect_identical(val("4.1"), 5)
  ## 7.1's row has no stamp and the date 20261021, before the next file's modified time: all are read
  drive$downloads <- character()
  latestSpreadFits("folder", d, polygonIDs = "7.1")
  expect_identical(drive$downloads, c(old, new, legacy))
})

test_that("ledgerFileOrder(): the last 8-digit date of the suffix, undated last, then modified time", {
  n <- paste0("fireSenseParams_1986-2024_linearFuel_esc50_20260925", c("", "_a_20261101", "_x20261024", "_rep3", "_D05_20261102"), ".rds")
  m <- c("2026-11-03", "2026-11-01", "2026-11-04", "2026-11-05", "2026-11-02")
  expect_identical(n[ledgerFileOrder(n, m, "_linearFuel_esc50_20260925")], n[c(5, 2, 3, 4, 1)])
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
