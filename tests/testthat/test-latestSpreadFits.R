## A fake Drive folder: `files` maps a file name to list(modified, ledger). drive_download writes the
## ledger to `path` and counts the call.
fakeDrive <- function(files, env = parent.frame()) {
  store <- new.env()
  store$downloads <- character()
  ## the bytes saveRDS writes are what the MD5 is taken of, so keep them
  bytes <- lapply(files, function(f) {
    tf <- tempfile(fileext = ".rds"); on.exit(unlink(tf))
    saveRDS(f$ledger, tf); readBin(tf, "raw", file.size(tf))
  })
  ls <- data.frame(name = names(files))
  ls$drive_resource <- lapply(names(files), function(n) {
    tf <- tempfile(); on.exit(unlink(tf)); writeBin(bytes[[n]], tf)
    list(modifiedTime = files[[n]]$modified, md5Checksum = unname(tools::md5sum(tf)))
  })
  testthat::local_mocked_bindings(
    drive_ls = function(path, ...) ls,
    drive_download = function(file, path, overwrite, ...) {
      store$downloads <- c(store$downloads, file$name)
      writeBin(bytes[[file$name]], path)
      invisible(file)
    },
    .package = "googledrive", .env = env)
  store
}

ledgerRows <- function(ids, value) {
  geom <- sf::st_sfc(lapply(seq_along(ids), function(i) sf::st_polygon(list(rbind(
    c(i, 0), c(i + 1, 0), c(i + 1, 1), c(i, 1), c(i, 0))))))
  df <- data.frame(polygonID = ids, objFunVal = value)
  df$geometry <- geom
  df
}

test_that("spreadFitFilenameFor() names a fit by its fire years and the model tag", {
  expect_identical(spreadFitFilenameFor(1985:2024), "fireSenseParams_1985-2024_linearFuel_esc50.rds")
  expect_identical(spreadFitFilenameFor(c(NA, 2001, 1990)), "fireSenseParams_1990-2001_linearFuel_esc50.rds")
  expect_error(spreadFitFilenameFor(NA), "no years")
})

test_that("each polygon comes from the newest file that has it; other models' files are ignored", {
  d <- withr::local_tempdir()
  drive <- fakeDrive(list(
    "fireSenseParams_1985-2024_linearFuel_esc50.rds" = list(modified = "2026-09-20T10:00:00Z",
                                                      ledger = ledgerRows(c("4.1", "4.3"), c(1, 1))),
    "fireSenseParams_1985-2025_linearFuel_esc50.rds" = list(modified = "2026-09-24T10:00:00Z",
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
  expect_identical(attr(out, "spreadFitFiles")[["4.3"]], "fireSenseParams_1985-2024_linearFuel_esc50.rds")
  expect_false("fireSenseParams_1985-2024.rds" %in% drive$downloads)
  expect_false("fireSenseParams_1985-2024_linearFuel.rds" %in% drive$downloads)
  expect_s3_class(sf::st_as_sf(out), "sf")                        # still a ledger CacheGeo can read
})

test_that("reading stops at the first file that has every polygon asked for, and reuses local copies", {
  d <- withr::local_tempdir()
  drive <- fakeDrive(list(
    "fireSenseParams_1985-2024_linearFuel_esc50.rds" = list(modified = "2026-09-20T10:00:00Z",
                                                      ledger = ledgerRows("4.3", 1)),
    "fireSenseParams_1985-2025_linearFuel_esc50.rds" = list(modified = "2026-09-24T10:00:00Z",
                                                      ledger = ledgerRows("4.1", 2))))
  out <- latestSpreadFits("folder", d, polygonIDs = "4.1")
  expect_identical(out$polygonID, "4.1")
  expect_identical(drive$downloads, "fireSenseParams_1985-2025_linearFuel_esc50.rds")
  ## second call: the same file is on disk with Drive's MD5, so nothing is downloaded
  latestSpreadFits("folder", d, polygonIDs = "4.1")
  expect_identical(drive$downloads, "fireSenseParams_1985-2025_linearFuel_esc50.rds")
  ## an ELF in no file: every file is read, and the result has only what exists
  out <- latestSpreadFits("folder", d, polygonIDs = "9.9")
  expect_setequal(out$polygonID, c("4.1", "4.3"))
})

test_that("no matching file gives NULL", {
  d <- withr::local_tempdir()
  fakeDrive(list("fireSenseParams.rds" = list(modified = "2026-09-20T10:00:00Z", ledger = ledgerRows("4.1", 1))))
  expect_null(suppressMessages(latestSpreadFits("folder", d)))
})

test_that("driveDownloadAtomic() never shows a reader a partial file, and ends with the new content", {
  d <- withr::local_tempdir()
  path <- finalPath <- file.path(d, "fireSenseParams_1985-2024_linearFuel_esc50.rds")
  ledgerOld <- ledgerRows("4.1", 1)
  ledgerNew <- ledgerRows(c("4.1", "4.3"), c(2, 2))
  bytesOf <- function(x) { tf <- tempfile(); on.exit(unlink(tf)); saveRDS(x, tf); readBin(tf, "raw", file.size(tf)) }
  newBytes <- bytesOf(ledgerNew)
  writeBin(bytesOf(ledgerOld), path)                       # the file other jobs are reading
  seen <- list()
  testthat::local_mocked_bindings(
    ## writes in two parts; a concurrent reader looks at `path` between them
    drive_download = function(file, path, overwrite, ...) {
      half <- length(newBytes) %/% 2
      con <- file(path, "wb"); on.exit(close(con))
      writeBin(newBytes[seq_len(half)], con); flush(con)
      seen[[length(seen) + 1]] <<- readRDS(finalPath)   # must be complete
      writeBin(newBytes[-seq_len(half)], con)
      invisible(file)
    },
    .package = "googledrive")
  tf <- tempfile(); writeBin(newBytes, tf)
  dribble <- data.frame(name = basename(path))
  dribble$drive_resource <- list(list(md5Checksum = unname(tools::md5sum(tf))))
  expect_true(driveDownloadAtomic(dribble, path))
  expect_length(seen, 1L)
  expect_identical(seen[[1]]$objFunVal, 1)                 # the reader saw the old, complete file
  expect_identical(readRDS(path)$objFunVal, c(2, 2))       # and the final file is the new one
  expect_identical(list.files(d, all.files = TRUE, no.. = TRUE), basename(path))  # no temp left
  ## already current: no second download
  expect_false(driveDownloadAtomic(dribble, path))
})

test_that("latestSpreadFits() replacing a stale local ledger never leaves it partial for another reader", {
  d <- withr::local_tempdir()
  fname <- "fireSenseParams_1985-2024_linearFuel_esc50.rds"
  finalPath <- file.path(d, fname)
  bytesOf <- function(x) { tf <- tempfile(); on.exit(unlink(tf)); saveRDS(x, tf); readBin(tf, "raw", file.size(tf)) }
  newBytes <- bytesOf(ledgerRows(c("4.1", "4.3"), c(2, 2)))
  writeBin(bytesOf(ledgerRows("4.1", 1)), finalPath)      # stale copy another job is reading
  tf <- tempfile(); writeBin(newBytes, tf)
  ls <- data.frame(name = fname)
  ls$drive_resource <- list(list(modifiedTime = "2026-09-20T10:00:00Z", md5Checksum = unname(tools::md5sum(tf))))
  seen <- list()
  testthat::local_mocked_bindings(
    drive_ls = function(path, ...) ls,
    drive_download = function(file, path, overwrite, ...) {
      half <- length(newBytes) %/% 2
      con <- file(path, "wb"); on.exit(close(con))
      writeBin(newBytes[seq_len(half)], con); flush(con)
      seen[[length(seen) + 1]] <<- tryCatch(readRDS(finalPath), error = function(e) e)
      writeBin(newBytes[-seq_len(half)], con)
      invisible(file)
    },
    .package = "googledrive")
  out <- latestSpreadFits("folder", d)
  expect_length(seen, 1L)
  expect_false(inherits(seen[[1]], "error"))
  expect_setequal(out$polygonID, c("4.1", "4.3"))
})
