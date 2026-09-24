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
  expect_identical(spreadFitFilenameFor(1985:2024), "fireSenseParams_1985-2024_linearFuel.rds")
  expect_identical(spreadFitFilenameFor(c(NA, 2001, 1990)), "fireSenseParams_1990-2001_linearFuel.rds")
  expect_error(spreadFitFilenameFor(NA), "no years")
})

test_that("each polygon comes from the newest file that has it; other models' files are ignored", {
  d <- withr::local_tempdir()
  drive <- fakeDrive(list(
    "fireSenseParams_1985-2024_linearFuel.rds" = list(modified = "2026-09-20T10:00:00Z",
                                                      ledger = ledgerRows(c("4.1", "4.3"), c(1, 1))),
    "fireSenseParams_1985-2025_linearFuel.rds" = list(modified = "2026-09-24T10:00:00Z",
                                                      ledger = ledgerRows("4.1", 2)),
    ## log fuel: newest of all, and has 5.1, but must never be read
    "fireSenseParams_1985-2024.rds" = list(modified = "2026-09-25T10:00:00Z",
                                           ledger = ledgerRows(c("4.1", "5.1"), c(9, 9)))))
  out <- latestSpreadFits("folder", d)
  expect_setequal(out$polygonID, c("4.1", "4.3"))
  expect_identical(out$objFunVal[out$polygonID == "4.1"], 2)      # the newer file's fit
  expect_identical(out$objFunVal[out$polygonID == "4.3"], 1)
  expect_identical(attr(out, "spreadFitFiles")[["4.3"]], "fireSenseParams_1985-2024_linearFuel.rds")
  expect_false("fireSenseParams_1985-2024.rds" %in% drive$downloads)
  expect_s3_class(sf::st_as_sf(out), "sf")                        # still a ledger CacheGeo can read
})

test_that("reading stops at the first file that has every polygon asked for, and reuses local copies", {
  d <- withr::local_tempdir()
  drive <- fakeDrive(list(
    "fireSenseParams_1985-2024_linearFuel.rds" = list(modified = "2026-09-20T10:00:00Z",
                                                      ledger = ledgerRows("4.3", 1)),
    "fireSenseParams_1985-2025_linearFuel.rds" = list(modified = "2026-09-24T10:00:00Z",
                                                      ledger = ledgerRows("4.1", 2))))
  out <- latestSpreadFits("folder", d, polygonIDs = "4.1")
  expect_identical(out$polygonID, "4.1")
  expect_identical(drive$downloads, "fireSenseParams_1985-2025_linearFuel.rds")
  ## second call: the same file is on disk with Drive's MD5, so nothing is downloaded
  latestSpreadFits("folder", d, polygonIDs = "4.1")
  expect_identical(drive$downloads, "fireSenseParams_1985-2025_linearFuel.rds")
  ## an ELF in no file: every file is read, and the result has only what exists
  out <- latestSpreadFits("folder", d, polygonIDs = "9.9")
  expect_setequal(out$polygonID, c("4.1", "4.3"))
})

test_that("no matching file gives NULL", {
  d <- withr::local_tempdir()
  fakeDrive(list("fireSenseParams.rds" = list(modified = "2026-09-20T10:00:00Z", ledger = ledgerRows("4.1", 1))))
  expect_null(suppressMessages(latestSpreadFits("folder", d)))
})
