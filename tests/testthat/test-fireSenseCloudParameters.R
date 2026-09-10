## fireSenseCloudParameters() read the shared fitted-parameter file with
## prepInputs(purge = 7, overwrite = TRUE), expecting a fresh copy. It never was:
## prepInputs() returns the copy already on disk whenever it matches CHECKSUMS.txt,
## `purge` only rebuilds CHECKSUMS.txt entries, and `overwrite` only governs the
## written output. A changed file on Drive was never seen. It must download on every
## call.
##
## Drive is stubbed: `drive_download()` writes whatever the "remote" currently holds.

test_that("fireSenseCloudParameters downloads the file on every call", {
  remoteContent <- data.frame(polygonID = "4.1")
  downloads <- 0L
  tmpUsed <- character()
  local_mocked_bindings(
    drive_get = function(id, ...) structure(list(id = id), class = "fakeDribble"),
    is_folder = function(x) FALSE,
    drive_download = function(file, path, ...) {
      downloads <<- downloads + 1L
      tmpUsed <<- c(tmpUsed, path)
      saveRDS(remoteContent, path)
      invisible(file)
    },
    .package = "googledrive"
  )
  dp <- withr::local_tempdir()

  first <- fireSenseCloudParameters(destinationPath = dp)
  expect_identical(first, remoteContent)

  remoteContent <- data.frame(polygonID = c("4.1", "14.3"))   # the file on Drive changes
  second <- fireSenseCloudParameters(destinationPath = dp)

  expect_identical(downloads, 2L)
  expect_identical(second, remoteContent)
  expect_identical(readRDS(file.path(dp, "fireSenseParams.rds")), remoteContent)
  expect_false(any(file.exists(tmpUsed)))
})

test_that("fireSenseCloudParameters finds targetFile when url is a folder", {
  got <- NULL
  local_mocked_bindings(
    drive_get = function(id, ...) "folder",
    is_folder = function(x) identical(x, "folder"),
    drive_ls = function(path, ...) data.frame(name = c("fireSenseParams.rds", "fireSenseParams2026-09.rds"),
                                              id = c("a", "b")),
    drive_download = function(file, path, ...) {
      got <<- file
      saveRDS(data.frame(polygonID = "5.4"), path)
      invisible(file)
    },
    .package = "googledrive"
  )
  dp <- withr::local_tempdir()
  out <- fireSenseCloudParameters(url = "https://drive.google.com/drive/folders/1X9-mRjyLMNpgkP_cfqhbr_AQEPOsVCHf",
                                  targetFile = "fireSenseParams2026-09.rds", destinationPath = dp)
  expect_identical(got$id, "b")
  expect_identical(out$polygonID, "5.4")
  expect_true(file.exists(file.path(dp, "fireSenseParams2026-09.rds")))
})
