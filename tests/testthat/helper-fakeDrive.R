## A fake Drive folder: `files` maps a file name to list(modified, ledger). The download is
## reproducible::preProcess(): the mock writes the ledger to `destinationPath` and records the call.
## googledrive::drive_download() must never be called: it writes in place, so a job reading the same
## file at that moment can see it missing or half-written.
fakeDrive <- function(files, env = parent.frame()) {
  store <- new.env()
  store$downloads <- character()
  store$calls <- list()
  ## the bytes saveRDS writes are what the MD5 is taken of, so keep them
  bytes <- lapply(files, function(f) {
    tf <- tempfile(fileext = ".rds"); on.exit(unlink(tf))
    saveRDS(f$ledger, tf); readBin(tf, "raw", file.size(tf))
  })
  ls <- data.frame(name = names(files), id = paste0("id_", seq_along(files)))
  ls$drive_resource <- lapply(names(files), function(n) {
    tf <- tempfile(); on.exit(unlink(tf)); writeBin(bytes[[n]], tf)
    list(modifiedTime = files[[n]]$modified, md5Checksum = unname(tools::md5sum(tf)))
  })
  testthat::local_mocked_bindings(
    drive_ls = function(path, ...) ls,
    drive_download = function(...) stop("googledrive::drive_download() must not be called"),
    .package = "googledrive", .env = env)
  testthat::local_mocked_bindings(
    preProcess = function(targetFile = NULL, url = NULL, destinationPath = ".", purge = FALSE, ...) {
      store$downloads <- c(store$downloads, targetFile)
      store$calls[[length(store$calls) + 1L]] <- list(targetFile = targetFile, url = url,
                                                     destinationPath = destinationPath, purge = purge)
      writeBin(bytes[[targetFile]], file.path(destinationPath, targetFile))
      list(targetFilePath = file.path(destinationPath, targetFile))
    },
    .package = "reproducible", .env = env)
  store
}

ledgerRows <- function(ids, value) {
  geom <- sf::st_sfc(lapply(seq_along(ids), function(i) sf::st_polygon(list(rbind(
    c(i, 0), c(i + 1, 0), c(i + 1, 1), c(i, 1), c(i, 0))))))
  df <- data.frame(polygonID = ids, objFunVal = value)
  df$geometry <- geom
  df
}
