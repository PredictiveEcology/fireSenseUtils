#' The shared IgnitionFit ledger files
#'
#' fireSense_IgnitionFit writes each fit (ignition and escape together) as a row of an `rds`
#' file in a Google Drive folder (its `ignitionFitGoogleDriveFolder`); fireSense_ELFs reads them
#' for a multi-ELF study area. With `ignitionFitFilename = "latest"` those modules use these
#' functions instead of one named file. Parallel to the spread-fit ledger helpers
#' (`spreadFitFilenameFor()`, `latestSpreadFits()`).
#'
#' `ignitionFitFileTag` marks the files whose fits use the current model: xgboost Poisson
#' ignition and escape fits (`fireSense_IgnitionFit`'s `modelAlgorithm = "xgboost"` with objective
#' `count:poisson`). The date is the model version. A fit made under a new model needs a new tag,
#' so the old files drop out of `"latest"`; files with an older tag, e.g. `"_xgboost"` (Tweedie
#' fits), are never read.
#' The tag names the objective ([ignitionObjectiveTxt]); a test fails if the objective changes and the tag does not.
#'
#' @param fireYears The fire years the fit used; only the first and last matter.
#' @param suffix Text put after the model tag in the file name, e.g. `"_bufIndep20261024"`; end it with
#'   the date of the fits (`yyyymmdd`), which orders the files for `"latest"`. A run with a suffix writes,
#'   and `"latest"` reads, only files with that suffix. With `""` (the default) `"latest"` reads every file
#'   of the model tag in the folder, newest date first. Experiment fits go in the folder's `"testing"`
#'   subfolder (`spreadFitLedgerFolder(testing = TRUE)`), which `"latest"` over the main folder never reads. The same `fitFileSuffix` parameter
#'   of `fireSense_spreadFit` and `fireSense_ignitionFit` sets it for both ledgers.
#'
#' @return `ignitionFitFilenameFor()`: the file name a fit over `fireYears` is written to,
#'   e.g. `"fireSenseIgnitionParams_1985-2024_xgboostPoisson_20261013.rds"` (then `suffix`).
#' @export
#' @rdname ignitionFitLedger
ignitionFitFilenameFor <- function(fireYears, suffix = "") {
  fireYears <- as.integer(fireYears[is.finite(fireYears)])
  if (!length(fireYears))
    stop("ignitionFitFilenameFor(): `fireYears` has no years")
  paste0(ignitionFitFilePrefix, min(fireYears), "-", max(fireYears), ignitionFitFileTag, suffix, ".rds")
}

ignitionFitFilePrefix <- "fireSenseIgnitionParams_"

#' @export
#' @rdname ignitionFitLedger
ignitionFitFileTag <- "_xgboostPoisson_20261013"

#' @description
#' `latestIgnitionFits()` returns, for every polygon, its newest row in the ledger files: the one with
#' the latest `writtenAt` (stamped by [reproducible::CacheGeoWrite()]). A row written before that stamp
#' existed has its file's date instead: the date in the file's suffix ([ledgerFileOrder()]), else older
#' than any dated file; ties go by the file's modified time. No row is newer than its file was last
#' modified, so files are read newest-modified first; with `polygonIDs`, reading stops once each has a
#' row newer than the next file. Files are read, also with no internet, as [latestSpreadFits()] reads them.
#'
#' @param cloudFolderID The Google Drive folder (url or id) holding the ledger files; default [ignitionFitLedgerFolder()].
#' @param destinationPath Local folder for the downloaded files.
#' @param polygonIDs Optional. The polygons (ELF ids) wanted; `NULL` reads every file.
#' @param suffix A suffix reads only the files named with it; `""` (default) reads every file of the
#'   model tag, whatever its suffix. See `ignitionFitFilenameFor()`.
#'
#' @return `latestIgnitionFits()`: the ledger rows, in the form fireSense_IgnitionFit writes them
#'   (a `data.frame` with a `geometry` column), or `NULL` when no file matches. Attribute
#'   `"ignitionFitFiles"` names the file each polygon came from.
#' @export
#' @rdname ignitionFitLedger
latestIgnitionFits <- function(cloudFolderID = ignitionFitLedgerFolder(), destinationPath, polygonIDs = NULL,
                              suffix = "") {
  .latestLedgerFits(cloudFolderID, destinationPath, polygonIDs, filePrefix = ignitionFitFilePrefix,
                    fileTag = ignitionFitFileTag, suffix = suffix, caller = "latestIgnitionFits", attrName = "ignitionFitFiles")
}
