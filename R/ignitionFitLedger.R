#' The shared IgnitionFit ledger files
#'
#' fireSense_IgnitionFit writes each fit (ignition and escape together) as a row of an `rds`
#' file in a Google Drive folder (its `ignitionFitGoogleDriveFolder`); fireSense_ELFs reads them
#' for a multi-ELF study area. With `ignitionFitFilename = "latest"` those modules use these
#' functions instead of one named file. Parallel to the spread-fit ledger helpers
#' (`spreadFitFilenameFor()`, `latestSpreadFits()`).
#'
#' `ignitionFitFileTag` marks the files whose fits use the current model: xgboost ignition and
#' escape fits (`fireSense_IgnitionFit`'s `modelAlgorithm = "xgboost"`). A fit made under a new
#' model needs a new tag, so the old files drop out of `"latest"`.
#'
#' @param fireYears The fire years the fit used; only the first and last matter.
#'
#' @return `ignitionFitFilenameFor()`: the file name a fit over `fireYears` is written to,
#'   e.g. `"fireSenseIgnitionParams_1985-2024_xgboost.rds"`.
#' @export
#' @rdname ignitionFitLedger
ignitionFitFilenameFor <- function(fireYears) {
  fireYears <- as.integer(fireYears[is.finite(fireYears)])
  if (!length(fireYears))
    stop("ignitionFitFilenameFor(): `fireYears` has no years")
  paste0(ignitionFitFilePrefix, min(fireYears), "-", max(fireYears), ignitionFitFileTag, ".rds")
}

ignitionFitFilePrefix <- "fireSenseIgnitionParams_"

#' @export
#' @rdname ignitionFitLedger
ignitionFitFileTag <- "_xgboost"

#' @description
#' `latestIgnitionFits()` returns, for every polygon, its rows from the most recently modified
#' ledger file that has it. Files are read newest first; with `polygonIDs`, reading stops once
#' all of them are found. A file already in `destinationPath` with the same MD5 as on Drive is not
#' downloaded again; one that differs is downloaded again with [reproducible::preProcess()].
#'
#' @param cloudFolderID The Google Drive folder (url or id) holding the ledger files.
#' @param destinationPath Local folder for the downloaded files.
#' @param polygonIDs Optional. The polygons (ELF ids) wanted; `NULL` reads every file.
#'
#' @return `latestIgnitionFits()`: the ledger rows, in the form fireSense_IgnitionFit writes them
#'   (a `data.frame` with a `geometry` column), or `NULL` when no file matches. Attribute
#'   `"ignitionFitFiles"` names the file each polygon came from.
#' @export
#' @rdname ignitionFitLedger
latestIgnitionFits <- function(cloudFolderID, destinationPath, polygonIDs = NULL) {
  .latestLedgerFits(cloudFolderID, destinationPath, polygonIDs, filePrefix = ignitionFitFilePrefix,
                    fileTag = ignitionFitFileTag, caller = "latestIgnitionFits", attrName = "ignitionFitFiles")
}
