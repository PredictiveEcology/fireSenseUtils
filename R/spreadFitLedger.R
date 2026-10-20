#' The shared SpreadFit ledger files
#'
#' fireSense_spreadFit writes each fit as a row of an `rds` file in a Google Drive folder
#' (its `spreadFitGoogleDriveFolder`); fireSense_ELFs and fireSense_dataPrepFit read them.
#' With `spreadFitFilename = "latest"` those modules use these functions instead of one named file.
#'
#' `spreadFitFileTag` marks the files whose fits use the current model: fuel biomass on the
#' linear scale ([fuelLogToLinear()]), escaped fires starting at 50 ha (`escapeSizeHa`), and the
#' objective with the annual-area and area-distribution terms. The date (2026-09-25) is when
#' this model started, i.e. the model version. Files with an older tag (`"_linearFuel_esc50"`:
#' the same model before it was dated; `"_linearFuel"`: 1-pixel escape; none: log fuel) hold fits of an earlier model, which
#' fireSense_spreadPredict would apply wrongly, so `"latest"` never reads them.
#' A fit made under a new model needs a new tag, so the old files drop out of `"latest"`.
#'
#' @param fireYears The fire years the fit used; only the first and last matter.
#' @param suffix Text put after the model tag in the file name, e.g. `"_boundfix_rep81"`; `""` (the
#'   default) is the production file. A test run with its own suffix writes, and `"latest"` reads, only
#'   files with that suffix, so test and production fits never mix. The same `fitFileSuffix` parameter
#'   of `fireSense_spreadFit` and `fireSense_ignitionFit` sets it for both ledgers.
#'
#' @return `spreadFitFilenameFor()`: the file name a fit over `fireYears` is written to,
#'   e.g. `"fireSenseParams_1985-2024_linearFuel_esc50_20260925.rds"` (then `suffix`).
#' @export
#' @rdname spreadFitLedger
spreadFitFilenameFor <- function(fireYears, suffix = "") {
  fireYears <- as.integer(fireYears[is.finite(fireYears)])
  if (!length(fireYears))
    stop("spreadFitFilenameFor(): `fireYears` has no years")
  paste0("fireSenseParams_", min(fireYears), "-", max(fireYears), spreadFitFileTag, suffix, ".rds")
}

#' @export
#' @rdname spreadFitLedger
spreadFitFileTag <- "_linearFuel_esc50_20260925"

#' @description
#' `latestSpreadFits()` returns, for every polygon, its rows from the most recently modified
#' ledger file that has it. Files are read newest first; with `polygonIDs`, reading stops once
#' all of them are found. A file already in `destinationPath` with the same MD5 as on Drive is not
#' downloaded again; one that differs is downloaded again with [reproducible::preProcess()], into a
#' temporary folder first.
#' Each listing of the folder is saved in `destinationPath`. When Google Drive cannot be reached
#' (a compute node with no internet), that saved listing and the local copies are used instead,
#' with any local ledger file written since then read first.
#'
#' @param cloudFolderID The Google Drive folder (url or id) holding the ledger files; default [spreadFitLedgerFolder()].
#' @param destinationPath Local folder for the downloaded files.
#' @param polygonIDs Optional. The polygons (ELF ids) wanted; `NULL` reads every file.
#' @param suffix Only files named with this suffix are read; see `spreadFitFilenameFor()`.
#'
#' @return `latestSpreadFits()`: the ledger rows, in the form fireSense_spreadFit writes them
#'   (a `data.frame` with a `geometry` column), or `NULL` when no file matches. Attribute
#'   `"spreadFitFiles"` names the file each polygon came from.
#' @export
#' @rdname spreadFitLedger
latestSpreadFits <- function(cloudFolderID = spreadFitLedgerFolder(), destinationPath, polygonIDs = NULL,
                             suffix = "") {
  .latestLedgerFits(cloudFolderID, destinationPath, polygonIDs, filePrefix = "fireSenseParams_",
                    fileTag = spreadFitFileTag, suffix = suffix, caller = "latestSpreadFits", attrName = "spreadFitFiles")
}
