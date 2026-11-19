#' The Google Drive folders of the shared fit ledgers and the ELF map cache
#'
#' The spread-fit and ignition-fit ledger files (see [spreadFitFilenameFor()] and
#' [ignitionFitFilenameFor()]) each have their own folder in the fireSense Drive folder
#' `FireSenseData`. The modules that write and read the ledgers (fireSense_spreadFit,
#' fireSense_ignitionFit, fireSense_dataPrepFit, fireSense_ELFs) use these as the defaults of their
#' `spreadFitGoogleDriveFolder` and `ignitionFitGoogleDriveFolder` parameters, and
#' [latestSpreadFits()] and [latestIgnitionFits()] use them as the default `cloudFolderID`.
#' `elfCacheFolder()` is where fireSense_ELFs shares the national ELF map through the cloud cache.
#'
#' @section Options:
#' To use another folder, set an option; every module and function then follows it.
#' \describe{
#'   \item{`fireSenseUtils.spreadFitLedgerFolder`}{the spread-fit ledger folder (url or id).
#'     Default: `FireSenseData/SpreadFitParams`.}
#'   \item{`fireSenseUtils.ignitionFitLedgerFolder`}{the ignition-fit ledger folder (url or id).
#'     Default: `FireSenseData/IgnitionFitParams`.}
#'   \item{`fireSenseUtils.elfCacheFolder`}{the ELF map cache folder (url or id).
#'     Default: `FireSenseData/ELFs`.}
#' }
#'
#' @section Experiments:
#' Fits made to try something (a new objective, weights, a data change) go in each ledger folder's
#' subfolder `"testing"`: `testing = TRUE` gives it (options
#' `fireSenseUtils.spreadFitLedgerFolderTesting`, `fireSenseUtils.ignitionFitLedgerFolderTesting`).
#' An experiment sets them as the modules' `spreadFitGoogleDriveFolder` and
#' `ignitionFitGoogleDriveFolder`. [latestSpreadFits()] and [latestIgnitionFits()] list a folder's
#' own files only (`googledrive::drive_ls()` is not recursive), so the default `"latest"` over the
#' main folder never reads an experiment's fits. Each ledger folder, and each `"testing"` folder,
#' also has a subfolder `"archived"` for ledgers that no longer hold the newest fit of any ELF;
#' `"latest"` never reads those either.
#'
#' @param folder The folder (url or id). An explicit value wins over the option.
#' @param testing Logical. `TRUE`: the ledger folder's `"testing"` subfolder, for experiment fits.
#'
#' @return A character string, the folder's url or id.
#' @name fitLedgerFolders
#' @rdname fitLedgerFolders
#' @export
spreadFitLedgerFolder <- function(
    folder = if (isTRUE(testing))
      getOption("fireSenseUtils.spreadFitLedgerFolderTesting",
                "https://drive.google.com/drive/folders/1nqSVuGRYxIvcgZMn7TIbn5HfRiUHR2xq")
    else getOption("fireSenseUtils.spreadFitLedgerFolder",
                   "https://drive.google.com/drive/folders/1X9-mRjyLMNpgkP_cfqhbr_AQEPOsVCHf"),
    testing = FALSE) {
  folder
}

#' @rdname fitLedgerFolders
#' @export
ignitionFitLedgerFolder <- function(
    folder = if (isTRUE(testing))
      getOption("fireSenseUtils.ignitionFitLedgerFolderTesting",
                "https://drive.google.com/drive/folders/1nwJL5nz6E_7_9Ao7asNsNlmReCOIuOXl")
    else getOption("fireSenseUtils.ignitionFitLedgerFolder",
                   "https://drive.google.com/drive/folders/1Hhdtinm5pWt7P6A-p1MnzFHIJP_0ELq8"),
    testing = FALSE) {
  folder
}

#' @rdname fitLedgerFolders
#' @export
elfCacheFolder <- function(
    folder = getOption("fireSenseUtils.elfCacheFolder",
                       "https://drive.google.com/drive/folders/1Gl9CWW9o36Eqkm54x_fZdHPMB56Qoefd")) {
  folder
}
