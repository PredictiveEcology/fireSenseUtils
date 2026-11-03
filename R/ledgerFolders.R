#' The Google Drive folder of the shared fit ledgers
#'
#' The folder that holds the spread-fit and ignition-fit ledger files (see
#' [spreadFitFilenameFor()] and [ignitionFitFilenameFor()]). The modules that write and read
#' the ledgers (fireSense_spreadFit, fireSense_ignitionFit, fireSense_dataPrepFit,
#' fireSense_ELFs) use these as the defaults of their `spreadFitGoogleDriveFolder` and
#' `ignitionFitGoogleDriveFolder` parameters, and [latestSpreadFits()] and
#' [latestIgnitionFits()] use them as the default `cloudFolderID`.
#'
#' @section Options:
#' To use another folder, set an option; every module and function then follows it.
#' \describe{
#'   \item{`fireSenseUtils.spreadFitLedgerFolder`}{the spread-fit ledger folder (url or id).
#'     Default: the shared folder of the fireSense project.}
#'   \item{`fireSenseUtils.ignitionFitLedgerFolder`}{the ignition-fit ledger folder (url or id).
#'     Default: the value of the spread-fit ledger folder, so one option moves both.}
#' }
#'
#' @section Experiments:
#' Fits made to try something (a new objective, weights, a data change) go in the folder's
#' subfolder `"testing"`: `testing = TRUE` gives it (options
#' `fireSenseUtils.spreadFitLedgerFolderTesting`, `fireSenseUtils.ignitionFitLedgerFolderTesting`).
#' An experiment sets it as the modules' `spreadFitGoogleDriveFolder` and
#' `ignitionFitGoogleDriveFolder`. [latestSpreadFits()] and [latestIgnitionFits()] list a folder's
#' own files only (`googledrive::drive_ls()` is not recursive), so the default `"latest"` over the
#' main folder never reads an experiment's fits. The spread and ignition ledgers share one folder,
#' so they share one `"testing"` subfolder.
#'
#' @param folder The folder (url or id). An explicit value wins over the option.
#' @param testing Logical. `TRUE`: the folder's `"testing"` subfolder, for experiment fits.
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
      getOption("fireSenseUtils.ignitionFitLedgerFolderTesting", spreadFitLedgerFolder(testing = TRUE))
    else getOption("fireSenseUtils.ignitionFitLedgerFolder", spreadFitLedgerFolder()),
    testing = FALSE) {
  folder
}
