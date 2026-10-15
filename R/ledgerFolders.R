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
#' @param folder The folder (url or id). An explicit value wins over the option.
#'
#' @return A character string, the folder's url or id.
#' @name fitLedgerFolders
#' @rdname fitLedgerFolders
#' @export
spreadFitLedgerFolder <- function(
    folder = getOption("fireSenseUtils.spreadFitLedgerFolder",
                       "https://drive.google.com/drive/folders/1X9-mRjyLMNpgkP_cfqhbr_AQEPOsVCHf")) {
  folder
}

#' @rdname fitLedgerFolders
#' @export
ignitionFitLedgerFolder <- function(
    folder = getOption("fireSenseUtils.ignitionFitLedgerFolder", spreadFitLedgerFolder())) {
  folder
}
