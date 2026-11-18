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
#' fireSense_spreadPredict would apply wrongly; their rows record no fuel scale, so `"latest"` skips them
#' (see `spreadFuelScale`).
#'
#' @param fireYears The fire years the fit used; only the first and last matter.
#' @param suffix Text put after the model tag in the file name, e.g. `"_bufIndep20261024"`; end it with
#'   the date of the fits (`yyyymmdd`), which orders the files for `"latest"`. A run with a suffix writes,
#'   and `"latest"` reads, only files with that suffix. With `""` (the default) `"latest"` reads every file
#'   of the model tag in the folder, newest date first. Experiment fits go in the folder's `"testing"`
#'   subfolder (`spreadFitLedgerFolder(testing = TRUE)`), which `"latest"` over the main folder never reads. The same `fitFileSuffix` parameter
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
#' `spreadFuelScale` is the scale of the fuel covariates of the spread model ([fuelLogToLinear()]);
#' fireSense_spreadFit records it in each ledger row (column `spreadFitFuelScaleTxt`). With the escape
#' definition the row also records ([escapeDefinition()]), it is what prediction must apply as the fit did, so
#' [latestSpreadFits()] skips a row whose fuel scale or escape definition differs from this run's. A row written
#' before the fuel scale was recorded counts as `spreadFuelScale` when its file carries `spreadFitFileTag`, and is
#' skipped otherwise. A change of model is then recorded by the fit itself: `spreadFitFileTag` is no longer
#' changed by hand, and only dates the files written before.
#' @export
#' @rdname spreadFitLedger
spreadFuelScale <- "linear"

#' @export
#' @rdname spreadFitLedger
spreadFitFuelScaleTxt <- "fuelScale"

## Why this run cannot use a SpreadFit ledger row of `file`, or NULL. `escape`: the run's escapeDefinition(), or NULL.
.spreadRowProblem <- function(row, file, escape = NULL) {
  scale <- if (spreadFitFuelScaleTxt %in% names(row)) as.character(row[[spreadFitFuelScaleTxt]][[1]])
  if (is.null(scale) || is.na(scale)) {
    if (!grepl(spreadFitFileTag, file, fixed = TRUE))
      return(paste0("its row records no fuel scale and its file is not of the model '", spreadFitFileTag, "'"))
  } else if (!identical(scale, spreadFuelScale)) {
    return(paste0("its fuel covariates are on the ", scale, " scale, not ", spreadFuelScale))
  }
  .escapeProblem(row, escape)
}

## checkEscapeDefinition() as a reason (NULL when the row agrees or records none)
.escapeProblem <- function(row, escape) {
  if (is.null(escape)) return(NULL)
  tryCatch({checkEscapeDefinition(row, escape); NULL}, error = function(e) conditionMessage(e))
}

#' @description
#' `latestSpreadFits()` returns, for every polygon, its newest usable row in the ledger files (any file of the
#' folder; rows of another fuel scale or escape definition are skipped, see `spreadFuelScale`): the one with
#' the latest `writtenAt` (stamped by [reproducible::CacheGeoWrite()]). A row written before that stamp
#' existed has its file's date instead: the date in the file's suffix ([ledgerFileOrder()]), else older
#' than any dated file; ties go by the file's modified time. No row is newer than its file was last
#' modified, so files are read newest-modified first; with `polygonIDs`, reading stops once each has a
#' row newer than the next file. Each file is read with [reproducible::CacheGeoRead()] into
#' `destinationPath`: downloaded only if it changed on Drive, and from the local copy when Drive
#' cannot be reached. The folder's file list is saved there too, and used when Drive cannot be
#' reached (a compute node with no internet).
#'
#' @param cloudFolderID The Google Drive folder (url or id) holding the ledger files; default [spreadFitLedgerFolder()].
#' @param destinationPath Local folder for the downloaded files.
#' @param polygonIDs Optional. The polygons (ELF ids) wanted; `NULL` reads every file.
#' @param suffix A suffix reads only the files named with it; `""` (default) reads every ledger file.
#'   See `spreadFitFilenameFor()`.
#' @param escapeDefinition This run's [escapeDefinition()] (e.g. [escapeDefinitionOf()]): rows fitted under
#'   another are skipped. `NULL` (default): not checked.
#'
#' @return `latestSpreadFits()`: the ledger rows, in the form fireSense_spreadFit writes them
#'   (a `data.frame` with a `geometry` column), or `NULL` when no file matches. Attribute
#'   `"spreadFitFiles"` names the file each polygon came from.
#' @export
#' @rdname spreadFitLedger
latestSpreadFits <- function(cloudFolderID = spreadFitLedgerFolder(), destinationPath, polygonIDs = NULL,
                             suffix = "", escapeDefinition = NULL) {
  .latestLedgerFits(cloudFolderID, destinationPath, polygonIDs, filePrefix = "fireSenseParams_",
                    fileTag = spreadFitFileTag, suffix = suffix, caller = "latestSpreadFits", attrName = "spreadFitFiles",
                    rowProblem = function(row, file) .spreadRowProblem(row, file, escapeDefinition))
}
