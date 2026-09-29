#' The shared SpreadFit ledger files
#'
#' fireSense_spreadFit writes each fit as a row of an `rds` file in a Google Drive folder
#' (its `spreadFitGoogleDriveFolder`); fireSense_ELFs and fireSense_dataPrepFit read them.
#' With `spreadFitFilename = "latest"` those modules use these functions instead of one named file.
#'
#' `spreadFitFileTag` marks the files whose fits use the current model: fuel biomass on the
#' linear scale ([fuelLogToLinear()]), escaped fires starting at 50 ha (`escapeSizeHa`), and the
#' objective with the annual-area and area-distribution terms. Files with an older tag
#' (`"_linearFuel"`: 1-pixel escape; none: log fuel) hold fits of an earlier model, which
#' fireSense_spreadPredict would apply wrongly, so `"latest"` never reads them.
#' A fit made under a new model needs a new tag, so the old files drop out of `"latest"`.
#'
#' @param fireYears The fire years the fit used; only the first and last matter.
#'
#' @return `spreadFitFilenameFor()`: the file name a fit over `fireYears` is written to,
#'   e.g. `"fireSenseParams_1985-2024_linearFuel_esc50.rds"`.
#' @export
#' @rdname spreadFitLedger
spreadFitFilenameFor <- function(fireYears) {
  fireYears <- as.integer(fireYears[is.finite(fireYears)])
  if (!length(fireYears))
    stop("spreadFitFilenameFor(): `fireYears` has no years")
  paste0("fireSenseParams_", min(fireYears), "-", max(fireYears), spreadFitFileTag, ".rds")
}

#' @export
#' @rdname spreadFitLedger
spreadFitFileTag <- "_linearFuel_esc50"

#' @description
#' `driveDownloadAtomic()` downloads one Drive file to `path` without ever leaving a partial file
#' there. Other processes read the ledger files in `destinationPath` (the two held-out folds of an ELF
#' share one), so `googledrive::drive_download()` straight onto `path` lets a reader see a truncated
#' file. The download goes to a temporary file in the same folder, which then replaces `path` with
#' `file.rename()` (atomic on POSIX). On Windows, where `file.rename()` cannot be relied on to replace
#' a file, it falls back to `file.copy(overwrite = TRUE)`, which is not atomic. A `path` that already
#' has the MD5 Drive reports for `file` is left alone.
#'
#' @param file A one-row `dribble` (as from `googledrive::drive_ls()`) with the Drive file.
#' @param path The local file to create or replace.
#'
#' @return `driveDownloadAtomic()`: `TRUE` if it downloaded, `FALSE` if `path` was already current,
#'   invisibly.
#' @export
#' @rdname spreadFitLedger
driveDownloadAtomic <- function(file, path) {
  remoteMd5 <- file$drive_resource[[1]]$md5Checksum
  if (file.exists(path) && identical(unname(tools::md5sum(path)), remoteMd5))
    return(invisible(FALSE))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(".", basename(path), "_"), tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  googledrive::drive_download(file, path = tmp, overwrite = TRUE)
  ok <- if (.Platform$OS.type == "windows") {
    file.copy(tmp, path, overwrite = TRUE)
  } else {
    file.rename(tmp, path)
  }
  if (!isTRUE(ok))
    stop("driveDownloadAtomic(): could not put the download at ", path)
  invisible(TRUE)
}

#' @description
#' `latestSpreadFits()` returns, for every polygon, its rows from the most recently modified
#' ledger file that has it. Files are read newest first; with `polygonIDs`, reading stops once
#' all of them are found. A file already in `destinationPath` with the same MD5 as on Drive is not
#' downloaded again.
#'
#' @param cloudFolderID The Google Drive folder (url or id) holding the ledger files.
#' @param destinationPath Local folder for the downloaded files.
#' @param polygonIDs Optional. The polygons (ELF ids) wanted; `NULL` reads every file.
#'
#' @return `latestSpreadFits()`: the ledger rows, in the form fireSense_spreadFit writes them
#'   (a `data.frame` with a `geometry` column), or `NULL` when no file matches. Attribute
#'   `"spreadFitFiles"` names the file each polygon came from.
#' @export
#' @rdname spreadFitLedger
latestSpreadFits <- function(cloudFolderID, destinationPath, polygonIDs = NULL) {
  files <- googledrive::drive_ls(cloudFolderID)
  files <- files[grepl(paste0("^fireSenseParams_.*", spreadFitFileTag, "\\.rds$"), files$name), ]
  if (!NROW(files)) {
    message("latestSpreadFits(): no fireSenseParams_*", spreadFitFileTag, ".rds file in ", cloudFolderID)
    return(NULL)
  }
  modified <- vapply(files$drive_resource, function(r) r$modifiedTime, character(1))
  files <- files[order(modified, decreasing = TRUE), ]

  dir.create(destinationPath, recursive = TRUE, showWarnings = FALSE)
  rows <- list()
  from <- character()                              # polygonID -> file
  for (i in seq_len(NROW(files))) {
    fname <- files$name[i]
    local <- file.path(destinationPath, fname)
    driveDownloadAtomic(files[i, ], local)
    ledger <- as.data.frame(readRDS(local))
    ids <- as.character(ledger[[polygonIDTxt]])
    new <- !ids %in% names(from)
    if (any(new)) {
      rows[[fname]] <- ledger[new, , drop = FALSE]
      from <- c(from, setNames(rep(fname, length(unique(ids[new]))), unique(ids[new])))
    }
    if (!is.null(polygonIDs) && all(as.character(polygonIDs) %in% names(from)))
      break
  }
  if (!length(rows))
    return(NULL)
  out <- if (length(rows) == 1L) {
    rows[[1]]
  } else {
    ## as reproducible::CacheGeo() appends rows to a ledger
    as.data.frame(data.table::rbindlist(lapply(rows, data.table::as.data.table),
                                        fill = TRUE, use.names = TRUE))
  }
  rownames(out) <- NULL
  for (f in unique(from))
    message("latestSpreadFits(): ", paste(names(from)[from == f], collapse = ", "), " from ", f)
  attr(out, "spreadFitFiles") <- from
  out
}
