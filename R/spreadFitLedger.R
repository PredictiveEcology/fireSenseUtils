#' The shared SpreadFit ledger files
#'
#' fireSense_SpreadFit writes each fit as a row of an `rds` file in a Google Drive folder
#' (its `spreadFitGoogleDriveFolder`); fireSense_ELFs and fireSense_dataPrepFit read them.
#' With `spreadFitFilename = "latest"` those modules use these functions instead of one named file.
#'
#' `spreadFitFileTag` marks the files whose fits use the current model: fuel biomass on the
#' linear scale ([fuelLogToLinear()]). Files without it hold fits on log fuel, which an older
#' fireSense_SpreadPredict would apply wrongly, so `"latest"` never reads them.
#' A fit made under a new model needs a new tag, so the old files drop out of `"latest"`.
#'
#' @param fireYears The fire years the fit used; only the first and last matter.
#'
#' @return `spreadFitFilenameFor()`: the file name a fit over `fireYears` is written to,
#'   e.g. `"fireSenseParams_1985-2024_linearFuel.rds"`.
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
spreadFitFileTag <- "_linearFuel"

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
#' @return `latestSpreadFits()`: the ledger rows, in the form fireSense_SpreadFit writes them
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
    remoteMd5 <- files$drive_resource[[i]]$md5Checksum
    if (!file.exists(local) || !identical(unname(tools::md5sum(local)), remoteMd5))
      googledrive::drive_download(files[i, ], path = local, overwrite = TRUE)
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
