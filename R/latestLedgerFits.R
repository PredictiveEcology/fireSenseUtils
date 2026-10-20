#' Read the newest ledger rows of each polygon from the `rds` files of a Drive folder
#'
#' The one implementation behind [latestSpreadFits()] and [latestIgnitionFits()].
#'
#' @param cloudFolderID,destinationPath,polygonIDs As in [latestSpreadFits()].
#' @param filePrefix,fileTag,suffix The ledger files are those named `<filePrefix>*<fileTag><suffix>.rds`.
#' @param caller Name of the calling function, for messages.
#' @param attrName Name of the attribute that records the file each polygon came from.
#' @return As [latestSpreadFits()].
#' @keywords internal
#' @noRd
## The link of a Google Drive file, as recorded in reproducible's download ledger and used to download it.
.driveFileUrl <- function(id) paste0("https://drive.google.com/file/d/", id)

.latestLedgerFits <- function(cloudFolderID, destinationPath, polygonIDs, filePrefix, fileTag, suffix = "", caller, attrName) {
  ## Which ledger files exist, newest first. Each listing is saved beside the local copies, so a run
  ## with no internet (a compute node) picks the same files; CacheGeoRead() then reads each one,
  ## from Drive if it changed, else (or offline) from the local copy.
  listing <- file.path(destinationPath, paste0(".ledgerListing_", googledrive::as_id(cloudFolderID), ".rds"))
  files <- tryCatch({
    drive <- googledrive::drive_ls(cloudFolderID)
    drive <- data.frame(name = drive$name, modifiedTime = vapply(drive$drive_resource, `[[`, character(1), "modifiedTime"))
    dir.create(destinationPath, recursive = TRUE, showWarnings = FALSE)
    saveRDS(drive, listing)
    drive
  }, error = function(e) {
    if (!reproducible::isConnectionError(e) || !file.exists(listing)) stop(e)
    message(caller, "(): Google Drive cannot be reached; using the file list saved ", file.mtime(listing))
    readRDS(listing)
  })
  ## the tag is followed directly by the suffix, so "" matches only unsuffixed files
  tagAndSuffix <- gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", paste0(fileTag, suffix))
  files <- files[grepl(paste0("^", filePrefix, ".*", tagAndSuffix, "\\.rds$"), files$name), ]
  if (!NROW(files)) {
    message(caller, "(): no ", filePrefix, "*", fileTag, suffix, ".rds file in ", cloudFolderID)
    return(NULL)
  }
  files <- files[order(files$modifiedTime, decreasing = TRUE), ]

  rows <- list()
  from <- character()                              # polygonID -> file
  for (fname in files$name) {
    ledger <- reproducible::CacheGeoLedger(fname, remote = cloudFolderID, destinationPath = destinationPath)
    ledger <- as.data.frame(reproducible::CacheGeoRead(ledger, area = NULL))
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
    ## as reproducible::CacheGeo() appends rows to a ledger. rbind() of data.frames, not rbindlist():
    ## that assigns into list columns, which fails for the ALTREP lists of xgboost models read by readRDS()
    cols <- unique(unlist(lapply(rows, names), use.names = FALSE))
    do.call(rbind, unname(lapply(rows, function(r) {
      r[setdiff(cols, names(r))] <- NA
      r[cols]
    })))
  }
  rownames(out) <- NULL
  for (f in unique(from))
    message(caller, "(): ", paste(names(from)[from == f], collapse = ", "), " from ", f)
  attr(out, attrName) <- from
  out
}
