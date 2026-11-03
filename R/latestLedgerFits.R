#' Read the newest ledger rows of each polygon from the `rds` files of a Drive folder
#'
#' The one implementation behind [latestSpreadFits()] and [latestIgnitionFits()].
#'
#' @param cloudFolderID,destinationPath,polygonIDs As in [latestSpreadFits()].
#' @param filePrefix,fileTag,suffix The ledger files are those named `<filePrefix>*<fileTag><suffix>.rds`;
#'   with `suffix = ""`, every `<filePrefix>*<fileTag>*.rds`, whatever its suffix ([ledgerFileOrder()]).
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
  ## The tag is followed directly by the suffix. A suffix reads only its own files; "" (the default) reads
  ## every file of the model tag, whatever its suffix, newest by the date in its suffix (ledgerFileOrder()).
  ## drive_ls() is not recursive: files in a subfolder (the "testing" folder) are never listed.
  esc <- function(x) gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
  files <- files[grepl(paste0("^", filePrefix, ".*", esc(fileTag), if (nzchar(suffix)) esc(suffix) else ".*",
                              "\\.rds$"), files$name), ]
  if (!NROW(files)) {
    message(caller, "(): no ", filePrefix, "*", fileTag, if (nzchar(suffix)) suffix else "*", ".rds file in ",
            cloudFolderID)
    return(NULL)
  }
  files <- files[ledgerFileOrder(files$name, files$modifiedTime, fileTag), ]

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

#' The order to read ledger files in: newest first, by the date in the suffix
#'
#' A ledger file is named `<prefix>*<fileTag><suffix>.rds`. Its date is the last 8-digit `yyyymmdd` in
#' the suffix (`"_gpscCause_20261029"` is 20261029, `"_bufIndep20261024"` 20261024): when its fits were
#' made, which the file's modified time is not (a file is modified whenever a fit is added). A file
#' without a date in its suffix (e.g. no suffix) is older than any dated one; files of the same date, or
#' both undated, are ordered by modified time.
#'
#' @param names File names.
#' @param modifiedTime Their modified times (ISO 8601 strings, as Google Drive gives them).
#' @param fileTag The model tag ([spreadFitFileTag], [ignitionFitFileTag]).
#' @return Integer: the order of `names`, newest first.
#' @export
ledgerFileOrder <- function(names, modifiedTime, fileTag) {
  suffix <- sub("\\.rds$", "", vapply(strsplit(names, fileTag, fixed = TRUE), function(x)
    if (length(x) > 1L) x[length(x)] else "", character(1)))
  date <- vapply(regmatches(suffix, gregexpr("[0-9]{8}", suffix)), function(d)
    if (length(d)) d[length(d)] else "00000000", character(1))
  order(date, as.character(modifiedTime), decreasing = TRUE)
}
