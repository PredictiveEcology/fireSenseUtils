#' Read the newest ledger rows of each polygon from the `rds` files of a Drive folder
#'
#' The one implementation behind [latestSpreadFits()] and [latestIgnitionFits()].
#'
#' @param cloudFolderID,destinationPath,polygonIDs As in [latestSpreadFits()].
#' @param filePrefix,fileTag The ledger files are those named `<filePrefix>*<fileTag>.rds`.
#' @param caller Name of the calling function, for messages.
#' @param attrName Name of the attribute that records the file each polygon came from.
#' @return As [latestSpreadFits()].
#' @keywords internal
#' @noRd
.latestLedgerFits <- function(cloudFolderID, destinationPath, polygonIDs, filePrefix, fileTag, caller, attrName) {
  files <- googledrive::drive_ls(cloudFolderID)
  files <- files[grepl(paste0("^", filePrefix, ".*", fileTag, "\\.rds$"), files$name), ]
  if (!NROW(files)) {
    message(caller, "(): no ", filePrefix, "*", fileTag, ".rds file in ", cloudFolderID)
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
    if (!file.exists(local) || !identical(unname(tools::md5sum(local)), remoteMd5)) {
      ## reproducible downloads into a temporary folder and then replaces `local`, so a job
      ## reading the same file (two jobs on one ELF share `destinationPath`) never sees it half-written.
      ## `purge = 7` because CHECKSUMS.txt still matches the old local file.
      local <- reproducible::preProcess(
        url = paste0("https://drive.google.com/file/d/", files$id[i]),
        targetFile = fname, destinationPath = destinationPath, fun = NA,
        purge = if (file.exists(local)) 7 else FALSE)$targetFilePath
    }
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
    message(caller, "(): ", paste(names(from)[from == f], collapse = ", "), " from ", f)
  attr(out, attrName) <- from
  out
}
