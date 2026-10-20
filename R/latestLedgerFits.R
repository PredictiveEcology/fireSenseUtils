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

## Where the listing of a ledger folder is kept, beside the local copies of its files
.ledgerListingPath <- function(cloudFolderID, destinationPath) {
  file.path(destinationPath, paste0(".ledgerListing_", googledrive::as_id(cloudFolderID), ".rds"))
}

## The files of the ledger folder (name, id, md5, modifiedTime). Each listing is saved in
## `destinationPath`; when Drive cannot be reached (a compute node with no internet) the saved one
## is used instead, with attribute "offline" TRUE.
.ledgerListing <- function(cloudFolderID, destinationPath, caller) {
  path <- .ledgerListingPath(cloudFolderID, destinationPath)
  files <- tryCatch(googledrive::drive_ls(cloudFolderID), error = function(e) {
    if (!reproducible::isConnectionError(e)) stop(e)
    NULL
  })
  if (is.null(files)) {
    if (!file.exists(path))
      stop(caller, "(): Google Drive cannot be reached, and ", destinationPath, " has no saved listing of ",
           cloudFolderID, ". Run once with internet (it is saved at every run) and copy ", path, " with the inputs.")
    listing <- readRDS(path)
    message(caller, "(): Google Drive cannot be reached; using the listing saved ",
            format(file.mtime(path), "%Y-%m-%d %H:%M"), " and the local copies in ", destinationPath)
    attr(listing, "offline") <- TRUE
    return(listing)
  }
  field <- function(f) vapply(files$drive_resource, function(r) if (is.null(r[[f]])) NA_character_ else r[[f]], character(1))
  listing <- data.frame(name = files$name, id = files$id, md5 = field("md5Checksum"), modifiedTime = field("modifiedTime"))
  dir.create(destinationPath, recursive = TRUE, showWarnings = FALSE)
  ## via a temporary file: two jobs on one ELF share `destinationPath`
  tmp <- tempfile(".ledgerListing_", tmpdir = destinationPath, fileext = ".rds")
  saveRDS(listing, tmp)
  file.rename(tmp, path)
  listing
}

.latestLedgerFits <- function(cloudFolderID, destinationPath, polygonIDs, filePrefix, fileTag, suffix = "", caller, attrName) {
  files <- .ledgerListing(cloudFolderID, destinationPath, caller)
  offline <- isTRUE(attr(files, "offline"))
  ## the tag is followed directly by the suffix, so "" matches only unsuffixed files
  tagAndSuffix <- gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", paste0(fileTag, suffix))
  pattern <- paste0("^", filePrefix, ".*", tagAndSuffix, "\\.rds$")
  files <- files[grepl(pattern, files$name), ]
  files <- files[order(files$modifiedTime, decreasing = TRUE), ]
  if (offline) {
    ## a local file not in the listing was written since, i.e. offline: newer than any listed file
    written <- setdiff(list.files(destinationPath, pattern = pattern), files$name)
    written <- written[order(file.mtime(file.path(destinationPath, written)), decreasing = TRUE)]
    files <- rbind(data.frame(name = written, id = rep(NA_character_, length(written)),
                              md5 = rep(NA_character_, length(written)),
                              modifiedTime = rep(NA_character_, length(written))), files)
  }
  if (!NROW(files)) {
    message(caller, "(): no ", filePrefix, "*", fileTag, suffix, ".rds file in ", cloudFolderID)
    return(NULL)
  }

  dir.create(destinationPath, recursive = TRUE, showWarnings = FALSE)
  rows <- list()
  from <- character()                              # polygonID -> file
  for (i in seq_len(NROW(files))) {
    fname <- files$name[i]
    local <- file.path(destinationPath, fname)
    remoteMd5 <- files$md5[i]
    if (!is.na(files$id[i])) {
      fileUrl <- .driveFileUrl(files$id[i])
      ## read from the local copy or downloaded: either way this file's data came from Drive
      reproducible::logUrlAccess(fileUrl, fn = caller, targetFile = fname, destinationPath = destinationPath)
    }
    if (offline) {
      ## the local copy as it is: it may hold rows written offline that Drive does not have yet
      if (!file.exists(local))
        stop(caller, "(): Google Drive cannot be reached and ", fname, " is not in ", destinationPath,
             ". Run once with internet and copy it with the inputs.")
    } else if (!file.exists(local) || !identical(unname(tools::md5sum(local)), remoteMd5)) {
      ## reproducible downloads into a temporary folder and then replaces `local`, so a job
      ## reading the same file (two jobs on one ELF share `destinationPath`) never sees it half-written.
      ## `purge = 7` because CHECKSUMS.txt still matches the old local file.
      local <- reproducible::preProcess(
        url = fileUrl,
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
