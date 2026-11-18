#' Read the newest ledger rows of each polygon from the `rds` files of a Drive folder
#'
#' The one implementation behind [latestSpreadFits()] and [latestIgnitionFits()].
#'
#' @param cloudFolderID,destinationPath,polygonIDs As in [latestSpreadFits()].
#' @param filePrefix,fileTag,suffix The ledger files are those named `<filePrefix>*.rds`; with a `suffix`, only
#'   `<filePrefix>*<suffix>.rds`. `fileTag` is the model tag in the file names, used to date them
#'   ([ledgerFileOrder()]) and by `rowProblem` for rows that record no model.
#' @param caller Name of the calling function, for messages.
#' @param attrName Name of the attribute that records the file each polygon came from.
#' @param rowProblem `NULL`, or `function(row, file)` returning `NULL` when a ledger row (one-row `data.frame`) of
#'   `file` can be used by this run, else why not (a string). Unusable rows are skipped: each polygon gets its
#'   newest usable row, and a message names the newer rows skipped and why.
#' @return As [latestSpreadFits()].
#' @keywords internal
#' @noRd
## The link of a Google Drive file, as recorded in reproducible's download ledger and used to download it.
.driveFileUrl <- function(id) paste0("https://drive.google.com/file/d/", id)

.latestLedgerFits <- function(cloudFolderID, destinationPath, polygonIDs, filePrefix, fileTag, suffix = "", caller, attrName,
                              rowProblem = NULL) {
  if (is.null(cloudFolderID))
    stop(caller, "(): `cloudFolderID` is NULL. The newest fits are found by listing a Google Drive folder; ",
         "to read a local ledger, name its file (the module's ledger filename parameter).", call. = FALSE)
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
  ## Every ledger file of the prefix, whatever its model: each row says what model it is, and `rowProblem` keeps
  ## the ones this run can use. A suffix reads only its own files. drive_ls() is not recursive: files in a
  ## subfolder (the "testing" folder) are never listed.
  esc <- function(x) gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
  files <- files[grepl(paste0("^", filePrefix, ".*", if (nzchar(suffix)) esc(suffix), "\\.rds$"), files$name), ]
  if (!NROW(files)) {
    message(caller, "(): no ", filePrefix, "*", if (nzchar(suffix)) suffix, ".rds file in ", cloudFolderID)
    return(NULL)
  }
  ## The newest row of each polygon: by its `writtenAt` (reproducible::CacheGeoWrite() stamps it); a row
  ## from before the stamp by its file's date as ledgerFileOrder() gives it (the date in the suffix, else
  ## older than any dated file; then the file's modified time). No row is newer than its file's ceiling,
  ## the later of its modified time and its suffix date, so files are read by ceiling, newest first, and
  ## reading stops once every polygon asked for has a row newer than the next file's ceiling.
  modified <- .ledgerTime(files$modifiedTime)
  nameDate <- .ledgerNameDate(files$name, fileTag)
  ceiling <- pmax(modified, ifelse(is.na(nameDate), -Inf, nameDate))
  read <- order(ceiling, modified, decreasing = TRUE)

  best <- list()                                   # polygonID -> list(file, row, key, tie)
  skipped <- list()                                # newer rows rowProblem() rejected
  for (k in read) {
    if (!is.null(polygonIDs) && all(as.character(polygonIDs) %in% names(best)) &&
        all(vapply(best[as.character(polygonIDs)], `[[`, numeric(1), "key") > ceiling[k]))
      break
    fname <- files$name[k]
    ledger <- reproducible::CacheGeoLedger(fname, remote = cloudFolderID, destinationPath = destinationPath)
    ledger <- as.data.frame(reproducible::CacheGeoRead(ledger, area = NULL))
    stamp <- if ("writtenAt" %in% names(ledger)) as.numeric(ledger$writtenAt) else rep(NA_real_, NROW(ledger))
    key <- ifelse(is.na(stamp), if (is.na(nameDate[k])) -Inf else nameDate[k], stamp)
    ids <- as.character(ledger[[polygonIDTxt]])
    for (i in seq_along(ids)) {
      b <- best[[ids[i]]]
      if (is.null(b) || key[i] > b$key || (key[i] == b$key && modified[k] > b$tie)) {
        why <- if (!is.null(rowProblem)) rowProblem(ledger[i, , drop = FALSE], fname)
        if (is.null(why))
          best[[ids[i]]] <- list(file = fname, row = ledger[i, , drop = FALSE], key = key[i], tie = modified[k])
        else
          skipped[[length(skipped) + 1L]] <- list(id = ids[i], file = fname, key = key[i], why = why)
      }
    }
  }
  ## the rows skipped that are newer than the one used (or with none used): a newer fit this run cannot use
  for (sk in skipped) {
    b <- best[[sk$id]]
    if (is.null(b) || sk$key > b$key)
      message(caller, "(): skipped the fit of ", sk$id, " in ", sk$file, ": ", sk$why)
  }
  from <- vapply(best, `[[`, character(1), "file")
  rows <- lapply(split(best, factor(from, levels = unique(from))), function(b)
    do.call(rbind, unname(lapply(b, `[[`, "row"))))
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
  date <- .ledgerNameDate(names, fileTag)
  order(ifelse(is.na(date), -Inf, date), as.character(modifiedTime), decreasing = TRUE)
}

## The last 8-digit `yyyymmdd` in the suffix after `fileTag`, as seconds (midnight UTC); NA if none
.ledgerNameDate <- function(names, fileTag) {
  suffix <- sub("\\.rds$", "", vapply(strsplit(names, fileTag, fixed = TRUE), function(x)
    if (length(x) > 1L) x[length(x)] else "", character(1)))
  date <- vapply(regmatches(suffix, gregexpr("[0-9]{8}", suffix)), function(d)
    if (length(d)) d[length(d)] else NA_character_, character(1))
  as.numeric(as.POSIXct(date, format = "%Y%m%d", tz = "UTC"))
}

## Drive's modifiedTime ("2026-10-29T10:00:00.123Z", or a date) as seconds
.ledgerTime <- function(x) {
  t <- as.POSIXct(sub("Z$", "", x), format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC")
  day <- is.na(t)
  t[day] <- as.POSIXct(x[day], format = "%Y-%m-%d", tz = "UTC")
  as.numeric(t)
}
