#' Per-run ELF metrics, for the national map
#'
#' `fireSense_summary` writes one file per run, named by `ELFmetricsFilename()`: a one-row table with
#' `runName` (the run's `.runName`), `ELF`, `spreadFitFile` and `ignitionFitFile` (the ledger files whose fits the
#' run used), then the run's metrics. A rerun of a run overwrites its own file, and no file is shared between
#' runs, so runs on different machines never write over each other.
#'
#' `readELFmetrics()` reads the files under `path`. With `current = TRUE`, it keeps only the runs made with the
#' current fits: a run made with an older fit describes that fit, not the current one. The current fits of an
#' ELF are, as with `"latest"`, the newest: among the ELF's runs, those whose fits were written last, by
#' `spreadFitWrittenAt` and `ignitionFitWrittenAt` (when the run's fits were written; a run over several ELFs is
#' as old as its oldest fit), or, for a metrics file without them, by the date in the ledger file names
#' ([ledgerFileOrder()]; none is oldest). With `files`, only runs whose ledger files are all among them count, so
#' one experiment's runs can be read. This relies on a refit always being a new ledger row, never replacing one.
#'
#' @param runName The run's `.runName`.
#' @param path Folder(s) searched, with their subfolders, for the files.
#' @param files `NULL` (default), or the ledger file names (spread and ignition) whose runs are current.
#' @param current If `TRUE` (default), drop the runs made with other fits.
#'
#' @return `ELFmetricsFilename()`: the file name, e.g. `"ELFmetrics_13.1_NRV_370_1.csv"`.
#'   `readELFmetrics()`: a `data.table` with one row per run, or `NULL` when there is no file.
#'   Pass it, with one row per ELF, to [plotELFs()]'s `metrics`.
#' @export
#' @rdname ELFmetrics
ELFmetricsFilename <- function(runName) paste0("ELFmetrics_", runName, ".csv")

#' @export
#' @rdname ELFmetrics
#' @examples
#' \dontrun{
#' m <- readELFmetrics("/mnt/fast/outputs")
#' ## one row per ELF: the mean over its runs
#' m <- m[, lapply(.SD, mean), by = "ELF", .SDcols = is.numeric]
#' plotELFs(metrics = m, metric = "areaBurnedHaPerYear")
#' }
readELFmetrics <- function(path, files = NULL, current = TRUE) {
  csvs <- list.files(path, pattern = "^ELFmetrics_.*\\.csv$", recursive = TRUE, full.names = TRUE)
  if (!length(csvs)) {
    message("readELFmetrics(): no ELFmetrics_*.csv under ", paste(path, collapse = ", "))
    return(NULL)
  }
  keys <- c("runName", "ELF", "spreadFitFile", "ignitionFitFile")
  ## ELF codes are text: "4.10" is not "4.1"
  m <- data.table::rbindlist(lapply(csvs, data.table::fread, colClasses = list(character = keys)),
                             use.names = TRUE, fill = TRUE)
  if (isTRUE(current)) {
    keep <- .newestFits(m$spreadFitFile, m$ELF, spreadFitFileTag, files, m[["spreadFitWrittenAt"]]) &
      .newestFits(m$ignitionFitFile, m$ELF, ignitionFitFileTag, files, m[["ignitionFitWrittenAt"]])
    if (!all(keep))
      message("readELFmetrics(): dropped ", sum(!keep), " of ", NROW(m), " runs made with other fits")
    m <- m[keep]
  }
  m
}

## TRUE for each run whose ledger files (", "-joined) are known (and among `only`, if given) and whose fits are the
## newest of its ELF's such runs: by `writtenAt` (when the fits were written), else by the date in the file names
## (.ledgerNameDate() with `fileTag`; none is oldest)
.newestFits <- function(files, ELF, fileTag, only = NULL, writtenAt = NULL) {
  parts <- strsplit(ifelse(is.na(files), "", files), ", ", fixed = TRUE)
  ok <- vapply(parts, function(f) length(f) > 0L && (is.null(only) || all(f %in% only)), logical(1))
  w <- if (is.null(writtenAt)) rep(NA_real_, length(parts)) else as.numeric(as.POSIXct(writtenAt, tz = "UTC"))
  date <- vapply(seq_along(parts), function(i) {
    d <- .ledgerRowDate(rep(w[i], length(parts[[i]])), parts[[i]], fileTag)   # one date per file, all `writtenAt` when it is known
    if (length(d)) min(ifelse(is.na(d), -Inf, d)) else -Inf
  }, numeric(1))
  date[!ok] <- NA
  newest <- stats::ave(date, as.character(ELF), FUN = function(x) if (all(is.na(x))) NA else max(x, na.rm = TRUE))
  ok & !is.na(newest) & date == newest
}
