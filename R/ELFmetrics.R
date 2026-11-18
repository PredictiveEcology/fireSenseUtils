#' Per-run ELF metrics, for the national map
#'
#' `fireSense_summary` writes one file per run, named by `ELFmetricsFilename()`: a one-row table with
#' `runName` (the run's `.runName`), `ELF`, `spreadFitFile` and `ignitionFitFile` (the ledger files whose fits the
#' run used), then the run's metrics. A rerun of a run overwrites its own file, and no file is shared between
#' runs, so runs on different machines never write over each other.
#'
#' `readELFmetrics()` reads the files under `path`. With `current = TRUE`, it keeps only the runs made with the
#' current fits: a run made with an older fit describes that fit, not the current one. The current fits of an
#' ELF are, as with `"latest"`, the newest: among the ELF's runs whose files are of the current model
#' ([spreadFitFileTag], [ignitionFitFileTag]), those whose files are newest by the date in their suffix
#' ([ledgerFileOrder()]; a run over several ELFs, with several files, by its oldest). With a `suffix`, only the
#' files [spreadFitFilenameFor()] and [ignitionFitFilenameFor()] give for `fireYears` and `suffix` count, so one
#' experiment's runs can be read. This relies on a refit always going to a new ledger file, never replacing a row.
#'
#' @param runName The run's `.runName`.
#' @param path Folder(s) searched, with their subfolders, for the files.
#' @param fireYears,suffix With a `suffix`, only the ledger files of these fire years and this suffix are current;
#'   see [spreadFitFilenameFor()]. `""` (default): the newest files of the current model, whatever their suffix.
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
readELFmetrics <- function(path, fireYears = defaultFireYears(), suffix = "", current = TRUE) {
  files <- list.files(path, pattern = "^ELFmetrics_.*\\.csv$", recursive = TRUE, full.names = TRUE)
  if (!length(files)) {
    message("readELFmetrics(): no ELFmetrics_*.csv under ", paste(path, collapse = ", "))
    return(NULL)
  }
  keys <- c("runName", "ELF", "spreadFitFile", "ignitionFitFile")
  ## ELF codes are text: "4.10" is not "4.1"
  m <- data.table::rbindlist(lapply(files, data.table::fread, colClasses = list(character = keys)),
                             use.names = TRUE, fill = TRUE)
  if (isTRUE(current)) {
    keep <- .newestFits(m$spreadFitFile, m$ELF, spreadFitFileTag,
                        if (nzchar(suffix)) spreadFitFilenameFor(fireYears, suffix)) &
      .newestFits(m$ignitionFitFile, m$ELF, ignitionFitFileTag,
                  if (nzchar(suffix)) ignitionFitFilenameFor(fireYears, suffix))
    if (!all(keep))
      message("readELFmetrics(): dropped ", sum(!keep), " of ", NROW(m), " runs made with other fits")
    m <- m[keep]
  }
  m
}

## TRUE for each run whose ledger files (", "-joined) are of the model `fileTag` (and among `only`, if given)
## and the newest of its ELF's such runs, by the date in their suffix (.ledgerNameDate(); none is oldest)
.newestFits <- function(files, ELF, fileTag, only = NULL) {
  parts <- strsplit(ifelse(is.na(files), "", files), ", ", fixed = TRUE)
  ok <- vapply(parts, function(f) length(f) > 0L && all(grepl(fileTag, f, fixed = TRUE)) &&
                 (is.null(only) || all(f %in% only)), logical(1))
  date <- vapply(parts, function(f) {
    d <- .ledgerNameDate(f, fileTag)
    if (length(d)) min(ifelse(is.na(d), -Inf, d)) else -Inf
  }, numeric(1))
  date[!ok] <- NA
  newest <- stats::ave(date, as.character(ELF), FUN = function(x) if (all(is.na(x))) NA else max(x, na.rm = TRUE))
  ok & !is.na(newest) & date == newest
}
