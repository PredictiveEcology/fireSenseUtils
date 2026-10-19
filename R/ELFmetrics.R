#' Per-run ELF metrics, for the national map
#'
#' `fireSense_summary` writes one file per run, named by `ELFmetricsFilename()`: a one-row table with
#' `runName` (the run's `.runName`), `ELF`, `spreadFitFile` and `ignitionFitFile` (the ledger files whose fits the
#' run used), then the run's metrics. A rerun of a run overwrites its own file, and no file is shared between
#' runs, so runs on different machines never write over each other.
#'
#' `readELFmetrics()` reads the files under `path`. With `current = TRUE`, it keeps only the runs made with the
#' current fits, those in [spreadFitFilenameFor()] and [ignitionFitFilenameFor()] of `fireYears` and `suffix`:
#' a run made with an older fit describes that fit, not the current one. This relies on a refit always going to
#' a new ledger file (a new `suffix`), never replacing a row of the current one.
#'
#' @param runName The run's `.runName`.
#' @param path Folder(s) searched, with their subfolders, for the files.
#' @param fireYears,suffix The fire years and ledger suffix of the current fits; see [spreadFitFilenameFor()].
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
    keep <- m$spreadFitFile %in% spreadFitFilenameFor(fireYears, suffix) &
      m$ignitionFitFile %in% ignitionFitFilenameFor(fireYears, suffix)
    if (!all(keep))
      message("readELFmetrics(): dropped ", sum(!keep), " of ", NROW(m), " runs made with other fits")
    m <- m[keep]
  }
  m
}
