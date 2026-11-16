## The columns of a year-effects table are those of fireSense_ignitionFit's fittedYearEffects():
## year, observed, expected, priorSD, multiplier, shape.
yearEffectColumns <- c("year", "observed", "expected", "priorSD", "multiplier", "shape")

#' Read the fitted year effects of national ignition ledgers
#'
#' The ignition fit stores, per ELF and process, the year effects a hindcast may use
#' (`modelList$yearEffects`, `fireSense_ignitionFit`'s `fittedYearEffects()`). This reads them into one long table.
#'
#' @param ledgerFiles Paths of ledgers `fireSenseIgnitionParams_*.rds`: a data.frame with one row per ELF
#'   (`polygonID`) and the list columns of [ignitionFitAdditionalColNamesTxt] (ignition, then escape), each element
#'   holding `$modelList$yearEffects`.
#' @return A data.frame, one row per ledger, ELF, process (`ignition`, `escape`) and year: `ledger`
#'   (file name without `.rds`), `ELF`, `process` and `year`, `observed`, `expected`, `priorSD` (the sd of the
#'   prior of the year effect), `multiplier` and `shape`. Ledger rows whose `yearEffects` is `NULL` (older fits)
#'   are skipped with a message. Stops when a file is missing or a ledger has two rows for one ELF.
#' @export
#' @examples
#' x <- data.frame(polygonID = "1.1")
#' ye <- data.frame(year = 2001:2003, observed = 5, expected = 5, priorSD = 0.5, multiplier = 1, shape = 10)
#' x$fireSense_IgnitionFitted <- list(list(modelList = list(yearEffects = ye)))
#' x$fireSense_EscapeFitted <- list(list(modelList = list(yearEffects = ye)))
#' f <- tempfile(fileext = ".rds"); saveRDS(x, f)
#' readYearEffects(f)
readYearEffects <- function(ledgerFiles) {
  if (anyNA(ledgerFiles) || !length(ledgerFiles)) stop("ledgerFiles must name at least one ledger file")
  missingFiles <- ledgerFiles[!file.exists(ledgerFiles)]
  if (length(missingFiles)) stop("ledger file(s) not found: ", paste(missingFiles, collapse = ", "))
  slots <- stats::setNames(ignitionFitAdditionalColNamesTxt, c("ignition", "escape"))
  pieces <- list()
  for (f in ledgerFiles) {
    x <- readRDS(f); ledger <- sub("\\.rds$", "", basename(f))
    if (anyDuplicated(x$polygonID)) stop(ledger, " has more than one row for ELF ", paste(unique(x$polygonID[duplicated(x$polygonID)]), collapse = ", "))
    for (i in seq_len(nrow(x))) for (p in names(slots)) {
      ye <- x[[slots[[p]]]][[i]]$modelList$yearEffects
      if (is.null(ye)) {
        message("no yearEffects for ", p, " of ELF ", x$polygonID[i], " in ", ledger, ": skipped")
        next
      }
      pieces[[length(pieces) + 1L]] <- data.frame(ledger = ledger, ELF = as.character(x$polygonID[i]), process = p,
                                                  ye[, yearEffectColumns, drop = FALSE], row.names = NULL)
    }
  }
  if (!length(pieces)) stop("none of the ledgers has any yearEffects")
  do.call(rbind, pieces)
}
