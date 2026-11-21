#' The date of each ELF's latest fit in the spread and ignition ledgers
#'
#' For each ELF, when its newest fit was written, from the ledger rows [latestSpreadFits()] and
#' [latestIgnitionFits()] return (every ledger file of each folder, not one suffix): the row's `writtenAt`, else
#' the date in its ledger file's name ([ledgerFileOrder()]), else `NA`. Subfolders (`"testing"`, `"archived"`) are
#' never listed, so the default folders give the production fits.
#'
#' @param spreadFolder,ignitionFolder The ledger folders (url or id); defaults [spreadFitLedgerFolder()] and
#'   [ignitionFitLedgerFolder()] (production, not `testing = TRUE`).
#' @param destinationPath Local folder for the downloaded ledger files, as in [latestSpreadFits()].
#' @param ELFs Optional. ELF names to list even when they have no fit (rows of `NA`), e.g. every ELF of
#'   the map; they come first, in this order.
#'
#' @return A `data.frame` with columns `ELF`, `spreadFitDate`, `spreadLedger`, `ignitionFitDate` and
#'   `ignitionLedger` (dates are `POSIXct`, UTC; ledgers are file names), one row per ELF.
#' @export
#' @seealso [plotELFs()] draws it with `metric = "spreadFitDate"`.
elfFitDates <- function(spreadFolder = spreadFitLedgerFolder(), ignitionFolder = ignitionFitLedgerFolder(),
                        destinationPath, ELFs = NULL) {
  spread <- .ledgerFitDates(latestSpreadFits(spreadFolder, destinationPath), "spreadFitFiles", spreadFitFileTag)
  ignition <- .ledgerFitDates(latestIgnitionFits(ignitionFolder, destinationPath), "ignitionFitFiles", ignitionFitFileTag)
  ids <- unique(c(as.character(ELFs), spread$ELF, ignition$ELF))
  i <- match(ids, spread$ELF)
  j <- match(ids, ignition$ELF)
  data.frame(ELF = ids,
             spreadFitDate = spread$date[i], spreadLedger = spread$ledger[i],
             ignitionFitDate = ignition$date[j], ignitionLedger = ignition$ledger[j],
             stringsAsFactors = FALSE)
}

## The date and ledger file of each polygon of `rows` (a result of .latestLedgerFits(), whose attribute
## `attrName` names each polygon's file), dated by .ledgerRowDate()
.ledgerFitDates <- function(rows, attrName, fileTag) {
  if (is.null(rows))
    return(data.frame(ELF = character(), date = as.POSIXct(character(), tz = "UTC"), ledger = character()))
  ids <- as.character(rows[[polygonIDTxt]])
  ledger <- unname(attr(rows, attrName)[ids])
  stamp <- if ("writtenAt" %in% names(rows)) rows$writtenAt else rep(NA_real_, length(ids))
  date <- .ledgerRowDate(stamp, ledger, fileTag)
  data.frame(ELF = ids, date = as.POSIXct(date, origin = "1970-01-01", tz = "UTC"), ledger = ledger,
             stringsAsFactors = FALSE)
}
