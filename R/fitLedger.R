#' The shared fit ledgers: one constructor per ledger, one read, one row rule
#'
#' The spread-fit and the ignition-fit ledgers are read and written only through
#' [reproducible::CacheGeoLedger()], [reproducible::CacheGeoRead()] and
#' [reproducible::CacheGeoWrite()]. These functions are the one place that says which file a
#' ledger is, how its rows are selected for a study area, and which row is a polygon's own fit.
#'
#' * `spreadFitLedger()` and `ignitionFitLedger()` return a `CacheGeoLedger` keyed on
#'   `polygonID` ([polygonIDTxt]). A named file is used as given. `"latest"` writes to the file
#'   named for the fire years and `fitFileSuffix` ([spreadFitFilenameFor()],
#'   [ignitionFitFilenameFor()]). A spread ledger read with `"latest"` holds, for each polygon,
#'   the newest row of any file with this suffix ([latestSpreadFits()]); those rows are saved to
#'   a local file named by their content, so the ledger has no remote.
#' * `readFitLedger()` selects the rows for a study area with the rule of the old `CacheGeo()`
#'   call form (`bufferOK = TRUE`): every row that intersects the study area, only
#'   if together they cover it, gaps closed by a 10 km buffer. A study area that spans several
#'   ELFs returns all their rows.
#' * `writeFitLedger()` writes one row, replacing the ledger's row with the same `polygonID`.
#' * `ledgerRowIndex()` and `ledgerFitForPolygon()` take rows already selected by area and
#'   find the polygon's own row by `polygonID`. The ID is never used to select rows.
#'
#' @param spreadFitFilename,ignitionFitFilename A file name, or `"latest"`.
#' @param cloudFolderID The Google Drive folder (url or id) of the ledger files; `NULL` for a
#'   local ledger (ignition only).
#' @param destinationPath The local folder of the ledger files.
#' @param fireYears The fire years the fit used; only the first and last matter.
#' @param annualNames Spread only: names of the annual covariates (`"year2001"`, ...), for the
#'   fire years when `fireYears` is `NULL` or has no year.
#' @param fitFileSuffix The `fitFileSuffix` parameter of the fit modules.
#' @param use `"read"` or `"write"`; matters for a spread ledger named `"latest"` only.
#' @param ledger A `CacheGeoLedger`.
#' @param row One ledger row to write: an `sf` object with `polygonID` ([polygonIDTxt]) and the fit columns.
#' @param area The study area (`sf`, `sfc` or `SpatVector`).
#' @param rows Ledger rows already selected for the study area, as [readFitLedger()] returns.
#' @param polygonID The polygon (ELF) whose fit is wanted, `sim$.ELFind`.
#' @param doFit `NA`, `TRUE` or `FALSE`, the `doFit` parameter. `TRUE` refits, so no fit is reused.
#'
#' @return `spreadFitLedger()`, `ignitionFitLedger()`: a `CacheGeoLedger`, or `NULL` when a
#'   spread `"latest"` read finds no file. `readFitLedger()`: an `sf` object, or `NULL` when no
#'   row is selected. `writeFitLedger()`: the whole ledger after the write. `ledgerRowIndex()`: the index in `rows` of the polygon's row, or `NA`.
#'   `ledgerFitForPolygon()`: that row (one-row data frame), or `NULL`.
#' @export
#' @rdname fitLedger
spreadFitLedger <- function(spreadFitFilename = "latest", cloudFolderID = spreadFitLedgerFolder(),
                            destinationPath, fireYears = NULL, annualNames = NULL,
                            fitFileSuffix = "", use = c("read", "write")) {
  use <- match.arg(use)
  if (!identical(spreadFitFilename, "latest"))
    return(.fitLedger(spreadFitFilename, cloudFolderID, destinationPath))
  if (identical(use, "write")) {
    if (!length(fireYears) || all(is.na(fireYears)))
      fireYears <- suppressWarnings(as.integer(gsub("\\D", "", annualNames)))
    return(.fitLedger(spreadFitFilenameFor(fireYears, fitFileSuffix), cloudFolderID, destinationPath))
  }
  latest <- latestSpreadFits(cloudFolderID, destinationPath, suffix = fitFileSuffix)
  if (is.null(latest))
    return(NULL)
  ## a local file named by its content: nothing to download, and a changed file has a new name
  tf <- tempfile(fileext = ".rds")
  on.exit(unlink(tf), add = TRUE)
  saveRDS(latest, tf)
  file <- paste0("fireSenseParams_latest_", substr(tools::md5sum(tf), 1, 12), ".rds")
  dir.create(destinationPath, recursive = TRUE, showWarnings = FALSE)
  file.copy(tf, file.path(destinationPath, file), overwrite = TRUE)
  .fitLedger(file, NULL, destinationPath)
}

#' @export
#' @rdname fitLedger
ignitionFitLedger <- function(ignitionFitFilename = "latest", cloudFolderID = ignitionFitLedgerFolder(),
                              destinationPath, fireYears, fitFileSuffix = "") {
  file <- if (identical(ignitionFitFilename, "latest"))
    ignitionFitFilenameFor(fireYears, fitFileSuffix) else ignitionFitFilename
  .fitLedger(file, cloudFolderID, destinationPath)
}

.fitLedger <- function(file, cloudFolderID, destinationPath) {
  if (is.null(cloudFolderID))
    reproducible::CacheGeoLedger(file, destinationPath = destinationPath, key = polygonIDTxt)
  else
    reproducible::CacheGeoLedger(file, remote = cloudFolderID, remoteType = "drive",
                                 destinationPath = destinationPath, key = polygonIDTxt)
}

#' @export
#' @rdname fitLedger
readFitLedger <- function(ledger, area) {
  if (is.null(ledger))
    return(NULL)
  if (inherits(area, "SpatVector")) area <- sf::st_as_sf(area)
  ## the old CacheGeo() call form is the rule of the selection: intersecting rows, used only if
  ## together they cover `area`, gaps closed by a 10 km buffer; NULL if none (see ?CacheGeoLedger)
  rows <- reproducible::CacheGeo(targetFile = ledger$file, cloudFolderID = ledger$remote,
                                 destinationPath = ledger$destinationPath, domain = area,
                                 action = "nothing", bufferOK = TRUE, useCache = FALSE)
  if (NROW(rows)) rows else NULL
}

#' @export
#' @rdname fitLedger
writeFitLedger <- function(ledger, row) {
  if (inherits(row, "SpatVector")) row <- sf::st_as_sf(row)
  domain <- sf::st_sf(geometry = sf::st_geometry(row))
  le <- function(x) x
  ## the old CacheGeo() call form with action = "update": the row replaces the ledger's row with the
  ## same polygonID (the key CacheGeo() sets when the rows have a polygonID column), or is added
  reproducible::CacheGeo(targetFile = ledger$file, cloudFolderID = ledger$remote,
                         destinationPath = ledger$destinationPath, domain = domain,
                         FUN = le(row), le = le, row = row, action = "update")
}

#' @export
#' @rdname fitLedger
ledgerRowIndex <- function(rows, polygonID) {
  if (is.null(polygonID) || !length(polygonID) || !is.data.frame(rows) || !NROW(rows) ||
      !polygonIDTxt %in% names(rows))
    return(NA_integer_)
  match(as.character(polygonID)[1], as.character(rows[[polygonIDTxt]]))
}

#' @export
#' @rdname fitLedger
ledgerFitForPolygon <- function(rows, polygonID, doFit = NA) {
  if (isTRUE(doFit))
    return(NULL)
  i <- ledgerRowIndex(rows, polygonID)
  if (is.na(i)) NULL else rows[i, ]
}

#' The ignition and escape model features, and the covariate names that give them
#'
#' `xgbFeatureNames()` is the one rule for which columns of the fit data `runXGBOOST()` uses as
#' features: every column except `pixelID`, `year`, `yearChar*` and the responses, in alphabetical order.
#'
#' `ignitionFeatureNames()` gives those features from what is known before any covariate is
#' built: the fuel classes (`sppEquiv[[fuelClassCol]]`, an `NA` class is not a fuel), the names of
#' the non-forest land-cover groups, `youngAge`, the ignition climate variables and the lightning
#' layers. It names every feature the covariate tables can have; a fuel class with no cohort in
#' any data year has no column in a table. The ignition prep and the check of a stored fit both
#' use it, and `checkIgnitionFeatureNames()` stops the prep when the table has a feature this
#' function does not give, so the two cannot drift apart.
#'
#' @param cols Column names of the fit data.
#' @param fuelClasses Fuel class of each species (`sppEquiv[[fuelClassCol]]`); `NA` and repeats are dropped.
#' @param nonForestedLCCGroups The non-forest land-cover groups (names are used) or their names.
#' @param missingLCCgroup The name given to forested pixels without cohorts; one of the groups.
#' @param climateVars Climate variable names, without underscores, as in `climateVariablesForFire$ignition`.
#' @param lightning Names of the lightning layers.
#' @param covariates A covariate `data.table` (or its column names).
#' @param expected Feature names, as `ignitionFeatureNames()` returns.
#' @return Character vectors; `checkIgnitionFeatureNames()` returns `covariates` invisibly.
#' @export
#' @rdname ignitionFeatureNames
xgbFeatureNames <- function(cols) {
  cols <- setdiff(cols, c("pixelID", yearTxt))
  cols <- grep("yearChar", cols, invert = TRUE, value = TRUE)
  sort(grep(paste0(ignitionsTxt, "|", escapesTxt), cols, value = TRUE, invert = TRUE))
}

#' @export
#' @rdname ignitionFeatureNames
ignitionFeatureNames <- function(fuelClasses, nonForestedLCCGroups, missingLCCgroup = NULL,
                                 climateVars, lightning = "lightningDays") {
  groups <- if (is.list(nonForestedLCCGroups)) names(nonForestedLCCGroups) else nonForestedLCCGroups
  xgbFeatureNames(unique(c(as.character(stats::na.omit(unique(fuelClasses))), groups, missingLCCgroup,
                           youngAgeTxt, climateVars, lightning)))
}

#' @export
#' @rdname ignitionFeatureNames
checkIgnitionFeatureNames <- function(covariates, expected) {
  cols <- if (is.character(covariates)) covariates else names(covariates)
  extra <- setdiff(xgbFeatureNames(cols), expected)
  if (length(extra))
    stop("The ignition covariates have feature(s) that ignitionFeatureNames() does not give: ",
         paste(extra, collapse = ", "), ". Update ignitionFeatureNames(), which decides whether a stored ",
         "fit can be used without building the covariates.", call. = FALSE)
  invisible(covariates)
}
