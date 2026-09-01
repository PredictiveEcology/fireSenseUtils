#' Get Fire `SpatialPoints` from Canadian Fire Database
#'
#' @param url Passed to `prepInputs`
#' @template studyArea
#' @template rasterToMatch
#' @param redownloadIn Numeric Time in YEARS that we tolerate the data to be "old" i.e.
#'   0.5 would mean "redownload data older than 6 months"
#' @param years Numeric vector of consecutive years to fetch.
#' @param fireSizeColName Character describing the name of the column containing fire size information.
#' @param NFDB_pointPath Passed to `destinationPath` in `prepInputs`
#'
#' @return A `sf` spatial points object, retaining all columns of the source data
#'   (notably `CAUSE`, `YEAR` and `SIZE_HA`). When `rasterToMatch` is supplied, an
#'   additional `size` column gives fire size in pixels.
#'
#' @export
#' @importFrom LandR asInteger
#' @importFrom reproducible Cache Checksums prepInputs
#' @importFrom sf st_read
#' @importFrom terra res
#' @importFrom utils getFromNamespace
getFirePoints_NFDB <- function(url = NULL,
                               studyArea = NULL, rasterToMatch = NULL,
                               redownloadIn = 1,
                               years = 1991:2017,
                               fireSizeColName = "SIZE_HA",
                               NFDB_pointPath) {
  if (is.null(url)) {
    url <- "http://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_pnt/current_version/NFDB_point.zip"
  }

  check <- Checksums(NFDB_pointPath,
    checksumFile = file.path(NFDB_pointPath, "CHECKSUMS.txt"),
    write = TRUE
  )
  whRowIsShp <- grep("NFDB_point.*shp$", check$expectedFile)
  whIsOK <- which(check$result[whRowIsShp] == "OK")
  needNewDownload <- TRUE
  if (any(whIsOK)) {
    filesToCheck <- unlist(lapply(check[whRowIsShp[whIsOK], "expectedFile"], as.character))
    filesToCheck <- getFromNamespace("filePathSansExt", "reproducible")(filesToCheck)
    dateOfFile <- substr(
      x = filesToCheck,
      start = nchar(filesToCheck) - 8 + 1, nchar(filesToCheck)
    )
    if (any((as.Date(dateOfFile, format = "%Y%m%d") + 365 * redownloadIn) > Sys.Date())) {
      # e.g., redownloadIn = 0.5 tolerates data up to 6 months old
      needNewDownload <- FALSE
    }
  }
  if (needNewDownload) {
    print("downloading NFDB")
    firePoints <- Cache(prepInputs,
      url = url,
      studyArea = studyArea,
      fun = "sf::st_read",
      destinationPath = NFDB_pointPath,
      useCache = "overwrite",
      useSAcrs = TRUE,
      omitArgs = c("NFDB_pointPath", "overwrite")
    )
  } else {
    NFDBs <- grep(list.files(NFDB_pointPath), pattern = "^NFDB", value = TRUE)
    shps <- grep(list.files(NFDB_pointPath), pattern = ".shp$", value = TRUE)
    aFile <- NFDBs[NFDBs %in% shps][1] # in case there are multiple files
    firePoints <- Cache(sf::st_read, file.path(NFDB_pointPath, aFile))
    oldOpts <- options("reproducible.cacheSaveFormat" = "rds")
    on.exit(options(oldOpts), add = TRUE)
    firePoints <- Cache(postProcess,
      firePoints,
      studyArea = studyArea,
      useSAcrs = TRUE,
      userTags = c("cacheTags", "NFDB")
    )
  }
  firePoints <- firePoints[firePoints$YEAR <= max(years) &
    firePoints$YEAR >= min(years), ]
  if (!is.null(rasterToMatch)) {
    firePoints$size <- asInteger(firePoints[[fireSizeColName]] / prod(res(rasterToMatch)) * 1e4)
  }

  return(firePoints)
}

#' Get Fire `SpatialPoints` from Canadian Fire Database
#'
#' @param url Passed to `prepInputs`
#' @template studyArea
#' @param redownloadIn Numeric Time in YEARS that we tolerate the data to be "old" i.e.
#'   0.5 would mean "re-download data older than 6 months"
#' @param years Numeric vector of consecutive years to fetch.
#' @param fireSizeColName Character describing the name of the column containing fire size information.
#' @param NFDB_pointPath Passed to `destinationPath` in `prepInputs`
#' @param fun Character. The function (as a string, e.g., `"sf::st_read"`)
#'   used to load the downloaded data. Forwarded to `prepInputs(fun = ...)`, and
#'   used to load already-downloaded data, so the class returned does not depend
#'   on whether a download was needed.
#' @param plot logical indicating whether to produce plot of fire points. Default FALSE.
#'
#' @return A `sf` spatial points object (or whatever class `fun` returns), retaining
#'   all columns of the source data (notably `CAUSE`, `YEAR` and `SIZE_HA`).
#'
#' @export
#' @importFrom crayon green yellow
#' @importFrom sf st_read
#' @importFrom reproducible Cache Checksums prepInputs projectInputs postProcess
getFirePoints_NFDB_V2 <- function(url = NULL,
                                  studyArea = NULL,
                                  redownloadIn = 1,
                                  years = 1991:2017,
                                  fireSizeColName = "SIZE_HA",
                                  NFDB_pointPath = NULL,
                                  fun = "sf::st_read",
                                  plot = FALSE) {
  if (is.null(NFDB_pointPath)) {
    stop("NFDB_pointPath cannot be NULL. Specify a file path.")
  }
  if (is.null(url)) {
    url <- "http://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_pnt/current_version/NFDB_point.zip"
  }
  check <- Checksums(NFDB_pointPath,
    checksumFile = file.path(NFDB_pointPath, "CHECKSUMS.txt"),
    write = TRUE
  )
  whRowIsShp <- grep("NFDB_point.*shp$", check$expectedFile)
  whIsOK <- which(check$result[whRowIsShp] == "OK")
  needNewDownload <- TRUE
  if (any(whIsOK)) {
    filesToCheck <- unlist(lapply(check[whRowIsShp[whIsOK], "expectedFile"], as.character))
    filesToCheck <- getFromNamespace("filePathSansExt", "reproducible")(filesToCheck)
    dateOfFile <- substr(
      x = filesToCheck,
      start = nchar(filesToCheck) - 8 + 1, nchar(filesToCheck)
    )
    if (any((as.Date(dateOfFile, format = "%Y%m%d") + 365 * redownloadIn) > Sys.Date())) {
      needNewDownload <- FALSE
    }
  }
  if (needNewDownload) {
    print("downloading NFDB...") # put prepInputs here
    firePoints <- Cache(prepInputs,
      url = url,
      fun = fun,
      studyArea = studyArea,
      destinationPath = NFDB_pointPath,
      useSAcrs = TRUE,
      omitArgs = c("NFDB_pointPath", "overwrite")
    )
  } else {
    print("NFDB present. Loading...") # put prepInputs here
    NFDBs <- grep(list.files(NFDB_pointPath), pattern = "^NFDB", value = TRUE)
    shps <- grep(list.files(NFDB_pointPath), pattern = ".shp$", value = TRUE)
    aFile <- NFDBs[NFDBs %in% shps][1] # in case there are multiple files
    ## use `fun` here too, so the returned class is the same in both branches
    loadFun <- if (grepl("::", fun, fixed = TRUE)) {
      getExportedValue(sub("::.*$", "", fun), sub("^.*::", "", fun))
    } else {
      match.fun(fun)
    }
    firePoints <- loadFun(file.path(NFDB_pointPath, aFile))
    firePoints <- postProcess(firePoints, studyArea = studyArea, useSAcrs = TRUE)
  }

  firePoints <- firePoints[firePoints$YEAR <= max(years) & firePoints$YEAR >= min(years), ]

  return(firePoints)
}
