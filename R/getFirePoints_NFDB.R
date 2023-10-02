#' Get Fire \code{SpatialPoints} from Canadian Fire Database
#'
#' @param url Passed to \code{prepInputs}
#' @template studyArea
#' @template rasterToMatch
#' @param redownloadIn Numeric Time in YEARS that we tolerate the data to be "old" i.e.
#'   0.5 would mean "redownload data older than 6 months"
#' @param years Numeric vector of consecutive years to fetch.
#' @param fireSizeColName Character describing the name of the column containing fire size information.
#' @param NFDB_pointPath Passed to \code{destinationPath} in \code{prepInputs}
#'
#' @return A \code{sf} spatial points object.
#'
#' @export
#' @importFrom terra res
#' @importFrom sf st_read
#' @importFrom reproducible Cache Checksums prepInputs
#' @importFrom utils getFromNamespace
getFirePoints_NFDB <- function(url = NULL,
                               studyArea = NULL, rasterToMatch = NULL,
                               redownloadIn = 1,
                               years = 1991:2017,
                               fireSizeColName = "SIZE_HA",
                               NFDB_pointPath) {
  if (!requireNamespace("SpaDES.core", quietly = TRUE)) {
    stop("Please install.packaes('SpaDES.core').")
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
    if (as.Date(dateOfFile, format = "%Y%m%d") + 365/redownloadIn > Sys.Date()) {
      # can change dyear(...) to whatever... e.g., dyear(0.5) would be 6 months
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
    # firePoints <- Cache(shapefile, file.path(NFDB_pointPath, aFile))
    firePoints <- Cache(sf::st_read, file.path(NFDB_pointPath, aFile))
    # firePoints1 <- as(firePoints, "Spatial")
    options("reproducible.cacheSaveFormat" = "rds")
    on.exit({
      options("reproducible.cacheSaveFormat" = "rds")
    })
    a <- Sys.time()
    firePoints <- Cache(postProcess,
                        firePoints,
                        studyArea = studyArea,
                        userTags = c("cacheTags", "NFDB")
    )
  }
  firePoints <- firePoints[firePoints$YEAR <= max(years) &
                             firePoints$YEAR >= min(years), ]
  firePoints <- firePoints[, c("YEAR", fireSizeColName)]
  if (!is.null(rasterToMatch)){
    firePoints$size <- asInteger(firePoints[[fireSizeColName]] / prod(res(rasterToMatch)) * 1e4)
  }
  names(firePoints)[1:2] <- c("date", "size_ha")

  return(firePoints)
}

#' Get Fire \code{SpatialPoints} from Canadian Fire Database
#'
#' @param url Passed to \code{prepInputs}
#' @template studyArea
#' @param redownloadIn Numeric Time in YEARS that we tolerate the data to be "old" i.e.
#'   0.5 would mean "re-download data older than 6 months"
#' @param years Numeric vector of consecutive years to fetch.
#' @param fireSizeColName Character describing the name of the column containing fire size information.
#' @param NFDB_pointPath Passed to \code{destinationPath} in \code{prepInputs}
#' @param plot logical indicating whether to produce plot of fire points. Default FALSE.
#'
#' @return A \code{sf} spatial points object.
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
                                  plot = FALSE) {
  if (is.null(url)) {
    url <- "http://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_pnt/current_version/NFDB_point.zip"
  }
  check <- Checksums(NFDB_pointPath,
                     checksumFile = file.path(NFDB_pointPath, "CHECKSUMS.txt"),
                     write = TRUE)
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
    if (as.Date(dateOfFile, format = "%Y%m%d") + 365/redownloadIn > Sys.Date()) {
      needNewDownload <- FALSE
    }
  }
  if (needNewDownload) {
    print("downloading NFDB...") # put prepInputs here
    firePoints <- Cache(prepInputs,
                        url = url,
                        fun = "st_read",
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
    firePoints <- st_read(file.path(NFDB_pointPath, aFile))
    firePoints <- postProcess(firePoints, studyArea = studyArea, useSAcrs = TRUE)
  }

  firePoints <- firePoints[firePoints$YEAR <= max(years) & firePoints$YEAR >= min(years), ]

  return(firePoints)
}
