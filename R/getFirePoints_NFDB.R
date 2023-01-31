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
#' @return A \code{SpatialPointsDataFrame}.
#'
#' @export
#' @importFrom raster crs crs<- res
#' @importFrom reproducible Cache Checksums prepInputs
#' @importFrom utils getFromNamespace
getFirePoints_NFDB <- function(url = NULL,
                               studyArea = NULL, rasterToMatch = NULL,
                               redownloadIn = 1,
                               years = 1991:2017,
                               fireSizeColName = "SIZE_HA",
                               NFDB_pointPath = NULL) {
  if (!requireNamespace("SpaDES.core", quietly = TRUE)) {
    stop("Please install.packaes('SpaDES.core').")
  }

  if (is.null(NFDB_pointPath)) {
    stop("NFDB_pointPath must be specified and non-NULL.")
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
    filesToCheck <- getFromNamespace("filePathSansExt", "reproducible")(unlist(lapply( # don't use tools::file_path_sans_ext to keep package number down
      check[whRowIsShp[whIsOK], "expectedFile"], as.character
    )))
    dateOfFile <- substr(
      x = filesToCheck,
      start = nchar(filesToCheck) - 8 + 1, nchar(filesToCheck)
    )
    if ((as.Date(dateOfFile, format = "%Y%m%d") + SpaDES.core::dyear(redownloadIn)) > Sys.Date()) {
      # can change dyear(...) to whatever... e.g., dyear(0.5) would be 6 months
      needNewDownload <- FALSE
    }
  }
  if (needNewDownload) {
    print("downloading NFDB")
    firePoints <- Cache(prepInputs,
      url = url,
      studyArea = studyArea,
      fun = "shapefile",
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
    firePoints <- Cache(sf::read_sf, file.path(NFDB_pointPath, aFile))
    # firePoints1 <- as(firePoints, "Spatial")
    options("reproducible.cacheSaveFormat" = "rds")
    on.exit({
      options("reproducible.cacheSaveFormat" = "rds")
    })
    a <- Sys.time()
    firePoints <- Cache(prepInputs,
      targetFile = file.path(NFDB_pointPath, aFile),
      destinationPath = NFDB_pointPath,
      # x = firePoints, fun = sf::read_sf,
      studyArea = studyArea, filename2 = NULL,
      rasterToMatch = rasterToMatch,
      userTags = c("cacheTags", "NFDB")
    )
  }
  firePoints <- firePoints[firePoints$YEAR <= max(years) &
    firePoints$YEAR >= min(years), ]
  firePoints <- firePoints[, c("YEAR", fireSizeColName)]
  firePoints$fireSize <- asInteger(firePoints[[fireSizeColName]] / prod(res(rasterToMatch)) * 1e4)
  names(firePoints) <- c("date", "size_ha", "size")

  #    rasterTemp <- setValues(pixelGroupMap2001, values = 1:ncell(pixelGroupMap2001))
  crs(firePoints) <- crs(studyArea)
  return(firePoints)
}

#' Get Fire \code{SpatialPoints} from Canadian Fire Database
#'
#' @param url Passed to \code{prepInputs}
#' @template studyArea
#' @template rasterToMatch
#' @param redownloadIn Numeric Time in YEARS that we tolerate the data to be "old" i.e.
#'   0.5 would mean "re-download data older than 6 months"
#' @param years Numeric vector of consecutive years to fetch.
#' @param fireSizeColName Character describing the name of the column containing fire size information.
#' @param NFDB_pointPath Passed to \code{destinationPath} in \code{prepInputs}
#' @param plot logical indicating whether to produce plot of fire points. Default FALSE.
#'
#' @return A \code{SpatialPointsDataFrame}.
#'
#' @export
#' @importFrom crayon green yellow
#' @importFrom raster crs crs<- plot res shapefile
#' @importFrom reproducible Cache Checksums prepInputs projectInputs postProcess
#' @importFrom sp coordinates<-
getFirePoints_NFDB_V2 <- function(url = NULL,
                                  studyArea = NULL,
                                  rasterToMatch = NULL,
                                  redownloadIn = 1,
                                  years = 1991:2017,
                                  fireSizeColName = "SIZE_HA",
                                  NFDB_pointPath = NULL,
                                  plot = FALSE) {
  if (!requireNamespace("SpaDES.core", quietly = TRUE)) {
    stop("Please install.packaes('SpaDES.core').")
  }

  if (is.null(NFDB_pointPath)) {
    stop("NFDB_pointPath must be specified and non-NULL.")
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
    filesToCheck <- getFromNamespace("filePathSansExt", "reproducible")(unlist(lapply( # don't use tools::file_path_sans_ext to keep package number down
      check[whRowIsShp[whIsOK], "expectedFile"], as.character
    )))
    dateOfFile <- substr(
      x = filesToCheck,
      start = nchar(filesToCheck) - 8 + 1, nchar(filesToCheck)
    )
    if ((as.Date(dateOfFile, format = "%Y%m%d") + SpaDES.core::dyear(redownloadIn)) > Sys.Date()) {
      # can change dyear(...) to whatever... e.g., dyear(0.5) would be 6 months
      needNewDownload <- FALSE
    }
  }
  if (needNewDownload) {
    print("downloading NFDB...") # put prepInputs here
    firePoints <- Cache(prepInputs,
      url = url,
      fun = "shapefile",
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
    firePoints <- raster::shapefile(file.path(NFDB_pointPath, aFile))
  }

  # Fix for messed up bbox
  message(crayon::yellow("Correcting original data problem..."))
  DT <- as.data.frame(firePoints@data)
  coordinates(DT) <- cbind(firePoints$LONGITUDE, firePoints$LATITUDE)
  correctCRS <- "+init=epsg:4269 +proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs +towgs84=0,0,0"
  crs(DT) <- correctCRS
  firePointsReady <- projectInputs(DT,
    destinationPath = NFDB_pointPath,
    filename2 = "NFDBpointsProjected",
    targetCRS = crs(rasterToMatch)
  )
  firePoints <- postProcess(firePointsReady, studyArea = studyArea, rasterToMatch = rasterToMatch)
  message(crayon::green("Fire points corrected"))
  if (isTRUE(plot)) {
    raster::plot(firePoints, col = "red")
    raster::plot(studyArea, add = TRUE)
  }
  firePoints <- firePoints[firePoints$YEAR <= max(years) & firePoints$YEAR >= min(years), ]
  firePoints$fireSize <- asInteger(firePoints[[fireSizeColName]] / prod(res(rasterToMatch)) * 1e4)
  return(firePoints)
}
