#' Get Fire SpatialPoints from Canadian Fire Database
#' @param url Passed to \code{prepInputs}
#' @param studyArea Passed to \code{prepInputs}
#' @param rasterToMatch Passed to \code{prepInputs}
#' @param redownloadIn Numeric Time in YEARS that we tolerate the data to be "old" i.e.
#'   0.5 would mean "redownload data older than 6 months"
#' @param NFDB_pointPath Passed to \code{destinationPath} in \code{prepInputs}
#' @export
#' @return
#' A \code{SpatialPointsDataFrame}.
getFirePoints_NFDB <- function(url = "http://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_pnt/current_version/NFDB_point.zip",
                               studyArea = NULL, rasterToMatch = NULL,
                               redownloadIn = 1,
                               years = 1991:2017,
                               fireSizeColName = "SIZE_HA",
                               NFDB_pointPath # Can't be NULL. Needs to be an existing location for the fire points
                               ){

  check <- Checksums(NFDB_pointPath, checksumFile = file.path(NFDB_pointPath, "CHECKSUMS.txt"), write = TRUE)
  whRowIsShp <- grep("NFDB_point.*shp$", check$expectedFile)
  whIsOK <- which(check$result[whRowIsShp] == "OK")
  needNewDownload <- TRUE
  if (any(whIsOK)) {
    filesToCheck <- tools::file_path_sans_ext(unlist(lapply(check[whRowIsShp[whIsOK], "expectedFile"], as.character)))
    dateOfFile <- substr(x = filesToCheck, start = nchar(filesToCheck) - 8 +
                           1, nchar(filesToCheck))
    if ((as.Date(dateOfFile, format = "%Y%m%d") + dyear(redownloadIn)) > Sys.Date()) {
      # can change dyear(...) to whatever... e.g., dyear(0.5) would be 6 months
      needNewDownload <- FALSE
    }
  }
  if (needNewDownload) {
    print("downloading NFDB")# put prepInputs here
    firePoints <- Cache(prepInputs, url = url,
                        studyArea = studyArea, fun = "shapefile",
                        destinationPath = NFDB_pointPath, useCache = "overwrite",
                        useSAcrs = TRUE, omitArgs = c("NFDB_pointPath", "overwrite"))
  } else {
    NFDBs <- grep(list.files(NFDB_pointPath), pattern = "^NFDB", value = TRUE)
    shps <- grep(list.files(NFDB_pointPath), pattern = ".shp$", value = TRUE)
    aFile <- NFDBs[NFDBs %in% shps][1] #in case there are multiple files
    #firePoints <- Cache(shapefile, file.path(NFDB_pointPath, aFile))
    firePoints <- Cache(sf::read_sf, file.path(NFDB_pointPath, aFile))
    #firePoints1 <- as(firePoints, "Spatial")
    options('reproducible.cacheSaveFormat' = 'rds')
    on.exit({
      options('reproducible.cacheSaveFormat' = 'rds')
    })
    a <- Sys.time()
    firePoints <- Cache(prepInputs, targetFile = file.path(NFDB_pointPath, aFile),
                         destinationPath = NFDB_pointPath,
                         #x = firePoints, fun = sf::read_sf,
                         studyArea = studyArea, filename2 = NULL,
                         rasterToMatch = rasterToMatch,
                         userTags = c("cacheTags", "NFDB"))
  }
  firePoints <- firePoints[firePoints$YEAR <= max(years) &
                                               firePoints$YEAR >= min(years),]
  firePoints <- firePoints[, c("YEAR", fireSizeColName)]
  firePoints$fireSize <- asInteger(firePoints[[fireSizeColName]] / prod(res(rasterToMatch)) * 1e4)
  names(firePoints) <- c("date", "size_ha", "size")

  #    rasterTemp <- setValues(pixelGroupMap2001, values = 1:ncell(pixelGroupMap2001))
  crs(firePoints) <- crs(studyArea)
  return(firePoints)
}