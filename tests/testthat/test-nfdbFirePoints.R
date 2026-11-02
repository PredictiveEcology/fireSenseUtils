## NFDB fire points (moved from fireSense_dataPrepFit with the functions). fireRecordShapefile() reads a
## local archive here: `url = NULL` with `archive`.

## A zipped shapefile beside a non-spatial file, as in the CFS archives.
fireArchive <- function(v, name, dir) {
  src <- file.path(dir, "src")
  dir.create(src, showWarnings = FALSE)
  terra::writeVector(v, file.path(src, paste0(name, ".shp")), overwrite = TRUE)
  writeLines("[fires]", file.path(src, "schema.ini"))
  zipFile <- file.path(dir, paste0(name, ".zip"))
  withr::with_dir(src, utils::zip(zipFile, dir(src), flags = "-q"))
  zipFile
}

crs <- "EPSG:3978"
studyArea <- terra::vect(terra::ext(0, 10000, 0, 10000), crs = crs)

test_that("nfdbFirePoints keeps fires of every size and cause in the study area and years", {
  pts <- terra::vect(cbind(c(1000, 2000, 3000, 4000, 20000), c(1000, 2000, 3000, 4000, 20000)), crs = crs)
  pts$FIRE_ID <- c("a", "b", "c", "d", "e")
  pts$YEAR <- c(1995, 2002, 2003, 2010, 2002)
  pts$SIZE_HA <- c(0, 0.1, 50, 1, 5)
  pts$CAUSE <- c("L", "N", "H", "L", "L")
  dir <- withr::local_tempdir()
  shp <- fireRecordShapefile(url = NULL, destinationPath = dir, archive = fireArchive(pts, "NFDB_point_test", dir))

  out <- nfdbFirePoints(shp, years = 1995:2005, studyArea = studyArea)

  ## "d" is from 2010, "e" is outside the study area
  expect_identical(sort(out$FIRE_ID), c("a", "b", "c"))
  expect_equal(out$SIZE_HA[order(out$FIRE_ID)], c(0, 0.1, 50))
  expect_identical(out$CAUSE[order(out$FIRE_ID)], c("L", "N", "H"))
})

## NFDB has rows that repeat another row's date and location, under another FIRE_ID (554 natural
## fires in release 20260811, 319 in BC); each was counted as a separate ignition.
test_that("dropRepeatedFirePoints keeps one row, the largest, per date and location", {
  pts <- terra::vect(cbind(c(1000, 1000, 1000, 2000, 1000), c(1000, 1000, 1000, 2000, 1000)), crs = crs)
  pts$YEAR <- c(2002, 2002, 2002, 2002, 2003)
  pts$MONTH <- c(7, 7, 7, 7, 7)
  pts$DAY <- c(16, 16, 17, 16, 16)
  pts$SIZE_HA <- c(0.1, 80, 0.1, 0.1, 0.1)
  pts$FIRE_ID <- c("VQ029", "VQ030", "VQ031", "VQ032", "VQ001")

  out <- dropRepeatedFirePoints(pts)

  expect_identical(out$FIRE_ID, c("VQ030", "VQ031", "VQ032", "VQ001"))
})

test_that("pointsOfCauses: natural keeps L and N; adding unknown keeps U, blank and missing", {
  pts <- terra::vect(cbind(1:7, 1), crs = "EPSG:3978")
  pts$CAUSE <- c("L", "N", "H", "U", "", NA, "H-PB")
  expect_identical(pointsOfCauses(pts, "natural")$CAUSE, c("L", "N"))
  expect_equal(NROW(pointsOfCauses(pts, c("natural", "unknown"))), 5)
  expect_identical(pointsOfCauses(pts, c("natural", "unknown"))$CAUSE[1:4], c("L", "N", "U", ""))
})

fires <- data.frame(SIZE_HA = c(0.1, 3, 6, 12, 49.9, 50, 51, 800), id = 1:8)

test_that("the sizes of the fires that did not escape are kept, in ha", {
  expect_equal(nonEscapedFireSizes(fires, 50), c(0.1, 3, 6, 12, 49.9))
  expect_equal(nonEscapedFireSizes(data.frame(SIZE_HA = c(NA, 0, 5)), 50), 5)
})

test_that("nfdbFirePointsOfCauses: the study area's fires of the years and causes, one row per fire", {
  pts <- terra::vect(cbind(c(1000, 1000, 3000, 4000, 20000), c(1000, 1000, 3000, 4000, 20000)), crs = crs)
  pts$FIRE_ID <- c("a", "a2", "c", "d", "e")
  pts$YEAR <- c(2002, 2002, 2003, 2003, 2002)
  pts$MONTH <- 7; pts$DAY <- 16
  pts$SIZE_HA <- c(0.1, 30, 50, 1, 5)
  pts$CAUSE <- c("L", "L", "H", "N", "L")
  dir <- withr::local_tempdir()
  shp <- fireRecordShapefile(url = NULL, destinationPath = dir, archive = fireArchive(pts, "NFDB_point_test", dir))
  out <- nfdbFirePointsOfCauses(shp, studyArea = studyArea, years = 1995:2005, causes = "natural")
  ## "a" repeats "a2" (kept: larger), "c" is human-caused, "e" is outside the study area
  expect_setequal(out$FIRE_ID, c("a2", "d"))
  expect_equal(nonEscapedFireSizes(out, 50)[order(out$FIRE_ID)], c(30, 1))
})
