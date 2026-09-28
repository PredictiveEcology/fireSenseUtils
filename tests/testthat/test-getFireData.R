## Fire records: NFDB points (downloaded, or reused when a recent copy is on disk), NBAC polygons by
## year, and makeLociList(). Downloads are stubbed.

nfdbPoints <- function() {
  sf::st_as_sf(data.frame(FIRE_ID = 1:4, YEAR = c(1990, 1995, 2000, 2020), SIZE_HA = c(1, 10, 100, 1000),
                          x = c(100, 300, 500, 700), y = 500),
               coords = c("x", "y"), crs = 3978)
}

## a recent NFDB shapefile, named as the NFDB names them: NFDB_point_<yyyymmdd>
writeNFDB <- function(dir, date = Sys.Date() - 10) {
  shp <- file.path(dir, paste0("NFDB_point_", format(date, "%Y%m%d"), ".shp"))
  sf::st_write(nfdbPoints(), shp, quiet = TRUE)
  shp
}

localNFDBoptions <- function(env = parent.frame()) {
  withr::local_options(reproducible.cachePath = withr::local_tempdir(.local_envir = env),
                       reproducible.useCache = FALSE, reproducible.verbose = -2, .local_envir = env)
}

test_that("getFirePoints_NFDB downloads when there is no recent copy, and keeps the asked years", {
  localNFDBoptions()
  dir <- withr::local_tempdir()
  got <- NULL
  local_mocked_bindings(prepInputs = function(url, ...) {
    got <<- url
    nfdbPoints()
  })
  rtm <- terra::rast(nrows = 10, ncols = 10, extent = c(0, 1000, 0, 1000), crs = "EPSG:3978")

  out <- suppressWarnings(capture.output(
    pts <- getFirePoints_NFDB(url = "https://x/NFDB.zip", years = 1991:2017, rasterToMatch = rtm,
                              NFDB_pointPath = dir)))

  expect_identical(got, "https://x/NFDB.zip")
  expect_identical(pts$YEAR, c(1995, 2000))
  ## fire size in cells: 10 ha and 100 ha on 100 x 100 m cells (1 ha)
  expect_identical(pts$size, c(10L, 100L))
})

test_that("getFirePoints_NFDB reuses a copy on disk that is newer than redownloadIn", {
  localNFDBoptions()
  dir <- withr::local_tempdir()
  writeNFDB(dir)
  local_mocked_bindings(prepInputs = function(...) stop("should not download"))

  pts <- getFirePoints_NFDB(url = "https://x/NFDB.zip", years = 1991:2017, NFDB_pointPath = dir)

  expect_identical(sort(pts$YEAR), c(1995, 2000))
})

test_that("getFirePoints_NFDB_V2 needs a path, and reads a local copy with `fun`", {
  expect_error(getFirePoints_NFDB_V2(), "NFDB_pointPath cannot be NULL")
  localNFDBoptions()
  dir <- withr::local_tempdir()
  writeNFDB(dir)
  local_mocked_bindings(prepInputs = function(...) stop("should not download"))

  out <- capture.output(
    pts <- getFirePoints_NFDB_V2(url = "https://x/NFDB.zip", years = 1991:2017, NFDB_pointPath = dir,
                                 fun = "terra::vect"))

  expect_s4_class(pts, "SpatVector")
  expect_identical(sort(pts$YEAR), c(1995, 2000))
})

test_that("getFirePoints_NFDB_V2 downloads when the copy on disk is too old", {
  localNFDBoptions()
  dir <- withr::local_tempdir()
  writeNFDB(dir, date = Sys.Date() - 800)
  local_mocked_bindings(prepInputs = function(url, fun, ...) nfdbPoints())

  out <- capture.output(
    pts <- getFirePoints_NFDB_V2(url = "https://x/NFDB.zip", years = 2000:2030, NFDB_pointPath = dir))

  expect_identical(pts$YEAR, c(2000, 2020))
})

test_that("getFirePolygons splits polygons by year, with their area in ha", {
  polys <- sf::st_buffer(nfdbPoints(), 100)
  polys$YEAR <- as.character(c(2001, 2001, 2003, 2003))
  local_mocked_bindings(prepInputs = function(url, useCache, ...) polys)

  out <- getFirePolygons(url = "https://x/NBAC.zip", years = 2001:2003)

  expect_named(out, c("year2001", "year2002", "year2003"))
  expect_null(out$year2002)
  expect_identical(out$year2001$FIRE_ID, 1:2)
  ## sf gave m2: st_area() has no `unit` argument
  expect_equal(out$year2003$POLY_HA, rep(3.14, 2), tolerance = 0.01)

  local_mocked_bindings(prepInputs = function(url, useCache, ...) terra::vect(polys))
  vout <- getFirePolygons(url = "https://x/NBAC.zip", years = 2003)
  expect_equal(vout$year2003$POLY_HA, rep(3.14, 2), tolerance = 0.01)
})

test_that("latestNBACUrl falls back to the known file when the listing cannot be read", {
  expect_identical(suppressWarnings(latestNBACUrl(base = withr::local_tempdir(), fallback = "fb")), "fb")
})

test_that("makeLociList gives each year's fires as cells, sizes in cells and ids", {
  ras <- terra::rast(nrows = 10, ncols = 10, extent = c(0, 1000, 0, 1000), crs = "EPSG:3978", vals = 1)
  pts <- nfdbPoints()
  names(pts)[names(pts) == "SIZE_HA"] <- "POLY_HA"
  byYear <- split(pts, pts$YEAR)

  out <- makeLociList(ras, byYear)

  expect_named(out, paste0("year", c(1990, 1995, 2000, 2020)))
  expect_identical(out$year2000$ids, 3L)
  expect_identical(out$year2000$size, 100L)                       # 100 ha on 1-ha cells
  expect_equal(out$year2000$cells, terra::cellFromXY(ras, cbind(500, 500)))
  expect_named(makeLociList(ras, lapply(byYear, terra::vect)), names(out))
  expect_error(makeLociList(ras, byYear, sizeColUnits = "acres"), "either ha or m2")
})
