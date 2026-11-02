## spreadFitFilename = "latest" (latestSpreadFits(), mocked here: no Drive) and a named file. Moved from
## fireSense_dataPrepFit with readSpreadFitLedger().

square <- function(x0, crs = 3978) {
  sf::st_sfc(sf::st_polygon(list(rbind(c(x0, 0), c(x0 + 1e4, 0), c(x0 + 1e4, 1e4), c(x0, 1e4), c(x0, 0)))),
             crs = crs)
}
## a ledger as fireSense_spreadFit writes it: one row per polygon, geometry, crs, parameters
ledger <- function(ids, value) {
  geoms <- do.call(c, lapply(seq_along(ids), function(i) square(i * 2e4)))
  df <- data.frame(polygonID = ids, objFunVal = value, crs = I(rep(sf::st_crs(3978)$wkt, length(ids))))
  df$geometry <- geoms
  df
}

test_that("'latest' reads the combined ledger and keeps only this study area's rows", {
  d <- withr::local_tempdir()
  testthat::local_mocked_bindings(latestSpreadFits = function(cloudFolderID, destinationPath, polygonIDs = NULL, suffix = "")
    ledger(c("4.1", "4.3"), c(1, 2)), .package = "fireSenseUtils")
  domain <- sf::st_as_sf(square(4e4))                       # the second polygon, 4.3
  out <- suppressMessages(readSpreadFitLedger("latest", "folder", domain, d))
  expect_identical(as.character(out$polygonID), "4.3")
  expect_equal(out$objFunVal, 2)
})

test_that("a changed ledger is read again, not returned from an earlier read", {
  d <- withr::local_tempdir()
  domain <- sf::st_as_sf(square(4e4))
  value <- 2
  testthat::local_mocked_bindings(latestSpreadFits = function(cloudFolderID, destinationPath, polygonIDs = NULL, suffix = "")
    ledger(c("4.1", "4.3"), c(1, value)), .package = "fireSenseUtils")
  first <- suppressMessages(readSpreadFitLedger("latest", "folder", domain, d))
  value <- 5                                                 # a newer fit of 4.3 arrives
  second <- suppressMessages(readSpreadFitLedger("latest", "folder", domain, d))
  expect_equal(c(first$objFunVal, second$objFunVal), c(2, 5))
})

test_that("'latest' with no ledger file is no fit", {
  d <- withr::local_tempdir()
  testthat::local_mocked_bindings(latestSpreadFits = function(...) NULL, .package = "fireSenseUtils")
  expect_null(readSpreadFitLedger("latest", "folder", sf::st_as_sf(square(0)), d))
})

test_that("a named file is read from the Drive folder, as before", {
  rec <- new.env()
  local_mocked_bindings(CacheGeo = function(...) { rec$args <- list(...); NULL }, .package = "reproducible")
  readSpreadFitLedger("fireSenseParams_1985-2024_linearFuel.rds", "folderURL", "domain", "dp")
  expect_identical(rec$args$targetFile, "fireSenseParams_1985-2024_linearFuel.rds")
  expect_identical(rec$args$cloudFolderID, "folderURL")
  expect_identical(rec$args$action, "nothing")
})

test_that("'latest' reads only the files with this run's fitFileSuffix", {
  d <- withr::local_tempdir()
  asked <- NULL
  testthat::local_mocked_bindings(latestSpreadFits = function(cloudFolderID, destinationPath, polygonIDs = NULL, suffix = "") {
    asked <<- suffix
    ledger("4.3", 2)
  }, .package = "fireSenseUtils")
  domain <- sf::st_as_sf(square(2e4 * 1))
  suppressMessages(readSpreadFitLedger("latest", "folder", domain, d, fitFileSuffix = "_boundfix_rep81"))
  expect_identical(asked, "_boundfix_rep81")
  suppressMessages(readSpreadFitLedger("latest", "folder", domain, d))
  expect_identical(asked, "")
})
