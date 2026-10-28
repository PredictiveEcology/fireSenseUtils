## Toy ledger: three ELFs as 10 km squares; E1 and E2 share an edge, E3 is apart (gap 20-30 km)
crsToy <- "EPSG:32618"
sq <- function(id, x0, x1, y0 = 0, y1 = 10) {
  box <- sf::st_polygon(list(rbind(c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1), c(x0, y0))))
  out <- sf::st_sf(polygonID = id, geometry = sf::st_sfc(box * 1000, crs = crsToy))
  out$params <- I(list(list(p = id)))
  out
}
toyLedger <- function(dir) {
  led <- ignitionFitLedger("toy.rds", cloudFolderID = NULL, destinationPath = dir, fireYears = 1986:2024)
  reproducible::CacheGeoWrite(led, rbind(sq("E1", 0, 10), sq("E2", 10, 20), sq("E3", 30, 40)), verbose = 0)
  led
}
area <- function(x0, x1) sq("area", x0, x1, 2, 8)[, "polygonID"]
ids <- function(x) if (is.null(x)) NULL else sort(x$polygonID)

test_that("ledger file names: named as given, latest from the fire years and suffix", {
  d <- withr::local_tempdir()
  expect_identical(ignitionFitLedger("mine.rds", NULL, d, 1986:2024)$file, "mine.rds")
  expect_identical(ignitionFitLedger("latest", NULL, d, c(1986, 2024), "_x")$file,
                   ignitionFitFilenameFor(1986:2024, "_x"))
  l <- spreadFitLedger("latest", NULL, d, fireYears = c(1986, 2024), fitFileSuffix = "_x", use = "write")
  expect_identical(l$file, spreadFitFilenameFor(1986:2024, "_x"))
  ## the window from the annual covariate names when there are no fire years
  l <- spreadFitLedger("latest", NULL, d, fireYears = NULL, annualNames = c("year2001", "year2003"), use = "write")
  expect_identical(l$file, spreadFitFilenameFor(2001:2003))
  expect_identical(spreadFitLedger("mine.rds", NULL, d, use = "write")$file, "mine.rds")
  expect_identical(l$key, polygonIDTxt)
  expect_identical(ignitionFitLedger("x.rds", NULL, d, 1986)$remote, NULL)
  r <- ignitionFitLedger("x.rds", "https://drive.google.com/drive/folders/abc", d, 1986)
  expect_identical(r$remoteType, "drive")
})

test_that("readFitLedger: the legacy selection, several ELFs together", {
  testthat::skip_if_not_installed("reproducible", "3.2.1.9075")
  d <- withr::local_tempdir()
  led <- toyLedger(d)
  suppressMessages({
    ## a study area spanning two ELFs returns both rows
    expect_identical(ids(readFitLedger(led, area(5, 15))), c("E1", "E2"))
    ## inside one ELF: that row
    expect_identical(ids(readFitLedger(led, area(2, 8))), "E1")
    ## the 10 km buffer closes the 20-30 km gap; beyond the last row + 10 km it is not covered: NULL
    expect_identical(ids(readFitLedger(led, area(5, 45))), c("E1", "E2", "E3"))
    expect_null(readFitLedger(led, area(5, 60)))
    ## a SpatVector area is accepted
    expect_identical(ids(readFitLedger(led, terra::vect(area(5, 15)))), c("E1", "E2"))
  })
  ## a missing ledger file: NULL
  none <- ignitionFitLedger("none.rds", NULL, d, 1986)
  expect_null(suppressMessages(readFitLedger(none, area(5, 15))))
  expect_null(readFitLedger(NULL, area(5, 15)))
})

test_that("ledgerFitForPolygon takes the own row from rows already selected by area", {
  testthat::skip_if_not_installed("reproducible", "3.2.1.9075")
  d <- withr::local_tempdir()
  led <- toyLedger(d)
  rows <- suppressMessages(readFitLedger(led, area(5, 15)))
  expect_identical(ledgerFitForPolygon(rows, "E2")$polygonID, "E2")
  expect_identical(ledgerRowIndex(rows, "E2"), match("E2", rows$polygonID))
  ## doFit = TRUE reuses nothing; FALSE and NA reuse
  expect_null(ledgerFitForPolygon(rows, "E2", doFit = TRUE))
  expect_identical(ledgerFitForPolygon(rows, "E2", doFit = FALSE)$polygonID, "E2")
  ## a polygon with no row, no rows, no id
  expect_null(ledgerFitForPolygon(rows, "E9"))
  expect_null(ledgerFitForPolygon(NULL, "E1"))
  expect_true(is.na(ledgerRowIndex(rows, NULL)))
  ## the ID never selects: E3 is in the ledger but not selected for this area, so it is not found
  expect_null(ledgerFitForPolygon(rows, "E3"))
})

test_that("a row with the right polygonID but a geometry elsewhere is not used", {
  testthat::skip_if_not_installed("reproducible", "3.2.1.9075")
  d <- withr::local_tempdir()
  led <- toyLedger(d)
  ## E1's id on a polygon far from the study area
  reproducible::CacheGeoWrite(led, sq("E1", 100, 110), verbose = 0)
  suppressMessages(rows <- readFitLedger(led, area(2, 8)))
  expect_null(ledgerFitForPolygon(rows, "E1"))
})

test_that("a row within the 10 km buffer of the study area is used", {
  testthat::skip_if_not_installed("reproducible", "3.2.1.9075")
  d <- withr::local_tempdir()
  led <- ignitionFitLedger("b.rds", NULL, d, 1986)
  reproducible::CacheGeoWrite(led, sq("E1", 0, 10), verbose = 0)
  ## the study area runs 8 km beyond the ledger polygon: covered with the buffer
  suppressMessages(rows <- readFitLedger(led, sq("a", 2, 18, 2, 8)[, "polygonID"]))
  expect_identical(ledgerFitForPolygon(rows, "E1")$polygonID, "E1")
  ## 20 km beyond: not covered
  suppressMessages(rows <- readFitLedger(led, sq("a", 2, 30, 2, 8)[, "polygonID"]))
  expect_null(rows)
})

test_that("ignitionFeatureNames gives the features the covariate table can have", {
  nf <- list(nfLCC_a = c(1, 2), nfLCC_b = 3)
  nm <- ignitionFeatureNames(c("Pice_mar", "Pice_gla", NA, "Pice_mar"), nf, "nfLCC_a",
                             c("CMD", "cumMDC"))
  expect_identical(nm, sort(c("Pice_gla", "Pice_mar", "nfLCC_a", "nfLCC_b", "youngAge", "CMD", "cumMDC",
                              "lightningDays")))
  tab <- data.table::data.table(pixelID = 1, ignitions = 0, year = 2001, youngAge = 0, CMD = 1,
                                nfLCC_a = 0, Pice_mar = 1, lightningDays = 3)
  expect_identical(checkIgnitionFeatureNames(tab, nm), tab)
  tab[, Surprise := 1]
  expect_error(checkIgnitionFeatureNames(tab, nm), "Surprise")
  ## responses, pixelID, year and yearChar are not features
  expect_identical(xgbFeatureNames(c("pixelID", "year", "yearChar", "ignitions", "escapes", "CMD")), "CMD")
})

test_that("writeFitLedger replaces the row of the same polygonID and adds new ones", {
  testthat::skip_if_not_installed("reproducible", "3.2.1.9075")
  d <- withr::local_tempdir()
  led <- ignitionFitLedger("w.rds", NULL, d, 1986)
  suppressMessages({
    writeFitLedger(led, sq("E1", 0, 10))
    all <- writeFitLedger(led, sq("E2", 10, 20))
    expect_identical(ids(all), c("E1", "E2"))
    r1 <- sq("E1", 0, 10); r1$params <- I(list(list(p = "new")))
    all <- writeFitLedger(led, r1)
    expect_identical(NROW(all), 2L)
    expect_identical(all$params[all$polygonID == "E1"][[1]], list(p = "new"))
    expect_identical(ids(readFitLedger(led, area(2, 8))), "E1")
  })
})
