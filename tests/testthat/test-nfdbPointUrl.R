## getFirePoints_NFDB() and getFirePoints_NFDB_V2() used `.../current_version/NFDB_point.zip`, which CFS
## renamed to `NFDB_point_shp.zip`; the old name returns HTTP 404, so no release after the copy on
## disk (fires to 2024) could be downloaded.
test_that("the NFDB point default URL is the full shapefile archive", {
  url <- fireSenseUtils:::nfdbPointUrl()
  expect_match(url, "/current_version/NFDB_point_shp\\.zip$")
  expect_no_match(url, "large_fires")
})

test_that("the NFDB point default URL is reachable", {
  skip_on_cran()
  skip_if_offline("cwfis.cfs.nrcan.gc.ca")
  skip_if_not_installed("httr2")
  resp <- httr2::request(fireSenseUtils:::nfdbPointUrl()) |>
    httr2::req_method("HEAD") |>
    httr2::req_error(is_error = function(resp) FALSE) |>
    httr2::req_timeout(30) |>
    httr2::req_perform()
  expect_identical(httr2::resp_status(resp), 200L)
})
