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
  expect_identical(attr(curlGetHeaders(fireSenseUtils:::nfdbPointUrl(), timeout = 30), "status"), 200L)
})
