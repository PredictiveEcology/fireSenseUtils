## fireRecordShapefile() returns the .shp among the files reproducible::preProcess() extracts.
## preProcess() downloads, so it is replaced and its arguments recorded.

test_that("the shapefile is picked out of the extracted files", {
  seen <- NULL
  local_mocked_bindings(preProcess = function(...) {
    seen <<- list(...)
    list(targetFilePath = c("/in/NFDB_point_20250519.dbf", "/in/NFDB_point_20250519.shp",
                            "/in/NFDB_point_20250519.shp.xml", "/in/NFDB_point_20250519.shx"))
  }, .package = "reproducible")

  out <- fireRecordShapefile("https://example.org/NFDB_point.zip", destinationPath = "/in")
  expect_identical(out, "/in/NFDB_point_20250519.shp")        # not .shp.xml, not .shx
  expect_identical(seen$url, "https://example.org/NFDB_point.zip")
  expect_identical(seen$destinationPath, "/in")
  expect_identical(seen$fun, NA)                               # extract only; do not load
})

test_that("extra arguments reach preProcess()", {
  seen <- NULL
  local_mocked_bindings(preProcess = function(...) {
    seen <<- list(...)
    list(targetFilePath = "a.shp")
  }, .package = "reproducible")
  fireRecordShapefile("u", destinationPath = "d", overwrite = TRUE, purge = 7)
  expect_identical(seen$overwrite, TRUE)
  expect_identical(seen$purge, 7)
})
