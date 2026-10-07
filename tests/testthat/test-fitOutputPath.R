test_that("fitOutputPath() names a folder for the polygon and the ledger file, and creates nothing", {
  d <- withr::local_tempdir()
  p <- fitOutputPath(d, "14.3", spreadFitFilenameFor(1985:2024))
  expect_identical(p, file.path(d, "fits", "14.3_fireSenseParams_1985-2024_linearFuel_esc50"))
  expect_false(dir.exists(file.path(d, "fits")))
  expect_identical(basename(fitOutputPath(d, "14.3", ignitionFitFilenameFor(1985:2024))),
                   "14.3_fireSenseIgnitionParams_1985-2024_xgboost")
  ## the scenario does not enter: the same polygon and fit give the same folder
  expect_identical(fitOutputPath(d, "14.3", "a/b/fit.rds"), file.path(d, "fits", "14.3_fit"))
  expect_false(identical(fitOutputPath(d, "14.3", "x.rds"), fitOutputPath(d, "14.4", "x.rds")))
})
