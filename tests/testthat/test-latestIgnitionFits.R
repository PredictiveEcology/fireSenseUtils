test_that("ignitionFitFilenameFor() names a fit by its fire years and the model tag", {
  expect_identical(ignitionFitFilenameFor(1985:2024), "fireSenseIgnitionParams_1985-2024_xgboostPoisson_20261013.rds")
  expect_identical(ignitionFitFilenameFor(c(NA, 2001, 1990)), "fireSenseIgnitionParams_1990-2001_xgboostPoisson_20261013.rds")
  expect_error(ignitionFitFilenameFor(NA), "no years")
})

test_that("a file with the old Tweedie tag \"_xgboost\" is not a ledger file; the Poisson tag is", {
  d <- withr::local_tempdir()
  drive <- fakeDrive(list(
    "fireSenseIgnitionParams_1985-2025_xgboost.rds" = list(modified = "2026-10-12T10:00:00Z",
                                                      ledger = ledgerRows(c("4.1", "5.1"), c(9, 9))),
    "fireSenseIgnitionParams_1985-2024_xgboostPoisson_20261013.rds" = list(modified = "2026-09-01T10:00:00Z",
                                                      ledger = ledgerRows("4.1", 1))))
  out <- latestIgnitionFits("folder", d)
  expect_identical(out$polygonID, "4.1")
  expect_identical(out$objFunVal, 1)
  expect_identical(drive$downloads, "fireSenseIgnitionParams_1985-2024_xgboostPoisson_20261013.rds")
  ## and the old tag's pattern does not match the new name
  oldPattern <- paste0("^", fireSenseUtils:::ignitionFitFilePrefix, ".*_xgboost\\.rds$")
  expect_false(grepl(oldPattern, ignitionFitFilenameFor(1985:2024)))
  expect_true(grepl(oldPattern, "fireSenseIgnitionParams_1985-2025_xgboost.rds"))
})

test_that("each polygon comes from the newest file that has it; other models' files are ignored", {
  d <- withr::local_tempdir()
  drive <- fakeDrive(list(
    "fireSenseIgnitionParams_1985-2024_xgboostPoisson_20261013.rds" = list(modified = "2026-09-20T10:00:00Z",
                                                      ledger = ledgerRows(c("4.1", "4.3"), c(1, 1))),
    "fireSenseIgnitionParams_1985-2025_xgboostPoisson_20261013.rds" = list(modified = "2026-09-24T10:00:00Z",
                                                      ledger = ledgerRows("4.1", 2)),
    ## an older model's file: has 5.1, but must never be read
    "fireSenseIgnitionParams_1985-2024.rds" = list(modified = "2026-09-25T10:00:00Z",
                                           ledger = ledgerRows(c("4.1", "5.1"), c(9, 9)))))
  out <- latestIgnitionFits("folder", d)
  expect_setequal(out$polygonID, c("4.1", "4.3"))
  expect_identical(out$objFunVal[out$polygonID == "4.1"], 2)      # the newer file's fit
  expect_identical(out$objFunVal[out$polygonID == "4.3"], 1)
  expect_identical(attr(out, "ignitionFitFiles")[["4.3"]], "fireSenseIgnitionParams_1985-2024_xgboostPoisson_20261013.rds")
  expect_false("fireSenseIgnitionParams_1985-2024.rds" %in% drive$downloads)
  expect_s3_class(sf::st_as_sf(out), "sf")                        # still a ledger CacheGeo can read
})

test_that("reading stops at the first file that has every polygon asked for, and reuses local copies", {
  d <- withr::local_tempdir()
  drive <- fakeDrive(list(
    "fireSenseIgnitionParams_1985-2024_xgboostPoisson_20261013.rds" = list(modified = "2026-09-20T10:00:00Z",
                                                      ledger = ledgerRows("4.3", 1)),
    "fireSenseIgnitionParams_1985-2025_xgboostPoisson_20261013.rds" = list(modified = "2026-09-24T10:00:00Z",
                                                      ledger = ledgerRows("4.1", 2))))
  out <- latestIgnitionFits("folder", d, polygonIDs = "4.1")
  expect_identical(out$polygonID, "4.1")
  expect_identical(drive$downloads, "fireSenseIgnitionParams_1985-2025_xgboostPoisson_20261013.rds")
  ## second call: the same file is on disk with Drive's MD5, so nothing is downloaded
  latestIgnitionFits("folder", d, polygonIDs = "4.1")
  expect_identical(drive$downloads, "fireSenseIgnitionParams_1985-2025_xgboostPoisson_20261013.rds")
  ## an ELF in no file: every file is read, and the result has only what exists
  out <- latestIgnitionFits("folder", d, polygonIDs = "9.9")
  expect_setequal(out$polygonID, c("4.1", "4.3"))
})

test_that("no matching file gives NULL", {
  d <- withr::local_tempdir()
  fakeDrive(list("fireSenseIgnitionParams.rds" = list(modified = "2026-09-20T10:00:00Z", ledger = ledgerRows("4.1", 1))))
  expect_null(suppressMessages(latestIgnitionFits("folder", d)))
})

test_that("rows of several ledger files are combined without touching the xgboost models in them", {
  skip_if_not_installed("xgboost")
  set.seed(1)
  newdata <- matrix(rnorm(20), ncol = 2)
  ## a ledger row holds AsIs lists of fitted objects whose boosters, once read back with readRDS(),
  ## are backed by ALTREP lists that cannot be assigned into
  fitOne <- function() {
    x <- matrix(rnorm(200), ncol = 2)
    y <- rpois(100, exp(0.5 * x[, 1]))
    booster <- xgboost::xgboost(x = x, y = y, objective = "count:poisson", nrounds = 3, nthreads = 1)
    list(modelList = list(model = list(Fold1 = booster), rescales = NULL, fittingRes = 1,
                          lambdaRescaleFactor = 1),
         scaleData = list(dim = dim(x)))
  }
  fits <- list(fitOne(), fitOne(), fitOne())
  expected <- lapply(fits, function(f) predict(f$modelList$model$Fold1, newdata))
  ledger <- function(ids, which) {
    df <- ledgerRows(ids, seq_along(ids))
    df$fireSense_IgnitionFitted <- I(fits[which])
    df$fireSense_EscapeFitted <- I(fits[which])
    df
  }
  d <- withr::local_tempdir()
  fakeDrive(list(
    "fireSenseIgnitionParams_1985-2024_xgboostPoisson_20261013.rds" = list(modified = "2026-09-20T10:00:00Z",
                                                      ledger = ledger(c("4.3", "5.1"), 1:2)),
    "fireSenseIgnitionParams_1985-2025_xgboostPoisson_20261013.rds" = list(modified = "2026-09-24T10:00:00Z",
                                                      ledger = ledger("4.1", 3))))
  out <- latestIgnitionFits("folder", d)
  expect_identical(out$polygonID, c("4.1", "4.3", "5.1"))
  expect_s3_class(out$fireSense_IgnitionFitted, "AsIs")
  expect_s3_class(sf::st_as_sf(out), "sf")
  expect_identical(.row_names_info(out), -3L)     # automatic row names
  for (col in c("fireSense_IgnitionFitted", "fireSense_EscapeFitted")) {
    got <- lapply(out[[col]], function(f) predict(f$modelList$model$Fold1, newdata))
    expect_equal(got, expected[c(3, 1, 2)])
  }
})
