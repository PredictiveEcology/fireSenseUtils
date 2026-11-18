test_that("ignitionFitFilenameFor() names a fit by its fire years and the model tag", {
  expect_identical(ignitionFitFilenameFor(1985:2024), "fireSenseIgnitionParams_1985-2024_xgboostPoisson_20261013.rds")
  expect_identical(ignitionFitFilenameFor(c(NA, 2001, 1990)), "fireSenseIgnitionParams_1990-2001_xgboostPoisson_20261013.rds")
  expect_error(ignitionFitFilenameFor(NA), "no years")
})

test_that("a row fitted with another objective is skipped, whatever its file's name", {
  ## the old Tweedie fits ("_xgboost" files) are found by their models, not by the tag in the file name
  skip_if_not_installed("xgboost")
  withModels <- function(ids, value, objective, writtenAt) {
    r <- ledgerRows(ids, value, writtenAt)
    fit <- storedFit(objective)
    r$fireSense_IgnitionFitted <- I(rep(list(fit), length(ids)))
    r$fireSense_EscapeFitted <- I(rep(list(fit), length(ids)))
    r
  }
  d <- withr::local_tempdir()
  fakeDrive(list(
    "fireSenseIgnitionParams_1985-2025_xgboost.rds" = list(modified = "2026-10-12T10:00:00Z",
                                                      ledger = withModels(c("4.1", "5.1"), c(9, 9), "reg:tweedie", "2026-10-12")),
    "fireSenseIgnitionParams_1985-2024_xgboostPoisson_20261013.rds" = list(modified = "2026-09-01T10:00:00Z",
                                                      ledger = withModels("4.1", 1, "count:poisson", "2026-09-01")),
    "fireSenseIgnitionParams_myLabel.rds" = list(modified = "2026-10-13T10:00:00Z",
                                                      ledger = withModels("6.1", 2, "count:poisson", "2026-10-13"))))
  expect_message(out <- latestIgnitionFits("folder", d),
                 "skipped the fit of 4.1 in fireSenseIgnitionParams_1985-2025_xgboost.rds: its IgnitionFitted model was fitted with reg:tweedie")
  expect_setequal(out$polygonID, c("4.1", "6.1"))                 # 5.1 has only a Tweedie fit
  expect_identical(out$objFunVal[out$polygonID == "4.1"], 1)
  expect_identical(attr(out, "ignitionFitFiles")[["6.1"]], "fireSenseIgnitionParams_myLabel.rds")
})

test_that("each polygon comes from the newest file that has it", {
  d <- withr::local_tempdir()
  drive <- fakeDrive(list(
    "fireSenseIgnitionParams_1985-2024_xgboostPoisson_20261013.rds" = list(modified = "2026-09-20T10:00:00Z",
                                                      ledger = ledgerRows(c("4.1", "4.3"), c(1, 1))),
    "fireSenseIgnitionParams_1985-2025_xgboostPoisson_20261013.rds" = list(modified = "2026-09-24T10:00:00Z",
                                                      ledger = ledgerRows("4.1", 2))))
  out <- latestIgnitionFits("folder", d)
  expect_setequal(out$polygonID, c("4.1", "4.3"))
  expect_identical(out$objFunVal[out$polygonID == "4.1"], 2)      # the newer file's fit
  expect_identical(out$objFunVal[out$polygonID == "4.3"], 1)
  expect_identical(attr(out, "ignitionFitFiles")[["4.3"]], "fireSenseIgnitionParams_1985-2024_xgboostPoisson_20261013.rds")
  expect_s3_class(sf::st_as_sf(out), "sf")                        # still a ledger CacheGeo can read
})

test_that("a row fitted under another escape definition is skipped when the run gives its own", {
  d <- withr::local_tempdir()
  esc <- function(ids, value, size, writtenAt) cbind(ledgerRows(ids, value, writtenAt), as.data.frame(escapeDefinition(size)))
  fakeDrive(list(
    "fireSenseIgnitionParams_a.rds" = list(modified = "2026-10-01T10:00:00Z", ledger = esc("4.1", 1, 50, "2026-10-01")),
    "fireSenseIgnitionParams_b.rds" = list(modified = "2026-11-01T10:00:00Z", ledger = esc("4.1", 2, 100, "2026-11-01"))))
  expect_message(out <- latestIgnitionFits("folder", d, escapeDefinition = escapeDefinition(50)), "the escape definition")
  expect_identical(out$objFunVal, 1)
  expect_identical(suppressMessages(latestIgnitionFits("folder", d))$objFunVal, 2)
})

test_that("reading stops once each polygon asked for has a row newer than the next file", {
  d <- withr::local_tempdir()
  drive <- fakeDrive(list(
    "fireSenseIgnitionParams_1985-2024_xgboostPoisson_20261013.rds" = list(modified = "2026-09-20T10:00:00Z",
                                                      ledger = ledgerRows("4.3", 1, "2026-09-20 09:00")),
    "fireSenseIgnitionParams_1985-2025_xgboostPoisson_20261013.rds" = list(modified = "2026-09-24T10:00:00Z",
                                                      ledger = ledgerRows("4.1", 2, "2026-09-24 09:00"))))
  out <- latestIgnitionFits("folder", d, polygonIDs = "4.1")
  expect_identical(out$polygonID, "4.1")
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

test_that("a suffix goes after the model tag; \"\" is the production name", {
  expect_identical(ignitionFitFilenameFor(1986:2024, "_boundfix_rep81"),
                   "fireSenseIgnitionParams_1986-2024_xgboostPoisson_20261013_boundfix_rep81.rds")
  expect_identical(ignitionFitFilenameFor(1986:2024, ""), ignitionFitFilenameFor(1986:2024))
})

test_that("\"latest\" reads every suffix, newest date in the suffix first; a suffix reads only its own files", {
  d <- withr::local_tempdir()
  prod <- ignitionFitFilenameFor(1986:2024)                                  # no suffix: no date, oldest
  test81 <- ignitionFitFilenameFor(1986:2024, "_boundfix_rep81")             # no date in the suffix either
  gpsc <- ignitionFitFilenameFor(1986:2024, "_gpscCause_20261029")
  bufI <- ignitionFitFilenameFor(1986:2024, "_bufIndep20261024")
  drive <- fakeDrive(setNames(list(
    list(modified = "2026-11-01T10:00:00Z", ledger = ledgerRows(c("4.1", "7.1"), c(1, 1))),
    list(modified = "2026-10-15T10:00:00Z", ledger = ledgerRows(c("4.1", "5.1", "7.1"), c(81, 81, 81))),
    list(modified = "2026-10-29T10:00:00Z", ledger = ledgerRows("4.1", 29)),
    ## modified after the gpsc file (a fit added later), but its suffix date is older
    list(modified = "2026-11-02T10:00:00Z", ledger = ledgerRows(c("4.1", "5.1"), c(24, 24)))),
    c(prod, test81, gpsc, bufI)))
  out <- latestIgnitionFits("folder", d)
  val <- function(id) out$objFunVal[out$polygonID == id]
  expect_identical(val("4.1"), 29)      # 20261029 beats 20261024, whatever the modified times
  expect_identical(val("5.1"), 24)
  expect_identical(val("7.1"), 1)       # undated files: the later modified first
  ## files are read by modified time (no row is newer than its file), and the date in the suffix
  ## when that is later; no polygonIDs: every file
  expect_identical(drive$downloads, c(bufI, prod, gpsc, test81))
  ## reading stops once every ELF asked for is found
  drive$downloads <- character()
  out <- latestIgnitionFits("folder", d, polygonIDs = c("4.1", "5.1"))
  ## gpsc's ceiling (its 2026-10-29 modified time) is after 4.1's row in bufI (20261024), so it is
  ## read; test81's (2026-10-15) is before every row kept, so reading stops there
  expect_identical(drive$downloads, c(bufI, prod, gpsc))
  ## a suffix reads only its own files
  out <- latestIgnitionFits("folder", d, suffix = "_boundfix_rep81")
  expect_setequal(out$polygonID, c("4.1", "5.1", "7.1"))
  expect_identical(unique(out$objFunVal), 81)
})

test_that("ledgerFileOrder(): the last 8-digit date of the suffix, undated last, then modified time", {
  n <- paste0("fireSenseIgnitionParams_1986-2024_xgboostPoisson_20261013", c("", "_a_20261101", "_x20261024", "_rep3", "_D05_20261102"), ".rds")
  m <- c("2026-11-03", "2026-11-01", "2026-11-04", "2026-11-05", "2026-11-02")
  expect_identical(n[ledgerFileOrder(n, m, "_xgboostPoisson_20261013")], n[c(5, 2, 3, 4, 1)])
})

test_that("latestIgnitionFits() and latestSpreadFits() stop clearly without a Drive folder", {
  expect_error(latestIgnitionFits(NULL, destinationPath = tempdir()), "latestIgnitionFits\\(\\): `cloudFolderID` is NULL")
  expect_error(latestSpreadFits(NULL, destinationPath = tempdir()), "latestSpreadFits\\(\\): `cloudFolderID` is NULL")
})
