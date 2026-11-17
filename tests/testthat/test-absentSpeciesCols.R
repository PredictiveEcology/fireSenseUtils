## A stored model trained with a species column that this run's vegetation lacks: the species has zero
## biomass there, so the column is 0, not "incompatible".
speciesFit <- function() {
  set.seed(3)
  x <- matrix(rnorm(300), ncol = 3, dimnames = list(NULL, c("CMD", "Pinu_con", "youngAge")))
  y <- rpois(100, exp(0.3 * x[, 1] + 0.5 * x[, 2]))
  b <- xgboost::xgboost(x, y, objective = "count:poisson", nrounds = 3, nthreads = 1, max_depth = 2)
  list(fit = list(modelList = list(model = list(Fold1 = b, rocs = list()))), x = x)
}

test_that("foldFeatureNames() are the features the fold models were trained with", {
  skip_if_not_installed("xgboost")
  expect_setequal(foldFeatureNames(speciesFit()$fit), c("CMD", "Pinu_con", "youngAge"))
})

test_that("speciesCovCols() leaves out climate, indicators, interactions and bookkeeping columns", {
  expect_identical(speciesCovCols(c("pixelID", "year", "CMD", "CMDsm", "MDC", "Tmax_sm", "lightning", "lightningDays", "cumMDC", "youngAge", "nfLCC_1_2",
                                    "treedWetland", "Pinu_con", "Abie_ama", "dom_agb_Pice_mar", "treedWetland_agb",
                                    "CMD:dom_agb_Pice_mar")),
                   c("Pinu_con", "Abie_ama", "dom_agb_Pice_mar", "treedWetland_agb"))
})

test_that("addAbsentSpeciesCols() adds a missing species column as the value, with a message, and nothing else", {
  dt <- data.table::data.table(pixelID = 1:3, CMD = c(1, 2, 3), youngAge = 0)
  expect_message(out <- addAbsentSpeciesCols(dt, c("CMD", "Pinu_con", "youngAge")), "Pinu_con")
  expect_identical(out$Pinu_con, c(0, 0, 0))
  expect_false("Pinu_con" %in% names(dt))                  # the input is not changed
  out <- addAbsentSpeciesCols(dt, c("CMD", "Pinu_con"), value = 3.6)
  expect_identical(out$Pinu_con, rep(3.6, 3))
  expect_identical(addAbsentSpeciesCols(dt, c("CMD", "youngAge")), dt)       # nothing missing: unchanged
  ## a missing non-species feature is not invented
  out <- addAbsentSpeciesCols(dt, c("CMD", "CMDsp", "Pinu_con"))
  expect_false("CMDsp" %in% names(out))
  expect_true("Pinu_con" %in% names(out))
  expect_identical(absentSpeciesCols(c("CMD", "CMDsp", "Pinu_con"), c("CMD", "youngAge")), "Pinu_con")
})

test_that("predictions with an absent species added as 0 equal those with an explicit 0 column", {
  skip_if_not_installed("xgboost")
  f <- speciesFit()
  cur <- data.table::data.table(CMD = f$x[, "CMD"], youngAge = f$x[, "youngAge"])
  explicit <- data.table::copy(cur)
  data.table::set(explicit, NULL, "Pinu_con", 0)
  added <- suppressMessages(addAbsentSpeciesCols(cur, foldFeatureNames(f$fit)))
  expect_equal(foldMeanPrediction(f$fit, as.matrix(added[, c("CMD", "Pinu_con", "youngAge")])),
               foldMeanPrediction(f$fit, as.matrix(explicit[, c("CMD", "Pinu_con", "youngAge")])))
})

test_that("a renamed or merged species class is not an absent species: unseen species column means no zero-fill", {
  expect_identical(absentSpeciesCols(c("CMD", "Pinu_con"), c("CMD", "youngAge")), "Pinu_con")
  expect_identical(absentSpeciesCols(c("CMD", "Pinu_con"), c("CMD", "Pinu_new")), character())
  dt <- data.table::data.table(pixelID = 1:2, CMD = 1, Pinu_new = 2)
  expect_identical(addAbsentSpeciesCols(dt, c("CMD", "Pinu_con")), dt)
})
