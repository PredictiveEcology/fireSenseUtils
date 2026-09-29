## 2026-09-29: a FireSense spread fit logged "Using a 3 parameter logistic equation / There will be 3
## logit terms & 11 terms in all: logit1, logit2, logit3, CMD, ...". hillSlope is fixed at 1; the
## three non-formula parameters were maxAsymptote, inflectionPoint1 and the per-year random effect
## yearSpreadSD. The names are on `lower`; runDEoptim() must print those.

test_that("runDEoptim() names the parameters from names(lower), grouped", {
  form <- ~ 0 + CMD + youngAge + dom_agb_Pice_mar
  covs <- c("CMD", "youngAge", "dom_agb_Pice_mar")
  lower <- stats::setNames(c(0.25, 0.1, 0, 0, 0, 0),
                           c("maxAsymptote", "inflectionPoint1", covs, "yearSpreadSD"))
  testthat::local_mocked_bindings(
    clusterSetup = function(...) list(itermax = 5, trace = FALSE, strategy = 2L, NP = 40L),
    DEoptimIterative = function(fn, lower, upper, control, ...) list(),
    .package = "clusters")
  withr::local_options(reproducible.useCache = FALSE)
  msgs <- testthat::capture_messages(suppressWarnings(runDEoptim(
    landscape = NULL, annualDTx1000 = NULL, nonAnnualDTx1000 = NULL,
    fireBufferedListDT = NULL, historicalFires = NULL,
    itermax = 5, trace = FALSE, strategy = 2L, cores = c("hostA", "hostB"),
    paths = list(cachePath = withr::local_tempdir()),
    lower = lower, upper = lower + 1, mutuallyExclusive = NULL, formulaToFit = form,
    objFunCoresInternal = 1L, covMinMax = NULL, maxFireSpread = 0.3, Nreps = 1L, thresh = 512,
    logPath = file.path(withr::local_tempdir(), "fit.log"), .verbose = FALSE)))
  msgs <- paste(msgs, collapse = "")
  expect_match(msgs, paste0("Fitting 6 parameters: logistic: maxAsymptote, inflectionPoint1; ",
                            "covariates: CMD, youngAge, dom_agb_Pice_mar; year effect: yearSpreadSD"),
               fixed = TRUE)
  expect_match(msgs, "objectiveFunction threshold SNLL to run all years after first 2 years: 512",
               fixed = TRUE)
  expect_false(grepl("logit", msgs))
})
