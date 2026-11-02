## The fire models' climate variables are set by fireSense_dataPrepFit unless a user supplies them, and canClimateData's
## `climateVariables` follows from them (Eliot, 2026-09-23: "I want the user to NOT need to select").
## Ignition (xgboost) takes several correlated variables; spread takes one.

test_that("the default gives ignition several variables and lets spread choose among them", {
  expect_identical(defaultClimateVariablesForFire$ignition, c("CMD", "cumMDC", "CMD_sm", "CMD_sp"))
  expect_identical(defaultClimateVariablesForFire$spread, "auto")
})

test_that("climate variable names are accepted with or without underscores", {
  expect_identical(climateNAnames(c("CMD_sm", "CMDsm", "cumMDC", "CMDsp")), c("CMD_sm", "CMD_sm", "cumMDC", "CMD_sp"))
  expect_error(climateNAnames(c("CMDsm", "notAVariable")), "notAVariable")
})

test_that("the climate layers cover every variable once, for the fire years", {
  cl <- fireClimateLayers(defaultClimateVariablesForFire, historicalYears = 1985:2024, projectedYears = NULL)
  expect_setequal(names(cl), c("historical_CMD", "historical_CMDsm", "historical_cumMDC", "historical_CMDsp"))
  expect_identical(cl$historical_CMDsm$.dots$historical_years, 1985:2024)
  expect_identical(cl$historical_CMDsm$fun, quote(calcAsIs))
  ## cumMDC is derived from monthly variables and starts 5 years early (climateData::climateLayers)
  expect_identical(cl$historical_cumMDC$fun, quote(calcCumMDC))
  expect_identical(cl$historical_cumMDC$.dots$historical_years, 1980:2024)
  ## no-underscore names from an older project give the same layers
  cl2 <- fireClimateLayers(list(ignition = c("CMD", "cumMDC", "CMDsm", "CMDsp"), spread = "auto"), 1985:2024)
  expect_identical(cl2, cl)
})

## Projected layers follow canClimateData's `projectedClimateYears` and nothing else: NRV runs set
## that empty (and a real climateGCM, which canClimateData validates), so they request no projections.
test_that("no projected years requests no projected layers; some request exactly those years", {
  for (none in list(NULL, integer(0))) {
    cl <- fireClimateLayers(defaultClimateVariablesForFire, historicalYears = 1985:2024, projectedYears = none)
    expect_false(any(grepl("^projected_", names(cl))))
  }
  yrs <- 2025:2044
  cl <- fireClimateLayers(defaultClimateVariablesForFire, 1985:2024, projectedYears = yrs)
  proj <- cl[grepl("^projected_", names(cl))]
  expect_setequal(names(proj), c("projected_CMD", "projected_CMDsm", "projected_cumMDC", "projected_CMDsp"))
  expect_identical(proj$projected_CMDsm$.dots$future_years, yrs)
  ## the fit variables added later follow the same rule
  have <- fireClimateLayers(list(ignition = "CMD_sm"), 1985:2024)
  expect_false(any(grepl("^projected_", names(addFitClimateVariables(have, "cumMDC", 1985:2024)))))
  expect_true(any(grepl("^projected_", names(addFitClimateVariables(have, "cumMDC", 1985:2024, yrs)))))
})

test_that("canClimateProjectedYears reads canClimateData's parameter, and is empty without it", {
  mp <- withr::local_tempdir()
  SpaDES.core::newModule("canClimateData", path = mp, open = FALSE)
  yearsIn <- function(yrs) {
    sim <- SpaDES.core::simInit(modules = "canClimateData", paths = list(modulePath = mp),
                                params = list(canClimateData = list(projectedClimateYears = yrs)))
    canClimateProjectedYears(sim)
  }
  expect_identical(yearsIn(2030:2040), 2030:2040)
  expect_length(yearsIn(integer(0)), 0)
  expect_length(canClimateProjectedYears(SpaDES.core::simInit(paths = list(modulePath = mp))), 0)
})

test_that("prediction prepares every fit's variables, keeping those already defined", {
  have <- fireClimateLayers(list(ignition = "CMD_sm"), historicalYears = 1985:2024)
  ## two ELFs whose fits chose CMD_sm and cumMDC
  out <- addFitClimateVariables(have, c("CMD_sm", "cumMDC"), historicalYears = 1990:2000, projectedYears = NULL)
  expect_setequal(names(out), c("historical_CMDsm", "historical_cumMDC"))
  expect_identical(out$historical_CMDsm, have$historical_CMDsm)            # kept, years untouched
  expect_identical(out$historical_cumMDC$fun, quote(calcCumMDC))
  ## nothing to add
  expect_identical(addFitClimateVariables(have, "CMDsm", 1990:2000), have)
})

test_that("spreadFitClimateVariables() gives every ledger row's climate variables, once", {
  ## as spreadFitLedgerInputs() returns them: two ELFs whose fits chose CMD and CMD_sm, a third CMD
  ins <- list(list(theseClimVars = "CMD"), list(theseClimVars = "CMD_sm"), list(theseClimVars = "CMD"))
  expect_identical(spreadFitClimateVariables(ins), c("CMD", "CMD_sm"))
  expect_null(spreadFitClimateVariables(list()))
})
