## fireSense_dataPrepFit builds the per-ELF species, fuel and climate inputs from the SpreadFit ledger
## before the fit; fireSense_spreadFit rebuilds them after it writes its row, so prediction in the
## same run finds the new ELF. Both use these two functions.

ledgerRow <- function(id, sppEquiv, groups, missing, climVar) {
  pars <- data.frame(maxAsymptote = 0.3, hillSlope1 = 1, inflectionPoint1 = 1)
  for (nm in c(climVar, "youngAge", unique(sppEquiv$FuelClass), names(groups))) pars[[nm]] <- 0.1
  row <- data.frame(numIterations = 10L)
  row[[polygonIDTxt]] <- id
  row$params <- list(pars)
  row$sppEquiv <- list(sppEquiv)
  row$nonForestedLCCGroups <- list(groups)
  row$missingLCCgroup <- missing
  row
}

threeRowLedger <- function() {
  rbind(
    ledgerRow("4.2.2", data.table::data.table(LandR = c("Pice_mar", "Popu_tre"), FuelClass = c("class1", "class2")),
              list(nf_A = c(40, 50)), "nf_A", "CMD"),
    ## a row whose sppEquiv has a species without a LandR name: harmonizing drops it
    ledgerRow("4.2.3", data.table::data.table(LandR = c("Pice_gla", NA), FuelClass = c("class1", "class3")),
              list(nf_B = c(40, 50, 100)), "nf_B", "CMDsm"),
    ledgerRow("4.2.1", data.table::data.table(LandR = "Pinu_ban", FuelClass = "class4"),
              list(nf_C = 50), "nf_C", "PPTsm")
  )
}

test_that("spreadFitLedgerInputs() gives each ledger row its harmonized species, groups and climate variables", {
  ledger <- threeRowLedger()
  out <- suppressWarnings(suppressMessages(
    spreadFitLedgerInputs(ledger, sppEquivCol = "LandR", fuelClassCol = "FuelClass", studyArea = NULL)))
  expect_identical(names(out), c("4.2.2", "4.2.3", "4.2.1"))
  expect_identical(names(out[["4.2.2"]]),
                   c("missingLCCgroup", "nonForestedLCCGroups", "sppColorVect", "sppEquiv", "sppNameVector",
                     "theseClimVars", "theseClimVarsNoUnderscore"))
  expect_identical(out[["4.2.2"]]$sppEquiv$LandR, c("Pice_mar", "Popu_tre"))
  expect_identical(out[["4.2.3"]]$sppEquiv$LandR, "Pice_gla")          # the NA species is dropped
  expect_identical(names(out[["4.2.1"]]$sppColorVect), c("Pinu_ban", "Mixed"))
  expect_identical(out[["4.2.3"]]$nonForestedLCCGroups, list(nf_B = c(40, 50, 100)))
  expect_identical(out[["4.2.1"]]$missingLCCgroup, "nf_C")
  expect_identical(out[["4.2.2"]]$theseClimVars, "CMD")
  expect_identical(out[["4.2.3"]]$theseClimVars, "CMD_sm")
  expect_identical(out[["4.2.3"]]$theseClimVarsNoUnderscore, "CMDsm")
  expect_identical(out[["4.2.1"]]$theseClimVars, "PPT_sm")
})

test_that("spreadFitLedgerInputs() reads an sf ledger, and one without a column in colNames", {
  ledger <- threeRowLedger()
  ledger$nonForestedLCCGroups <- NULL                              # a ledger written before the column
  ledger <- sf::st_sf(ledger, geometry = sf::st_sfc(rep(list(sf::st_point(c(0, 0))), 3)))
  out <- suppressWarnings(suppressMessages(
    spreadFitLedgerInputs(ledger, sppEquivCol = "LandR", fuelClassCol = "FuelClass", studyArea = NULL)))
  expect_length(out, 3L)
  expect_false("nonForestedLCCGroups" %in% names(out[[1]]))
  expect_identical(out[["4.2.3"]]$missingLCCgroup, "nf_B")
  ## colNames limits the list-columns read, as fireSense_dataPrepFit's spreadFitAdditionalColNames does
  out2 <- suppressWarnings(suppressMessages(
    spreadFitLedgerInputs(ledger, sppEquivCol = "LandR", fuelClassCol = "FuelClass", studyArea = NULL,
                          colNames = "sppEquiv")))
  expect_false("missingLCCgroup" %in% names(out2[[1]]))
})

test_that("spreadFitLedgerInputLists() gives the five per-ELF lists, named by polygonID", {
  out <- suppressWarnings(suppressMessages(
    spreadFitLedgerInputs(threeRowLedger(), sppEquivCol = "LandR", fuelClassCol = "FuelClass", studyArea = NULL)))
  lists <- spreadFitLedgerInputLists(out)
  expect_identical(names(lists), c("sppEquivs", "sppColorVects", "sppNameVectors",
                                   "nonForestedLCCGroupsList", "missingLCCgroupList"))
  for (l in lists) expect_identical(names(l), c("4.2.2", "4.2.3", "4.2.1"))
  expect_identical(lists$sppEquivs[["4.2.3"]], out[["4.2.3"]]$sppEquiv)
  expect_identical(lists$missingLCCgroupList, list(`4.2.2` = "nf_A", `4.2.3` = "nf_B", `4.2.1` = "nf_C"))
})
