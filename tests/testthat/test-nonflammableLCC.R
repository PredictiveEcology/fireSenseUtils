## Root cause: fireSense_dataPrepFit.R:81 and fireSense_dataPrepPredict.R:46 hard-code
## `nonflammableLCC = c(0, 20, 31, 32, 33)`, the NTEMS non-flammable codes, and never
## included SCANFI's combined rock/exposed code. `makeFireSenseLCC()`'s SCANFI path
## (`LandR::prepInputs_SCANFI_LCC_FAO()`) produces that code, so rock entered ELF fits
## as flammable non-forest (20% of ELF 14.4's "flammable" pixels, 2026-09 land cover).
## `fireSenseNonflammableLCC` is the single source of truth both the modules and
## `makeFireSenseLCC()`'s own dominant-flammable-class step take their default from.

test_that("fireSenseNonflammableLCC includes LandR's SCANFI rock/exposed code", {
  ## LandR::convert_SCANFI_LCC_codes() recodes SCANFI's raw 1-8 landcover classes with
  ## `oldVals <- 1:8; newVals <- c(...)`, in class order "Bryoids, herbs, rock/exposed,
  ## shrubs, broadleaf, conifer, mixedwood, water". Read those two vectors out of the
  ## installed LandR function itself, not a copy, so a LandR change is caught here too.
  exprs <- as.list(body(LandR::convert_SCANFI_LCC_codes))
  getAssigned <- function(varName) {
    for (e in exprs) {
      if (is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name(varName))) {
        return(eval(e[[3]]))
      }
    }
    stop("could not find assignment to ", varName, " in LandR::convert_SCANFI_LCC_codes")
  }
  oldVals <- getAssigned("oldVals")
  newVals <- getAssigned("newVals")
  labels <- c("Bryoids", "herbs", "rock/exposed", "shrubs", "broadleaf", "conifer", "mixedwood", "water")
  stopifnot(identical(oldVals, seq_along(labels)))

  rockCode  <- newVals[match("rock/exposed", labels)]
  waterCode <- newVals[match("water", labels)]

  expect_identical(rockCode, 30)
  expect_true(rockCode %in% fireSenseNonflammableLCC)
  expect_true(waterCode %in% fireSenseNonflammableLCC)
})

test_that("fireSenseNonflammableLCC also has the NTEMS non-flammable codes", {
  ## 0 = no data, 20 = water, 31 = snow/ice, 32 = rock/rubble, 33 = exposed/barren land
  ## (LandR::prepInputs_NTEMS_LCC_FAO(), LandR::prepInputs_SCANFI_LCC_FAO()).
  expect_true(all(c(0, 20, 31, 32, 33) %in% fireSenseNonflammableLCC))
})

test_that("ELFflammableArea() and makeFireSenseLCC() default to fireSenseNonflammableLCC", {
  expect_identical(eval(formals(ELFflammableArea)$nonflammableLCC), fireSenseNonflammableLCC)
  expect_identical(eval(formals(makeFireSenseLCC)$nonflammableLCC), fireSenseNonflammableLCC)
})
