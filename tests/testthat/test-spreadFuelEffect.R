## spreadFuelEffect(): the change in spread probability from 0 to a fuel's median non-zero value, everything else at its centre.
## The rows are small synthetic ones, as latestSpreadFits() returns them.

lowerP <- 0.13
fuelRow <- function(pars, medianNonZero = list(dom_agb_Pice_mar = 0.3), centre = list(dom_agb_Pice_mar = 0.1), ...) {
  params <- data.frame(maxAsymptote = 0.25, hillSlope1 = 1, inflectionPoint1 = 1, yearSpreadSD = 1, pars, check.names = FALSE)
  row <- data.frame(polygonID = "4.1", ...)
  row$params <- list(params)
  if (!is.null(medianNonZero)) row[[spreadFitCovMedianNonZeroTxt]] <- list(medianNonZero)
  if (!is.null(centre)) row[[spreadFitCovCentreTxt]] <- list(centre)
  row
}
## the default logistic link: p = f + (pc - f) * plogis(x)  (logistic3p with hillSlope1 = inflectionPoint1 = 1)
pLogis <- function(x, pc = 0.25) lowerP + (pc - lowerP) * plogis(x)

test_that("spreadFuelEffect() is the change in p from 0 to the median non-zero value, the others at their centres", {
  row <- fuelRow(list(`(Intercept)` = 1, dom_agb_Pice_mar = 2))
  ## x = intercept + beta * (value - centre): fuel absent -> 1 + 2 * (0 - 0.1); at its median -> 1 + 2 * (0.3 - 0.1)
  expect_equal(unname(spreadFuelEffect(row)), pLogis(1 + 2 * 0.2) - pLogis(1 + 2 * -0.1))
  expect_named(spreadFuelEffect(row), "dom_agb_Pice_mar")
})

test_that("spreadFuelEffect() averages over the parameter sets and has one value per fuel term, other covariates at their centre", {
  params <- data.frame(maxAsymptote = 0.25, hillSlope1 = 1, inflectionPoint1 = 1, `(Intercept)` = c(0.5, 1.5),
                       dom_agb_Pice_mar = c(1, 3), sec_agb_Pinu_con = 2, CMD = 7, check.names = FALSE)
  row <- data.frame(polygonID = "4.1"); row$params <- list(params)
  row[[spreadFitCovMedianNonZeroTxt]] <- list(list(dom_agb_Pice_mar = 0.3, sec_agb_Pinu_con = 0.05))
  row[[spreadFitCovCentreTxt]] <- list(list(dom_agb_Pice_mar = 0.1, sec_agb_Pinu_con = 0.02, CMD = 5))
  out <- spreadFuelEffect(row)
  expect_named(out, c("dom_agb_Pice_mar", "sec_agb_Pinu_con"))
  ## CMD sits at its centre (0 centred), so its coefficient does not matter
  dom <- mean(pLogis(c(0.5, 1.5) + c(1, 3) * 0.2) - pLogis(c(0.5, 1.5) + c(1, 3) * -0.1))
  sec <- mean(pLogis(c(0.5, 1.5) + 2 * 0.03) - pLogis(c(0.5, 1.5) + 2 * -0.02))
  expect_equal(unname(out), c(dom, sec))
})

test_that("a fit at its ceiling has an effect near 0, whatever the coefficient", {
  big <- spreadFuelEffect(fuelRow(list(`(Intercept)` = 30, dom_agb_Pice_mar = 80)))
  expect_lt(abs(big), 1e-6)
  ## the same coefficient where the intercept leaves room moves p
  expect_gt(spreadFuelEffect(fuelRow(list(`(Intercept)` = 0, dom_agb_Pice_mar = 80))), 0.05)
})

test_that("spreadFuelEffect() applies the row's link: the power law is evaluated by the one link function", {
  row <- fuelRow(list(`(Intercept)` = 1, dom_agb_Pice_mar = 2), spreadLink = "powerLaw", powerLawGamma_spread = 1.9,
                 powerLawShape_spread = 1)
  pl <- function(x) 0.25 - (0.25 - lowerP) * exp(-sign(x) * abs(x / 1.9))
  expect_equal(unname(spreadFuelEffect(row)), pl(1 + 2 * 0.2) - max(pl(1 + 2 * -0.1), 0))
  ## a power law is not the logistic
  expect_false(isTRUE(all.equal(unname(spreadFuelEffect(row)), pLogis(1.4) - pLogis(0.8))))
})

test_that("a row without the fit-data summary, or a fuel never non-zero, has NA", {
  row <- fuelRow(list(`(Intercept)` = 1, dom_agb_Pice_mar = 2), medianNonZero = NULL)
  expect_identical(spreadFuelEffect(row), c(dom_agb_Pice_mar = NA_real_))
  expect_identical(spreadFuelEffect(fuelRow(list(`(Intercept)` = 1, dom_agb_Pice_mar = 2), medianNonZero = list(dom_agb_Pice_mar = NA))),
                   c(dom_agb_Pice_mar = NA_real_))
  expect_length(spreadFuelEffect(fuelRow(list(`(Intercept)` = 1, youngAge = -2))), 0L)         # no fuel term
})

test_that("a fit without an intercept has no centres: the covariate is at 0 and at its median", {
  row <- fuelRow(list(dom_agb_Pice_mar = 2), centre = NULL)
  expect_equal(unname(spreadFuelEffect(row)), pLogis(2 * 0.3) - pLogis(0))
})

test_that("spreadCovStats() gives each fuel's median non-zero value and sd on the scale of the centre", {
  ## biomass 0, 2000, 4000, 6000 over 1e4 (covMinMax c(0, 1e4)): 0 and 0.2, 0.4, 0.6; and a second year
  annual <- list(`2001` = data.table::data.table(pixelID = 1:4, CMD = c(0L, 0L, 0L, 0L)),
                 `2002` = data.table::data.table(pixelID = 1:4, CMD = c(0L, 0L, 0L, 0L)))
  nonAnnual <- list(`2000` = data.table::data.table(pixelID = 1:4, dom_agb_Pice_mar = c(0L, 2000L, 4000L, 6000L) * 1000L))  # x1000 integers, as the fit has them
  s <- spreadCovStats(annual, nonAnnual, "~ 1 + dom_agb_Pice_mar + CMD", covMinMax = list(dom_agb_Pice_mar = fuelLinearRange) |> data.table::as.data.table(),
                         mutuallyExclusive = NULL)
  expect_equal(s$medianNonZero$dom_agb_Pice_mar, 0.4)
  expect_equal(s$sd$dom_agb_Pice_mar, sd(rep(c(0, 0.2, 0.4, 0.6), 2)))
  expect_named(s$medianNonZero, "dom_agb_Pice_mar")                          # CMD is not a fuel
  expect_null(spreadCovStats(annual, nonAnnual, "~ 1 + CMD")$medianNonZero)
  ## the centre comes from the same pass, and is spreadCovCentre()'s
  cmm <- data.table::as.data.table(list(dom_agb_Pice_mar = fuelLinearRange))
  expect_equal(spreadCovStats(annual, nonAnnual, "~ 1 + dom_agb_Pice_mar", cmm, NULL)$centre,
               spreadCovCentre(annual, nonAnnual, "~ 1 + dom_agb_Pice_mar", cmm, NULL))
  expect_equal(spreadCovStats(annual, nonAnnual, "~ 1 + dom_agb_Pice_mar", cmm, NULL)$centre$dom_agb_Pice_mar, mean(c(0, 0.2, 0.4, 0.6)))
  expect_null(spreadCovStats(annual, nonAnnual, "~ 0 + dom_agb_Pice_mar", cmm, NULL)$centre)   # no intercept: no centre, still a summary
})

test_that("elfFuelEffects() and elfFeatures() give each member species of a grouped fuel its effect, in the 'Spread fuel effect' group", {
  sppEquiv <- data.frame(LandR = c("Popu_tre", "Betu_pap", "Pice_mar"), FuelClass = c("Pp_tr.Bt_pa", "Pp_tr.Bt_pa", "Pice_mar"))
  row <- fuelRow(list(`(Intercept)` = 1, `dom_agb_Pp_tr.Bt_pa` = 2, dom_agb_Pice_mar = 3),
                 medianNonZero = list(`dom_agb_Pp_tr.Bt_pa` = 0.3, dom_agb_Pice_mar = 0.3),
                 centre = list(`dom_agb_Pp_tr.Bt_pa` = 0.1, dom_agb_Pice_mar = 0.1), spreadLink = "ceiling", spreadClimate = "predictor",
                 fuelScale = "linear")
  row$sppEquiv <- list(sppEquiv)
  old <- fuelRow(list(`(Intercept)` = 1, dom_agb_Pice_mar = 3), medianNonZero = NULL, spreadLink = "ceiling", spreadClimate = "predictor",
                 fuelScale = "linear")
  old$polygonID <- "4.3"; old$sppEquiv <- list(sppEquiv)
  ## the ledger reader fills the columns a row lacks with NA
  for (col in setdiff(names(row), names(old))) old[[col]] <- NA
  rows <- rbind(row, old)
  eff <- elfFuelEffects(rows)
  e <- function(elf, co) eff$value[eff$ELF == elf & eff$coefficient == co]
  expect_equal(e("4.1", "dom_agb_Popu_tre"), e("4.1", "dom_agb_Betu_pap"))
  expect_equal(e("4.1", "dom_agb_Popu_tre"), pLogis(1 + 2 * 0.2) - pLogis(1 + 2 * -0.1))
  expect_identical(eff$note[eff$ELF == "4.1" & eff$coefficient == "dom_agb_Popu_tre"], "shared with Betu_pap")
  expect_true(is.na(e("4.3", "dom_agb_Pice_mar")))
  expect_identical(eff$note[eff$ELF == "4.3"], "no fit-data summary")
})

test_that("elfFeatures() emits spreadEffect_<fuel> in the group 'Spread fuel effect', with a label and units, and none where NA", {
  local_mocked_bindings(
    ELFfireRates = function(...) data.table::data.table(ELF = "5.1", naturalFiresPerYear = 1, naturalFiresPer100km2PerYear = 1,
                                                        haBurnedPerYear = 1, proportionBurnedPerYear = 0.01),
    latestSpreadFits = function(...) {
      row <- fuelRow(list(`(Intercept)` = 1, dom_agb_Pice_mar = 2), spreadLink = "ceiling", spreadClimate = "predictor", fuelScale = "linear")
      row$polygonID <- "5.1"
      old <- row; old$polygonID <- "3.1.1"; old[[spreadFitCovMedianNonZeroTxt]] <- NA
      rbind(row, old)
    })
  feat <- elfFeatures(toyELFs(), NULL, fireYears = 2001:2003, pixelAreaHa = 2500, destinationPath = tempdir(),
                      dates = data.frame(ELF = "5.1", spreadFitDate = as.POSIXct(NA), ignitionFitDate = as.POSIXct(NA)))
  eff <- feat[feat$group == "Spread fuel effect", ]
  expect_identical(eff$ELF, "5.1")                                        # the ELF without a summary has no row
  expect_identical(eff$feature, "spreadEffect_dom_agb_Pice_mar")
  expect_identical(eff$label, "Dominant AGB: Black spruce (Pice_mar) \u2014 change in spread probability")
  expect_identical(eff$units, "probability")
  expect_equal(eff$value, pLogis(1 + 2 * 0.2) - pLogis(1 + 2 * -0.1))
})

test_that("youngAge is 0 in both scenarios (a fuelled pixel is not young), the other covariates at their centres", {
  ## x = intercept + beta_fuel * (value - centre) + beta_young * (0 - centreYoung): the young term is the same in both scenarios
  params <- data.frame(maxAsymptote = 0.25, hillSlope1 = 1, inflectionPoint1 = 1, `(Intercept)` = 1, dom_agb_Pice_mar = 2,
                       youngAge = -3, check.names = FALSE)
  row <- data.frame(polygonID = "4.1"); row$params <- list(params)
  row[[spreadFitCovMedianNonZeroTxt]] <- list(list(dom_agb_Pice_mar = 0.3))
  row[[spreadFitCovCentreTxt]] <- list(list(dom_agb_Pice_mar = 0.1, youngAge = 0.2))
  young <- -3 * (0 - 0.2)
  expect_equal(unname(spreadFuelEffect(row)), pLogis(1 + young + 2 * 0.2) - pLogis(1 + young + 2 * -0.1))
  ## not the young term at its centre (0): that would be another value
  expect_false(isTRUE(all.equal(unname(spreadFuelEffect(row)), pLogis(1 + 2 * 0.2) - pLogis(1 + 2 * -0.1))))
})

test_that("a fuel x climate interaction is built as the fit builds it, from its centred covariates", {
  params <- data.frame(maxAsymptote = 0.25, hillSlope1 = 1, inflectionPoint1 = 1, `(Intercept)` = 1, CMD = 1,
                       dom_agb_Pice_mar = 2, `CMD:dom_agb_Pice_mar` = 5, check.names = FALSE)
  row <- data.frame(polygonID = "4.1"); row$params <- list(params)
  row[[spreadFitCovMedianNonZeroTxt]] <- list(list(dom_agb_Pice_mar = 0.3))
  row[[spreadFitCovCentreTxt]] <- list(list(dom_agb_Pice_mar = 0.1, CMD = 0.5))
  ## CMD at its centre (0), so the interaction is 0 in both scenarios
  expect_equal(unname(spreadFuelEffect(row)), pLogis(1 + 2 * 0.2) - pLogis(1 + 2 * -0.1))
})

test_that("with climate in the ceiling (climateRef) the climate is at its reference, so it adds nothing to the fuel's effect", {
  params <- data.frame(maxAsymptote = 0.2, hillSlope1 = 1, inflectionPoint1 = 1, `(Intercept)` = 1, CMD = 5, dom_agb_Pice_mar = 2,
                       check.names = FALSE)
  mk <- function(mode) {
    row <- data.frame(polygonID = "4.1", spreadLink = "powerLawSmooth", spreadClimate = mode, lowerAsymptote_spread = 0.1,
                      powerLawGamma_spread = 1.9, powerLawShape_spread = 1, powerLawPc50_spread = 0.25)
    row$params <- list(params)
    row$climateRef_spread <- list(list(CMD = 0.4))
    row[[spreadFitCovMedianNonZeroTxt]] <- list(list(dom_agb_Pice_mar = 0.3))
    row[[spreadFitCovCentreTxt]] <- list(list(dom_agb_Pice_mar = 0.1, CMD = 0.4))
    row
  }
  ceiling <- spreadFuelEffect(mk("ceiling"))
  expect_false(is.na(ceiling))
  expect_gt(ceiling, 0)
  ## the same params, climate under predictor (CMD at its centre, 0): the same x for the fuel, so the same effect
  expect_equal(unname(ceiling), unname(spreadFuelEffect(mk("predictor"))))
})

test_that("elfFuelEffects() gives NA, excluded as 'different model', to an ELF whose model differs; spreadCovStats() has nothing to find without intercept or fuel", {
  mk <- function(id, link) {
    row <- fuelRow(list(`(Intercept)` = 1, dom_agb_Pice_mar = 2), spreadLink = link, spreadClimate = "predictor", fuelScale = "linear")
    row$polygonID <- id
    row
  }
  eff <- elfFuelEffects(rbind(mk("a", "ceiling"), mk("b", "ceiling"), mk("c", "powerLaw")))
  expect_identical(eff$excluded[eff$ELF == "c"], "different model")
  expect_true(is.na(eff$value[eff$ELF == "c"]))
  expect_false(is.na(eff$value[eff$ELF == "a"]))
  annual <- list(`2001` = data.table::data.table(pixelID = 1:2, CMD = c(0L, 0L)))
  nonAnnual <- list(`2000` = data.table::data.table(pixelID = 1:2, nf = c(0L, 0L)))
  expect_identical(spreadCovStats(annual, nonAnnual, "~ 0 + CMD"), list(centre = NULL, medianNonZero = NULL, sd = NULL))
})
