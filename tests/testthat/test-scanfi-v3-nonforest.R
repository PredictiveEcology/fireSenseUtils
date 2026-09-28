## LandR::prepInputs_SCANFI_LCC_FAO(dataVersion = "V3") (SCANFI v3, annual 1985-2025, not merged
## in LandR yet) adds one Canada-LCC-space code for "burn scar" that has no counterpart in the
## NTEMS/SCANFI-v2 code table `fireSenseNonflammableLCC` documents. The value is not settled
## upstream ("60 is likely"), so every reference to it here goes through this one constant.
burnScarLCC <- 60L

test_that("makeFireSenseLCC(scanfiVersion = 'V3') stops clearly when LandR lacks SCANFI V3 support", {
  ## R/makeFireSenseLCC.R: LandR::prepInputs_SCANFI_LCC_FAO()'s own dataVersion checks (as of
  ## LandR 1.2.0.9035) branch on "V1", and fall through to V2 behaviour for anything else --
  ## including "V3" -- without erroring. So this cannot be tested by asking LandR to fail; the
  ## guard reads the installed LandR's own source instead.
  installedSupportsV3 <- any(grepl("\"V3\"", deparse(body(LandR::prepInputs_SCANFI_LCC_FAO)),
                                    fixed = TRUE))
  skip_if(installedSupportsV3, "installed LandR already supports SCANFI V3")

  ## the check runs before any download, so this is a fast, offline call
  expect_error(
    makeFireSenseLCC(neededYear = 2020, to = terra::rast(nrows = 2, ncols = 2, vals = 1),
                     destinationPath = tempdir(), scanfiVersion = "V3"),
    "SCANFI V3"
  )

  ## other dataVersions are untouched by the guard
  expect_true(.checkScanfiVersionSupported("V2"))
  expect_true(.checkScanfiVersionSupported("V1"))
})

test_that("a new flammable non-forest code (burn scar) is kept and grouped, like 40/50/100", {
  withr::local_package("terra")
  withr::local_package("data.table")
  set.seed(1)

  ## fuelClassPrep()/assessFuelClasses(): a landscape of non-forest pixels only (B is NA
  ## everywhere), three land covers that burn at different rates so k-means on the glm
  ## coefficients has something to cluster -- burnScarLCC is one of them, exactly like the
  ## existing 40/50/80 test, just with an LCC code fireSenseUtils has never seen before.
  n <- 300L
  landscape <- data.table(
    cell = seq_len(n), speciesCode = NA_character_,
    lcc = rep(c(50L, burnScarLCC, 100L), each = n / 3L),
    B = NA_integer_, totalBiomass = NA_integer_, year = 2020L
  )
  burnRates <- c(0.3, 0.55, 0.05)
  names(burnRates) <- as.character(c(50L, burnScarLCC, 100L))
  landscape[, burned := rbinom(.N, 1, burnRates[as.character(lcc)])]
  noSpp <- data.table(LandR = character(0), FuelClass = character(0))

  out <- assessFuelClasses(landscape = landscape, fuelCol = "FuelClass", sppEquiv = noSpp,
                           sppEquivCol = "LandR", nonforestLCC = c(50L, burnScarLCC, 100L))

  ## kept: burnScarLCC is not dropped, not NA'd, and not silently merged into "missingForest"
  ## (there is nothing to lump it with here -- nonforestLCC covers every code present)
  expect_true(burnScarLCC %in% as.numeric(unlist(out$nonForestedLCCGroups)))
  burnScarGroup <- names(out$nonForestedLCCGroups)[
    vapply(out$nonForestedLCCGroups, function(x) burnScarLCC %in% x, logical(1))
  ]
  expect_length(burnScarGroup, 1L)

  ## makeLandcoverDT(): a synthetic land-cover raster carrying burnScarLCC alongside 50, 100 and
  ## a forest code. The burn-scar pixels must land in their k-means group's column, not vanish as
  ## an unrecognized `lcc_<code>` column (which is what a hard-coded nonForestedLCCGroups/
  ## forestedLCC list would do).
  rstLCC <- terra::rast(nrows = 10, ncols = 10,
                        vals = rep(c(50L, burnScarLCC, 100L, 210L), 25))
  flammableRTM <- terra::rast(rstLCC, vals = 1)
  landcoverDT <- makeLandcoverDT(rstLCC = rstLCC, flammableRTM = flammableRTM,
                                 forestedLCC = c(210, 220, 230, 240),
                                 nonForestedLCCGroups = out$nonForestedLCCGroups)

  expect_true(sum(landcoverDT[[burnScarGroup]]) > 0)
})
