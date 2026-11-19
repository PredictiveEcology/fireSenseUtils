## Climate in the power-law ceiling (spreadLinkInputs(), spreadProbFromCovs(climateRef =, powerLawPc50 =)): climate
## moves pc toward the landscape's runaway threshold p50, never past it:
## pc_i = max(p50 - (p50 - pc0) * exp(-shift_i), floor), shift_i = sum(beta_clim * (clim - climateRef)).

p50 <- 0.3
pc0 <- 0.27
b <- c(clim = 1.5, fuel = 3, youngAge = -1)
ref <- list(clim = 0.3)

test_that("without climateRef the power-law links give what they gave, bit for bit", {
  mat <- matOf("2004")
  for (link in spreadPowerLawLinks) {
    old <- logisticAll(lp3, mat, b, f, link = link, gamma = 1.9, lowerAsymptote = 0.05)
    expect_identical(spreadProbFromCovs(mat, lp3, b, f, link = link, powerLawGamma = 1.9, lowerAsymptote = 0.05), old)
    expect_identical(spreadProbFromCovs(mat, lp3, b, f, link = link, powerLawGamma = 1.9, lowerAsymptote = 0.05,
                                        powerLawPc50 = p50), old)       # p50 is used only with climateRef
    inputs <- spreadLinkInputs(mat, lp3, b, f, link, powerLawGamma = 1.9, lowerAsymptote = 0.05)
    expect_identical(linearLinkProb(inputs$x, link, inputs$pars), old)
    expect_identical(inputs$pars, linearLinkPars(link, lp3[1], 1.9, f, spreadPowerLawShape, 0.05))
  }
})

test_that("in ceiling mode p never passes p50, and pc is pc0 at the reference climate", {
  for (link in spreadPowerLawLinks) {
    mat <- matOf("2005")
    for (shiftBy in c(-5, -1, 0, 1, 5, 50)) {
      m <- mat; m[, "clim"] <- 0.3 + shiftBy
      p <- spreadProbFromCovs(m, lp3, b, f, link = link, climateRef = ref, powerLawPc50 = p50)
      expect_lte(max(p), p50)
    }
    m <- mat; m[, "clim"] <- 0.3
    inputs <- spreadLinkInputs(m, lp3, b, f, link, ref, p50)
    expect_equal(unname(inputs$pars[, "pc"]), rep(pc0, nrow(m)))
    ## and p there is the fuel-only p with pc = pc0 (climateAtRef gives it for any climate)
    p0 <- spreadProbFromCovs(mat, lp3, b, f, link = link, climateRef = ref, powerLawPc50 = p50, climateAtRef = TRUE)
    expect_equal(as.numeric(spreadProbFromCovs(m, lp3, b, f, link = link, climateRef = ref, powerLawPc50 = p50)),
                 as.numeric(p0))
    expect_equal(as.numeric(p0), as.numeric(spreadProbFromCovs(mat[, c("fuel", "youngAge")], lp3, b[-1], f, link = link)))
  }
})

test_that("log(p50 - pc) is linear in climate with slope -beta, and pc is floored for a very wet year", {
  mat <- matOf("2005")
  clim <- seq(-0.2, 0.9, length.out = nrow(mat))
  mat[, "clim"] <- clim
  pc <- spreadLinkInputs(mat, lp3, b, f, "powerLaw", ref, p50)$pars[, "pc"]
  ok <- pc > f
  expect_equal(unname(coef(lm(log(p50 - pc[ok]) ~ clim[ok]))[2]), -b[["clim"]])
  expect_equal(unname(coef(lm(log(p50 - pc[ok]) ~ clim[ok]))[1]), log(p50 - pc0) + b[["clim"]] * 0.3)
  wet <- mat; wet[, "clim"] <- -20
  expect_true(all(spreadLinkInputs(wet, lp3, b, f, "powerLaw", ref, p50)$pars[, "pc"] == f))
  ## the fuel part of x has no climate in it
  expect_equal(spreadLinkInputs(mat, lp3, b, f, "powerLaw", ref, p50)$x, unname(mat[, "fuel"] * 3 - mat[, "youngAge"]))
})

test_that("spreadProbFromCovs() is linearLinkProb() of spreadLinkInputs(), with several climate columns and an intercept", {
  mat <- cbind(`(Intercept)` = 1, matOf("2004", c("fuel", "clim", "clim2", "youngAge")))
  bb <- c(`(Intercept)` = -0.5, fuel = 3, clim = 2, clim2 = 1, youngAge = -1)
  rr <- c(clim2 = 0.5, clim = 0.4)
  for (link in spreadPowerLawLinks) {
    inputs <- spreadLinkInputs(mat, lp3, bb, f, link, rr, p50, lowerAsymptote = 0.05)
    expect_identical(dim(inputs$pars), c(nrow(mat), length(spreadLinkParNames[[link]])))
    expect_identical(colnames(inputs$pars), spreadLinkParNames[[link]])
    expect_identical(spreadProbFromCovs(mat, lp3, bb, f, link = link, climateRef = rr, powerLawPc50 = p50, lowerAsymptote = 0.05),
                     linearLinkProb(inputs$x, link, inputs$pars))
    shift <- 2 * (mat[, "clim"] - 0.4) + 1 * (mat[, "clim2"] - 0.5)
    expect_equal(unname(inputs$pars[, "pc"]), pmax(p50 - (p50 - pc0) * exp(-shift), f))
  }
})

test_that("ceiling mode stops without a link it can use or without p50", {
  mat <- matOf("2004")
  expect_error(spreadProbFromCovs(mat, lp3, b, f, link = "richards", climateRef = ref, powerLawPc50 = p50), "power-law link")
  expect_error(spreadLinkInputs(mat, lp3, b, f, "richards", ref, p50), "power-law link")
  expect_error(spreadProbFromCovs(mat, lp3, b, f, link = "powerLaw", climateRef = ref), "powerLawPc50")
  expect_error(spreadLinkInputs(mat, lp3, b, f, "powerLaw", ref, NA_real_), "powerLawPc50")
  expect_error(spreadLinkInputs(mat, lp3, b, f, "powerLaw", c(CMD = 0.3), p50), "must be a column of `mat`")
  expect_error(spreadLinkInputs(mat, lp3, b, f, "logistic3pUpper"), "must be one of")
})

test_that("powerLawPcBounds(climateCeiling = TRUE) has p50 as its upper bound", {
  expect_equal(powerLawPcBounds(0.3), c(spreadPowerLawPcBounds[1], 0.305))
  expect_equal(powerLawPcBounds(0.3, climateCeiling = TRUE), c(spreadPowerLawPcBounds[1], 0.3))
})

test_that("the objective gives the same p on the whole landscape and in its buffers, and climate raises it", {
  skip_if_not_installed("SpaDES.tools")
  args <- list(clim = function(pix) 0.5 + 0.4 * cos(pix), link = "powerLaw", powerLawGamma = 1.9,
               climateRef = c(clim = 0.5), powerLawPc50 = 0.5, returnSims = TRUE, returnBurned = TRUE, seed = 3)
  for (more in list(list(), list(yearSpreadSD = 0.5))) {
    full <- do.call(wlObjective, c(list(fx = closedBuffer(), annPix = seq_len(3600L)), args, more))
    lazy <- do.call(wlObjective, c(list(whole = TRUE), args, more))
    expect_identical(lazy$burned, full$burned)
    expect_identical(lazy$sims, full$sims)
  }
  ## climate above its reference gives bigger fires than climate at it
  hot <- do.call(wlObjective, utils::modifyList(c(list(whole = TRUE), args), list(clim = function(pix) rep(1.2, length(pix)))))
  ref0 <- do.call(wlObjective, utils::modifyList(c(list(whole = TRUE), args), list(clim = function(pix) rep(0.5, length(pix)))))
  expect_gt(mean(hot$sims$sim), mean(ref0$sims$sim))
})

test_that("rows with no climate cannot burn, and the parameters of such cells come from one helper", {
  mat <- matOf("2005")
  mat[c(3, 7), "clim"] <- NA
  for (link in spreadPowerLawLinks) {
    inputs <- spreadLinkInputs(mat, lp3, b, f, link, ref, p50)
    expect_false(anyNA(inputs$pars))
    expect_true(all(is.na(inputs$x[c(3, 7)])))
    expect_identical(inputs$pars[c(3, 7), , drop = FALSE], unburnableLinkPars(link, 2L))
    expect_false(anyNA(inputs$x[-c(3, 7)]))
    expect_true(all(is.na(spreadProbFromCovs(mat, lp3, b, f, link = link, climateRef = ref, powerLawPc50 = p50)[c(3, 7)])))
  }
  u <- unburnableLinkPars("powerLawSmooth", 3L)
  expect_identical(dim(u), c(3L, length(spreadLinkParNames$powerLawSmooth)))
  expect_identical(colnames(u), spreadLinkParNames$powerLawSmooth)
  expect_true(all(u[, c("pc", "floor", "lower")] == 0))
  expect_equal(linearLinkProb(c(1, -1, 5), "powerLawSmooth", u), c(0, 0, 0))
})

test_that("pc0 above p50 stops in ceiling mode", {
  expect_error(spreadLinkInputs(matOf("2004"), c(0.35, 1, 1), b, f, "powerLaw", ref, p50), "powerLawPcBounds\\(climateCeiling = TRUE\\)")
})

test_that("the response plot's dashed ceiling is p50 in ceiling mode, with climate curves that stay below it", {
  par <- c(0.2, clim = 1.5, fuel = 3, youngAge = -1)
  r <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20, vals = 1)
  fb <- lapply(climFix$historicalFires, function(h) dtC(pixelID = 1:nPixC, buffer = as.integer(1:nPixC <= 30), ids = 1L))
  d <- spreadFitValidationData(par, r, climFix$annualDTx1000, climFix$nonAnnualDTx1000, formC,
                               climFix$historicalFires, fb, mutuallyExclusive = NULL, Nreps = 2L, link = "powerLaw",
                               powerLawGamma = 1.5, climateRef = c(clim = 0.3), powerLawPc50 = p50, seed = 1)
  expect_identical(attr(d, "spreadFitValidation")$powerLawPc50, p50)
  expect_lte(max(d$p), p50)
  g <- plotSpreadFitResponse(d)
  hl <- Filter(function(l) inherits(l$geom, "GeomHline"), g$layers)[[1]]
  expect_equal(hl$data$yintercept[2], qlogis(p50))
})
