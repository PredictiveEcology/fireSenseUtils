# fireSenseUtils (development version)

* New `ignitionFitFilenameFor()` and `latestIgnitionFits()`, parallel to `spreadFitFilenameFor()`/
  `latestSpreadFits()`: name and read the shared, geo-keyed ledger of `fireSense_IgnitionFit`
  fits (ignition and escape together), so `fireSense_ELFs` can assemble per-ELF fit lists for a
  multi-ELF study area the same way it does for spread fits.
* New constant `ignitionFitAdditionalColNamesTxt` (`"fireSense_IgnitionFitted"`,
  `"fireSense_EscapeFitted"`): the list-column names of an ignition-fit ledger row.

# fireSenseUtils 0.2.3.9076

* New exported `spreadProbGates()` tells, without running any spread, whether a parameter set (or a list of them) passes the objective's spreadProb gates ("Too burny a landscape", "Not spread out enough", median out of range) in the first block of years. The gate test itself is now `spreadProbGateTest()`, which `.objfunSpreadFit()` calls, so the objective and the screen cannot disagree; the objective's results are unchanged. fireSense_spreadFit uses it to draw its threshold-calibration trials from the logistic's active range.

# fireSenseUtils 0.2.3.9075

* `.objfunSpreadFit()` no longer caps a simulated fire at `multiplier()` of its observed size; spread is bounded only by the year's buffers. A fire is a runaway, and censored as the cap-hit rule did, when it burns any pixel of the edge ring of its own buffer: a buffer pixel with a queen neighbour outside that fire's buffer (another fire's buffer, a removed non-flammable pixel, `NA`, the raster edge). `penaliseRunaways = TRUE` (default) replaces `penaliseCapHits`; `capSizes` and `penaliseCapHits` are deprecated and ignored with a warning, and `runDEoptim()` takes `penaliseRunaways`. The SNLL log line reports the runaway share. Without the cap the first-block SNLL has a different scale, so `thresh` values calibrated before this change will differ, and a calibration must use this objective. `spreadFitValidationData(capSizes =)` is removed (validation never capped).
* New `bufferEdge()` and `addBufferEdge()` compute the edge ring. `runDEoptim()` adds an `edge` column to `fireBufferedListDT` once, before the tables go to the workers (about 0.09 s per 90,000 buffer pixels; the check per evaluation is about 8 ms for 200,000 burned pixels); a table without the column gets it computed on every call.
* `.objfunSpreadFit(returnTerms = TRUE)` also returns `firstBlockSNLL`, the first block's SNLL per year as `thresh` is compared with it, and `bailed` (1/0), whether a year of that block was refused as "too burny" or "not spread out enough" or the block was above the threshold. Both are available with `thresh = Inf`.

# fireSenseUtils 0.2.3.9074

* `cohortsToFuelClasses()` builds the fuel classes on cell vectors instead of one `SpatRaster` per class, and has `asTable = TRUE` to return the `pixelID` table that `fireSenseCovariatesCreate()` makes from the stack, so that function no longer builds the rasters and converts them back (about 6 s of its 10 s on 6.7M cells). The raster and the table are identical to before.
* `fireSenseCovariatesCreate()` takes `fuelClassTable`, the result of `cohortsToFuelClasses(..., asTable = TRUE)`, so a caller making the `"species"` and `"domSecWetland"` covariates from the same cohorts builds the fuel classes once. The table is copied, not modified.

# fireSenseUtils 0.2.3.9073

* New exported `isYoungAge(age, cutoffForYoungAge)` is the one place the `youngAge` rule lives (`age <= cutoffForYoungAge`, `NA` not young); `makeTSD()`, `castCohortData()`, `cohortsToFuelClasses()`, `youngAgeAtYear()`, `calcYoungAge()`, `calcNonForestYoungAge()` and `fireSenseCovariatesCreate()` all call it. `calcNonForestYoungAge()` used `<`, so prediction differed from fitting for non-forest pixels whose age equals the cutoff. `castCohortData()` now gives `youngAge = 0` (was `NA`) for an `NA` stand age, and the deprecated `calcYoungAge()` no longer counts an `NA` age as young.
* `cohortsToFuelClasses()` no longer fails when terra is not attached: `as.int` is now imported from terra.

# fireSenseUtils 0.2.3.9072

* `.objfunSpreadFit()` treats a simulated fire that reaches its cap (`multiplier()` of the observed size) as censored, a runaway, instead of a fire of the capped size. Its size becomes `runawaySize` (default: the landscape's non-`NA` pixels) in the `adTest`, annual-area and `mad` terms, and in the per-fire size likelihood it has no density at the observed size (the likelihood is that of the other replicates times their share). `penaliseCapHits = FALSE` restores the old scoring; `capSizes = FALSE` and `returnSims` are unchanged. The per-fire likelihood is still floored at `minLik`, and `thresh` values calibrated before this change will differ. The SNLL log line now reports the capHit share. `runDEoptim()` takes `penaliseCapHits` (default `TRUE`) and passes it to the fit and the re-score; it changes the DEoptim cache key.

# fireSenseUtils 0.2.3.9071

* `youngAge` can be resolved per fire year again. New `youngAgeAtYear()` gives it for any set of pixels: time since disturbance at the data year, aged to the fire year and reset by every fire since (all fires, not only the fitted buffers); young is `<= cutoffForYoungAge` and `NA` time since disturbance is never young (the old `calcYoungAge()` called it young). `firePixelsByYear()` makes its per-year fire pixel lists from polygons or a fire-year raster. `calcYoungAge()` is deprecated.
* `fireSenseCovariatesCreate()` has `youngAge = TRUE`. With `FALSE` it builds no `youngAge` column and zeroes no fuel, non-forest land cover or treed wetland; the default keeps the previous behaviour.
* New `prepare_FuelCovsCoarseByYear()` builds the coarse ignition fuel rasters with `youngAge` (and the clearing it implies) per fire year; `stackAndExtract()` and `mergePreparedCovs()` accept these per-year lists. Fits made with `youngAge` fixed at the data year need refitting.

# fireSenseUtils 0.2.3.9070

* The pooled `other_agb` spread covariate is removed, and `fuelCovariates = "domSecOther"` is renamed `"domSecWetland"`: `fireSenseCovariatesCreate()` now builds only `dom_agb_<class>`, `sec_agb_<class>` and (with `rstLCC`) `treedWetland_agb`. `treedWetland_agb` still holds all tree AGB on treed-wetland pixels, including classes that are neither dom nor sec. Spread fits made with `other_agb` need refitting. `collapseFuelClassesToDomSecOther()` is now `collapseFuelClassesToDomSec()`.

# fireSenseUtils 0.2.3.9069

* `latestSpreadFits()` downloads a ledger file with `reproducible::preProcess()` (into a temporary folder, then placed in `destinationPath`) instead of `googledrive::drive_download()`, which wrote the file in place: a second job on the same ELF, sharing `destinationPath`, could find it missing or half-written. A local file whose MD5 differs from Drive's is fetched again with `purge = 7`.


# fireSenseUtils 0.2.3.9068

* `runELFs()` finds `fireSense_ELFs` among the project's modules or the children of any listed parent module
  (e.g. `PredictiveEcology/fireSense@development`), through the new internal `.findModuleInProject()`. It used
  `grep("ELFs", modules)`, which matched nothing once a parent was listed, so no module ran and the Drive upload
  failed with "missing value where TRUE/FALSE needed". `runELFs()` now stops naming the modules searched when there
  is no ELFs module, and stops with "No ELF output file" before `drive_update()` when the run saved none. Needs
  SpaDES.core (>= 3.2.1.9027) for `moduleMetadata(defineModuleListItems = "childModules")`.

# fireSenseUtils 0.2.3.9067

* `runDEoptim()`'s default `logPath` is named `runDEoptim_<time>_*.log` (was `fireSense_SpreadFit_<time>_*.log`),
  so no package code names a module; the module passes its own `logPath`. Documentation and the tutorial use the
  new module names (`fireSense_spreadFit`, `fireSense_ignitionFit`, `fireSense_spreadPredict`,
  `fireSense_ignitionPredict`).

# fireSenseUtils 0.2.3.9066

* `runDEoptim()` names the fitted parameters from `names(lower)`, grouped as logistic parameters,
  formula covariates and the year effect, e.g. "Fitting 11 parameters: logistic: maxAsymptote,
  inflectionPoint1; covariates: CMD, youngAge, ...; year effect: yearSpreadSD". It called
  `termsInDEoptim()`, which counted every non-formula parameter as a "logit" term and so reported
  yearSpreadSD as a third logistic parameter. Requires clusters >= 0.0.54, whose progress figures are
  labelled from `names(lower)` too.
* `termsInDEoptim()` is deprecated and will be removed in the next release.

# fireSenseUtils 0.2.3.9065

* `runDEoptim()`'s cache of the whole fit now includes the data the objective runs on (the digest
  `clusterSetup()` makes of the objects it ships to the workers). Before, two fits that differed only in
  that data, such as the two held-out folds of an ELF, shared one fit; clusters 0.0.53 fixes the same gap in
  the per-generation cache. Calls `clusters::DEoptimIterative()`, the new name of `DEoptimIterative2()`.

# fireSenseUtils 0.2.3.9064

* `runDEoptim()` gains `plotEvery` (default 25), passed to `clusters::DEoptimIterative2()` (exported
  since clusters 0.0.51; it was called with `:::`). The DEoptim progress figures are drawn every
  `plotEvery` generations and at the end, not after every generation, which took 16% of a fit's wall
  time. `plotEvery` is left out of the fit's cache key, so changing it does not refit. Needs clusters
  >= 0.0.52.

# fireSenseUtils 0.2.3.9063

* Added: `spreadFitValidationData()`, `plotSpreadFitValidation()` and `plotSpreadFitResponse()`
  (`?spreadFitValidation`). A response curve of `logit(p)` against a covariate looks right by
  construction, because its points are the model's own predictions. These compare the fit with the
  data instead: the share of pixel-years that burned against the share that burned in the fit's own
  simulations, binned by each covariate. The response curves are kept as a second figure, titled as
  the model's response. Both plot functions take the data first and return a `ggplot`, for
  `SpaDES.core::Plots()`.
* Added: `.objfunSpreadFit(returnBurned = TRUE)` records the pixels burned in each simulated
  replicate (`attr(, "burned")`). Off by default; the objective's value is unchanged.

# fireSenseUtils 0.2.3.9062

* Added: exported defaults shared by `fireSense_dataPrepFit` and `fireSense_dataPrepPredict`
  (`fireSenseForestedLCC`, `fireSenseYoungAgeCutoff`, `fireSenseNonForestCanBeYoungAge`,
  `fireSenseFlammabilityThreshold`, `fireSenseFuelClassCol`, `fireSenseIgAggFactor`,
  `fireSenseSCANFIVersion`; see `?fireSenseSharedDefaults`). The modules take their parameter
  defaults from them, so a fit and its predictions cannot silently use different values. The
  functions here with the same arguments (`makeTSD()`, `calcYoungAge()`, `castCohortData()`,
  `cohortsToFuelClasses()`, `makeFireSenseLCC()`) default to them too. Values are unchanged.

# fireSenseUtils 0.2.3.9061

* Added: `makeFireSenseLCC()` gains `scanfiVersion` (default `"V3"`, SCANFI's annual 1985-2025
  land cover), forwarded to `LandR::prepInputs_SCANFI_LCC_FAO(dataVersion = )`. LandR's own
  `dataVersion` check silently treats anything other than `"V1"` as `"V2"`, so `makeFireSenseLCC()`
  now checks the installed LandR for `"V3"` support itself and stops with a clear message instead
  of building V2 land cover under a `"V3"` request.
* Verified: SCANFI v3's new "burn scar" land-cover code (value not yet settled upstream) needs no
  change to non-forest grouping -- `assessFuelClasses()`/`makeLandcoverDT()` already group whatever
  non-forest codes are handed to them, so an unrecognized flammable non-forest code is kept and
  clustered like any other, not dropped. Added a test confirming this for a synthetic code.

# fireSenseUtils 0.2.3.9060

* Fixed: SCANFI's rock/exposed land-cover code (`30`, from `LandR::convert_SCANFI_LCC_codes()`)
  was missing from `makeFireSenseLCC()`'s and `ELFflammableArea()`'s `nonflammableLCC` default,
  so rock entered fits as flammable non-forest. New exported constant `fireSenseNonflammableLCC`
  is the single source of truth for the non-flammable codes of the FireSense land cover;
  `makeFireSenseLCC()` and `ELFflammableArea()` now default to it.

# fireSenseUtils 0.2.3.9059

## Removed

* Unused functions, with no callers found in PredictiveEcology code: `compareClimate()`, `compareMDC()`, `burnClassGenerator()`, `burnClassSummary()`, `burnClassPredict()`, `burnProbFromClass()`, `bufferIgnitionPoints()`, and the old ignition-fit objective functions `pw()`, `oom()`, `extractSpecial()`, `.objFunIgnition()`, `.objFunIgnitionPW()`, `objNlminb()`. `mclust` and `pROC` are no longer imported. To use one, install fireSenseUtils 0.2.3.9058.

# fireSenseUtils 0.2.3.9058

* Fixed: `getFirePolygons()` gave `POLY_HA` in m2, not ha, for `sf` polygons (`sf::st_area()` has no `unit` argument).
* Fixed: `rescaleCovariates()` failed for non-xgboost models (`inRange()` was not imported).
* Fixed: `stackAndExtract()` gave `NA` ignitions, not 0, when `fires` is `NULL`.

# fireSenseUtils 0.2.3.9057

* New spread-fuel representation: `fireSenseCovariatesCreate(fuelCovariates = "domSecOther")`
  collapses the per-species fuel-class columns to exactly four AGB terms, `dom_agb_<class>` and
  `sec_agb_<class>` (the ELF's two fuel classes with the most total treed AGB), `other_agb` (the
  rest, pooled) and `treedWetland_agb` (all tree AGB on treed-wetland pixels, removed from the
  other three there). New `chooseDomSecFuelClasses()` picks `domClass`/`secClass` once per ELF, so
  a prediction can be told the same classes the fit used instead of re-choosing them from its own
  area. The default remains `fuelCovariates = "species"` (the previous per-species columns),
  unchanged.
* New constant `treedWetlandAgbTxt` (`"treedWetland_agb"`).

# fireSenseUtils 0.2.3.9056

## Removed

* `DEoptimIterative()`: unused since `runDEoptim()` switched to `clusters:::DEoptimIterative2()`. Its plotting helpers go with it; `DEoptim` moves to Suggests.

# fireSenseUtils 0.2.3.9054

* Fixed: `ELFsInStudyArea()` failed when the study area overlapped a single ELF.
* Changed: `ELFsInStudyArea()` labels a buffer cell shared by several ELFs with the ELF whose core is nearest, rather than the one listed last.

# fireSenseUtils 0.2.3.9053

* Fixed: `harmonizeFireData()` lost the fire polygons of later years when a whole year was dropped for lying outside the study area.
* New `harmonizeFireDataDeps()`: the functions `harmonizeFireData()` calls, for a cached call's `.cacheExtra`.

# fireSenseUtils 0.2.3.9052

* Fixed: `cleanUpSpreadFirePoints()` matched points by row number instead of fire ID, so ignition points on non-flammable pixels were never moved and fires with no flammable pixel were never dropped.

# fireSenseUtils 0.2.3.9049

* `hillSlope1` (the spread link's slope) is fixed at 1, not fitted. With the linear predictor
  `x = covariates %*% beta`, `hillSlope1` enters `logistic3p()`/`logistic3pUpper()` only as
  `hillSlope1 * x`, so scaling every covariate coefficient by `k` and dividing `hillSlope1` by `k`
  leaves every prediction unchanged: it was never identifiable, and letting DEoptim fit it let every
  coefficient drift along that ridge. `.objfunSpreadFit()` now calls the new internal
  `fixHillSlope1()` on every `par` it receives, reinserting `hillSlope1 = 1` as the 2nd logistic
  parameter, so callers (`fireSense_SpreadFit`) no longer include it in `lower`/`upper`.

# fireSenseUtils 0.2.3.9048

* Fixed: `makeMutuallyExclusive()` recomputed the rows to zero from the key column after each
  pattern, so once the key column itself matched one of its own patterns (as `youngAge` could) it
  was zeroed and every later pattern zeroed nothing. Rows are now fixed once, before any zeroing,
  and a pattern can no longer match the key column itself.
* Fixed: `fireSenseCovariatesCreate()` applied `makeMutuallyExclusive()` before `youngAge` was
  finalized for non-forest pixels, so those pixels kept their `nfLCC_*` land-cover value. The
  exclusivity step now runs after `youngAge` is final, and also zeroes `treedWetland`: a burned bog
  is no longer treated as "still wet" once it is young.
* New: `youngAgeExclusiveCols()`, exported, builds the `youngAge` mutually-exclusive column list
  from a covariate name vector and a caller-supplied fuel-column list (`fireSense_SpreadPredict`
  uses this; `fireSense_SpreadFit` derives its own list from its non-annual covariate table).

# fireSenseUtils 0.2.3.9047

* `plotELFs()` gains `which` and `fill`, to fill the named ELFs, e.g. `plotELFs("inputs", which = c("13.1", "4.1"), fill = "green")`.

# fireSenseUtils 0.2.3.9046

* `spreadFitFileTag` is now `"_linearFuel_esc50"`: fits made with the 50-ha escape and the new objective terms go
  to a new ledger file (`fireSenseParams_<years>_linearFuel_esc50.rds`), and `"latest"` no longer reads the
  `_linearFuel` files of the previous model. An ELF with only an old fit counts as unfitted and is fitted again.

# fireSenseUtils 0.2.3.9045

* `.objfunSpreadFit()` and `runDEoptim()` gain options, off by default:
  - `jumpTries`, `jumpMeanDist`: with `escapeSizeHa`, passed to `SpaDES.tools::spreadCpp()` so a fire stuck
    below the escape size can jump to nearby burnable land.
  - `yearAreaWeight`: scores each fit year's observed area burned against its simulated totals with the per-fire
    size likelihood (`yearAreaNLL()`); a year the objective refused to simulate gets `minLik`. `"auto"` weights
    it equally with the per-fire term.
  - `areaDistWeight`: compares simulated and observed sizes by the share of area burned they make up
    (`areaWeightedCvM()`), which weights the large fires. `"auto"` uses the AD term's weight.
  - `returnTerms = TRUE` returns each term separately.
* With `escapeSizeHa`, a fire's size cap is raised to the escape size where it was smaller, as `spreadCpp()`
  needs `minSize <= maxSize`.
* Changed: the Anderson-Darling term now uses only simulated fires larger than 1 pixel, as the per-fire
  likelihood already did; fires that never spread were counted as a size mismatch. Objective values of fits
  that use `adTest` change.

# fireSenseUtils 0.2.3.9044

* New argument `escapeSizeHa` in `.objfunSpreadFit()` and `runDEoptim()` (default `NULL`, the old behaviour). When
  set, the spread model is fitted as "given the fire escaped": only observed fires of at least
  `escapeSizePixels(escapeSizeHa, landscape)` pixels are fitted, and every simulated fire burns its first that
  many cells regardless of spread probability (`SpaDES.tools::spreadCpp(minSize =)`) before spreading normally.
  Before, an escaped fire was any fire over 1 pixel, and 11-16% of simulated fires never left their first pixel.
  New exported helper `escapeSizePixels()`, so fireSense uses the same threshold.
  Needs SpaDES.tools >= 2.1.3.9009 (`spreadCpp(minSize =)`).

# fireSenseUtils 0.2.3.9043

* `latestSpreadFits()` reads the shared SpreadFit ledger when a module's `spreadFitFilename` is `"latest"`: for
  each polygon, its rows from the most recently modified `fireSenseParams_*_linearFuel.rds` file in the Drive folder
  that has it, reading newest first and stopping once the polygons asked for are found. Files without the
  `spreadFitFileTag` (`"_linearFuel"`) hold log-fuel fits and are never read. `spreadFitFilenameFor(fireYears)`
  names the file a fit is written to.

# fireSenseUtils 0.2.3.9042

* Needs `LandR >= 1.2.0.9022` (was `>= 1.1.0.9066`). The default land cover, SCANFI with wetland classes,
  calls `LandR::prepInputs_CWIM()` and `LandR::wetlandToLCC()`, which arrived in 1.2.0.9022
  (PredictiveEcology/LandR#228). With an older LandR the package installed and the run stopped at the land
  cover instead.

# fireSenseUtils 0.2.3.9041

* The spread objective's random effect is now per YEAR, not per fire, and is named `yearSpreadSD`
  (was `fireSpreadSD`; `fitYearSpreadSD` was `fitFireSpreadSD`). It is a seasonal departure: each replicate draws
  one `eps ~ N(0, yearSpreadSD^2)` for the year, and every pixel of that year spreads with
  `plogis(qlogis(p) + eps)`, so all of a year's fires are bigger or smaller together. fireSense does the same in a
  forecast. `yearSpreadSD = 0` is the model without it, exactly.

# fireSenseUtils 0.2.3.9040

* Per-fire random effect in the spread objective. Weather during a fire is not in the model, so every fire of a year
  saw identical spread probabilities and the simulated sizes were too alike: medians too large and the largest fires
  too small at once; weighting big fires (`weighted = "sqrt"`) only shifted the whole distribution. With
  `fireSpreadSD` as the last parameter, each fire draws `eps ~ N(0, fireSpreadSD^2)` per replicate and every pixel
  of its buffer spreads with `plogis(qlogis(p) + eps)`; `spreadCpp()` is unchanged. `.objfunSpreadFit()` uses it
  when `par` names it (`fitFireSpreadSD = NULL`); `runDEoptim()` when `lower` names it. `fireSpreadSD = 0` is the
  previous model exactly.

# fireSenseUtils 0.2.3.9039

* Fit diagnostics (`?fitDiagnostics`), the checks made by hand on the phase-2 fits, as functions:
  `simulateFireSizes()` (the observed fires simulated without the size cap or the "too burny" gate),
  `scoreFireSizes()` (per-fire and per-year error, quantiles, AD), `linkSaturation()` (share of pixel-years at the
  spread-probability ceiling), `coefIdentifiability()`, `profileCoefficients()` and `identifiedInIsolation()` (is
  each covariate's coefficient pinned by the population, and does dropping it worsen the fit), `fitConvergence()`.
* `.objfunSpreadFit(returnSims = TRUE)` returns the simulated fires instead of the objective; `capSizes = FALSE`
  lifts the size cap. These replace the `trace()` used for held-out validation.
* `runDEoptim()` gains `profileReps` and `simulateMembers` (both off by default), which run the profile and the
  uncapped simulations on the fit's workers after the re-score, as `attr(DE, "profile")` and `attr(DE, "fitSims")`.

# fireSenseUtils 0.2.3.9038

* New spread link `logistic3pUpper()`: `logistic3p()` with Stukel's (1988) generalized-logistic upper tail, one
  extra parameter `upperTail1`. In `logistic3p()` the curve approaches its ceiling at a rate set by the slope alone,
  so no parameter could change the upper end without moving everything else; in fitted models most spreadable
  pixels sat pressed against the ceiling (72% of pixel-years in ELF 5.2.1). `upperTail1 < 0` slows that approach,
  `> 0` speeds it, and `0` is `logistic3p()` exactly. `logisticAll()` chooses it when the logistic parameters
  include `upperTail1` or when `link = "logistic3pUpper"`; `.objfunSpreadFit()` and `runDEoptim()` pass `link`
  through to the fit and its re-score. Checked against `sirt::pgenlogis()`.

# fireSenseUtils 0.2.3.9037

* `runDEoptim()` gains `sizeLik`, `sizeLikDf`, `weighted` and `adWeight`, passed to `.objfunSpreadFit()` both in
  the fit and in the re-score of the final population. Before this, every fit and every re-score used the
  objective's defaults whatever the caller wanted -- `weighted = TRUE` (a log(size) weight) with `sizeLik = "kde"`.
  The defaults here are the objective's, so a caller who sets nothing sees no change.

# fireSenseUtils 0.2.3.9036

* The spread objective uses `SpaDES.tools::spreadCpp()` instead of `spread()`. It follows the same rules
  but makes its own random draws, so **results change** and fits made before this are not comparable
  with fits made after. A full objective evaluation is 2.3-2.5x faster; end to end through DEoptim it saves
  a median ~9.5 s per evaluation on both ELFs tested. Requires `SpaDES.tools (>= 2.1.3.9008)`.

# fireSenseUtils 0.2.3.9035

## New features

* `.objfunSpreadFit()` gains `adWeight`: what the `"adTest"` statistic is multiplied by before it is
  added to the fire-size SNLL. It defaults to `"auto"`, the new `adWeightAuto()` = `c * sqrt(nFires)`,
  which holds the AD term's share of the influence on the objective to 0.36-0.68 across six ELFs and
  the four `sizeLik` x `weighted` combinations. **This changes the default objective**: it was a fixed
  50, which gave 0.16-0.94 over the same range, balancing the two terms only for `sizeLik = "kde"`
  with `weighted = FALSE`. Pass `adWeight = 50` for the old behaviour.

# fireSenseUtils 0.2.3.9034

## Bug fixes

* `adStatistic()` (0.2.3.9031) multiplied integer counts, which overflow past about 46,000 simulated
  fires (e.g. 1,800 fires at `Nreps = 50`): the adTest term became `NA` with a warning, and DEoptim
  stops on a non-finite objective. It now computes in doubles. (9032 and 9033 are #73 and #74.)

# fireSenseUtils 0.2.3.9033

## Bug fixes

* `weighted = TRUE` did not weight. It computed `log(lik * log(size))`, which is
  `log(lik) + log(log(size))`: an offset, under which every fire's likelihood moved the SNLL by the
  same amount whatever its size. The weight now multiplies the log-likelihood. `weighted` takes
  `FALSE`, `TRUE` or `"log"` (`log(size)`), or `"sqrt"` (`sqrt(size)`); weights are divided by their
  mean over the fitted fires, so the SNLL keeps its scale. Objective values change when
  `weighted` is not `FALSE`.

# fireSenseUtils 0.2.3.9032

## New features

* `.objfunSpreadFit()` and `objFunInner()` gain `sizeLik` and `sizeLikDf`. `sizeLik = "t"` takes a
  fire's size likelihood from a Student-t on the square-root scale (mean and standard deviation of
  its simulated sizes, `sizeLikDf` degrees of freedom, default 5) instead of their kernel density.
  The kernel density is zero away from the simulated sizes, so a fire they never reach scores the
  `minLik` floor however far off it is; the t has no floor and penalises a miss by its distance.
  The default, `"kde"`, is unchanged. (0.2.3.9031 is the adTest fix, #72.)

# fireSenseUtils 0.2.3.9031

## Bug fixes

* The `"adTest"` term of `.objfunSpreadFit()` compared observed fires with each simulated fire's
  *mean* size over `Nreps`, which has a much shorter tail than single fires: on ELF 4.3 a model
  compared with one of its own replicates scored `50 * AD` = 3260. It now gets every replicate's
  fire (42 on the same comparison). Objective values change.
* New internal `adStatistic()` computes the Anderson-Darling statistic directly. `kSamples::ad.test()`
  also standardises it, in time quadratic in the sample: 4 s per evaluation with the larger sample,
  against 0.003 s. `kSamples` moves to Suggests.
  
# fireSenseUtils 0.2.3.9030

## Bug fixes

* `objFunInner()`: the fire-size likelihood compared the observed size with the square root of the
  simulated sizes (`demp(x = size, obs = sqrt(N))`, since 2021-02-10). Simulated fires are capped at
  `multiplier(size)`, so any fire above ~12 pixels scored the `minLik` floor whatever the parameters;
  only the smallest fires informed a fit. Both are now square-rooted. Objective values change, so
  the early-bail `thresh` of an existing fit no longer applies.

# fireSenseUtils 0.2.3.9029

## New features

* `fuelLogToLinear()`, `fuelLinearRange` and `isLinearFuelRange()`: the spread model now takes fuel
  biomass on the linear scale, divided by a fixed 1e4. `fireSenseCovariatesCreate()` still returns it
  as `logMinB()`, because it also builds the ignition covariates and its output is cached for every
  fitted polygon; `fireSense_SpreadFit` and `fireSense_SpreadPredict` both undo the log with
  `fuelLogToLinear()`. `fuelLinearRange` (`c(0, 1e4)`) is the `covMinMax` that performs the division,
  and `isLinearFuelRange()` recognises a linear fit from its stored `covMinMax_spread`, so parameters
  fitted on the log scale keep predicting as before. (0.2.3.9028 is `covCentre`, #67.)
# fireSenseUtils 0.2.3.9028

## Enhancements

* `.objfunSpreadFit()`, `objFunInner()` and `spreadProbFromIntegerCovs()` gain `covCentre`: values
  subtracted from the rescaled covariates. `NULL` (default) changes nothing. Centring is applied after
  the mutual-exclusivity step and the covariate range assertion.

# fireSenseUtils 0.2.3.9027

## Performance

* `objFunInner()` runs `SpaDES.tools::spread()` on the bounding box of each fire year's pixels
  instead of the whole landscape (new internal `cropToCells()`). `spread()` allocates
  landscape-length state on every call, and it is called `Nreps` times per fire year. With a
  one-cell margin the crop draws the same random numbers, so results do not change: identical
  objective values with identical seeds on ELFs 5.3.1, 5.3.2 and 13.1, and an evaluation 1.6-2.7x,
  2.6-3.8x and 1.4-1.7x faster.

# fireSenseUtils 0.2.3.9026

## Performance

* `objFunInner()` no longer scans the whole landscape once per fire year. The spreadProb values its
  bail tests summarise were recovered with `cells[cells > a | cells > b]`, four passes over a
  landscape-length vector that is zero except at that year's pixels (6.0M cells for ~41k pixels on
  ELF 5.3.1); they are now read from the spreadProb column. The vector is filled only when
  `spread()` will run, and is allocated numeric so filling it does not coerce it. Identical
  objective values with identical seeds; an evaluation is 1.1-1.4x faster on ELF 5.3.1.

# fireSenseUtils 0.2.3.9024

## Bug fixes

* new `makeFireSenseLCCDeps()`: the functions whose code affects what `makeFireSenseLCC()` returns,
  to pass as a cached call's `.cacheExtra`. `reproducible::Cache()` digests only the called
  function's own code, so callers have to name the functions it calls -- and which those are depends
  on `lccSource`, a run-time option. A caller that writes the list out is therefore wrong for the
  other source and goes stale when the default moves. That already happened: `fireSense_dataPrepFit`
  pinned `LandR::prepInputs_NTEMS_LCC_FAO()` and kept it after the default became SCANFI, so every
  SCANFI and CWIM change was invisible to the cache while NTEMS changes invalidated it for nothing.

# fireSenseUtils 0.2.3.9023

## Enhancements

* `makeFireSenseLCC()` gains `lccSource`, `"SCANFI"` by default (`options(fireSense.lccSource = )`), with
  `"NTEMS"` still available and unchanged. Biomass_borealDataPrep's default land cover is moving from NTEMS
  to SCANFI, and the fire models must be fitted on the land cover the simulation predicts with. SCANFI has
  no wetland classes, so they are added exactly as that module does: `LandR::prepInputs_CWIM()` for the
  wetland site layer and `LandR::wetlandToLCC()` to recode it, wet treed pixels to 81 and other wet
  flammable pixels to 80. Non-flammable (0) and water are never recoded. Checked on ELF 5.3.2's fit grid
  (8.1 M cells): the NTEMS path gives byte-identical output to 0.2.3.9022.

  The two LandR wetland functions are looked up at run time (PredictiveEcology/LandR#228), so this package
  installs and checks against a LandR without them; `lccSource = "SCANFI"` then stops with a message naming
  what is missing.

* `fireSenseCovariatesCreate()` gains `rstLCC` (and `treedWetlandLCC`, default 81). Given it, the covariates
  carry `treedWetland` (new constant `treedWetlandTxt`), 1 on treed wetland. Class 81 is a forested class, so
  those pixels otherwise reach the fire models only through their fuel biomass and look like upland forest.
  It is a site attribute, not a fuel state, so it is added after the `youngAge` exclusivity -- a burned bog is
  still wet -- and `mergePreparedCovs()` leaves it out of the ignition data's all-cover-is-zero filter.
  Without `rstLCC` nothing changes, so existing callers and their caches are untouched.
# fireSenseUtils 0.2.3.9022

## Enhancements

* `runDEoptim()` re-scores the final population after the fit: each member is evaluated `rescoreReps`
  (default 10) more times, on the same workers, with the early stop off, and the result is attached as
  `attr(DE, "finalRescore")`. New `rescorePopulation()` and `bestByReplicatedMean()` do the work, the
  second picking the `n` distinct members with the lowest replicated means. DEoptim evaluates a
  surviving member only once, so the values it holds are partly luck: in eight FireSense fits its best
  member ranked 1st to 6th of 60 by replicated mean, and the fitted-parameter ledger's "5 best" were five
  copies of that one member. The re-score costs about ten generations' worth of evaluations.

* `spreadFitAdditionalColNamesTxt` gains `covMinMax_spread`: prediction cannot rescale covariates as the
  fit did without it, and the ledger never stored it. Ledger readers must tolerate rows that predate it
  (see `fireSense_dataPrepFit`).

# fireSenseUtils 0.2.3.9020

## Bug fixes

* `rasterFireBufferDT()` drew its per-year seeds from the session's random stream, so two calls with
  identical inputs built different buffers unless the caller happened to have set the same seed --
  and nothing upstream did. `fireSense_dataPrepFit` caches this call, so its spread-fit buffers,
  spread points and covariates changed between runs of the same study area, and cached entries holding
  different data were kept apart only by an accident of their cache keys. #57 made the buffers the same
  across forked workers; this makes them the same across sessions. The new `seed` argument defaults to
  a digest of the inputs, the draws run under `withr::with_seed()`, and the session's random stream is
  no longer advanced. Buffers differ from those built by earlier versions, so cached buffers are rebuilt.

# fireSenseUtils 0.2.3.9019

## Enhancements

* `.objfunSpreadFit()` gains `pruneAbove` (default `Inf`), which makes its existing early-bail
  test adaptive. After the first of its two batches of fire years it compares the accumulated
  SNLL against `thresh * numYrsDone` -- a bound fixed for the whole fit, which therefore cannot
  tighten as the population improves. The test is now `min(thresh * numYrsDone, pruneAbove)`, so a
  caller running one DEoptim generation per call can pass the worst value the current population
  would accept. This is exact rather than heuristic: the second batch contributes a non-negative
  SNLL, so the first batch's value is a lower bound on the total, and any trial above `pruneAbove`
  would have been rejected by selection anyway -- the search trajectory is unchanged. It is worth
  doing because a DEoptim generation is synchronous: its wall time is the slowest of its `NP`
  evaluations, so the tail sets the clock.

# fireSenseUtils 0.2.3.9018

## Bug fixes

* `runDEoptim()` asks for about 10 workers per estimated parameter (`nCoresNeeded`, default
  `10 * length(lower)`) instead of a fixed 100, and DEoptim's `NP` is the number of workers the
  cluster was built with (clusters >= 0.0.33). It passes only the settings it chooses to
  `clusters:::DEoptimIterative2()`, whose defaults fill the rest; a complete `DEoptim.control()`
  list would have overridden them.
* `runDEoptim()` passes DEoptim settings to DEoptim. `.c` was sent to the objective function, which
  ignores it, so DEoptim always used its default `c`; it is now DEoptim's `c`. The new
  `DEoptimControl` list carries any other `DEoptim.control()` setting (`CR`, `F`, `p`, `reltol`, ...)
  through `clusters::clusterSetup(controlArgs = )`.
* `runDEoptim()` no longer fails with "missing value where TRUE/FALSE needed" where R has no
  OpenMP, as in the CRAN macOS builds: `RhpcBLASctl::omp_get_max_threads()` returns `NA` there.

# fireSenseUtils 0.2.3.9017

## Bug fixes

* `ELFmergePlan()` no longer fails with "non-character argument" when an ELF with too few fires has no
  neighbour at all, such as an ELF whose core touches no other ELF. That ELF is not fitted, as when its
  neighbours share no base with it.

# fireSenseUtils 0.2.3.9016

## New features

* ELFs with too few fires can be merged with a neighbour. `ELFneighbours()` measures the core border each
  pair of ELFs shares. `ELFmergePlan()` takes every ELF that `ELFfitStatus()` calls `"zero"` or `"few"`
  and pairs it with the neighbouring ELF of the same base and depth (another piece of the same split
  ecoprovince, or another whole ecoprovince of the same ecozone) that shares the longest border. If the two
  together reach both thresholds they merge; otherwise neither is fitted. `mergeELFs()` applies the plan to
  the ELF maps, `ELFmergedName()` names a merged ELF by its shared base and members' last parts (3.2.1 with
  3.2.4 is `"3.2.1_4"`), `ELFsSkipped()` lists the ELFs not fitted and `ELFrunName()` maps a merged member
  to its merged ELF.

## Bug fixes

* `ELFfireCounts()`, `ELFfitStatus()`, `ELFsExcluded()` and `ELFflammableArea()` are exported. They were
  documented as exported in 0.2.3.9014, but `NAMESPACE` was not regenerated, so no other package or module
  could call them.

# fireSenseUtils 0.2.3.9015

## Bug fixes

* `bufferToArea()` and `rasterFireBufferDT()` give the same buffers for the same seed. They pick buffer
  pixels at random, in forked workers when `cores > 1`, and each forked child seeded itself
  independently: the same call gave different buffers on every run, and different ones again with one
  core (61 of 5,000 buffer pixels differed between two runs of one year of ELF 11.2). Each polygon set,
  or year, now runs under a seed drawn from the session's stream, forked or not. Spread-fit data built
  with these functions change once as a result.

# fireSenseUtils 0.2.3.9007

## Bug fixes

* `getFirePoints_NFDB()` and `getFirePoints_NFDB_V2()` download the National Fire Database
  points from `.../current_version/NFDB_point_shp.zip`. CFS renamed the archive from
  `NFDB_point.zip`, which now returns HTTP 404, so no release after the copy already on disk
  (fires to 2024) could be fetched. The URL is in the internal `nfdbPointUrl()`, with a test.

# fireSenseUtils 0.2.3.9006

## Maintenance

* `parallel` is now declared in `Imports` (it was imported in `NAMESPACE` only). `covr` records
  coverage in forked children only for packages that declare `parallel`, so the forked code in
  `bufferToArea()` and `rasterFireBufferDT()` was reported as never run.

# fireSenseUtils 0.2.3.9005

## Bug fixes

* `bufferToArea()` and `rasterFireBufferDT()` with `cores > 1` no longer hang. GDAL keeps one
  worker-thread pool per process, created at the first multi-threaded raster write. A forked
  worker inherited that pool but none of its threads, and waited for them forever with 0 CPU.
  Forked workers now write rasters single-threaded (`terraOptions(threads = 1)`).

# fireSenseUtils 0.2.3.9004

## Bug fixes

* `fireSenseCloudParameters()` now downloads the shared parameter file from Google
  Drive on every call. It used `prepInputs(purge = 7, overwrite = TRUE)`, which never
  downloads again once a copy on disk matches CHECKSUMS.txt (`purge` only rebuilds
  those entries; `overwrite` only affects the written output), so a changed file on
  Drive was not seen. `url` may now also be the folder containing `targetFile`;
  `useCache` is ignored.

# fireSenseUtils 0.2.3.9003

## Performance

* `objFunInner()` no longer wraps its `spread()` replicate loop in
  `system.time()`. `system.time()` defaults to `gcFirst = TRUE`, so every fire year
  forced a full garbage collection, and the result was assigned to a variable that
  was never read. Measured on ELF 12.4 (1,737,132 cells, 18 fire years, `Nreps = 25`,
  in a 20-40 GB process), three consecutive objective-function evaluations went from
  82.7 / 84.0 / 83.9 s to 9.8 / 8.7 / 8.8 s -- about 9.5x -- with warm caches and no
  other change. In the profile the forced collections were 109.96 s of a 134 s
  evaluation, 82% of self time, against 6.44 s for all the `spread()` calls they were
  timing. It also explains why replicate count used to make no difference: the
  collection is once per fire year, outside the replicate loop. With it gone,
  `Nreps = 5` takes 3.2 s against 8.8 s at 25, so roughly 80% of an evaluation now
  scales with replicates as it always should have.

# fireSenseUtils 0.2.3.9002

## Bug fixes

* `objFunInner()` returned `ret` unconditionally but built it only inside
  `if (isTRUE(doFitting))`. With no test selected -- which is how
  `fireSense_SpreadFit`'s `mode = "debug"` calls the chain, passing `tests = ""`
  -- the branch is skipped and the call failed with `object 'ret' not found`,
  making debug mode unusable. `ret` is now initialised before the branch, so the
  function returns an empty list when nothing was asked of it.

# fireSenseUtils 0.2.3.9001

* `objFunInner()` now asks `SpaDES.tools::spread()` to skip its input checks by
  the name that function actually uses. It passed `skipChecks = TRUE`, which
  belongs to `spread3()`; `spread()` has no such argument, so it landed in `...`
  and did nothing. Every call therefore re-validated the whole per-cell
  `spreadProb` vector -- `na.omit()` copies it, `inRange()` scans it -- which is
  work proportional to the number of cells in the landscape, independent of how
  much actually burns, and this function calls `spread()` `Nreps` times per fire
  year. Measured on synthetic landscapes, with identical output:

  | landscape | before | after | speedup |
  | --- | --- | --- | --- |
  | 1M cells | 0.030 s | 0.016 s | 1.9x |
  | 4M cells | 0.074 s | 0.022 s | 3.4x |
  | 9M cells | 0.183 s | 0.036 s | 5.1x |

# fireSenseUtils 0.2.3

* `objFunSpread()`: the Anderson-Darling test now compares fire size distributions
  drawn from the *same* years on both sides. Years are simulated in two batches --
  the 2 with the largest total area burned, then the rest -- and the `adTest` block
  was gated on `ii == 2`, reading `results$fireSizes` after `results` had been
  overwritten by the second batch. So the simulated sample held only the
  small-fire years, while the observed sample used `historicalFiresAboveMin`
  unsubset, i.e. *every* year. Since the omitted years are selected as
  largest-area-burned, the observed sample retained an upper tail the simulated
  sample structurally could not have; `ad.test` is tail-sensitive and its result is
  multiplied by 50, so the only way `DEoptim` could shrink it was to inflate
  simulated fire sizes in ordinary years, biasing the fit toward over-prediction.
  Simulated sizes are now pooled across both batches (they were already being
  computed and discarded) and the test runs once after the loop, via the new
  internal `pooledFireSizes()`. Behaviour on the early-bail path is unchanged: a
  parameter set that fails the first batch still skips `adTest` entirely.
* `makeELFs()`: only run the NA-hole-filling `focal()` step when `x` is a
  `SpatRaster`. When `x` is (or defaults to) an `sf` of fire regime polygons there
  are no NA slivers to fill and `focal()` has no method for it, failing with
  `no method found for signature sf, data.frame`.

* `mergeAndSplitRas()`: write each call's per-ecoprovince `.tif` files into their own
  `ELFs_<digest>` subdirectory of `destinationPath`, keyed on `ecopRseg`, `ecopLCC`,
  `maxArea` and `field`. They previously went straight into `destinationPath` named
  only by the province code (`4.1.tif`), scattering ~36 anonymous files through a
  shared inputs directory and letting successive calls overwrite one another. Note
  this protects `destinationPath` only: `reproducible`'s cache-restore path collapses
  these to `cachePath/<basename>.tif`, dropping the cacheId, so distinct calls still
  collide there -- that is a `reproducible` issue this cannot fix.

* `fireSenseCloudParametersMap()` and `plotELFs()`: use `ELFs$poly` rather than
  treating the `makeELFs()` return value as a `SpatVector`. `makeELFs()` returns a
  list of rasters plus a `poly` element (`ELFsInStudyArea()` already unwraps it this
  way), so `terra::plot(ELFs)`, `ELFs$buffer`, `ELFs$ID` and `ELFs[keep, ]` were all
  operating on the list. Both call sites are fixed.

* Documentation regenerated with roxygen2 8.1.0 (was 8.0.0). Mostly formatting:
  8.1.0 emits one multi-symbol `importFrom()` per package instead of one line per
  symbol, so NAMESPACE shrinks considerably with no change to what is imported.
  One **API change** falls out of the upgrade: `spreadFitAdditionalColNames`
  (`R/objDefaults.R`) carries `@export` in its source but was not being exported --
  roxygen2 8.0.0 silently dropped it, and 8.1.0 honours it. It is now exported and
  documented, alongside the existing `spreadFitAdditionalColNamesTxt`.

* `getFirePoints_NFDB()` no longer drops columns. It previously subset to
  `c("YEAR", fireSizeColName)` and renamed those to `date`/`size_ha`, discarding
  `CAUSE` among everything else -- the reason the scfm modules explicitly avoid
  fireSenseUtils (`## NOTE: do not use fireSenseUtils - it removes the cause
  column`). All source columns are now retained, matching
  `getFirePoints_NFDB_V2()` and `scfmutils::getFirePoints_NFDB_scfm()`. The
  derived `size` column (fire size in pixels) is still added when `rasterToMatch`
  is supplied. **Breaking:** the `date` and `size_ha` output names are gone; use
  `YEAR` and `SIZE_HA` (#32).

* `getFirePoints_NFDB_V2()` now uses `fun` to load already-downloaded data, not
  just downloaded data. The cached branch hardcoded `st_read()`, so a caller
  passing `fun = "terra::vect"` (as `fireSense_dataPrepFit` does) got a
  `SpatVector` on the first run and an `sf` on every run after. The return class
  is now the same in both branches (#32).

* `getFirePoints_NFDB()` and `getFirePoints_NFDB_V2()`: fix inverted
  `redownloadIn`. The staleness threshold was `365 / redownloadIn`, so
  `redownloadIn = 0.5` gave a 730-day tolerance rather than the documented
  "redownload data older than 6 months". Now `365 * redownloadIn`. Only the
  default of `1` was unaffected (#32).

* `getFirePoints_NFDB()` and `getFirePoints_NFDB_V2()`: the staleness check
  errored with "the condition has length > 1" when `NFDB_pointPath` held more
  than one `NFDB_point*.shp`; now wrapped in `any()` (#32).

* `getFirePoints_NFDB_V2()`: restore the `NFDB_pointPath` non-NULL check, so the
  `NULL` default fails with a clear message instead of inside `Checksums()` (#32).

* `getFirePoints_NFDB()`: pass `useSAcrs = TRUE` to `postProcess()` in the
  already-downloaded branch, so the CRS no longer depends on whether a download
  happened; restore `reproducible.cacheSaveFormat` to its previous value on exit
  instead of overwriting it; drop the dead `SpaDES.core` requirement check (the
  package is only in Suggests and nothing in the function uses it) (#32).
* Replace the single `ggpubr::ggarrange()` call in `visualizeDE()` with
  `cowplot::plot_grid()`, and move `cowplot` from Suggests to Imports (it was
  already used, behind a `requireNamespace()` guard, in `plot_summaries.R`).
  `ggpubr` is dropped from Imports. This removes 14 packages from the
  dependency tree -- ggpubr, car, carData, pbkrtest, doBy, Deriv, lme4,
  quantreg, SparseM, MatrixModels, nloptr, minqa, abind and Formula.

  It also fixes R-CMD-check on the `oldrel-2` leg. `Deriv 4.3.0` (published
  2026-07-23) uses `R_ClosureFormals`, an R >= 4.5 C-API entry point, but
  declares no minimum R version, so it fails to compile on R 4.4.x with
  `'R_ClosureFormals' was not declared in this scope`. The dependency install
  then aborted and `rcmdcheck` died at `loadVignetteBuilder()`, which made the
  failure look like a vignette problem.

* `makeTSD()` is now general purpose: the pixels to age from fire history and the
  flammable mask can be supplied directly via the new `pixToUpdate` and
  `flammablePixels` arguments, instead of only through `lcc` (the `landcoverDT`
  from `fireSense_dataPrepFit`). `lcc` is now optional and, when supplied, still
  derives both vectors (the `fireSense_dataPrepFit`-specific behaviour is
  preserved behind an `if`); explicit arguments take precedence. Fully backwards
  compatible (#18).

* `Firesense_LCC_flammability` vignette: probe the LCC (`ftp.maps.canada.ca`) and
  fire-polygon (`cwfis.cfs.nrcan.gc.ca`) data servers up front and skip the live
  download/analysis chunks when either is unreachable, degrading to
  documentation-only (same as the existing macOS/`archive` path). This stops
  R-CMD-check from failing when the servers are unreachable from CI runners.
  Adds `curl` to Suggests.

* `makeFireSenseLCC()`: remove commented-out dead code (an unused majority-NA
  block aggregation path) from the aggregation step.

# fireSenseUtils 0.1.5

* Fix namespace conflict warning: remove blanket `import(data.table)` in favour of
  explicit `importFrom` declarations; remove `importFrom(purrr, transpose)` which
  was silently overriding `data.table::transpose` (all call sites already use
  `purrr::transpose()` with explicit namespace).
* Add missing `importFrom(data.table, setorder)` used in `DEoptimIterative`.

# fireSenseUtils 0.0.5

* TODO
