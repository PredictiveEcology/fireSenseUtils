#' Hard-coded name/text constants
#'
#' Character (or character-vector) constants used across the package and
#' downstream modules to keep column names, identifiers, and other free-text
#' tokens consistent. All such constants end in `Txt`.
#'
#' @format
#'   - `polygonIDTxt`: `character(1)`. Name of the polygon-ID column
#'     (currently `"polygonID"`).
#'   - `nonNFColNamesTxt`: `character` vector of bookkeeping column names
#'     that are *not* non-forest landcover classes (currently `"pixelID"`
#'     and `polygonIDTxt`); used with [base::setdiff()] to select the
#'     non-forest landcover columns when summing/aggregating across rows
#'     (see [makeTSD()]).
#'   - `yearTxt`: `character(1)`. Token used as a year-column name or as a
#'     prefix on year-suffixed layer/column names (currently `"year"`).
#'   - `youngAgeTxt`: `character(1)`. Name of the "young-age" cohort class
#'     (currently `"youngAge"`).
#'   - `burnableHaTxt`: `character(1)`. Name of the column of `spreadFirePoints`
#'     holding each fire's burned area (ha) on the pixels the spread fit can burn,
#'     the fit's target size (currently `"BURNABLE_HA"`).
#'   - `ignitionsTxt`: `character(1)`. Name of the ignitions column
#'     (currently `"ignitions"`).
#'   - `escapesTxt`: `character(1)`. Name of the escapes column
#'     (currently `"escapes"`).
#'   - `ignitionObjectiveTxt`: `character(1)`. The xgboost objective the ignition and escape models
#'     are fitted with (currently `"count:poisson"`); read it back from a fitted model with
#'     [xgbObjective()].
#'   - `spreadFitAdditionalColNamesTxt`: `character` vector of extra
#'     simList-slot/column names attached to spread-fit outputs
#'     (`"numIterations"`, `"objFunVal"`, `"params"`, `"sppEquiv"`,
#'     `"nonForestedLCCGroups"`, `"missingLCCgroup"`, `"covMinMax_spread"`). `covMinMax_spread` is
#'     what prediction needs to rescale covariates exactly as the fit did.
#'   - `ignitionFitAdditionalColNamesTxt`: `character` vector of the list-column names attached
#'     to ignition-fit ledger rows (`"fireSense_IgnitionFitted"`, `"fireSense_EscapeFitted"`).
#'   - `spreadInterceptTxt`: `character(1)`. Name of the spread model's intercept, in a formula's
#'     design (`spreadDesignCols()`), in `lower`/`upper` and in a ledger row's `params`
#'     (currently `"(Intercept)"`, what `stats::lm()` calls it).
#'   - `spreadFitCovCentreTxt`: `character(1)`. Name of the ledger column that holds the covariate
#'     centres of a fit made with an intercept (currently `"covCentre_spread"`). It is not in
#'     `spreadFitAdditionalColNamesTxt`: a fit without an intercept writes no such column.
#'   - `spreadLinkTxt`, `spreadFitClimateRefTxt`, `spreadFitPowerLawGammaTxt`: `character(1)` each. Names of the
#'     ledger columns a fit made with a non-default `spreadLink` writes: `"spreadLink"` (holding
#'     `"climateOutside"` or `"powerLaw"`), `"climateRef_spread"` (the fit's climate reference from
#'     [spreadClimateRef()], with `"climateOutside"`) and `"powerLawGamma_spread"` (the fit's `gamma`, with
#'     `"powerLaw"`), so prediction applies the same link ([spreadProbFromCovs()]). A fit with the default
#'     `"ceiling"` link writes none of them; a row without them, every row written before, predicts as it
#'     always did.
#'   - `spreadLinkOptions`: `character` vector of the values `fireSense_spreadFit`'s `spreadLink` can
#'     take: `"ceiling"` (the default: climate inside the logistic, under `maxAsymptote`),
#'     `"climateOutside"` and `"powerLaw"` ([powerLawLink()]).
#'   - `spreadPowerLawGamma`: `numeric(1)`. The default `gamma` of the `"powerLaw"` link, 1.9, from a percolation
#'     calibration on a uniform landscape. `fireSense_spreadFit` calibrates `gamma` on each ELF's own flammable map
#'     instead ([calibratePowerLaw()]), unless a number is given.
#'   - `spreadFitPowerLawPc50Txt`, `spreadFitPowerLawR2Txt`: `character(1)` each. Names of the ledger columns that
#'     hold, with `"powerLaw"` and a calibrated `gamma`, the calibration's `p50` (the spread probability at which
#'     half the fires run away, to compare with the fitted `p_c`) and the `R2` of its `gamma` regression
#'     (`"powerLawPc50_spread"`, `"powerLawR2_spread"`).
#'   - `spreadFitPowerLawPcMarginTxt`: `character(1)`. Name of the ledger column that holds the margin added to `p50`
#'     for the upper bound of `p_c` ([powerLawPcBounds()]; `"powerLawPcMargin_spread"`), so the bounds a fit used can be
#'     recovered.
#'   - `spreadPowerLawPGrid`: `numeric`. The spread probabilities [calibratePowerLaw()] simulates by default.
#'
#'   - `treedWetlandAgbTxt`: `character(1)`. Name of the pooled treed-wetland biomass column
#'     built by `fireSenseCovariatesCreate(fuelCovariates = "domSecWetland")` (currently
#'     `"treedWetland_agb"`).
#'
#' @name fireSenseUtils-constants
#' @aliases polygonIDTxt nonNFColNamesTxt yearTxt youngAgeTxt burnableHaTxt treedWetlandTxt treedWetlandAgbTxt ignitionsTxt escapesTxt ignitionObjectiveTxt spreadFitAdditionalColNamesTxt ignitionFitAdditionalColNamesTxt spreadInterceptTxt spreadFitCovCentreTxt spreadLinkTxt spreadFitClimateRefTxt spreadFitPowerLawGammaTxt spreadLinkOptions spreadPowerLawGamma spreadFitPowerLawPc50Txt spreadFitPowerLawR2Txt spreadFitPowerLawPcMarginTxt spreadPowerLawPGrid
NULL

#' @export
polygonIDTxt <- "polygonID"

#' @export
nonNFColNamesTxt <- c("pixelID", polygonIDTxt)

#' @export
yearTxt <- "year"

#' @export
youngAgeTxt <- "youngAge"

#' @export
burnableHaTxt <- "BURNABLE_HA"

#' @export
treedWetlandTxt <- "treedWetland"

#' @export
treedWetlandAgbTxt <- "treedWetland_agb"

#' @export
ignitionsTxt <- "ignitions"

#' @export
escapesTxt <- "escapes"

#' @export
ignitionObjectiveTxt <- "count:poisson"

#' @export
spreadFitAdditionalColNamesTxt <- c(
  "numIterations", "objFunVal", "params",
  "sppEquiv", "nonForestedLCCGroups", "missingLCCgroup",
  "covMinMax_spread"
)

#' @export
ignitionFitAdditionalColNamesTxt <- c("fireSense_IgnitionFitted", "fireSense_EscapeFitted")

#' @export
spreadInterceptTxt <- "(Intercept)"

#' @export
spreadFitCovCentreTxt <- "covCentre_spread"

#' @export
spreadLinkTxt <- "spreadLink"

#' @export
spreadFitClimateRefTxt <- "climateRef_spread"

#' @export
spreadFitPowerLawGammaTxt <- "powerLawGamma_spread"

#' @export
spreadLinkOptions <- c("ceiling", "climateOutside", "powerLaw")

#' @export
spreadPowerLawGamma <- 1.9

#' @export
spreadFitPowerLawPc50Txt <- "powerLawPc50_spread"

#' @export
spreadFitPowerLawR2Txt <- "powerLawR2_spread"

#' @export
spreadFitPowerLawPcMarginTxt <- "powerLawPcMargin_spread"

#' @export
spreadPowerLawPGrid <- sort(unique(round(c(seq(0.15, 0.21, by = 0.02), seq(0.22, 0.32, by = 0.005), seq(0.34, 0.38, by = 0.02)), 4)))

#' Memory (GB) assumed per DEoptim worker before a fit has a memory record
#'
#' `runDEoptim()` sets `options(clusters.workerMemoryGB)` to this, unless the user has set it, so
#' `clusters` (>= 0.0.75) caps workers per host by free memory from a fit's first cluster build.
#' Measured 2026-10-09 over 457 spread-fit workers on 15 hosts: peak resident memory median 4.5 GB,
#' 90th percentile about 6.5 GB, maximum 14.1 GB. Once a fit has run a chunk, its own record is used.
#' @export
spreadFitWorkerMemoryGB <- 14
