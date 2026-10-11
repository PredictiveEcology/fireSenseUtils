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
#'   - `ignitionsTxt`: `character(1)`. Name of the ignitions column
#'     (currently `"ignitions"`).
#'   - `escapesTxt`: `character(1)`. Name of the escapes column
#'     (currently `"escapes"`).
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
#'
#'   - `treedWetlandAgbTxt`: `character(1)`. Name of the pooled treed-wetland biomass column
#'     built by `fireSenseCovariatesCreate(fuelCovariates = "domSecWetland")` (currently
#'     `"treedWetland_agb"`).
#'
#' @name fireSenseUtils-constants
#' @aliases polygonIDTxt nonNFColNamesTxt yearTxt youngAgeTxt treedWetlandTxt treedWetlandAgbTxt ignitionsTxt escapesTxt spreadFitAdditionalColNamesTxt ignitionFitAdditionalColNamesTxt spreadInterceptTxt spreadFitCovCentreTxt
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
treedWetlandTxt <- "treedWetland"

#' @export
treedWetlandAgbTxt <- "treedWetland_agb"

#' @export
ignitionsTxt <- "ignitions"

#' @export
escapesTxt <- "escapes"

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

#' Memory (GB) assumed per DEoptim worker before a fit has a memory record
#'
#' `runDEoptim()` sets `options(clusters.workerMemoryGB)` to this, unless the user has set it, so
#' `clusters` (>= 0.0.75) caps workers per host by free memory from a fit's first cluster build.
#' Measured 2026-10-09 over 457 spread-fit workers on 15 hosts: peak resident memory median 4.5 GB,
#' 90th percentile about 6.5 GB, maximum 14.1 GB. Once a fit has run a chunk, its own record is used.
#' @export
spreadFitWorkerMemoryGB <- 14
