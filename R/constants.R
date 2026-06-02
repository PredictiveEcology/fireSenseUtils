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
#'     `"nonForestedLCCGroups"`, `"missingLCCgroup"`).
#'
#' @name fireSenseUtils-constants
#' @aliases polygonIDTxt nonNFColNamesTxt yearTxt youngAgeTxt ignitionsTxt escapesTxt spreadFitAdditionalColNamesTxt
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
ignitionsTxt <- "ignitions"

#' @export
escapesTxt <- "escapes"

#' @export
spreadFitAdditionalColNamesTxt <- c(
  "numIterations", "objFunVal", "params",
  "sppEquiv", "nonForestedLCCGroups", "missingLCCgroup"
)
