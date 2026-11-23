#' Cause groups and the NFDB / NBAC values in them
#'
#' The one place that says which `CAUSE` (NFDB points) and `FIRECAUS` (NBAC polygons) values are which cause, so
#' that fitting (`fireSense_dataPrepFit`) and validation ([hindcastObserved()], [ELFfireCounts()]) agree.
#' `natural` is NFDB `N` (and `L`, in older releases) and NBAC `Natural`. `unknown` is NFDB `U` and NBAC
#' `Undetermined` (NBAC 20260513; `Unknown` in case a release spells it so), and a blank or missing cause, as
#' neither source says what started the fire. `human` is NFDB `H` (and `H-PB`, prescribed burn) and NBAC `Human`;
#' `other` is NFDB `RE` (re-burn). NFDB 20260811 has only H, N and U; NBAC 20260513 only Human, Natural and
#' Undetermined.
#'
#' @param groups Cause groups: any of `natural`, `unknown`, `human`, `other`.
#' @param source `"nfdb"` (values of `CAUSE`) or `"nbac"` (values of `FIRECAUS`).
#' @return The values in `groups`; `NA` is among those of `unknown`.
#' @export
fireCauseCodes <- function(groups, source = c("nfdb", "nbac")) {
  source <- match.arg(source)
  codes <- list(
    natural = list(nfdb = c("L", "N"),    nbac = "Natural"),
    unknown = list(nfdb = c("U", "", NA), nbac = c("Undetermined", "Unknown", "", NA)),
    human   = list(nfdb = c("H", "H-PB"), nbac = "Human"),
    other   = list(nfdb = "RE",           nbac = character()))
  bad <- setdiff(groups, names(codes))
  if (length(bad)) stop("unknown fire cause group: ", paste(bad, collapse = ", "),
                        " (use ", paste(names(codes), collapse = ", "), ")")
  unique(unlist(lapply(codes[groups], `[[`, source), use.names = FALSE))
}

## TRUE for each row of the data.frame `d` that is a prescribed burn, by any flag the records carry: NBAC `PRESCRIBED`
## ("true" or missing), NFDB `PRESCRIBED` ("PB" or "No"; "yes" in the polygons), NFDB `FIRE_TYPE` ("PB", "Prescribed
## Burn") and NFDB `CAUSE` "H-PB". In NFDB 20260811, 38 natural and unknown fires are flagged only by `PRESCRIBED`.
prescribedFires <- function(d) {
  yes <- c("true", "t", "y", "yes", "1", "pb")
  prescribed <- rep(FALSE, nrow(d))
  if ("PRESCRIBED" %in% names(d)) prescribed <- prescribed | tolower(as.character(d$PRESCRIBED)) %in% yes
  if ("FIRE_TYPE" %in% names(d)) prescribed <- prescribed | d$FIRE_TYPE %in% c("PB", "Prescribed Burn")
  if ("CAUSE" %in% names(d)) prescribed <- prescribed | d$CAUSE %in% "H-PB"
  prescribed
}

#' The fires the models are fitted to and scored on: of the modelled causes, and not prescribed burns
#'
#' The one function that selects historic fire polygons (NBAC, NFDB) or points (NFDB) by cause. It drops every
#' prescribed burn (see Details), whatever its cause, because no fit, hindcast, score or observed record uses a
#' prescribed fire. With `fireCauses` (or `codes`) it also keeps only the rows whose cause is one of them.
#'
#' A row is a prescribed burn if its `PRESCRIBED` is `"true"`, `"yes"`, `"PB"` or `"1"` (NBAC: `"true"` or missing;
#' NFDB points: `"PB"` or `"No"`; NFDB polygons: `"yes"`), its `FIRE_TYPE` is `"PB"` or `"Prescribed Burn"`, or its
#' `CAUSE` is `"H-PB"`. A column the records lack is not looked at.
#'
#' @param x `SpatVector`, `sf` or `data.frame` of fire polygons or points; `NULL` is returned as `NULL`.
#' @param fireCauses `NULL` (default): no cut on cause, only prescribed burns are dropped. Cause groups, see
#'   [fireCauseCodes()] (e.g. [fireSenseFireCauses]).
#' @param source `"nbac"` (cause column [polyCauseColTxt]) or `"nfdb"` (`CAUSE`), as in [fireCauseCodes()].
#' @param causeCol The cause column; default by `source`.
#' @param codes The cause values to keep; default those of `fireCauses` in `source`. Give values of `causeCol` to
#'   select by them directly (e.g. `naturalPolyCauses` of [hindcastObserved()]); `NULL` keeps every cause.
#' @return `x` without the other rows. A missing cause column is an error when a cut on cause is asked for.
#' @export
modelledFires <- function(x, fireCauses = NULL, source = c("nbac", "nfdb"), causeCol = NULL,
                          codes = if (!is.null(fireCauses)) fireCauseCodes(fireCauses, source)) {
  if (is.null(x)) return(NULL)
  source <- match.arg(source)
  d <- as.data.frame(x)
  keep <- !prescribedFires(d)
  if (!is.null(codes)) {
    if (is.null(causeCol)) causeCol <- if (source == "nbac") polyCauseColTxt else "CAUSE"
    if (!causeCol %in% names(d)) stop("fire records have no '", causeCol, "' column: cannot select by cause")
    keep <- keep & d[[causeCol]] %in% codes
  }
  x[which(keep), ]
}
