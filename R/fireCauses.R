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
