#' @param simFiles Optional character vector of simulation output file paths
#'    (e.g. the `file` column of `SpaDES.core::outputs(sim)`), used to locate the
#'    per-replicate files instead of reconstructing them from `outputDir`.
#'    Needed when the replicate outputs do not live in `<outputDir>/repNN/`, e.g.
#'    when summarizing runs restored from another machine, or when the replicate
#'    directories are unpadded (`rep1` rather than `rep01`).
#'    When supplied, the replicates present in `simFiles` are used, and `Nreps`
#'    is ignored.
