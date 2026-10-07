#' The folder a fit writes its figures and held-out results to
#'
#' A fit depends only on its polygon (ELF), its fire years and its model, and the name of its
#' ledger file (`spreadFitFilenameFor()`, `ignitionFitFilenameFor()`) already holds the last two.
#' So a fit's figures and held-out results go in one folder named for the polygon and that file,
#' next to the ledger file's own folder, and not in the output folder of whichever
#' scenario and replicate happened to run the fit first. Nothing is created.
#'
#' @param inputPath The simulation's `inputPath(sim)`, where the ledger files are.
#' @param studyAreaName The polygon the fit is for, e.g. the ELF id `"14.3"`.
#' @param fitFilename The ledger file the fit writes, e.g. `spreadFitFilenameFor(1985:2024)`.
#'
#' @return A path: `<inputPath>/fits/<studyAreaName>_<fitFilename without extension>`.
#' @export
#' @examples
#' fitOutputPath("inputs", "14.3", spreadFitFilenameFor(1985:2024))
fitOutputPath <- function(inputPath, studyAreaName, fitFilename) {
  file.path(inputPath, "fits",
            paste0(studyAreaName, "_", tools::file_path_sans_ext(basename(fitFilename))))
}
