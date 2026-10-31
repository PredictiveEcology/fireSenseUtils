#' The phases a fireSense job can be split into
#'
#' `"all"` fits if the ELF has no fit, then predicts. `"inputs"` builds every input and cache and
#' stops before the spread threshold calibration (which forks). `"fit"` fits if the ELF has no fit,
#' then stops. See [phaseEvents()] and [runJob()].
#'
#' @export
fireSensePhases <- c("all", "inputs", "fit")

#' The `SpaDES.core` event barriers of a phase
#'
#' @param phase One of [fireSensePhases].
#' @param validates Whether the spread fit validates (its `mode` contains `"validate"`, or a
#'   `heldOutFold` is set): the `"fit"` phase then stops after `crossValidate`, not after `run`.
#'
#' @return `NULL` for `"all"`; otherwise the `events` list (`.stopBefore`, `.stopAfter`) that
#'   `SpaDES.core::simInit()` takes.
#'
#' @export
phaseEvents <- function(phase = fireSensePhases, validates = FALSE) {
  phase <- match.arg(phase)
  stopBeforeClimate <- list(fireSense_dataPrepPredict = "getClimateRasters")
  switch(phase,
         all = NULL,
         inputs = list(.stopBefore = c(list(fireSense_spreadFit = "estimateThreshold"), stopBeforeClimate)),
         fit = list(.stopAfter = list(fireSense_spreadFit = if (validates) "crossValidate" else "run"),
                    .stopBefore = stopBeforeClimate))
}

#' Run one fireSense job of an ELF
#'
#' Runs the simulation that `SpaDES.project::setupProject()` set up, up to the end of `phase`
#' ([phaseEvents()]), and returns it. The `sim` is copied first, so it is not changed.
#'
#' With `fitQueue`, a `"inputs"` job is skipped when `SpaDES.project::queueSkipForFits()` gives a
#' reason (the fits have claimed or will claim this ELF). With `afterQueue`, a `"fit"` job first calls
#' `SpaDES.project::queueWaitWhile()`, so it does not start while that queue is still building the ELF's
#' inputs. Both need `SpaDES.project`.
#'
#' @param sim The result of `SpaDES.project::setupProject()`.
#' @param elf The ELF's name, as the queues know it.
#' @param phase One of [fireSensePhases]; by default the `phase` parameter of `fireSense_spreadFit`,
#'   or `"all"` if it has none.
#' @param fitQueue The fit queue's `.rds` file, or `NA`.
#' @param afterQueue The inputs queue's `.rds` file, or `NA`.
#' @param diagFile A file to [source()] before the run, or `NA`.
#' @param rprofPrefix If not `NA`, the run is profiled with [Rprof()] to
#'   `<rprofPrefix>_<elf>_<pid>.out`.
#' @param uploadDir The Google Drive folder that a finished `"all"` job saves its simList to with
#'   `SpaDES.project::outSaveTarUpload()`, or `NULL` for no upload. A failed upload only warns.
#' @param hindcast Whether the job is a hindcast, which is not uploaded.
#'
#' @return The simList, invisibly; `NULL` if the job was skipped.
#'
#' @export
runJob <- function(sim, elf, phase = sim$params$fireSense_spreadFit$phase, fitQueue = NA,
                   afterQueue = NA, diagFile = NA, rprofPrefix = NA, uploadDir = NULL, hindcast = FALSE) {
  phase <- match.arg(phase, fireSensePhases)
  if (!is.na(fitQueue) || !is.na(afterQueue) || !is.null(uploadDir))
    if (!requireNamespace("SpaDES.project", quietly = TRUE))
      stop("Please install PredictiveEcology/SpaDES.project to use `fitQueue`, `afterQueue` or `uploadDir`")

  validates <- "validate" %in% sim$params$fireSense_spreadFit$mode ||
    any(!is.na(sim$params$.globals$heldOutFold))
  copy <- reproducible::Copy(sim)
  copy$events <- phaseEvents(phase, validates)

  skip <- NULL
  if (phase == "inputs" && !is.na(fitQueue))
    skip <- SpaDES.project::queueSkipForFits(fitQueue, elf)
  if (!is.null(skip)) {
    message("phase 'inputs' of ", elf, " skipped: ", skip)
    return(invisible(NULL))
  }
  if (phase == "fit" && !is.na(afterQueue))
    SpaDES.project::queueWaitWhile(afterQueue, elf)

  if (!is.na(diagFile)) source(diagFile)
  if (!is.na(rprofPrefix)) {
    utils::Rprof(paste0(rprofPrefix, "_", elf, "_", Sys.getpid(), ".out"), interval = 0.5, line.profiling = TRUE)
    on.exit(utils::Rprof(NULL))
  }

  simOut <- SpaDES.core::simInitAndSpades2(copy)

  if (phase == "all" && !hindcast && !is.null(uploadDir))
    tryCatch(SpaDES.project::outSaveTarUpload(runName = sim$runName, sim = simOut, gFolder = uploadDir),
             error = function(e) warning("outSaveTarUpload failed: ", conditionMessage(e), call. = FALSE))
  invisible(simOut)
}
