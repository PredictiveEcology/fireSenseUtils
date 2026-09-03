utils::globalVariables(c(
  "areaBurnedHa", "Nfires", "POLY_HA", "SIZE_HA", "sumAB", "sumBurn", "val", "var", "YEAR"
))

#' Locate per-replicate simulation output files
#'
#' @param simFiles character vector of simulation output file paths (e.g. the `file`
#'    column of `SpaDES.core::outputs(sim)`).
#'
#' @param filename character. The file name to locate within each replicate.
#'
#' @return character vector of file paths, named by replicate number, ordered by replicate.
#'
#' @noRd
.repOutputFiles <- function(simFiles, filename) {
  ## There is no sensible default. Reconstructing these paths from a `repNN`
  ## convention is what this function exists to avoid: it assumes zero-padded
  ## directories, which some projects do not use (`rep1`, not `rep01`), and
  ## silently yields paths that do not exist.
  if (is.null(simFiles) || length(simFiles) == 0L) {
    stop("`simFiles` is required: pass the output file paths, e.g. the `file` ",
         "column of SpaDES.core::outputs(sim).", call. = FALSE)
  }

  f <- simFiles[basename(simFiles) == filename]
  repNum <- suppressWarnings(as.integer(sub(".*[/\\\\]rep0*([0-9]+)[/\\\\].*", "\\1", f)))
  f <- f[!is.na(repNum)]
  repNum <- repNum[!is.na(repNum)]

  if (length(f) == 0L) {
    stop("no file named '", filename, "' found among `simFiles`.", call. = FALSE)
  }

  f <- f[order(repNum)]
  names(f) <- as.character(sort(repNum))
  f
}

#' Stop with the errors raised inside `mclapply()`
#'
#' `mclapply()` returns a `try-error` for each worker that failed, which otherwise
#' surfaces much later as an unrelated error (e.g. from `rbindlist()`).
#'
#' @param res list returned by `parallel::mclapply()`
#' @param ok logical vector: which elements are valid results
#' @param what character. What was being read, used in the error message.
#'
#' @return `NULL`, invisibly; called for the side effect of stopping on errors.
#'
#' @noRd
.stopOnMclapplyErrors <- function(res, ok, what) {
  if (!all(ok)) {
    stop("failed to read ", sum(!ok), " of ", length(ok), " ", what, ":\n",
         paste(vapply(res[!ok], function(x) paste(as.character(x), collapse = ""),
                      character(1)), collapse = ""))
  }
  invisible(NULL)
}

#' Extract the attribute table of a spatial object
#'
#' Handles `SpatVector`, `sf`, `Spatial*DataFrame`, and lists of any of these.
#'
#' @param x a spatial object, or a list of them
#'
#' @return a `data.frame` of the object's attributes
#'
#' @noRd
.attributes <- function(x) {
  if (is.list(x) && !inherits(x, c("sf", "data.frame"))) {
    return(do.call(rbind, lapply(x, .attributes)))
  }

  if (inherits(x, "SpatVector")) {
    terra::values(x)
  } else if (inherits(x, "sf")) {
    sf::st_drop_geometry(x)
  } else if (isTRUE(methods::.hasSlot(x, "data"))) {
    x@data
  } else {
    as.data.frame(x)
  }
}

#' Plot historic ignitions, escapes, and area burned
#'
#' @template summary_plots
#'
#' @param pixelSize raster pixel (Cell) size
#'
#' @param firePolys A `sf` spatial polygons of historic fire burn areas, from the Canadian
#'    National Fire Database.
#'
#' @param ignitionPoints A `sf` spatial points of historic fire ignitions, from the Canadian
#'    National Fire Database.
#'
#' @template simFiles
#'
#' @export
#' @importFrom data.table as.data.table setnames
plotHistoricFires <- function(climateScenario, studyAreaName, outputDir,
                              pixelSize, firePolys, ignitionPoints, simFiles) {
  if (requireNamespace("ggplot2", quietly = TRUE) &&
    requireNamespace("SpaDES.core", quietly = TRUE)) {
    ## A historical run has no scenario, and carries NA for it -- usually R's
    ## bare *logical* NA, which strsplit() rejects with "non-character
    ## argument". Split defensively and build the subtitle from whatever parts
    ## exist, so a historical run plots with no scenario rather than failing.
    scenLabel <- .scenarioLabel(climateScenario)
    scenParts <- .scenarioParts(climateScenario)
    gcm <- if (length(scenParts) >= 1L) scenParts[[1L]] else NA_character_
    ssp <- if (length(scenParts) >= 2L) scenParts[[2L]] else NA_character_
    scenSubtitle <- paste(scenParts, collapse = " ")

    historicalBurns <- as.data.table(.attributes(firePolys))

    ## restrict to escapes only, but sum poly_ha for burns
    historicalBurns <- historicalBurns[
      SIZE_HA > pixelSize,
      .(sumBurn = sum(as.numeric(POLY_HA)), nFires = .N), .(YEAR)
    ]
    setnames(historicalBurns, "YEAR", "year")
    historicalBurns[, stat := "observed"]

    ## only use the first rep
    f <- .repOutputFiles(simFiles, filename = "fireSense_burnSummary.csv")[1L]
    rep <- as.integer(names(f))

    burnDT <- data.table::fread(f)
    burnSummary <- data.table(
      year = burnDT[["year"]],
      N = burnDT[["N"]],
      areaBurnedHa = burnDT[["areaBurnedHa"]],
      rep = as.integer(rep)
    ) ## TODO: this is the BUFFERED studyArea, not the REPORTING one!!!!

    ## TODO: average by number of reps
    projectedEscapes <- burnSummary[areaBurnedHa > pixelSize, .(nFires = .N), .(year)]
    projectedBurns <- burnSummary[, .(sumBurn = sum(areaBurnedHa)), .(year)]
    projectedBurns <- projectedBurns[projectedEscapes, on = c("year")]
    projectedBurns[, stat := "projected"]
    dat <- rbind(projectedBurns, historicalBurns)

    trueHistoricalIgs <- as.data.table(.attributes(ignitionPoints)) %>%
      .[, .N, .(YEAR)] %>%
      setnames(., "YEAR", "year") %>%
      .[, stat := "observed"] %>%
      .[, year := as.numeric(year)]
    projectedIgs <- burnSummary[, .N, .(year)] %>%
      .[, stat := "projected"]
    dat2 <- rbind(trueHistoricalIgs, projectedIgs)

    gIgnitions <- ggplot2::ggplot(data = dat2, ggplot2::aes(x = year, y = N, col = stat)) +
      ggplot2::geom_point() +
      # ggplot2::geom_smooth() +
      ggplot2::ylim(0, max(dat2$N) * 1.2) +
      ggplot2::labs(
        y = "number of ignitions",
        title = studyAreaName,
        subtitle = scenSubtitle
      )

    gEscapes <- ggplot2::ggplot(data = dat, ggplot2::aes(x = year, y = nFires, col = stat)) +
      ggplot2::geom_point() +
      # ggplot2::geom_smooth() +
      ggplot2::ylim(0, max(dat$nFires) * 1.2) +
      ggplot2::labs(
        y = "number of escaped fires",
        title = studyAreaName,
        subtitle = scenSubtitle
      )

    gBurns <- ggplot2::ggplot(data = dat, ggplot2::aes(x = year, y = sumBurn, col = stat)) +
      ggplot2::geom_point() +
      # ggplot2::geom_smooth() +
      ggplot2::ylim(0, max(dat$sumBurn) * 1.1) +
      ggplot2::labs(
        y = "annual area burned (ha)",
        title = paste(studyAreaName, "rep", rep),
        subtitle = scenSubtitle
      )

    figs <- list(
      ignition = .figFile("simulated_Ignitions", outputDir, studyAreaName, scenLabel),
      escape   = .figFile("simulated_Escapes", outputDir, studyAreaName, scenLabel),
      spread   = .figFile("simulated_burnArea", outputDir, studyAreaName, scenLabel)
    )
    ggplot2::ggsave(plot = gIgnitions, filename = figs$ignition)
    ggplot2::ggsave(plot = gEscapes, filename = figs$escape)
    ggplot2::ggsave(plot = gBurns, filename = figs$spread)

    return(unlist(figs))
  }
}

#' Plot cumulative burn maps
#'
#' @template summary_plots
#'
#' @template Nreps
#'
#' @param years integer (length 2). start and end years for comparison.
#'
#' @template rasterToMatch
#'
#' @template simFiles
#'
#' @return a file path corresponding to the images and/or objects written to disk
#'
#' @export
#' @importFrom parallel mclapply
plotCumulativeBurns <- function(climateScenario, studyAreaName, outputDir,
                                Nreps, years, rasterToMatch, simFiles) {
  if (requireNamespace("ggplot2", quietly = TRUE) &&
      requireNamespace("raster", quietly = TRUE) &&
      requireNamespace("rasterVis", quietly = TRUE) &&
      requireNamespace("RColorBrewer", quietly = TRUE)) {
    scenLabel <- .scenarioLabel(climateScenario)
    fs <- .repOutputFiles(simFiles, filename = paste0("burnMap_year", years[2], ".tif"))

    burnMapAllReps <- parallel::mclapply(fs, function(f) raster::raster(f))
    .stopOnMclapplyErrors(burnMapAllReps,
                          vapply(burnMapAllReps, inherits, logical(1), "BasicRaster"),
                          "cumulative burn maps")

    ## `simFiles` determines which reps are actually available
    nreps <- length(burnMapAllReps)

    cumulBurnMap <- raster::calc(raster::stack(burnMapAllReps), fun = sum) / nreps

    ## TODO: drop this once the function is switched to terra/tidyterra
    if (inherits(rasterToMatch, "SpatRaster")) {
      rasterToMatch <- raster::raster(rasterToMatch)
    }
    cumulBurnMap <- raster::mask(raster::crop(cumulBurnMap, rasterToMatch), rasterToMatch)

    myPal <- RColorBrewer::brewer.pal("Reds", n = nreps + 1) ## include 0 ## TODO: max 9 cols!
    myTheme <- rasterVis::rasterTheme(region = myPal)

    fburnMap <- .figFile("cumulBurnMap", outputDir, studyAreaName, scenLabel)

    fig <- rasterVis::levelplot(cumulBurnMap,
      margin = list(FUN = "mean"), ## median?
      main = paste0("Cumulative burn map ", years[1], "-", years[2], " under ", scenLabel),
      colorkey = list(
        at = seq(0, raster::maxValue(cumulBurnMap), length.out = nreps + 1),
        space = "bottom",
        axis.line = list(col = "black"),
        width = 0.75
      ),
      par.settings = myTheme
    )

    ## levelplot (trellis graphics more generally) won't plot correctly inside loop w/o print()
    png(filename = fburnMap, height = 1000, width = 1000)
    print(fig)
    dev.off()

    return(fburnMap)
  }
}

#' Plot burn summary
#'
#' Create plot with subplots showing: a) area burned; b) number of fires; c) mean fire size.
#'
#' @template summary_plots
#'
#' @template Nreps
#'
#' @param years integer (length 2). start and end years for comparison.
#'
#' @param pixelSize raster pixel (Cell) size
#'
#' @template simFiles
#'
#' @export
#' @importFrom data.table data.table rbindlist
#' @importFrom parallel mclapply
#' @importFrom stats coefficients lm pf
plotBurnSummary <- function(climateScenario, studyAreaName, outputDir,
                            Nreps, years, pixelSize, simFiles) {
  if (requireNamespace("ggplot2", quietly = TRUE) &&
    requireNamespace("cowplot", quietly = TRUE)) {
    scenLabel <- .scenarioLabel(climateScenario)
    fs <- .repOutputFiles(simFiles, filename = "fireSense_burnSummary.csv")

    burnSummaryPerRep <- parallel::mclapply(names(fs), function(rep) {
      burnDT <- data.table::fread(fs[[rep]])
      burnSummary <- data.table(
        year = burnDT[["year"]],
        N = burnDT[["N"]],
        areaBurnedHa = burnDT[["areaBurnedHa"]],
        rep = as.integer(rep)
      )
      burnSummary ## TODO: this is the BUFFERED studyArea, not the REPORTING one!!!!
    })
    .stopOnMclapplyErrors(burnSummaryPerRep,
                          vapply(burnSummaryPerRep, is.data.frame, logical(1)),
                          "burn summaries")

    burnSummaryAllReps <- rbindlist(burnSummaryPerRep)

    # totAreaBurned <- burnSummaryAllReps[, lapply(.SD, sum), by = c("year", "rep"), .SDcols = "areaBurnedHa"]
    # totAreaBurend <- totAreaBurned[, lapply(.SD, mean), by = "year", .SDcols = "areaBurnedHa"]

    burnSummaryAllReps[, sumAB := sum(areaBurnedHa), by = c("year", "rep")]
    areaBurned <- unique(burnSummaryAllReps[, c("year", "rep", "sumAB")])

    tend <- lm(sumAB ~ year, data = areaBurned)
    coeff <- coefficients(tend)
    Fstats <- summary(tend)$fstatistic
    names(Fstats) <- NULL
    pValueA <- ifelse(pf(Fstats[1], Fstats[2], Fstats[3], lower.tail = FALSE) < 0.01,
      " \n(significant)", " \n(non-significant)"
    )

    areaBurned[, var := "area_burned"]
    areaBurned[, val := sumAB]

    # numberFires <- burnSummaryAllReps[, lapply(.SD, length), by = c("year", "rep"), .SDcols = "N"]
    # numberFires <- numberFires[, lapply(.SD, mean), by = "year", .SDcols = "N"]

    burnSummaryAllReps[, Nfires := length(N), by = c("year", "rep")]
    nFires <- unique(burnSummaryAllReps[, c("year", "rep", "Nfires")])

    tendF <- lm(Nfires ~ year, data = nFires)
    coeffF <- coefficients(tendF)
    Fstats <- summary(tendF)$fstatistic
    names(Fstats) <- NULL
    pValueF <- ifelse(pf(Fstats[1], Fstats[2], Fstats[3], lower.tail = FALSE) < 0.01,
      " \n(significant)", " \n(non-significant)"
    )
    nFires[, var := "number_fires"]
    nFires[, val := Nfires]

    # meanFireSize <- burnSummaryAllReps[, lapply(.SD, mean), by = c("year", "rep"), .SDcols = "areaBurnedHa"]
    # meanFireSize <- meanFireSize[, lapply(.SD, mean), by = "year", .SDcols = "areaBurnedHa"]

    pixelSizeHa <- (pixelSize^2) / 1e4
    burnSummaryAllReps[areaBurnedHa > pixelSizeHa, fireSize := mean(areaBurnedHa, na.rm = TRUE),
      by = c("year", "rep")
    ]
    fireSize <- na.omit(unique(burnSummaryAllReps[, c("year", "rep", "fireSize")]))

    tendS <- lm(fireSize ~ year, data = fireSize)
    coeffS <- coefficients(tendS)
    Fstats <- summary(tendS)$fstatistic
    names(Fstats) <- NULL
    pValueS <- ifelse(pf(Fstats[1], Fstats[2], Fstats[3], lower.tail = FALSE) < 0.01,
      " \n(significant)", " \n(non-significant)"
    )

    fireSize[, var := "fire_size"]
    fireSize[, val := fireSize]

    ## plotting
    coefXA <- round(coeff[2], 1)
    coefYA <- round(coeff[1], 1)
    coefXF <- round(coeffF[2], 1)
    coefYF <- round(coeffF[1], 1)
    coefXS <- round(coeffS[2], 1)
    coefYS <- round(coeffS[1], 1)

    replacementNames <- c(
      paste0(
        "Area burned:\n",
        "y = ", ifelse(coefXA < 10000, coefXA, formatC(coefXA, format = "e", digits = 2)),
        "x + ", ifelse(coefYA < 10000, coefYA, formatC(coefYA, format = "e", digits = 2)), pValueA
      ),
      paste0(
        "No fires:\n",
        "y = ", ifelse(coefXF < 10000, coefXF, formatC(coefXF, format = "e", digits = 2)),
        "x + ", ifelse(coefYF < 10000, coefYF, formatC(coefYF, format = "e", digits = 2)), pValueF
      ),
      paste0(
        "Mean fire size:\n",
        "y = ", ifelse(coefXS < 10000, coefXS, formatC(coefXS, format = "e", digits = 2)),
        "x + ", ifelse(coefYS < 10000, coefYS, formatC(coefYS, format = "e", digits = 2)), pValueS
      )
    )
    names(replacementNames) <- c("area_burned", "number_fires", "fire_size")

    dt <- rbind(areaBurned, nFires, fireSize, use.names = FALSE)
    # Now remove original variable. It uses the first item's nameL sumAB
    dt[, sumAB := NULL]

    p1 <- ggplot2::ggplot(data = dt[var == "area_burned", ], ggplot2::aes(x = year, y = val)) +
      ggplot2::geom_point(colour = "grey70") +
      ggplot2::stat_smooth(method = "lm", color = "darkred", fill = "red") +
      ggplot2::facet_grid(var ~ ., labeller = ggplot2::labeller(var = replacementNames)) +
      ggplot2::theme(
        legend.position = "none",
        strip.text.y = ggplot2::element_text(size = 9, face = "bold"),
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        plot.margin = ggplot2::unit(c(0.2, 0.2, -0.01, 0.2), "cm")
      ) +
      ggplot2::labs(y = "total area burned (ha)")
    p2 <- ggplot2::ggplot(data = dt[var == "number_fires", ], ggplot2::aes(x = year, y = val, colour = "blue")) +
      ggplot2::geom_point(colour = "grey70") +
      ggplot2::stat_smooth(method = "lm", fill = "blue", color = "darkblue") +
      ggplot2::facet_grid(var ~ ., labeller = ggplot2::labeller(var = replacementNames)) +
      ggplot2::theme(
        legend.position = "none",
        strip.text.y = ggplot2::element_text(size = 9, face = "bold"),
        plot.margin = ggplot2::unit(c(0.2, 0.2, -0.01, 0.2), "cm"),
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      ) +
      ggplot2::ylab(label = "no. of fires")
    p3 <- ggplot2::ggplot(data = dt[var == "fire_size", ], ggplot2::aes(x = year, y = val)) +
      ggplot2::geom_point(colour = "grey70") +
      ggplot2::stat_smooth(method = "lm", color = "orange", fill = "orange") +
      ggplot2::facet_grid(var ~ ., labeller = ggplot2::labeller(var = replacementNames)) +
      ggplot2::theme(
        legend.position = "none",
        strip.text.y = ggplot2::element_text(size = 9, face = "bold"),
        plot.margin = ggplot2::unit(c(-0.01, 0.2, 0.2, 0.2), "cm")
      ) +
      ggplot2::labs(y = "mean fire size (ha)")

    title <- cowplot::ggdraw() +
      cowplot::draw_label(paste("Fires in the", studyAreaName, "study area under",
                                scenLabel))

    p <- cowplot::plot_grid(p1, p2, p3, align = "h", nrow = 3, labels = "AUTO")

    fgg <- .figFile("burnSummary", outputDir, studyAreaName, scenLabel)
    gg <- cowplot::plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1))
    ggplot2::ggsave(gg, filename = fgg, height = 8, width = 11)

    return(fgg)
  }
}

## Label for a climate scenario, for use in titles and filenames.
##
## A run with no climate scenario carries NA -- typically R's bare *logical*
## NA. Pasting that into a title or filename gives literal "NA" text (e.g.
## "burnSummary_someArea_NA.png"), so it is rendered as "NRV" instead, and said
## out loud: an unset scenario and a deliberate NRV run look identical here, so
## the caller should be told which one we assumed.
##
## Anything that cannot be a scenario (NULL, length 0, length > 1) is an error,
## not an unset value -- relabelling those would hide a bug upstream behind a
## plausible-looking filename.
##
## Call once per function and reuse: it is used in several filenames and titles
## and should not repeat its message.
## @noRd
.scenarioLabel <- function(climateScenario) {
  if (is.null(climateScenario) || length(climateScenario) == 0L) {
    stop("`climateScenario` is missing; pass NA for a run with no scenario.",
         call. = FALSE)
  }
  if (length(climateScenario) != 1L) {
    stop("`climateScenario` must be length 1, got ", length(climateScenario), ".",
         call. = FALSE)
  }
  if (is.na(climateScenario) || !nzchar(as.character(climateScenario))) {
    message("climateScenario is NA, so declaring it to be NRV; ",
            "if this is incorrect, pass the name.")
    return("NRV")
  }
  as.character(climateScenario)
}

## Underscore-separated parts of a climate scenario ("CNRM-ESM2-1_ssp370" ->
## c("CNRM-ESM2-1", "ssp370")); character() for a historical run. Unlike a bare
## strsplit(), never errors when climateScenario is a logical NA.
## @noRd
.scenarioParts <- function(climateScenario) {
  if (length(climateScenario) != 1L || is.na(climateScenario)) return(character())
  parts <- strsplit(as.character(climateScenario), "_")[[1L]]
  parts[nzchar(parts)]
}

## Path of a figure file, and the directory to hold it.
##
## All three plot functions wrote to the same place under the same naming
## scheme -- <outputDir>/<studyAreaName>/figures/<prefix>_<studyAreaName>_
## <scenario>.png -- each building the path and creating the directory itself.
## One definition means the layout can change in one place, and no caller can
## forget the dir.create().
##
## @param prefix Figure kind, e.g. "burnSummary".
## @param create Create the containing directory. TRUE except when the caller
##   only wants the path.
## @noRd
.figFile <- function(prefix, outputDir, studyAreaName, scenarioLabel, create = TRUE) {
  f <- file.path(
    outputDir, studyAreaName, "figures",
    paste0(prefix, "_", studyAreaName, "_", scenarioLabel, ".png")
  )
  if (isTRUE(create)) {
    dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  }
  f
}
