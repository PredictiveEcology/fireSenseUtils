utils::globalVariables(c(
  "areaBurnedHa", "Nfires", "POLY_HA", "SIZE_HA", "sumAB", "sumBurn", "val", "var", "YEAR"
))

#' Plot historic ignitions, escapes, and area burned
#'
#' @template summary_plots
#' @param firePolys A `SpatialPolygonsDataFrame` of historic fire burn areas, from the Canadian
#'    National Fire Database.
#' @param ignitionPoints A `SpatialPointsDataFrame` of historic fire ignitions, from the Canadian
#'    National Fire Database.
#'
#' @export
#' @importFrom data.table as.data.table setnames
plotHistoricFires <- function(climateScenario, studyAreaName, outputDir, firePolys, ignitionPoints) {
  if (requireNamespace("ggplot2", quietly = TRUE) &&
      requireNamespace("SpaDES.core", quietly = TRUE)) {
    gcm <- strsplit(climateScenario, "_")[[1]][1]
    ssp <- strsplit(climateScenario, "_")[[1]][2]
    runName <- sprintf("%s_%s", studyAreaName, climateScenario) ## doesn't matter which run, all same
    run <- 1L
    sim <- SpaDES.core::loadSimList(file.path(outputDir, runName, "rep01", paste0(runName, "_rep01.qs")))
    burnSummary <- sim$burnSummary
    rm(sim)

    historicalBurns <- do.call(what = rbind, args = firePolys)
    historicalBurns <- as.data.table(historicalBurns@data)

    ## restrict to escapes only, but sum poly_ha for burns
    res <- ifelse(grepl("ROF", studyAreaName), 125, 250)
    historicalBurns <- historicalBurns[SIZE_HA > res,
                                       .(sumBurn = sum(as.numeric(POLY_HA)), nFires = .N), .(YEAR)]
    setnames(historicalBurns, "YEAR", "year")
    historicalBurns[, stat := "observed"]
    projectedEscapes <- burnSummary[areaBurnedHa > res, .(nFires = .N), .(year)]
    projectedBurns <- burnSummary[, .(sumBurn = sum(areaBurnedHa)), .(year)]
    projectedBurns <- projectedBurns[projectedEscapes, on = c("year")]
    projectedBurns[, stat := "projected"]
    dat <- rbind(projectedBurns, historicalBurns)

    trueHistoricalIgs <- as.data.table(ignitionPoints) %>%
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
      ggplot2::labs(y = "number of ignitions",
                    title = studyAreaName,
                    subtitle = paste(gcm, ssp))

    gEscapes <- ggplot2::ggplot(data = dat, ggplot2::aes(x = year, y = nFires, col = stat)) +
      ggplot2::geom_point() +
      # ggplot2::geom_smooth() +
      ggplot2::ylim(0, max(dat$nFires) * 1.2) +
      ggplot2::labs(y = "number of escaped fires",
                    title = studyAreaName,
                    subtitle = paste(gcm, ssp))

    gBurns <- ggplot2::ggplot(data = dat, ggplot2::aes(x = year, y = sumBurn, col = stat)) +
      ggplot2::geom_point() +
      # ggplot2::geom_smooth() +
      ggplot2::ylim(0, max(dat$sumBurn) * 1.1) +
      ggplot2::labs(y = "annual area burned (ha)",
                    title = paste(studyAreaName, "rep", run),
                    subtitle = paste(gcm, ssp))

    figDir <- file.path(outputDir, runName, "figures")
    figs <- list(
      ignition = file.path(figDir, paste0("simulated_Ignitions_", studyAreaName, "_", climateScenario, ".png")),
      escape = file.path(figDir, paste0("simulated_Escapes_", studyAreaName, "_", climateScenario, ".png")),
      spread = file.path(figDir, paste0("simulated_burnArea_", studyAreaName, "_", climateScenario, ".png"))
    )
    ggplot2::ggsave(plot = gIgnitions, filename = figs$ignition)
    ggplot2::ggsave(plot = gEscapes, filename = figs$escape)
    ggplot2::ggsave(plot = gBurns, filename = figs$spread)

    return(figs)
  }
}

#' Plot cumulative burn maps
#'
#' @template summary_plots
#' @template Nreps
#' @template rasterToMatch
#'
#' @return a file path corresponding to the images and/or objects written to disk
#'
#' @export
#' @importFrom parallel mclapply
#' @importFrom raster calc crop mask maxValue raster stack
plotCumulativeBurns <- function(studyAreaName, climateScenario, outputDir, Nreps, rasterToMatch) {
  if (requireNamespace("ggplot2", quietly = TRUE) &&
      requireNamespace("rasterVis", quietly = TRUE) &&
      requireNamespace("RColorBrewer", quietly = TRUE)) {
    burnMapAllReps <- parallel::mclapply(1:Nreps, function(rep) {
      runName <- sprintf("%s_%s", studyAreaName, climateScenario)
      resultsDir <- file.path(outputDir, runName, sprintf("rep%02d", rep))

      burnMap <- raster(file.path(resultsDir, "burnMap_2100_year2100.tif"))
    })

    cumulBurnMap <- calc(stack(burnMapAllReps), fun = sum) / Nreps
    cumulBurnMap <- mask(crop(cumulBurnMap, rasterToMatch), rasterToMatch)

    myPal <- RColorBrewer::brewer.pal("Reds", n = Nreps + 1) ## include 0 ## TODO: max 9 cols!
    myTheme <- rasterVis::rasterTheme(region = myPal)

    fburnMap <- file.path(outputDir, studyAreaName, "figures",
                          paste0("cumulBurnMap_", studyAreaName, "_", climateScenario, ".png"))

    fig <- rasterVis::levelplot(cumulBurnMap, margin = list(FUN = "mean"), ## median?
                                main = paste0("Cumulative burn map 2011-2100 under ", climateScenario),
                                colorkey = list(
                                  at = seq(0, maxValue(cumulBurnMap), length.out = Nreps + 1),
                                  space = "bottom",
                                  axis.line = list(col = "black"),
                                  width = 0.75
                                ),
                                par.settings = myTheme)

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
#' @template Nreps
#'
#' @export
#' @importFrom data.table data.table rbindlist
#' @importFrom parallel mclapply
#' @importFrom qs qread
#' @importFrom stats coefficients lm pf
plotBurnSummary <- function(studyAreaName, climateScenario, outputDir, Nreps) {
  if (requireNamespace("ggplot2", quietly = TRUE) &&
      requireNamespace("cowplot", quietly = TRUE)) {
    burnSummaryAllReps <- rbindlist(parallel::mclapply(1:Nreps, function(rep) {
      runName <- sprintf("%s_%s", studyAreaName, climateScenario)
      resultsDir <- file.path(outputDir, runName, sprintf("rep%02d", rep))

      burnDT <- qs::qread(file.path(resultsDir, "burnSummary_year2100.qs"))
      burnSummary <- data.table(year = burnDT[["year"]],
                                N = burnDT[["N"]],
                                areaBurnedHa = burnDT[["areaBurnedHa"]],
                                rep = as.integer(rep))
      burnSummary ## TODO: this is the BUFFERED studyArea, not the REPORTING one!!!!
    }))

    # totAreaBurned <- burnSummaryAllReps[, lapply(.SD, sum), by = c("year", "rep"), .SDcols = "areaBurnedHa"]
    # totAreaBurend <- totAreaBurned[, lapply(.SD, mean), by = "year", .SDcols = "areaBurnedHa"]

    burnSummaryAllReps[, sumAB := sum(areaBurnedHa), by = c("year", "rep")]
    areaBurned <- unique(burnSummaryAllReps[, c("year", "rep", "sumAB")])

    tend <- lm(sumAB ~ year, data = areaBurned)
    coeff <- coefficients(tend)
    Fstats <- summary(tend)$fstatistic
    names(Fstats) <- NULL
    pValueA <- ifelse(pf(Fstats[1], Fstats[2], Fstats[3], lower.tail = FALSE) < 0.01,
                      " \n(significant)", " \n(non-significant)")

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
                      " \n(significant)", " \n(non-significant)")
    nFires[, var := "number_fires"]
    nFires[, val := Nfires]

    # meanFireSize <- burnSummaryAllReps[, lapply(.SD, mean), by = c("year", "rep"), .SDcols = "areaBurnedHa"]
    # meanFireSize <- meanFireSize[, lapply(.SD, mean), by = "year", .SDcols = "areaBurnedHa"]

    burnSummaryAllReps[areaBurnedHa > 6.25, fireSize := mean(areaBurnedHa, na.rm = TRUE),
                       by = c("year", "rep")]
    fireSize <- na.omit(unique(burnSummaryAllReps[, c("year", "rep", "fireSize")]))

    tendS <- lm(fireSize ~ year, data = fireSize)
    coeffS <- coefficients(tendS)
    Fstats <- summary(tendS)$fstatistic
    names(Fstats) <- NULL
    pValueS <- ifelse(pf(Fstats[1], Fstats[2], Fstats[3], lower.tail = FALSE) < 0.01,
                      " \n(significant)", " \n(non-significant)")

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
      paste0("Area burned:\n",
             "y = ", ifelse(coefXA < 10000, coefXA, formatC(coefXA, format = "e", digits = 2)),
             "x + ", ifelse(coefYA < 10000, coefYA, formatC(coefYA, format = "e", digits = 2)), pValueA),
      paste0("No fires:\n",
             "y = ", ifelse(coefXF < 10000, coefXF, formatC(coefXF, format = "e", digits = 2)),
             "x + ", ifelse(coefYF < 10000, coefYF, formatC(coefYF, format = "e", digits = 2)), pValueF),
      paste0("Mean fire size:\n",
             "y = ", ifelse(coefXS < 10000, coefXS, formatC(coefXS, format = "e", digits = 2)),
             "x + ", ifelse(coefYS < 10000, coefYS, formatC(coefYS, format = "e", digits = 2)), pValueS)
    )
    names(replacementNames) <- c("area_burned", "number_fires", "fire_size")

    dt <- rbind(areaBurned, nFires, fireSize, use.names = FALSE)
    # Now remove original variable. It uses the first item's nameL sumAB
    dt[, sumAB := NULL]

    p1 <- ggplot2::ggplot(data = dt[var == "area_burned",], ggplot2::aes(x = year, y = val)) +
      ggplot2::geom_point(colour = "grey70") +
      ggplot2::stat_smooth(method = "lm", color = "darkred", fill = "red") +
      ggplot2::facet_grid(var ~ ., labeller = ggplot2::labeller(var = replacementNames)) +
      ggplot2::theme(legend.position = "none",
                     strip.text.y = ggplot2::element_text(size = 9, face = "bold"),
                     axis.title.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     plot.margin = ggplot2::unit(c(0.2, 0.2, -0.01, 0.2), "cm")) +
      ggplot2::labs(y = "total area burned (ha)")
    p2 <- ggplot2::ggplot(data = dt[var == "number_fires",], ggplot2::aes(x = year, y = val, colour = "blue")) +
      ggplot2::geom_point(colour = "grey70") +
      ggplot2::stat_smooth(method = "lm", fill = "blue", color = "darkblue") +
      ggplot2::facet_grid(var ~ ., labeller = ggplot2::labeller(var = replacementNames)) +
      ggplot2::theme(legend.position = "none",
                     strip.text.y = ggplot2::element_text(size = 9, face = "bold"),
                     plot.margin = ggplot2::unit(c(0.2, 0.2, -0.01, 0.2), "cm"),
                     axis.title.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank()) +
      ggplot2::ylab(label = "no. of fires")
    p3 <- ggplot2::ggplot(data = dt[var == "fire_size",], ggplot2::aes(x = year, y = val)) +
      ggplot2::geom_point(colour = "grey70") +
      ggplot2::stat_smooth(method = "lm", color = "orange", fill = "orange") +
      ggplot2::facet_grid(var ~ ., labeller = ggplot2::labeller(var = replacementNames)) +
      ggplot2::theme(legend.position = "none",
                     strip.text.y = ggplot2::element_text(size = 9, face = "bold"),
                     plot.margin = ggplot2::unit(c(-0.01, 0.2, 0.2, 0.2), "cm")) +
      ggplot2::labs(y = "mean fire size (ha)")

    title <- cowplot::ggdraw() +
      cowplot::draw_label(paste("Fires in the", studyAreaName, "study area under", climateScenario))

    p <- cowplot::plot_grid(p1, p2, p3, align = "h", nrow = 3, labels = "AUTO")

    fgg <- file.path(outputDir, studyAreaName, "figures",
                     paste0("burnSummary_", studyAreaName, "_", climateScenario, ".png"))
    gg <- cowplot::plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1))
    ggplot2::ggsave(gg, filename = fgg, height = 8, width = 11)
  }
}
