utils::globalVariables(c(
  "burned", "burnprob", "cell", "lcc", "pixelIndex", "speciesCode", "YEAR", "yearRange",
  "absCoef", "possGroups", "propB", "pvalue", "sig", "sim", "species", "above10PctRelB",
  "assignedFuelClass", "B_MgHa", "coef"
))

#' Calculate proportional burn of landcover and tree species
#'
#' @template pixelGroupMap
#'
#' @template cohortData
#'
#' @param rstLCC a landcover map
#'
#' @param nonflammableLCC nonflammable landcover in `rstLCC`
#'
#' @param fires a single `sf` or `SpatVector` object of fire polygons containing a `YEAR` column
#'
#' @param yearRange  the range of years represented by this landscape
#'
#' @return `data.table` with cell, biomass, tree species or lcc for non-forest, and year of fire
#'
#' @examples
#' # fuelClassPrep(
#' #   sim$pixelGroupMap2011, sim$cohortData2011, sim$rstLCC2011,
#' #   nonflammableLCC = P(sim)$nonflammableLCC,
#' #   fires = TODO,
#' #   yearRange = c(2012, 2020)
#' # )
#'
#' @export
#' @importFrom data.table as.data.table
#'
fuelClassPrep <- function(pixelGroupMap, cohortData, rstLCC, nonflammableLCC, fires, yearRange) {
  ## Filter fires by year range
  firesSubset <- fires[fires$YEAR >= yearRange[1] & fires$YEAR <= yearRange[2], ] |>
    rasterize(pixelGroupMap, field = "YEAR", fun = "min") |>
    as.data.frame(cells = TRUE) |> as.data.table()

  ## Add pixels to cohort data
  speciesData <- addPixels2CohortData(cohortData = cohortData, pixelGroupMap)

  ## Prepare landscape data excluding non-flammable land cover classes
  landscapeData <- as.data.table(as.data.frame(rstLCC, cells = TRUE)) |>
    setnames(c("cell", "lcc"))
  landscapeData <- landscapeData[!lcc %in% nonflammableLCC]

  ## Merge species data with landscape data
  landscapeData <- speciesData[landscapeData, on = c("pixelIndex" = "cell")]

  ## Merge fires data with landscape data
  landscapeData <- firesSubset[landscapeData, on = c("cell" = "pixelIndex")]
  landscapeData[, burned := ifelse(is.na(YEAR), 0, 1)]
  landscapeData[is.na(speciesCode), speciesCode := paste0("LCC", lcc)]
  landscapeData[, year := yearRange[1]]

  return(landscapeData)
}


#' wrapper on GLM
#'
#' @param species tree species in landscape$speciesCode
#' @param landscape created by `fuelClassPrep`
#' @importFrom stats glm binomial
makeGLM <- function(species, landscape) {
  landscape <- landscape[speciesCode == species,]
  out <- glm(data = landscape, formula = burned ~ B_MgHa, family = binomial(link = "logit"))
}

#' Semi-automated selection of fuel classed based on GLMs
#'
#' @param landscape data.table created by `prepFuelClasses`
#' @param fuelCol the column in `sppEquiv` with default fuel classes
#' @param rstLCC a landcover map
#' @template sppEquiv
#' @template sppEquivCol
#' @param yearRange  the range of years represented by this landscape
#' @importFrom data.table melt
#' @importFrom stats coefficients binomial glm
#' @export
#' @return a data.table with assigned fuel classes for each species

assessFuelClasses <- function(landscape, fuelCol, sppEquiv,
                              sppEquivCol, targetFuelClasses) {

  landscape$B_MgHa <- landscape$B/100


  actualSpecies <- unique(landscape[!is.na(B),]$speciesCode)
  fuelGLMs <- lapply(actualSpecies, makeGLM, landscape = landscape)
  names(fuelGLMs) <- actualSpecies
  coeffs <- sapply(fuelGLMs, coefficients)
  pvalues <- sapply(fuelGLMs, FUN = function(x){summary(x)$coefficients[2,4]})
  speciesStats <- data.table(species = names(fuelGLMs), coef = coeffs[2,], pvalue = pvalues)


  #diagnostic plots
  Bmeans <- landscape[speciesCode %in% actualSpecies, .(meanB = as.integer(mean(B_MgHa)),
                                                        maxB = as.integer(max(B_MgHa))),
                      .(speciesCode)]
  Bmeans <- melt(Bmeans, id.vars = "speciesCode", value.name = "B_MgHa")
  preds <- lapply(actualSpecies, function(spec, GLM = fuelGLMs,
                                          toJoin = Bmeans) {

    GLM <- GLM[[spec]]
    df <- data.table(B_MgHa = seq(0, max(toJoin$B_MgHa), 2),
                     speciesCode = spec)
    #make sure mean and max are in pred frame, for plotting
    toJoin <- toJoin[speciesCode %in% spec, .(B_MgHa, speciesCode)]
    df <- rbind(toJoin, df, fill = TRUE)
    df[, burnprob := predict(GLM, newdata = df)]
    return(df)
  })
  preds <- rbindlist(preds)
  predsOfNote <- preds[Bmeans, on = c("speciesCode", "B_MgHa")]
  #TODO: decide whether this is helpful
  # ggplot(preds, aes(x = B_MgHa, y = burnprob, col = speciesCode)) +
  #          geom_line() +
  #          labs(x = "B (Mg/ha)", y = "burn pred") +
  #   geom_point(data = predsOfNote, aes(x = B_MgHa, y = burnprob, shape = variable))

  ####TODO: figure out what to do next - probably plot these predictions by biomass.
  forest <- landscape[!is.na(B),]
  fuels <-  sppEquiv[, .SD, .SDcols = c(sppEquivCol, fuelCol)]

  setnames(fuels, old = sppEquivCol, new = "speciesCode")
  forest <- fuels[forest, on = c("speciesCode")]
  forest[, propB := B/totalBiomass]
  forestPix <- nrow(forest[, .N, .(cell, year)]) #this is the only way because cells change between years
  selectCols <- c("speciesCode", fuelCol)
  above10PctRelB <- forest[propB > 0.1, .(above10PctRelB = .N/forestPix), by = c("speciesCode", fuelCol)]
  speciesStats[, sig := pvalue < 0.001]
  speciesStats[, sign := ifelse(coef > 0, "positive", "negative")]
  speciesStats[sig == FALSE, sign := "neutral"]
  speciesStats <- speciesStats[, .(species, coef, pvalue, sign)][above10PctRelB, on = c("species" = "speciesCode")]
  #   1. start by identifying species that have <5% of pixels with >10% biomass for merging
  #   2. Only species that are named in same FBPS class can be merged (we need to look at that now)... e.g., Pinu_con and Pinu_ban
  #   3. Ensure these are fitting within "negatives", "positives", and "neutrals"
  #   So, totally unrelated ones can't be merged, but otherwise allows merging of e.g. pines
  #   These would only be merged if they fulfill Criterion 1 above. Otherwise merged.
  #   if there is a small amount of some species left over after all merging, then put it in one of:
  #   "other positive", "other negative",
  #   If the species is allowed to be merged based on the Potential Fuel Class column,
  #   Revisit this classification if there are >5 classes, starting with higher thresholds for "Criterion 1"

  # Function to assign NewFuelClass

  out <- combine_fuel_classes(df = speciesStats, targetFuelClasses = targetFuelClasses)
  return(out)
}


#' merging and assignment of fuel classes
#'
#' @param df created during `assessFuelClasses`
#' @param targetFuelClasses the number of classes at which to stop further merging
#' @return a data.table with assigned fuel classes for each species
combine_fuel_classes <- function(df, targetFuelClasses = 5) {
  # Sort by FuelClass, sign, and above10PctRelB in ascending order for consistent processing
  df <- copy(df)
  nSpecies <- nrow(df)
  neededJoins <- nSpecies - targetFuelClasses
  if (neededJoins > 0) {
    df[, possGroups := .N, .(FuelClass)]
    #if each fuel class is unique no merge is possible
    #of course may not be able to merge due to sign either
    if (sum(df$possGroups) - nSpecies <= neededJoins) {
      warning("will not be able to reduce fuel groups to 5 or fewer")
    }

    uniques <- df[possGroups == 1,]
    possMerge <- df[!species %in% uniques$species]
    setkey(possMerge, above10PctRelB)
    #least biomass is now at the top
    joined <- 0
    attemptedJoins <- 0
    toMerge <- possMerge[0]

    for (i in 1:nrow(possMerge) -1) {
      toMerge <- possMerge[1]
      #sort by closest coefficient
      possMerge[, absCoef := abs(coef - toMerge$coef)]
      mergeable <- possMerge[species != toMerge$species &
                               FuelClass == toMerge$FuelClass]
      if (toMerge$sign != "neutral"){
        mergeable <- mergeable[sign == "neutral" | sign == toMerge$sign]
      }
      #confirm this works with zero length (no more mergeable options)
      if (nrow(mergeable) > 0) {
        mergeable <- mergeable[absCoef == min(absCoef)] #reduce to 1 choice
        possMerge[species == mergeable$species, assignedFuelClass := FuelClass]
        toMerge[, assignedFuelClass := FuelClass]
        joined <- joined + 1
      }
      uniques <- rbind(uniques, toMerge, fill = TRUE) #keep reclassed species
      attemptedJoins <- attemptedJoins + 1
      possMerge <- possMerge[!species %in% toMerge$species,] #remove it from list
      if (joined == neededJoins) {
        break
      }
    }
    newClass <- rbind(uniques, possMerge, fill = TRUE)
    df <- newClass[, .(species, coef, sign, FuelClass, above10PctRelB, assignedFuelClass)]
    df[is.na(assignedFuelClass), assignedFuelClass := species]

  } else {
    df[, assignedFuelClass := species]
  }
  return(df[])
}
