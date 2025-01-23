utils::globalVariables(c(
  "above10PctRelB", "absCoef", "assignedFuelClass", "B_MgHa", "burned", "burnprob",
  "cell", "coef", "lcc", "pixelIndex", "newName", "possGroups", "propB",
  "pvalue", "sig", "sim", "species","speciesCode", "YEAR", "yearRange"
))

#' Calculate proportional burn of landcover and tree species
#'
#' @template pixelGroupMap
#' @template cohortData
#' @param rstLCC a landcover map
#' @param nonflammableLCC nonflammable landcover in `rstLCC`
#' @param nonforestLCC vector or list of vectors of flammable nonforest LCC
#' @param fires a single `sf` or `SpatVector` object of fire polygons containing a `YEAR` column
#' @param yearRange  the range of years represented by this landscape
#'
#' @return `data.table` with cell, biomass, tree species or lcc for non-forest, and year of fire
#'
#' @examples
#' # fuelClassPrep(
#' #   sim$pixelGroupMap2011, sim$cohortData2011, sim$rstLCC2011,
#' #   rstLCC = sim$rstLCC2011, nonForestLCC = sim$nonForestLCCGroups,
#' #   nonflammableLCC = P(sim)$nonflammableLCC,
#' #   fires = do.call(rbind, sim$spreadFirePolys),
#' #   yearRange = c(2012, 2020)
#' # )
#'
#' @export
#' @importFrom data.table as.data.table
#'
fuelClassPrep <- function(pixelGroupMap, cohortData, rstLCC,
                          nonflammableLCC, fires, nonforestLCC, yearRange) {
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

  #lcc must supersede
  landscapeData[lcc %in% unlist(nonforestLCC), speciesCode := NA]
  landscapeData[, year := yearRange[1]]

  return(landscapeData)
}


#' wrapper on GLM that creates separate GLMs for each item in fuel
#'
#' @param fuel vector of tree species or lcc in landscape by which to subset
#' @param landscape created by `fuelClassPrep`
#' @param form the formula to use - as a character
#' @importFrom stats glm binomial
makeGLM <- function(fuel, landscape, form) {

  landscape <- landscape[speciesCode == fuel,]
  out <- glm(data = landscape, formula = as.formula(form), family = binomial(link = "logit"))
}

#' Semi-automated selection of fuel classed based on GLMs
#'
#' @param landscape data.table created by `prepFuelClasses`
#' @param fuelCol the column in `sppEquiv` with default fuel classes
#' @template sppEquiv
#' @template sppEquivCol
#' @param targetNonForestClasses the number of non-forest fuel classes to generate
#' @param targetFuelClasses target number of treed fuel classes to generate
#' @param nonforestLCC vector or list of vectors of non-forest fuel classes
#' @param pValue for glm coef significance when deciding to merge fuel classes
#' @importFrom data.table melt
#' @importFrom stats coefficients binomial glm kmeans
#' @export
#' @return a list of three objects:
#' 1. modSppEquiv, a data.table with assigned fuel classes for each tree species
#' 2. nonForestedLCCGroups, a named list of non-forest landcovers grouped by fuel class
#' 3. missingLCCgroup, the class in nonForestLCCGroups to assign forested pixels missing from cohortData

assessFuelClasses <- function(landscape, fuelCol, sppEquiv, sppEquivCol,
                              targetNonForestClasses = 2,
                              targetFuelClasses = 5, nonforestLCC, pValue = 0.001) {

  ####sort non-forest classes###
  nfData <- landscape[is.na(B)]
  #this will assume there is always missing landcover
  #TODO: review this assumption when everything is NTEMS-ified. Does it matter?
  nfData[!lcc %in% nonforestLCC, .N, .(lcc)]
  nfData[, lcc := as.character(lcc)]
  nfData[!lcc %in% nonforestLCC, lcc := "missingForest"]
  nfData[, lcc := as.factor(lcc)]
  #do it this way in case some lcc is not present

  #the makeGLM function will subset by values in species column, so copy lcc over
  nf_glms <- glm(data = nfData, formula = burned ~ 0 + lcc)
  #kmeans
  nfData[, .N, .(lcc)]

  #for ease later, keep these as character
  origNames <- as.character(levels(nfData$lcc))
  coefsToSort <- coefficients(nf_glms)
  names(coefsToSort) <- origNames
  nf_classes <- kmeans(x = coefsToSort, centers = targetNonForestClasses)
  nf_classes <- nf_classes$cluster

  #remove missing forest as it must be different
  missingForest <- nf_classes["missingForest"]

  nf_classes <- nf_classes[!names(nf_classes) %in% "missingForest"]
  missingForestFriends <- names(nf_classes)[nf_classes == missingForest]
  #there must be another way to do this?....
  nf_groups <- lapply(unique(nf_classes), FUN = function(class) {
    names(nf_classes)[nf_classes == class]
  })
  nf_vals <- lapply(nf_groups, as.numeric) #coerce back to lcc classes
  nf_groups <- sapply(nf_groups, paste, collapse = "_") |>
    lapply(FUN = function(x){ paste0("nfLCC_", x)})
  names(nf_vals) <- nf_groups

  #sort missingForest - it either gets grouped with others, or is alone
  if (length(missingForestFriends) > 0) {
    whichName <- sapply(nf_vals, function(x){
      all(x %in% as.numeric(missingForestFriends))})
    missingForest <- names(nf_vals)[whichName]
  } else {
    #missingForest is in a class all of its own
    #TODO: should this be allowed? it seems it will frequently occur
    # but only because burned pixels are more likely to be absent from landR
    # due to the fact they have no biomass
    missingForest <- "missingForest"
    nf_vals["missingForest"] <- -1 #new class that is not possible in a sane LCC
    #this should in theory be possible
  }


  #####sort the forested fuel classes (i.e. tree species)####
  landscape$B_MgHa <- landscape$B/100
  treeSpecies <- unique(landscape[!is.na(B)]$speciesCode)
  #in the event forested wetland should remain non-forest,
  #this allows forested wetland biomass without it counting towards flammability
  landscape <- landscape[!lcc %in% nonforestLCC & !is.na(B)]

  fuelGLMs <- lapply(treeSpecies, makeGLM, landscape = landscape, form = "burned ~ B_MgHa")
  names(fuelGLMs) <- treeSpecies
  coeffs <- sapply(fuelGLMs, coefficients)
  pvalues <- sapply(fuelGLMs, FUN = function(x){summary(x)$coefficients[2,4]})
  speciesStats <- data.table(species = names(fuelGLMs), coef = coeffs[2,], pvalue = pvalues)

  ####TODO: figure out what to do next - probably plot these predictions by biomass.
  forest <- landscape[!is.na(B),]
  fuels <-  sppEquiv[, .SD, .SDcols = c(sppEquivCol, fuelCol)]

  setnames(fuels, old = sppEquivCol, new = "speciesCode")
  forest <- fuels[forest, on = c("speciesCode")]
  forest[, propB := B/totalBiomass]
  forestPix <- nrow(forest[, .N, .(cell, year)]) #this is the only way because cells change between years
  selectCols <- c("speciesCode", fuelCol)
  above10PctRelB <- forest[propB > 0.1, .(above10PctRelB = .N/forestPix), by = c("speciesCode", fuelCol)]
  speciesStats[, sig := pvalue < pValue]
  speciesStats[, sign := ifelse(coef > 0, "positive", "negative")]
  speciesStats[sig == FALSE, sign := "neutral"]
  speciesStats <- speciesStats[, .(species, coef, pvalue, sign)][above10PctRelB, on = c("species" = "speciesCode")]


  modSppEquiv <- combine_fuel_classes(df = speciesStats, targetFuelClasses = targetFuelClasses)

  return(list(modSppEquiv = modSppEquiv,
              nonForestedLCCGroups = nf_vals,
              missingLCCgroup = missingForest))
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

  #fix the names in the event the original fuel classes aren't accurate
  #(e.g. SprcFrLrch may not contain Fir and/or Larch)
  if (any(df$species != df$assignedFuelClass)) {
    needsNewNames <- df[assignedFuelClass != species, .(assignedFuelClass, species)]
    needsNewNames[, newName := abbreviate(species)]
    needsNewNames[, newName := paste(newName, collapse = "."), .(assignedFuelClass)]
    needsNewNames <- needsNewNames[, .(species, newName)]
    df <- needsNewNames[df, on = c("species")]
    df[!is.na(newName), assignedFuelClass := newName]
    df[, newName := NULL]
  }

  return(df[])
}
