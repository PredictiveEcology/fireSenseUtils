utils::globalVariables(c(
  "above10PctRelB", "absCoef", "assignedFuelClass", "B_MgHa", "burned",
  "burnprob", "cell", "coef", "dominantSign",  "genus", "lcc", "pixelIndex",
  "newName", "possGroups", "propB", "pvalue", "sig", "sim", "spec",
  "species","speciesCode", "tempSign", "YEAR", "yearRange"
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
#' @param nonforestLCC vector or list of vectors of flammable non-forest LCC
#'
#' @param fires a single `sf` or `SpatVector` object of fire polygons containing a `YEAR` column
#'
#' @param yearRange  the range of years represented by this landscape
#'
#' @return `data.table` with cell, biomass, tree species or LCC for non-forest, and year of fire
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


#' Wrapper on `glm` that creates separate GLMs for each item in fuel
#'
#' @param fuel vector of tree species or LCC in landscape by which to subset
#'
#' @param landscape created by `fuelClassPrep`
#'
#' @param form the formula to use - as a character
#'
#' @importFrom stats glm binomial
makeGLM <- function(fuel, landscape, form) {

  landscape <- landscape[speciesCode == fuel,]
  out <- glm(data = landscape, formula = as.formula(form), family = binomial(link = "logit"))
}

#' Semi-automated selection of fuel classed based on GLMs
#'
#' @param landscape data.table created by `prepFuelClasses`
#'
#' @param fuelCol the column in `sppEquiv` with default fuel classes
#'
#' @template sppEquiv
#'
#' @template sppEquivCol
#'
#' @param targetNonForestClasses the number of non-forest fuel classes to generate
#'
#' @param targetFuelClasses target number of treed fuel classes to generate
#'
#' @param nonforestLCC vector or list of vectors of non-forest fuel classes
#'
#' @param pValue for glm coefficient significance when deciding to merge fuel classes
#'
#' @return a list of three objects:
#' 1. `modSppEquiv`, a data.table with assigned fuel classes for each tree species;
#' 2. `nonForestedLCCGroups`, a named list of non-forest land cover classes grouped by fuel class;
#' 3. `missingLCCgroup`, the class in `nonForestLCCGroups` to assign forested pixels missing from `cohortData`;
#'
#' @export
#' @importFrom data.table melt
#' @importFrom stats coefficients binomial glm kmeans
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
  hasMissingForest <- FALSE
  if ("missingForest" %in% names(coefsToSort)) {
    hasMissingForest <- TRUE
    coefsToSort <- coefsToSort[setdiff(names(coefsToSort), "missingForest")]
  }
  nf_classes <- kmeans(x = coefsToSort, centers = targetNonForestClasses)
  nf_classes <- nf_classes$cluster

  #remove missing forest as it must be different
  #there must be another way to do this?....
  nf_groups <- lapply(unique(nf_classes), FUN = function(class) {
    names(nf_classes)[nf_classes == class]
  })
  nf_vals <- lapply(nf_groups, as.numeric) #coerce back to lcc classes
  nf_groups <- sapply(nf_groups, paste, collapse = "_") |>
    lapply(FUN = function(x){ paste0("nfLCC_", x)})
  names(nf_vals) <- nf_groups

  #sort missingForest - it either gets grouped with others, or is alone
  if (hasMissingForest) {
    origCoefs <- coefficients(nf_glms)
    names(origCoefs) <- origNames
    missingCoef <- origCoefs["missingForest"]
    other  <- origCoefs[names(origCoefs) != "missingForest"]
    other <- abs(other - missingCoef)
    coefMatch <- names(other)[which(other == min(other))]
    missingForest <- names(nf_vals)[grep(pattern = as.numeric(coefMatch), x = nf_vals )]
    #essentially - give it the new class name for whichever landcover it was closet to
  } else {
    #this doesn't matter as there is none but return something
    missingForest <- names(nf_vals)[1]
  }

  #####sort the forested fuel classes (i.e. tree species)####
  landscape$B_MgHa <- landscape$B/100
  #in the event forested wetland should remain non-forest,
  #this allows forested wetland biomass without it counting towards flammability
  landscape <- landscape[!lcc %in% nonforestLCC & !is.na(B)]
  treeSpecies <- unique(landscape[!is.na(B)]$speciesCode)

  fuelGLMs <- lapply(treeSpecies, makeGLM, landscape = landscape, form = "burned ~ B_MgHa")
  names(fuelGLMs) <- treeSpecies
  coeffs <- sapply(fuelGLMs, coefficients)
  pvalues <- sapply(fuelGLMs, FUN = function(x){summary(x)$coefficients[2,4]})
  speciesStats <- data.table(species = names(fuelGLMs), coef = coeffs[2,], pvalue = pvalues)


  forest <- landscape[!is.na(B),]
  fuels <-  unique(sppEquiv[, .SD, .SDcols = c(sppEquivCol, fuelCol)])

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

#' Merging and assignment of fuel classes
#'
#' @param df created during `assessFuelClasses`
#'
#' @param targetFuelClasses the number of classes at which to stop further merging
#'
#' @param lowThreshold species with `above10PctRelB` below this threshold will be merged with like
#' fuel classes, when available, regardless of coefficient
#'
#' @return a data.table with assigned fuel classes for each species
#'
#' @importFrom data.table copy setkey
combine_fuel_classes <- function(df, targetFuelClasses = 5, lowThreshold = 0.05) {
  ## Sort by FuelClass, sign, and above10PctRelB in ascending order for consistent processing
  df <- copy(df)
  df[, assignedFuelClass := character(0)]
  nSpecies <- nrow(df)

  df[, N := .N, .(FuelClass)]
  df[, dominantSign := sign]
  df_Sparse <- df[above10PctRelB < lowThreshold & N > 1]

  ## 1. merge these sparse ones first because DEoptim will struggle with small sample size
  if (nrow(df_Sparse) > 0) {
    for (i in 1:nrow(df_Sparse)) {
      df_Friends <- df[FuelClass %in% df_Sparse[i, ]$FuelClass,]
      df_Friends <- df_Friends[!species %in% df_Sparse$species,]
      ## if multiple of same class are in sparse these will be merged in 2
      ## keeping them in results in errors
      ## order by fuel class, sign, and abundance
      ## ensure the sign is prioritized
      df_Friends[, tempSign := sign]
      ## ensure setting key by sign yields like sign at top, alphabetically
      df_Friends[sign == df_Sparse[i,]$sign, sign := paste0("a-", sign)]
      setkey(df_Friends, FuelClass, sign, above10PctRelB)
      ## merge both
      df[species %in% df_Sparse[i,]$species, assignedFuelClass := FuelClass]
      ## merge with the 2nd least abundant
      ## TODO: this 2nd subset won't be correct if sign is different.
      df[species %in% df_Friends[1,]$species, assignedFuelClass := FuelClass]
      df[species %in% df_Sparse[i,]$species, dominantSign := df_Friends[1, ]$dominantSign]

      df_Friends[, sign := tempSign]
      df_Friends[, tempSign := NULL]
    }
  }

  neededJoins <- nSpecies - targetFuelClasses - nrow(df_Sparse) ## already joined
  ##. 2 join according to sign compatibility, then similar coefficient value,
  #     in order of increasing biomass
  if (neededJoins > 0) {

    df[, possGroups := .N, .(FuelClass)]
    ## if each fuel class is unique no merge is possible
    ## of course may not be able to merge due to sign either
    if (sum(df$possGroups) - nSpecies <= neededJoins) {
      warning("will not be able to reduce fuel groups to 5 or fewer")
    }

    uniques <- df[possGroups == 1,]
    possMerge <- df[!species %in% uniques$species]

    ## if none were merged above, create this column
    if (is.null(possMerge$assignedFuelClass)) {
      possMerge[, assignedFuelClass := character(0)]
    }
    setkey(possMerge, assignedFuelClass, above10PctRelB)
    ## least biomass is now at the top
    joined <- 0
    attemptedJoins <- 0
    toMerge <- possMerge[0]

    for (i in 1:nrow(possMerge) -1) {
      toMerge <- possMerge[1]
      ## sort by closest coefficient
      possMerge[, absCoef := abs(coef - toMerge$coef)]
      mergeable <- possMerge[species != toMerge$species &
                               FuelClass == toMerge$FuelClass]
      if (toMerge$sign != "neutral"){
        mergeable <- mergeable[sign == "neutral" | sign == toMerge$sign]
      }
      ## confirm this works with zero length (no more merge-able options)
      if (nrow(mergeable) > 0) {
        mergeable <- mergeable[absCoef == min(absCoef)] ## reduce to 1 choice
        possMerge[species == mergeable$species, assignedFuelClass := FuelClass]
        toMerge[, assignedFuelClass := FuelClass]
        joined <- joined + 1
      }
      uniques <- rbind(uniques, toMerge, fill = TRUE) ## keep reclassed species
      attemptedJoins <- attemptedJoins + 1
      possMerge <- possMerge[!species %in% toMerge$species,] ## remove it from list
      if (joined == neededJoins) {
        break
      }
    }
    newClass <- rbind(uniques, possMerge, fill = TRUE)
    df <- newClass[, .(species, coef, sign, FuelClass, above10PctRelB, assignedFuelClass)]
    df[is.na(assignedFuelClass), assignedFuelClass := species]

  } else {
    df[is.na(assignedFuelClass), assignedFuelClass := species]
  }

  ## fix the names in the event the original fuel classes aren't accurate
  ## (e.g. SprcFrLrch may not contain Fir and/or Larch)
  if (any(df$species != df$assignedFuelClass)) {
    df <- abbreviateSpNames(df)
    # needsNewNames <- df[assignedFuelClass != species, .(assignedFuelClass, species)]
    # #if there is an underscore assume it separates genus/species
    # hasUnderscore <- needsNewNames[grep("_", species)]
    #
    # myFun <- function(STRING, placement){
    #   out <- strsplit(STRING, split = "_")
    #   sapply(out, "[[", placement)
    # }
    # hasUnderscore[, genus := abbreviate(species)]
    # hasUnderscore[, genus := lapply(.SD, FUN = myFun, placement = 1), .SDcol = "genus"]
    # hasUnderscore[, spec := lapply(.SD, FUN = myFun, placement = 2), .SDcol = "species"]
    # hasUnderscore[, spec := substr(spec, start = 1, stop = 2)]
    # hasUnderscore[, newName := paste0(genus, "_", spec)]
    # hasUnderscore[, newName := paste(newName, collapse = "."), .(assignedFuelClass)]
    # hasUnderscore <- hasUnderscore[, .(assignedFuelClass, species, newName)]
    #
    # #else
    # needsNewNames <- needsNewNames[!species %in% hasUnderscore$species]
    # needsNewNames[, newName := abbreviate(species, minlength = 5, strict = TRUE)]
    # needsNewNames[, newName := paste(newName, collapse = "."), .(assignedFuelClass)]
    # #join
    # needsNewNames <- rbind(hasUnderscore, needsNewNames)
    # needsNewNames <- needsNewNames[, .(species, newName)]
    #
    #
    # df <- needsNewNames[df, on = c("species")]
    # df[!is.na(newName), assignedFuelClass := newName]
    # df[, newName := NULL]
  }
  df[, N := .N, .(assignedFuelClass)]

  return(df[])
}


#' Abbreviate species names for fuel classes
#'
#' @param df data.table
#' @return data.table
abbreviateSpNames <- function(df) {
  afcName <- "assignedFuelClass"
  if (is.null(df[[afcName]]))
    set(df, NULL, afcName, sapply(seq(NROW(df)), function(x) basename(tempfile())))
  needsNewNames <- df[assignedFuelClass != species, .(assignedFuelClass, species)]
  ## if there is an underscore assume it separates genus/species
  hasUnderscore <- needsNewNames[grep("_", species)]

  myFun <- function(STRING, placement){
    out <- strsplit(STRING, split = "_")
    sapply(out, "[[", placement)
  }
  hasUnderscore[, genus := abbreviate(species)]
  hasUnderscore[, genus := lapply(.SD, FUN = myFun, placement = 1), .SDcol = "genus"]
  hasUnderscore[, spec := lapply(.SD, FUN = myFun, placement = 2), .SDcol = "species"]
  hasUnderscore[, spec := substr(spec, start = 1, stop = 2)]
  hasUnderscore[, newName := paste0(genus, "_", spec)]
  hasUnderscore[, newName := paste(newName, collapse = "."), .(assignedFuelClass)]
  hasUnderscore <- hasUnderscore[, .(assignedFuelClass, species, newName)]

  ## else
  needsNewNames <- needsNewNames[!species %in% hasUnderscore$species]
  needsNewNames[, newName := abbreviate(species, minlength = 5, strict = TRUE)]
  needsNewNames[, newName := paste(newName, collapse = "."), .(assignedFuelClass)]

  ## join
  needsNewNames <- rbind(hasUnderscore, needsNewNames)
  needsNewNames <- needsNewNames[, .(species, newName)]

  df <- needsNewNames[df, on = c("species")]
  df[!is.na(newName), assignedFuelClass := newName]
  df[, newName := NULL]
}




#' Create FireSense Input Covariates
#'
#' Constructs pixel-level covariates required by FireSense (e.g., mutually
#' exclusive fuel-class indicators, transformed by \code{logMinB}, and a
#' \code{youngAge} indicator) by combining cohort-derived fuel classes with
#' land cover, flammability, and (optionally) non-forest time-since-disturbance.
#'
#' The function:
#' \enumerate{
#'   \item Calls \code{cohortsToFuelClasses()} to derive fuel-class rasters from cohorts.
#'   \item Converts these to a long data table keyed by \code{pixelID} and joins to \code{landcoverDT}.
#'   \item Verifies there are no \code{NA} values across covariate columns and raises an error otherwise.
#'   \item Marks pixels with all-zero covariates as the provided \code{missingLCCgroup}.
#'   \item Enforces mutual exclusivity among the specified fuel/landcover columns.
#'   \item Applies \code{logMinB} to all fuel-class columns (excluding \code{youngAge}).
#'   \item If \code{nonForestCanBeYoungAge = TRUE}, augments \code{youngAge} for non-forest pixels
#'         using \code{nonForest_timeSinceDisturbance} and the \code{cutoffForYoungAge}.
#' }
#'
#' @param cohortData \code{data.table} (or similar) of stand/age cohorts used to derive fuel classes.
#' @param pixelGroupMap A raster (SpatRaster) or mapping structure linking cohort groups to pixels.
#' @param flammableRTM A raster (SpatRaster) mask of flammable pixels (non-flammable typically \code{NA}/0).
#' @param sppEquiv Species equivalency table used to map species codes to fuel-class groupings.
#' @param landcoverDT \code{data.table} with at least \code{pixelID} and one or more land-cover
#'   indicators (e.g., logical or numeric columns for LCC groups) used in the covariate set.
#' @param fuelClassCol Character name of the fuel-class column produced/expected by \code{sppEquiv}.
#' @param sppEquivCol Character name of the species-equivalency column in \code{sppEquiv}.
#' @param missingLCCgroup A single unquoted column name (passed to \code{data.table:::=} via \code{eval})
#'   that will be set to 1 where all covariates are zero (forested LCC absent from \code{cohortData}).
#' @param nonForestedLCCGroups Named vector or list where names correspond to the non-forest land-cover
#'   columns in \code{landcoverDT} (used to detect non-forest pixels).
#' @param nonForest_timeSinceDisturbance A raster/vector of time-since-disturbance (TSD) for non-forest pixels;
#'   indexed by \code{pixelID} and used when \code{nonForestCanBeYoungAge = TRUE}.
#' @param cutoffForYoungAge Numeric threshold (years). Pixels with TSD \eqn{\le} this value are considered
#'   \code{youngAge}.
#' @param nonForestCanBeYoungAge Logical. If \code{TRUE}, non-forest pixels can be flagged as \code{youngAge}
#'   based on \code{nonForest_timeSinceDisturbance} and \code{cutoffForYoungAge}.
#' @param studyAreaName Character tag used for caching/user tags.
#' @param useCache Logical. If \code{TRUE}, enable caching in \code{Cache()} calls.
#'
#' @details
#' \itemize{
#'   \item All covariate columns must be non-\code{NA}; the function stops with an error if any \code{NA}
#'         are detected after joining \code{fuelClasses} and \code{landcoverDT}.
#'   \item Pixels where all covariates sum to zero are flagged by setting \code{missingLCCgroup := 1}.
#'   \item Mutual exclusivity is enforced via \code{makeMutuallyExclusive()} using all fuel-class columns
#'         and the land-cover columns (excluding \code{pixelID}); \code{youngAge} is made exclusive to them.
#'   \item Fuel-class columns (excluding \code{youngAge}) are transformed with \code{logMinB}.
#'   \item When \code{nonForestCanBeYoungAge = TRUE}, \code{youngAge} for non-forest is derived using
#'         \code{nonForest_timeSinceDisturbance} and \code{cutoffForYoungAge}. Any existing \code{youngAge}
#'         from fuel classes is combined with this non-forest contribution.
#' }
#'
#' @return A \code{data.table} named \code{spreadCovariates} containing one row per \code{pixelID}, with:
#' \itemize{
#'   \item Mutually exclusive fuel-class covariate columns (log-transformed by \code{logMinB}).
#'   \item Land-cover indicator columns from \code{landcoverDT}.
#'   \item A \code{youngAge} indicator (possibly augmented by non-forest TSD logic).
#' }
#'
#' @section Expectations/Assumptions:
#' \itemize{
#'   \item \code{fuelClasses} returned by \code{cohortsToFuelClasses()} includes one column per fuel class,
#'         plus optionally \code{youngAge}, and a \code{cell} (renamed to \code{pixelID} here).
#'   \item \code{landcoverDT} must contain a \code{pixelID} column aligned to rasters by cell index.
#'   \item \code{nonForestedLCCGroups} names correspond to column names in \code{landcoverDT}.
#'   \item \code{nonForest_timeSinceDisturbance} is indexable by \code{pixelID} (e.g., vector aligned to raster cells).
#' }
#'
#' @seealso
#'   \code{\link{cohortsToFuelClasses}}, \code{\link{makeMutuallyExclusive}},
#'   \code{\link{logMinB}}, \code{\link{putBackIntoRaster}}, \code{\link{calcNonForestYoungAge}}
#'
#' @examples
#' \dontrun{
#' covs <- fireSenseCovariatesCreate(
#'   cohortData = cohortsDT,
#'   pixelGroupMap = pgRas,
#'   flammableRTM = flammRas,
#'   sppEquiv = sppEquivDT,
#'   landcoverDT = landcoverDT,
#'   fuelClassCol = "fuelClass",
#'   sppEquivCol = "sppEquivCol",
#'   missingLCCgroup = "LCC_missing",
#'   nonForestedLCCGroups = c(NonForest = 1),
#'   nonForest_timeSinceDisturbance = NF_TSD,
#'   cutoffForYoungAge = 20,
#'   nonForestCanBeYoungAge = TRUE,
#'   studyAreaName = "MyArea",
#'   useCache = TRUE
#' )
#' }
#'
#' @export
fireSenseCovariatesCreate <- function(cohortData, 
                                   pixelGroupMap,
                                   flammableRTM,
                                   sppEquiv,
                                   landcoverDT,
                                   fuelClassCol,
                                   sppEquivCol,
                                   missingLCCgroup,
                                   nonForestedLCCGroups,
                                   nonForest_timeSinceDisturbance,
                                   cutoffForYoungAge,
                                   nonForestCanBeYoungAge,
                                   studyAreaName, useCache = TRUE) {
  
  fuelClasses <- cohortsToFuelClasses(
    cohortData = cohortData,
    pixelGroupMap = pixelGroupMap,
    flammableRTM = flammableRTM,
    landcoverDT = landcoverDT,
    sppEquiv = sppEquiv,
    fuelClassCol = fuelClassCol,
    sppEquivCol = sppEquivCol,
    cutoffForYoungAge = cutoffForYoungAge
  )
  
  
  ## make columns for each fuel class
  # fuelClasses <- terra::app(fuelClasses, fun = logMinB)
  # terra app is horrifically slow
  fcs <- setdiff(names(fuelClasses), "youngAge")
  fuelClasses <- as.data.table(as.data.frame(fuelClasses, cells = TRUE))
  setnames(fuelClasses, old = "cell", new = "pixelID")
  
  # make sure join is only landcoverDT -- this adds the nonForest that are in sim$landcoverDT
  spreadCovariates <- fuelClasses[landcoverDT, on = c("pixelID")]
  
  
  ## Nov 2023 - there should not be NA values - previously this used nafill
  ## if they return - use x <- as.data.table(nafill(vegData), 0) and setnames(x, names(vegData))
  spreadCovariates[, rowcheck := rowSums(.SD), .SD = setdiff(names(spreadCovariates), "pixelID")]
  if (any(is.na(spreadCovariates$rowcheck))) {
    stop("NA in vegData columns of fireSense_dataPrepPredict... please contact module developers")
  }
  # if all rows are 0, it must be a forested LCC absent from cohortData
  spreadCovariates[rowcheck == 0, (missingLCCgroup) := 1]
  set(spreadCovariates, NULL, "rowcheck", NULL)
  
  # Making exclusive has to be prior to logMinB, or else the 0 biomass become -0.59 or so
  #   --> they need to stay at the minimum of 3.605
  exclusiveCols <- c(fcs, names(landcoverDT))
  exclusiveCols <- setdiff(exclusiveCols, "pixelID")
  spreadCovariates <- makeMutuallyExclusive(dt = spreadCovariates,
                                            mutuallyExclusive = list("youngAge" = exclusiveCols))
  
  spreadCovariates <- spreadCovariates[, eval(fcs) := lapply(.SD, FUN = logMinB), .SDcols = fcs]
  
  
  if (nonForestCanBeYoungAge) {
    dig1 <- reproducible::.robustDigest(list(landcoverDT, flammableRTM))
    
    LCCras <- putBackIntoRaster(landcoverDT = landcoverDT, # list(sim$landcoverDT2010, sim$landcoverDT2020),
                                flammableMap = flammableRTM, # list(sim$flammableRTM2010, sim$flammableRTM2020),
                                lcc = names(nonForestedLCCGroups))
    # LCCras <- Map(
    #   f = putBackIntoRaster,
    #   landcoverDT = landcoverDT, # list(sim$landcoverDT2010, sim$landcoverDT2020),
    #   flammableMap = flammableRTM, # list(sim$flammableRTM2010, sim$flammableRTM2020),
    #   MoreArgs = list(lcc = names(nonForestedLCCGroups))
    # ) |>
    #   reproducible::Cache(.functionName = "putBackIntoRaster",
    #         .cacheExtra = dig1, omitArgs = c("landcoverDT", "flammableMap"),
    #         userTags = c("putBackIntoRaster", studyAreaName))
    # dig1b <- reproducible::.robustDigest(list(LCCras))
    
    # LCCras <- Map(
    #   f = calcNonForestYoungAge,
    #   landcoverDT = landcoverDT, # list(sim$landcoverDT2010, sim$landcoverDT2020),
    #   NFTSD = nonForest_timeSinceDisturbances, # list(sim$nonForest_timeSinceDisturbance2010,
    #   # sim$nonForest_timeSinceDisturbance2020),
    #   LCCras = LCCras, #  list(LCCras[[1]], LCCras[[2]]),
    #   MoreArgs = list(cutoffForYoungAge = cutoffForYoungAge)
    # ) |>
    #   reproducible::Cache(.cacheExtra = append(dig1), omitArgs = c("landcoverDT", "LCCras"))#, "flammableMap"),)
    
    LCCras <- calcNonForestYoungAge(
      landcoverDT = landcoverDT, # list(sim$landcoverDT2010, sim$landcoverDT2020),
      NFTSD = nonForest_timeSinceDisturbance, # list(sim$nonForest_timeSinceDisturbance2010,
      # sim$nonForest_timeSinceDisturbance2020),
      LCCras = LCCras, #  list(LCCras[[1]], LCCras[[2]]),
      cutoffForYoungAge = cutoffForYoungAge
    ) |>
      reproducible::Cache(.cacheExtra = dig1, useCache = useCache, 
                          omitArgs = c("landcoverDT", "LCCras"))#, "flammableMap"),)
    
    #for (i in names(fuelClasses)) {
    if (youngAgeName %in% names(fuelClasses)) {
      
      YA1 <- fuelClasses[[youngAgeName]]
      YA2 <- terra::values(LCCras[[youngAgeName]])[fuelClasses$pixelID]
      bothYA <- YA1 + YA2
      fuelClasses[[youngAgeName]] <- bothYA
    }  else {
      fuelClasses[[youngAgeName]] <- LCCras[[youngAgeName]]
    }
    # toKeep <- setdiff(names(LCCras), youngAgeName)
    # LCCras <- terra::subset(LCCras, toKeep) ## to avoid double-counting
    #}
    
    # this should only alter non-forest
    spreadCovariates[, isNonForest := rowSums(.SD) > 0, .SDcols = names(nonForestedLCCGroups)]
    spreadCovariates[, YA_NF := as.vector(nonForest_timeSinceDisturbance)[spreadCovariates$pixelID] <= cutoffForYoungAge &
                       isNonForest == TRUE]
    spreadCovariates[YA_NF == TRUE, youngAge := 1]
    spreadCovariates[, c("YA_NF", "isNonForest") := NULL]
  }
  spreadCovariates
}





#' Rebuild land cover layers back into a SpatRastser
#'
#' @description
#' Given a template raster (e.g., a flammability map) and a `data.table`
#' containing pixel values for one or more land cover classes, this function
#' reconstructs each class as a raster layer and returns a multi-layer
#' `SpatRaster` whose layer names match the class names supplied in `lcc`.
#'
#' @details
#' This function:
#' 1. Creates a working copy of the template raster's values (initialized to 0
#'    on all non-`NA` cells and keeping original `NA`s).
#' 2. For each class name in `lcc`, it pulls the corresponding vector from
#'    `landcoverDT`, assigns those values into the
#'    working vector at positions indicated by `landcoverDT$pixelID`, and then
#'    writes these values back into a fresh `terra::rast(flammableMap)` layer.
#' 3. Stacks all layers into a single `terra::rast()` and sets layer names to
#'    `lcc`.
#'
#' **Assumptions / requirements**
#' - `flammableMap` is a `terra::SpatRaster` providing the spatial geometry,
#'   extent, resolution, and `NA` mask.
#' - `landcoverDT` is a `data.table` with a column named `pixelID` that indexes
#'   cell positions in the template raster's value vector (i.e., values returned
#'   by `terra::values(flammableMap, mat = FALSE)`). Indices must be valid and
#'   1-based (R-style).
#' - Each element of `lcc` is the name of a numeric column in `landcoverDT`.
#' - All layers are created in the same CRS, resolution, and extent as
#'   `flammableMap`.
#'
#' **Notes**
#' - Non-`NA` cells in the template are initialized to 0 before assignment, so
#'   only positions referenced by `pixelID` will receive class values; other
#'   non-`NA` positions remain 0 in each layer.
#' - If you want to preserve existing values outside `pixelID`, modify the
#'   initialization step accordingly.
#'
#' @param lcc character vector.
#'   Names of the land cover class columns in `landcoverDT` to rebuild into
#'   raster layers (e.g., `c("Conifer", "Deciduous", "Shrub")`).
#' @param landcoverDT data.table.
#'   A `data.table` that must contain:
#'   - a column `pixelID` with integer indices into the template raster's value
#'     vector (as returned by `terra::values(flammableMap, mat = FALSE)`), and
#'   - one numeric column per class named in `lcc`.
#' @param flammableMap SpatRaster.
#'   A `terra::SpatRaster` used as the spatial template (geometry, extent,
#'   resolution, and `NA` mask) for the output layers.
#'
#' @return
#' A multi-layer `terra::SpatRaster` where each layer corresponds to one entry
#' in `lcc`, with layer names set to `lcc`.
#'
#' @examples
#' \dontrun{
#' library(terra)
#' library(data.table)
#'
#' # Template raster (5x5 example)
#' r <- rast(nrows = 5, ncols = 5, vals = NA)
#' # Mark some cells as valid (non-NA)
#' values(r) <- rep(NA_real_, ncell(r))
#' valid_idx <- c(3, 5, 7, 13, 17, 20)
#' tmp_vals <- values(r, mat = FALSE)
#' tmp_vals[valid_idx] <- 0
#' r <- setValues(r, tmp_vals)
#'
#' # Data table with pixel indices and two classes
#' dt <- data.table(
#'   pixelID = valid_idx,
#'   Conifer  = c(0.6, 0.2, 0.8, 0.1, 0.0, 0.4),
#'   Deciduous = c(0.1, 0.3, 0.0, 0.4, 0.5, 0.2)
#' )
#'
#' # Rebuild layers
#' out <- putBackIntoRaster(
#'   lcc = c("Conifer", "Deciduous"),
#'   landcoverDT = dt,
#'   flammableMap = r
#' )
#'
#' out
#' names(out)        # "Conifer" "Deciduous"
#' plot(out)
#' }
#'
#' @seealso [terra::rast()], [terra::values()], [terra::setValues()]
#'
#' @importFrom terra rast values setValues
#' @export
putBackIntoRaster <- function(lcc, landcoverDT, flammableMap) {
  
  stopifnot(is.character(lcc), length(lcc) > 0)
  stopifnot(all(lcc %in% colnames(landcoverDT)))
  stopifnot("pixelID" %in% colnames(landcoverDT))
  stopifnot(inherits(flammableMap, "SpatRaster"))
  
  lccRasters <- list()
  for (i in 1:length(lcc)) {
    flamMapVals <- terra::values(flammableMap, mat = FALSE)
    flamMapVals[!is.na(flamMapVals)] <- 0
    lccRasters[[i]] <- terra::rast(flammableMap)
    lccVals <- landcoverDT[, get(lcc[i])]
    flamMapVals[landcoverDT$pixelID] <- lccVals
    lccRasters[[i]] <- terra::setValues(lccRasters[[i]], flamMapVals)
  }
  lccRasters <- terra::rast(lccRasters)
  names(lccRasters) <- lcc
  
  return(lccRasters)
}




#' Zero-out young non-forest pixels across LCC layers and add a `youngAge` mask
#'
#' @description
#' Identifies pixels that have **any non-forest presence** (based on the
#' non-`pixelID` columns of `landcoverDT`) and whose **age** is **below**
#' `cutoffForYoungAge`. It then sets those pixels to 0 across **all layers** of
#' `LCCras` and appends a new single-layer raster named `youngAge` with 1 for
#' the affected pixels and 0 elsewhere (preserving `NA` where present in the
#' template).
#'
#' @details
#' The function proceeds as follows:
#' 1. Determines the set of non-forest columns as all columns in `landcoverDT`
#'    except `pixelID`, and computes a per-row sum (`sumRows`).
#' 2. Retrieves pixel ages from `NFTSD` using `pixelID` as indices.
#' 3. Selects pixels where `sumRows > 0` **and** `age < cutoffForYoungAge`;
#'    these pixels are considered "young non-forest".
#' 4. Sets those pixels to 0 in **every layer** of `LCCras` using cell indexing.
#' 5. Builds a new single-layer raster with 1 at those pixel locations and 0
#'    elsewhere (respecting `NA` from the template), and adds it to `LCCras`
#'    as a layer named `youngAge`.
#'
#' **Assumptions / requirements**
#' - `LCCras` is a `terra::SpatRaster` (possibly multi-layer). Cell indexing
#'   with `LCCras[pixToChange] <- 0` applies across all layers for the indexed
#'   cells.
#' - `landcoverDT` is a `data.table` with:
#'   - a `pixelID` column of integer cell indices (1-based) consistent with
#'     `LCCras` cell numbering, and
#'   - all other columns representing non-forest variables whose row-wise sum
#'     indicates presence (`sumRows > 0`).
#' - `NFTSD` is a numeric vector (or similar) indexable by `pixelID`, providing
#'   the age (in the same units as `cutoffForYoungAge`) for each cell.
#' - `cutoffForYoungAge` is a scalar numeric threshold; pixels with `age <
#'   cutoffForYoungAge` are flagged as "young".
#'
#' **Notes**
#' - The function modifies `LCCras` in place (by reference to the object passed
#'   in) and returns the modified raster with an extra `youngAge` layer.
#' - `NA` cells in `LCCras` remain `NA` in the created `youngAge` layer.
#' - The `landcoverDT` temporary columns `sumRows` and `age` are removed prior
#'   to return.
#'
#' @param landcoverDT data.table.
#'   A `data.table` containing a `pixelID` column of integer cell indices
#'   (1-based) and one or more non-forest columns (all columns except
#'   `pixelID` are interpreted as non-forest for the row-sum).
#' @param NFTSD numeric.
#'   A numeric vector (or similar) that can be indexed by `pixelID` to obtain
#'   per-cell age values.
#' @param LCCras SpatRaster.
#'   A `terra::SpatRaster` (single- or multi-layer) whose cells correspond to
#'   `pixelID`. Pixels identified as "young non-forest" are set to 0 in all
#'   layers. A new `youngAge` layer is appended.
#' @param cutoffForYoungAge numeric(1).
#'   Age threshold; pixels with `age < cutoffForYoungAge` and non-forest
#'   presence are flagged as young.
#'
#' @return
#' A `terra::SpatRaster` identical to `LCCras` but with:
#' - all "young non-forest" pixels set to 0 across its existing layers, and
#' - an additional layer named `youngAge` where flagged pixels are 1, others are
#'   0, and `NA`s follow the template mask.
#'
#' @examples
#' \dontrun{
#' library(terra)
#' library(data.table)
#'
#' # Build a small template raster
#' r <- rast(nrows = 3, ncols = 3)
#' values(r) <- c(NA, 0, 0,
#'                0, 0, NA,
#'                0, 0, 0)
#'
#' # Make a 2-layer LCC raster (copy template twice)
#' LCC <- c(r, r)
#' names(LCC) <- c("LCC1", "LCC2")
#'
#' # A data.table with pixel indices and two non-forest columns
#' dt <- data.table(
#'   pixelID = c(2, 3, 4, 5, 7, 8, 9),
#'   NF1 = c(1, 0, 0, 2, 0, 0, 1),
#'   NF2 = c(0, 1, 0, 0, 0, 3, 0)
#' )
#'
#' # Age vector indexed by cell ID (length ncell(LCC))
#' ages <- seq_len(ncell(LCC)) * 5  # toy ages
#'
#' out <- calcNonForestYoungAge(
#'   landcoverDT = dt,
#'   NFTSD = ages,
#'   LCCras = LCC,
#'   cutoffForYoungAge = 15
#' )
#'
#' out
#' names(out)     # "LCC1" "LCC2" "youngAge"
#' plot(out$youngAge)
#' }
#'
#' @seealso
#' - `terra::rast()`, `terra::values()`, `terra::setValues()` for raster creation
#'   and manipulation.
#' - `data.table` semantics for `.SD`, in-place updates (`:=`), and fast
#'   row/column operations.
#'
#' @references
#' - Terra package reference manual and vignettes for `SpatRaster` operations
#'   (Hijmans, R.J.).  
#' - data.table reference: `.SD`, `:=`, and efficient grouping/filtering
#'   (Dowle, M., & Srinivasan, A.).
#'
#' @importFrom terra rast values setValues
#' @export
calcNonForestYoungAge <- function(landcoverDT, NFTSD, LCCras, cutoffForYoungAge) {
  nfLCC <- setdiff(colnames(landcoverDT), "pixelID")
  landcoverDT[, sumRows := rowSums(.SD), .SDcol = nfLCC]
  landcoverDT[, age := NFTSD[pixelID]]
  # this need to be chagned in LCCras and also converted to a youngAge raster
  # as the rasters will be aggregated
  pixToChange <- landcoverDT[sumRows > 0 & age < cutoffForYoungAge]$pixelID
  
  landcoverDT[, c("sumRows", "age") := NULL]
  
  # check how terra works - this should change all pixels if length 2+ spat raster
  LCCras[pixToChange] <- 0
  youngAge <- rast(LCCras, nlyr = 1)
  
  temp <- values(LCCras[[1]])
  temp[!is.na(temp)] <- 0
  temp[pixToChange] <- 1
  
  youngAge <- setValues(youngAge, temp)
  
  LCCras$youngAge <- youngAge
  
  return(LCCras)
}
