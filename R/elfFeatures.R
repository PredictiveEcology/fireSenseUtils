#' Per-ELF features and outlines for an interactive map
#'
#' `elfSpreadCoefs()` gives each ELF's spread-fit coefficients, `elfFeatures()` puts them with the fire rates and the
#' fit dates in one long table, and `elfOutlines()` draws the ELFs' cores for a web map.
#'
#' `elfSpreadCoefs()` takes the rows of [latestSpreadFits()] (the newest production spread fit of each ELF). A
#' coefficient is the mean over the row's parameter sets (`params`, one row per set). Coefficients are comparable
#' only between fits of one model form, so only the ELFs whose row has the most common form get values: the same
#' `spreadLink`, `spreadClimate` (predictor, ceiling or activation) and `fuelScale`, and the same climate split
#' (a `_yr` and `_dev` term per climate variable, or one term; see [isClimateSplitCol()]). The climate variable itself
#' does not matter: each variable is its own coefficient (`CMD`, `cumMDC`, ...) and appears only for the ELFs that
#' fitted it, as do the fuel terms (`youngAge`, `treedWetland_agb`, `dom_agb_*`, `sec_agb_*`, `nfLCC_*`). A fuel term of a
#' grouped fuel class (e.g. `dom_agb_Pp_tr.Bt_pa.Pp_ba`) is one coefficient per member species (`dom_agb_Popu_tre`,
#' ...), each with the term's value and a `note` of the species it is shared with; the members come from the row's
#' `sppEquiv` (species to fuel class), not from the abbreviated name. The others get
#' one row with `NA` and `excluded = "different model"`. A non-forest term such as
#' `nfLCC_100_40_60` applies to each of its land-cover classes, so each class also gets a row `nfLCC_<class>`,
#' with the term's value.
#'
#' `elfFeatures()` is one table for the map's layers, with the groups `"Observed"` (the fire rates of
#' [ELFfireRates()] over the core: natural ignitions per year, hectares burned per year, natural ignitions per 100
#' km^2 per year, percent of the core burned per year), `"Fit"` (the dates of the latest production fits, from
#' [elfFitDates()]; ignition and escape are one ledger row, so one date) and `"Spread coefficients"` (those of
#' `elfSpreadCoefs()`, with a column per land-cover class instead of per non-forest term). ELFs are those of the map
#' after the merges in `plan`; a value that is missing gets no row.
#'
#' `elfOutlines()` is the core of each ELF (not its buffer), one polygon per ELF, simplified and in longitude and
#' latitude.
#'
#' @param rows The result of [latestSpreadFits()], or `NULL`.
#' @param ELFs,plan,fireYears,pixelAreaHa,destinationPath As in [ELFfireRates()]: the ELFs before the merges, and
#'   `plan` the merges (`NULL`: none). For `elfOutlines()`, `ELFs` is the list as it is mapped (after the merges).
#' @param spreadFolder,ignitionFolder The ledger folders, as in [elfFitDates()].
#' @param dates The result of [elfFitDates()] if the caller has it already, else `NULL` (default): it is read.
#' @param tolerance Simplification tolerance of `elfOutlines()`, in metres of the ELFs' projection.
#' @param coefficient Coefficient names, e.g. a `coefficient` column of `elfSpreadCoefs()`.
#'
#' @return `elfSpreadCoefs()`: a `data.frame` with `ELF`, `coefficient`, `value`, `note` and `excluded` (`NA`, or
#'   `"different model"`); attribute `model` describes the common model. `elfFeatures()`: a `data.frame` with
#'   `ELF`, `group`, `feature` (an id), `label`, `units`, `type` (`"numeric"` or `"date"`), `value` (a date is
#'   the number of days since 1970-01-01, UTC) and `note` (`""`, or for a species of a grouped fuel term the species
#'   that share its coefficient). `elfOutlines()`: a `SpatVector` (EPSG:4326) with `ELF` and `name`
#'   ([elfLabels()]). `elfCoefLabel()`: a character vector.
#' @export
#' @rdname elfFeatures
elfSpreadCoefs <- function(rows) {
  if (is.null(rows) || !NROW(rows))
    return(.coefsTable())
  rows <- as.data.frame(rows)
  ids <- as.character(rows[[polygonIDTxt]])
  means <- lapply(rows$params, function(p) vapply(as.list(p), mean, numeric(1), na.rm = TRUE))
  column <- function(x) if (x %in% names(rows)) as.character(rows[[x]]) else rep(NA_character_, length(ids))
  model <- paste(column(spreadLinkTxt), column(spreadClimateTxt), column(spreadFitFuelScaleTxt),
                 vapply(means, function(m) any(isClimateSplitCol(names(m))), logical(1)), sep = "|")
  common <- names(which.max(table(model)))
  same <- model == common
  out <- lapply(seq_along(ids), function(i) {
    if (!same[i])
      return(.coefsTable(ids[i], NA_character_, NA_real_, "", "different model"))
    m <- means[[i]]
    nf <- grep("^nfLCC_", names(m), value = TRUE)
    ## a term of several land-cover classes applies to each: one row per class (a one-class term already is one)
    perClass <- do.call(rbind, lapply(nf, function(term) {
      classes <- strsplit(sub("^nfLCC_", "", term), "_")[[1]]
      data.frame(coefficient = paste0("nfLCC_", classes), value = m[[term]])
    }))
    all <- rbind(data.frame(coefficient = names(m), value = unname(m)), perClass)
    all <- .fuelTermsBySpecies(all, rows$sppEquiv[[i]])
    all <- all[!duplicated(all$coefficient), ]
    .coefsTable(ids[i], all$coefficient, all$value, all$note, NA_character_)
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  attr(out, "model") <- strsplit(common, "|", fixed = TRUE)[[1]]
  out
}

## The table `elfSpreadCoefs()` returns, columns in this order (no arguments: no rows)
.coefsTable <- function(ELF = character(), coefficient = character(), value = numeric(), note = character(),
                        excluded = character())
  data.frame(ELF = ELF, coefficient = coefficient, value = value, note = note, excluded = excluded)

## A fuel term (`dom_agb_<class>`, `sec_agb_<class>`) of a grouped fuel class becomes one coefficient per member
## species, each with the term's value and a `note` naming the species that share it. The members are the species
## whose `fireSenseFuelClassCol` of the row's `sppEquiv` is the class; a class with no match (or no `sppEquiv`) is
## a species of its own. `coefs` has `coefficient` and `value`; the result adds `note`.
.fuelTermsBySpecies <- function(coefs, sppEquiv) {
  n <- NROW(coefs)
  coefs$note <- ""
  prefixed <- paste0("^(", domAgbTxt, "|", secAgbTxt, ")")
  fuel <- grepl(prefixed, coefs$coefficient)
  if (!any(fuel) || is.null(sppEquiv) || !fireSenseFuelClassCol %in% names(sppEquiv))
    return(coefs)
  byClass <- lapply(split(as.character(sppEquiv[["LandR"]]), as.character(sppEquiv[[fireSenseFuelClassCol]])), unique)
  members <- vector("list", n)                       # NULL: not a fuel term, or a class with no member
  members[fuel] <- byClass[sub(prefixed, "", coefs$coefficient[fuel])]
  hasMembers <- lengths(members) > 0L
  times <- pmax(lengths(members), 1L)
  out <- coefs[rep(seq_len(n), times), , drop = FALSE]
  at <- rep(hasMembers, times)
  prefix <- sub(paste0(prefixed, ".*"), "\\1", coefs$coefficient[hasMembers])
  species <- unlist(members[hasMembers])
  out$coefficient[at] <- paste0(rep(prefix, lengths(members[hasMembers])), species)
  group <- rep(seq_along(prefix), lengths(members[hasMembers]))
  out$note[at] <- vapply(seq_along(species), function(i) {
    others <- setdiff(members[hasMembers][[group[i]]], species[i])
    if (length(others)) paste("shared with", paste(others, collapse = ", ")) else ""
  }, character(1))
  out
}

#' @export
#' @rdname elfFeatures
elfCoefLabel <- function(coefficient) {
  species <- function(codes) {
    ## the species' English name from whichever column of the table holds the code; the code itself if it lacks one
    name <- LandR::equivalentName(codes, LandR::sppEquivalencies_CA, "EN_generic_full")
    ifelse(is.na(name), codes, paste0(name, " (", codes, ")"))
  }
  plain <- c("(Intercept)" = "Intercept", maxAsymptote = "Maximum spread probability", hillSlope1 = "Hill slope",
             inflectionPoint1 = "Inflection point", yearSpreadSD = "Year SD", youngAge = "Young stands",
             treedWetland_agb = "Treed wetland AGB")
  out <- coefficient
  isPlain <- coefficient %in% names(plain)
  out[isPlain] <- plain[coefficient[isPlain]]
  dom <- startsWith(coefficient, domAgbTxt)
  sec <- startsWith(coefficient, secAgbTxt)
  out[dom] <- paste0("Dominant AGB: ", species(sub(paste0("^", domAgbTxt), "", coefficient[dom])))
  out[sec] <- paste0("Secondary AGB: ", species(sub(paste0("^", secAgbTxt), "", coefficient[sec])))
  nf <- grepl("^nfLCC_", coefficient)
  out[nf] <- paste("Non-forest LCC", gsub("_", ", ", sub("^nfLCC_", "", coefficient[nf])))
  out
}

#' @export
#' @rdname elfFeatures
elfFeatures <- function(ELFs, plan, fireYears, pixelAreaHa, destinationPath,
                        spreadFolder = spreadFitLedgerFolder(), ignitionFolder = ignitionFitLedgerFolder(),
                        dates = NULL) {
  if (is.null(plan)) plan <- noELFPlan()
  onMap <- names(mergeELFs(ELFs, plan)$rasWhole)
  feature <- function(ELF, group, feature, label, units, type, value, note = "")
    data.frame(ELF = as.character(ELF), group = group, feature = feature, label = label, units = units, type = type,
               value = as.numeric(value), note = note, stringsAsFactors = FALSE)

  rates <- as.data.frame(ELFfireRates(ELFs, plan, fireYears = fireYears, pixelAreaHa = pixelAreaHa,
                                      destinationPath = destinationPath))
  rates$percentBurnedPerYear <- 100 * rates$proportionBurnedPerYear
  observed <- data.frame(
    feature = c("naturalFiresPerYear", "haBurnedPerYear", "naturalFiresPer100km2PerYear", "percentBurnedPerYear"),
    label = c("Natural ignitions per year", "Area burned per year", "Natural ignitions per 100 km2 of core per year",
              "Core burned per year"),
    units = c("ignitions/yr", "ha/yr", "ignitions/100 km2/yr", "%/yr"), stringsAsFactors = FALSE)
  out <- list(do.call(rbind, lapply(seq_len(nrow(observed)), function(i)
    feature(rates$ELF, "Observed", observed$feature[i], observed$label[i], observed$units[i], "numeric",
            rates[[observed$feature[i]]]))))

  if (is.null(dates)) dates <- elfFitDates(spreadFolder, ignitionFolder, destinationPath, ELFs = onMap)
  days <- function(x) as.numeric(as.Date(x))
  out <- c(out, list(
    feature(dates$ELF, "Fit", "spreadFitDate", "Spread fit date", "date", "date", days(dates$spreadFitDate)),
    feature(dates$ELF, "Fit", "ignitionFitDate", "Ignition and escape fit date", "date", "date",
            days(dates$ignitionFitDate))))

  coefs <- elfSpreadCoefs(latestSpreadFits(spreadFolder, destinationPath))
  ## a land-cover class has its own row; the term of several classes it came from has none
  coefs <- coefs[!is.na(coefs$coefficient) & !grepl("^nfLCC_[0-9]+_[0-9]", coefs$coefficient), ]
  if (NROW(coefs))     # (paste0() of no coefficients would still give "spread_")
    out <- c(out, list(feature(coefs$ELF, "Spread coefficients",
                               paste0("spread_", sub("^\\(Intercept\\)$", "Intercept", coefs$coefficient)),
                               elfCoefLabel(coefs$coefficient), "", "numeric", coefs$value, coefs$note)))

  out <- do.call(rbind, out)
  out <- out[out$ELF %in% onMap & !is.na(out$value), ]
  rownames(out) <- NULL
  out
}

#' A merge plan with no merges
#'
#' The `plan` of [mergeELFs()] and [plotELFs()] when no ELF is merged or left out.
#'
#' @return A `data.table` with `action`, `ELF` and `members`, and no rows.
#' @export
noELFPlan <- function() data.table::data.table(action = character(0), ELF = character(0), members = list())

#' @export
#' @rdname elfFeatures
#' @importFrom terra aggregate project simplifyGeom
elfOutlines <- function(ELFs, tolerance = 1000) {
  poly <- .ELFpoly(ELFs)
  core <- terra::aggregate(poly[poly$buffer %in% 2, "ID"], by = "ID")
  core <- terra::project(terra::simplifyGeom(core, tolerance = tolerance), "EPSG:4326")
  out <- core[, "ID"]
  names(out) <- "ELF"
  out$name <- elfLabels(out$ELF, "name")
  out
}
