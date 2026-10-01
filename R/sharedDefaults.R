#' Defaults shared by `fireSense_dataPrepFit` and `fireSense_dataPrepPredict`
#'
#' The two modules each have parameters with these names, and a fit and the predictions
#' made from it are only consistent if both modules use the same values. Both modules take
#' their parameter defaults from these constants, as they do `nonflammableLCC` from
#' [fireSenseNonflammableLCC], and the functions here that take the same arguments default
#' to them too.
#'
#' - `fireSenseForestedLCC`: land-cover codes treated as forest, in the land cover built by
#'   [makeFireSenseLCC()]: treed wetland (`81`), coniferous (`210`), broadleaf (`220`),
#'   mixedwood (`230`) and disturbed forest land (`240`).
#' - `fireSenseYoungAgeCutoff`: age (years) at and below which a pixel is `youngAge`.
#' - `fireSenseNonForestCanBeYoungAge`: whether burned non-forest is `youngAge` until it
#'   passes `fireSenseYoungAgeCutoff`.
#' - `fireSenseFlammabilityThreshold`: minimum proportion of flammable fine-resolution land
#'   cover for a coarser pixel to be flammable.
#' - `fireSenseFuelClassCol`: the column of `sppEquiv` that defines fuel classes.
#' - `fireSenseIgAggFactor`: aggregation factor (cells per side) for the ignition and
#'   escape covariates.
#' - `fireSenseSCANFIVersion`: the SCANFI land-cover version [makeFireSenseLCC()] reads for
#'   non-forest land cover (passed to [LandR::prepInputs_SCANFI_LCC_FAO()]'s `dataVersion`).
#'
#' @name fireSenseSharedDefaults
#' @aliases fireSenseForestedLCC fireSenseYoungAgeCutoff fireSenseNonForestCanBeYoungAge
#'   fireSenseFlammabilityThreshold fireSenseFuelClassCol fireSenseIgAggFactor
#'   fireSenseSCANFIVersion
#' @format Vectors of length 1, except `fireSenseForestedLCC` (`numeric`, length 5).
NULL

#' @rdname fireSenseSharedDefaults
#' @export
fireSenseForestedLCC <- c(81, 210, 220, 230, 240)

#' @rdname fireSenseSharedDefaults
#' @export
fireSenseYoungAgeCutoff <- 15

#' @rdname fireSenseSharedDefaults
#' @export
fireSenseNonForestCanBeYoungAge <- TRUE

#' @rdname fireSenseSharedDefaults
#' @export
fireSenseFlammabilityThreshold <- 0.1

#' @rdname fireSenseSharedDefaults
#' @export
fireSenseFuelClassCol <- "FuelClass"

#' @rdname fireSenseSharedDefaults
#' @export
fireSenseIgAggFactor <- 4

#' @rdname fireSenseSharedDefaults
#' @export
fireSenseSCANFIVersion <- "V3"
