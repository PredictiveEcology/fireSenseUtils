#' Defaults shared by the fireSense modules and functions
#'
#' The `fireSense_dataPrepFit` and `fireSense_dataPrepPredict` modules each have parameters
#' with these names, and a fit and the predictions made from it are only consistent if both
#' use the same values; `fireSense_spreadFit` takes its `runawayEdge*` defaults from here too.
#' The modules take their parameter defaults from these constants, as they do
#' `nonflammableLCC` from [fireSenseNonflammableLCC], and the functions here that take the
#' same arguments default to them too.
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
#' - `fireSenseRunawayEdgeFrac`, `fireSenseRunawayEdgeMin`: the defaults of `runawayEdgeFrac` and
#'   `runawayEdgeMin` of [.objfunSpreadFit()] and [runDEoptim()]: a simulated fire is a runaway when
#'   it burns at least `max(runawayEdgeMin, ceiling(runawayEdgeFrac * n))` of the `n` pixels of its
#'   buffer's edge ring.
#'
#' @name fireSenseSharedDefaults
#' @aliases fireSenseForestedLCC fireSenseYoungAgeCutoff fireSenseNonForestCanBeYoungAge
#'   fireSenseFlammabilityThreshold fireSenseFuelClassCol fireSenseIgAggFactor
#'   fireSenseSCANFIVersion fireSenseRunawayEdgeFrac fireSenseRunawayEdgeMin
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

#' @rdname fireSenseSharedDefaults
#' @export
fireSenseRunawayEdgeFrac <- 0.01

#' @rdname fireSenseSharedDefaults
#' @export
fireSenseRunawayEdgeMin <- 3L
