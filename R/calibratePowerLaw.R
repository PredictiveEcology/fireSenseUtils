#' Calibrate the `"powerLaw"` link's `gamma` on a flammable map
#'
#' The `"powerLaw"` link ([powerLawLink()]) is `p = p_c - (p_c - floor) * exp(-x / gamma)`. `gamma` cannot be
#' fitted: it is confounded with the scale of the coefficients. This sets it from the landscape itself, so that
#' one unit of `x` is the same percent change in fire size in every ELF. Fires are spread with a constant
#' spread probability `p` on the flammable pixels (0 elsewhere), each from a random start in the largest
#' (8-neighbour) connected flammable patch, with `SpaDES.tools::spreadCpp()` as `fireSense_burn` and the fit
#' call it (8 directions, an escaped fire burns its first [escapeSizePixels()] pixels, `jumpTries`,
#' `jumpMeanDist`). A fire runs away when it burns at least `runawayShare` of that patch.
#'
#' `p50` is where the share of runaway fires crosses 0.5: the share is made monotone in `p` with isotonic
#' regression ([stats::isoreg()], which needs no smoothing choice) and interpolated linearly between the two grid
#' points that straddle 0.5. `gamma` is the negative slope of `log(mean fire size, ha)` on `log(p50 - p)` for the
#' `p` in `pGrid` from `pFitMin` up to `p50`, with `p_c` fixed at `p50` (more stable per ELF than fitting `p_c`
#' as well). On five ELFs' flammable maps (240 m) `p50` was 0.258-0.283, `gamma` 2.2-3.3 and `R2` 0.92-0.985.
#'
#' The result depends on `seed` and the arguments only, not on `cores`: each fire has its own L'Ecuyer-CMRG
#' stream. The caller's random-number state is left as it was.
#'
#' @param flammable A `SpatRaster`: 1 on flammable pixels, 0 or `NA` elsewhere (e.g. `flammableRTM`).
#' @param pGrid Increasing spread probabilities to simulate. The default ([spreadPowerLawPGrid]) brackets
#'   thresholds of 0.22-0.32 with a fine step there.
#' @param nFires Fires per value of `p`.
#' @param cores Number of forked processes ([parallel::mclapply()]); `1` runs in this process, as does any
#'   value where forking is unavailable ([parallelly::supportsMulticore()], e.g. Windows).
#' @param seed Integer seed.
#' @param escapeSizeHa,jumpTries,jumpMeanDist As in [.objfunSpreadFit()] and `fireSense_burn`.
#'   `escapeSizeHa` `NULL` or `NA` spreads each fire from its start pixel only.
#' @param runawayShare A fire is a runaway if it burns at least this share of the largest flammable patch.
#' @param pFitMin Lowest `p` used for `gamma`; very small fires are dominated by the escape size.
#' @return A list: `p50`, `gamma`, `R2` (of the `gamma` regression), and `table`, a `data.frame` with one row
#'   per `p`: `p`, `nFires`, `meanHa`, `shareRunaway` and `shareRunawayMonotone` (the isotonic fit).
#' @export
#' @examples
#' r <- terra::rast(nrows = 61, ncols = 61, xmin = 0, xmax = 14640, ymin = 0, ymax = 14640, vals = 1)
#' cal <- calibratePowerLaw(r, pGrid = seq(0.15, 0.4, by = 0.025), nFires = 10)
#' cal$p50
calibratePowerLaw <- function(flammable, pGrid = spreadPowerLawPGrid, nFires = 400L, cores = 1L, seed = 1L,
                              escapeSizeHa = 50, jumpTries = 20L, jumpMeanDist = 3,
                              runawayShare = 0.25, pFitMin = 0.15) {
  stopifnot(is.unsorted(pGrid) == FALSE, all(pGrid > 0 & pGrid < 1), nFires >= 1)
  flam <- which(terra::values(flammable, mat = FALSE) %in% 1)
  if (!length(flam)) stop("calibratePowerLaw(): `flammable` has no flammable pixel")
  patches <- terra::values(terra::patches(flammable == 1, directions = 8, zeroAsNA = TRUE), mat = FALSE)
  big <- as.integer(names(which.max(table(patches))))
  bigPix <- which(patches %in% big)
  nBig <- length(bigPix)
  rm(patches)
  minPx <- if (length(escapeSizeHa) && !is.na(escapeSizeHa)) escapeSizePixels(escapeSizeHa, flammable) else 1L
  pixHa <- prod(terra::res(flammable)) / 1e4
  nCell <- terra::ncell(flammable)
  jobs <- expand.grid(rep = seq_len(nFires), p = pGrid)

  ## One L'Ecuyer-CMRG stream per fire, from `seed`: the same fires whatever `cores` is
  streams <- list()
  withr::with_seed(seed, .rng_kind = "L'Ecuyer-CMRG", {
    jobs$start <- bigPix[sample.int(nBig, nrow(jobs), replace = TRUE)]
    s <- .Random.seed
    for (i in seq_len(nrow(jobs))) { streams[[i]] <- s; s <- parallel::nextRNGStream(s) }
  })
  oneFire <- function(i) {
    assign(".Random.seed", streams[[i]], envir = globalenv())
    sp <- numeric(nCell)
    sp[flam] <- jobs$p[i]
    ix <- SpaDES.tools::spreadCpp(landscape = flammable, loci = jobs$start[i], spreadProb = sp, directions = 8L,
                                  minSize = minPx, jumpTries = as.integer(jumpTries),
                                  jumpMeanDist = jumpMeanDist)$indices
    c(px = length(ix), inBig = sum(ix %in% bigPix))
  }
  withr::local_preserve_seed() # oneFire() installs each fire's stream; the caller's state comes back on exit
  if (!parallelly::supportsMulticore()) cores <- 1L # forking is unavailable on Windows (and in RStudio)
  res <- if (cores > 1L) {
    out <- parallel::mclapply(seq_len(nrow(jobs)), oneFire, mc.cores = cores, mc.preschedule = FALSE, mc.set.seed = FALSE)
    bad <- vapply(out, function(x) inherits(x, "try-error") || is.null(x), logical(1))
    if (any(bad)) stop("calibratePowerLaw(): ", sum(bad), " of ", length(out), " fires failed: ",
                       paste(unique(vapply(out[bad], function(x) paste(as.character(x), collapse = ""), character(1))), collapse = " "))
    out
  } else lapply(seq_len(nrow(jobs)), oneFire)
  res <- do.call(rbind, res)
  jobs$ha <- res[, "px"] * pixHa
  jobs$runaway <- res[, "inBig"] >= runawayShare * nBig

  tab <- data.frame(p = pGrid,
                    nFires = as.integer(nFires),
                    meanHa = as.vector(tapply(jobs$ha, jobs$p, mean)),
                    shareRunaway = as.vector(tapply(jobs$runaway, jobs$p, mean)))
  tab$shareRunawayMonotone <- stats::isoreg(tab$p, tab$shareRunaway)$yf
  i <- which(tab$shareRunawayMonotone >= 0.5)[1L]
  if (is.na(i)) stop("calibratePowerLaw(): runaway fires never reach half of the fires within pGrid (up to ", max(pGrid), ")")
  if (i == 1L) stop("calibratePowerLaw(): half of the fires run away already at the lowest p in pGrid (", min(pGrid), ")")
  y <- tab$shareRunawayMonotone[c(i - 1L, i)]
  p50 <- tab$p[i - 1L] + (0.5 - y[1L]) / (y[2L] - y[1L]) * (tab$p[i] - tab$p[i - 1L])
  below <- tab[tab$p >= pFitMin & tab$p < p50 & tab$meanHa > 0, ]
  if (nrow(below) < 4L)
    stop("calibratePowerLaw(): fewer than 4 values of p in pGrid between ", pFitMin, " and p50 = ", round(p50, 3))
  fit <- stats::lm(log(meanHa) ~ log(p50 - p), below)
  list(p50 = p50, gamma = unname(-stats::coef(fit)[2L]), R2 = summary(fit)$r.squared, table = tab)
}

#' Per-ELF bounds of `p_c` for the `"powerLaw"` link
#'
#' A fit that tries `p_c` far above the ELF's own threshold ([calibratePowerLaw()]`$p50`) only simulates
#' runaway fires that are then rejected. So the upper bound is `p50 + margin`; the lower bound is that of
#' [spreadPowerLawPcBounds]. `p50` must lie inside those outer limits.
#'
#' @param p50 The calibrated threshold.
#' @param margin Added to `p50` for the upper bound.
#' @return `c(lower, upper)` for the `maxAsymptote` slot.
#' @export
powerLawPcBounds <- function(p50, margin = 0.005) {
  outer <- spreadPowerLawPcBounds
  if (!(length(p50) == 1L && is.finite(p50) && p50 > outer[1L] && p50 <= outer[2L]))
    stop("powerLawPcBounds(): the calibrated p50 (", format(p50), ") is outside spreadPowerLawPcBounds (",
         paste(outer, collapse = " to "), ")")
  c(outer[1L], min(p50 + margin, outer[2L]))
}
