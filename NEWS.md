# fireSenseUtils 0.2.3.9005

## Bug fixes

* `bufferToArea()` and `rasterFireBufferDT()` with `cores > 1` no longer hang. GDAL keeps one
  worker-thread pool per process, created at the first multi-threaded raster write. A forked
  worker inherited that pool but none of its threads, and waited for them forever with 0 CPU.
  Forked workers now write rasters single-threaded (`terraOptions(threads = 1)`).

# fireSenseUtils 0.2.3.9004

## Bug fixes

* `fireSenseCloudParameters()` now downloads the shared parameter file from Google
  Drive on every call. It used `prepInputs(purge = 7, overwrite = TRUE)`, which never
  downloads again once a copy on disk matches CHECKSUMS.txt (`purge` only rebuilds
  those entries; `overwrite` only affects the written output), so a changed file on
  Drive was not seen. `url` may now also be the folder containing `targetFile`;
  `useCache` is ignored.

# fireSenseUtils 0.2.3.9003

## Performance

* `objFunInner()` no longer wraps its `spread()` replicate loop in
  `system.time()`. `system.time()` defaults to `gcFirst = TRUE`, so every fire year
  forced a full garbage collection, and the result was assigned to a variable that
  was never read. Measured on ELF 12.4 (1,737,132 cells, 18 fire years, `Nreps = 25`,
  in a 20-40 GB process), three consecutive objective-function evaluations went from
  82.7 / 84.0 / 83.9 s to 9.8 / 8.7 / 8.8 s -- about 9.5x -- with warm caches and no
  other change. In the profile the forced collections were 109.96 s of a 134 s
  evaluation, 82% of self time, against 6.44 s for all the `spread()` calls they were
  timing. It also explains why replicate count used to make no difference: the
  collection is once per fire year, outside the replicate loop. With it gone,
  `Nreps = 5` takes 3.2 s against 8.8 s at 25, so roughly 80% of an evaluation now
  scales with replicates as it always should have.

# fireSenseUtils 0.2.3.9002

## Bug fixes

* `objFunInner()` returned `ret` unconditionally but built it only inside
  `if (isTRUE(doFitting))`. With no test selected -- which is how
  `fireSense_SpreadFit`'s `mode = "debug"` calls the chain, passing `tests = ""`
  -- the branch is skipped and the call failed with `object 'ret' not found`,
  making debug mode unusable. `ret` is now initialised before the branch, so the
  function returns an empty list when nothing was asked of it.

# fireSenseUtils 0.2.3.9001

* `objFunInner()` now asks `SpaDES.tools::spread()` to skip its input checks by
  the name that function actually uses. It passed `skipChecks = TRUE`, which
  belongs to `spread3()`; `spread()` has no such argument, so it landed in `...`
  and did nothing. Every call therefore re-validated the whole per-cell
  `spreadProb` vector -- `na.omit()` copies it, `inRange()` scans it -- which is
  work proportional to the number of cells in the landscape, independent of how
  much actually burns, and this function calls `spread()` `Nreps` times per fire
  year. Measured on synthetic landscapes, with identical output:

  | landscape | before | after | speedup |
  | --- | --- | --- | --- |
  | 1M cells | 0.030 s | 0.016 s | 1.9x |
  | 4M cells | 0.074 s | 0.022 s | 3.4x |
  | 9M cells | 0.183 s | 0.036 s | 5.1x |

# fireSenseUtils 0.2.3

* `objFunSpread()`: the Anderson-Darling test now compares fire size distributions
  drawn from the *same* years on both sides. Years are simulated in two batches --
  the 2 with the largest total area burned, then the rest -- and the `adTest` block
  was gated on `ii == 2`, reading `results$fireSizes` after `results` had been
  overwritten by the second batch. So the simulated sample held only the
  small-fire years, while the observed sample used `historicalFiresAboveMin`
  unsubset, i.e. *every* year. Since the omitted years are selected as
  largest-area-burned, the observed sample retained an upper tail the simulated
  sample structurally could not have; `ad.test` is tail-sensitive and its result is
  multiplied by 50, so the only way `DEoptim` could shrink it was to inflate
  simulated fire sizes in ordinary years, biasing the fit toward over-prediction.
  Simulated sizes are now pooled across both batches (they were already being
  computed and discarded) and the test runs once after the loop, via the new
  internal `pooledFireSizes()`. Behaviour on the early-bail path is unchanged: a
  parameter set that fails the first batch still skips `adTest` entirely.
* `makeELFs()`: only run the NA-hole-filling `focal()` step when `x` is a
  `SpatRaster`. When `x` is (or defaults to) an `sf` of fire regime polygons there
  are no NA slivers to fill and `focal()` has no method for it, failing with
  `no method found for signature sf, data.frame`.

* `mergeAndSplitRas()`: write each call's per-ecoprovince `.tif` files into their own
  `ELFs_<digest>` subdirectory of `destinationPath`, keyed on `ecopRseg`, `ecopLCC`,
  `maxArea` and `field`. They previously went straight into `destinationPath` named
  only by the province code (`4.1.tif`), scattering ~36 anonymous files through a
  shared inputs directory and letting successive calls overwrite one another. Note
  this protects `destinationPath` only: `reproducible`'s cache-restore path collapses
  these to `cachePath/<basename>.tif`, dropping the cacheId, so distinct calls still
  collide there -- that is a `reproducible` issue this cannot fix.

* `fireSenseCloudParametersMap()` and `plotELFs()`: use `ELFs$poly` rather than
  treating the `makeELFs()` return value as a `SpatVector`. `makeELFs()` returns a
  list of rasters plus a `poly` element (`ELFsInStudyArea()` already unwraps it this
  way), so `terra::plot(ELFs)`, `ELFs$buffer`, `ELFs$ID` and `ELFs[keep, ]` were all
  operating on the list. Both call sites are fixed.

* Documentation regenerated with roxygen2 8.1.0 (was 8.0.0). Mostly formatting:
  8.1.0 emits one multi-symbol `importFrom()` per package instead of one line per
  symbol, so NAMESPACE shrinks considerably with no change to what is imported.
  One **API change** falls out of the upgrade: `spreadFitAdditionalColNames`
  (`R/objDefaults.R`) carries `@export` in its source but was not being exported --
  roxygen2 8.0.0 silently dropped it, and 8.1.0 honours it. It is now exported and
  documented, alongside the existing `spreadFitAdditionalColNamesTxt`.

* `getFirePoints_NFDB()` no longer drops columns. It previously subset to
  `c("YEAR", fireSizeColName)` and renamed those to `date`/`size_ha`, discarding
  `CAUSE` among everything else -- the reason the scfm modules explicitly avoid
  fireSenseUtils (`## NOTE: do not use fireSenseUtils - it removes the cause
  column`). All source columns are now retained, matching
  `getFirePoints_NFDB_V2()` and `scfmutils::getFirePoints_NFDB_scfm()`. The
  derived `size` column (fire size in pixels) is still added when `rasterToMatch`
  is supplied. **Breaking:** the `date` and `size_ha` output names are gone; use
  `YEAR` and `SIZE_HA` (#32).

* `getFirePoints_NFDB_V2()` now uses `fun` to load already-downloaded data, not
  just downloaded data. The cached branch hardcoded `st_read()`, so a caller
  passing `fun = "terra::vect"` (as `fireSense_dataPrepFit` does) got a
  `SpatVector` on the first run and an `sf` on every run after. The return class
  is now the same in both branches (#32).

* `getFirePoints_NFDB()` and `getFirePoints_NFDB_V2()`: fix inverted
  `redownloadIn`. The staleness threshold was `365 / redownloadIn`, so
  `redownloadIn = 0.5` gave a 730-day tolerance rather than the documented
  "redownload data older than 6 months". Now `365 * redownloadIn`. Only the
  default of `1` was unaffected (#32).

* `getFirePoints_NFDB()` and `getFirePoints_NFDB_V2()`: the staleness check
  errored with "the condition has length > 1" when `NFDB_pointPath` held more
  than one `NFDB_point*.shp`; now wrapped in `any()` (#32).

* `getFirePoints_NFDB_V2()`: restore the `NFDB_pointPath` non-NULL check, so the
  `NULL` default fails with a clear message instead of inside `Checksums()` (#32).

* `getFirePoints_NFDB()`: pass `useSAcrs = TRUE` to `postProcess()` in the
  already-downloaded branch, so the CRS no longer depends on whether a download
  happened; restore `reproducible.cacheSaveFormat` to its previous value on exit
  instead of overwriting it; drop the dead `SpaDES.core` requirement check (the
  package is only in Suggests and nothing in the function uses it) (#32).
* Replace the single `ggpubr::ggarrange()` call in `visualizeDE()` with
  `cowplot::plot_grid()`, and move `cowplot` from Suggests to Imports (it was
  already used, behind a `requireNamespace()` guard, in `plot_summaries.R`).
  `ggpubr` is dropped from Imports. This removes 14 packages from the
  dependency tree -- ggpubr, car, carData, pbkrtest, doBy, Deriv, lme4,
  quantreg, SparseM, MatrixModels, nloptr, minqa, abind and Formula.

  It also fixes R-CMD-check on the `oldrel-2` leg. `Deriv 4.3.0` (published
  2026-07-23) uses `R_ClosureFormals`, an R >= 4.5 C-API entry point, but
  declares no minimum R version, so it fails to compile on R 4.4.x with
  `'R_ClosureFormals' was not declared in this scope`. The dependency install
  then aborted and `rcmdcheck` died at `loadVignetteBuilder()`, which made the
  failure look like a vignette problem.

* `makeTSD()` is now general purpose: the pixels to age from fire history and the
  flammable mask can be supplied directly via the new `pixToUpdate` and
  `flammablePixels` arguments, instead of only through `lcc` (the `landcoverDT`
  from `fireSense_dataPrepFit`). `lcc` is now optional and, when supplied, still
  derives both vectors (the `fireSense_dataPrepFit`-specific behaviour is
  preserved behind an `if`); explicit arguments take precedence. Fully backwards
  compatible (#18).

* `Firesense_LCC_flammability` vignette: probe the LCC (`ftp.maps.canada.ca`) and
  fire-polygon (`cwfis.cfs.nrcan.gc.ca`) data servers up front and skip the live
  download/analysis chunks when either is unreachable, degrading to
  documentation-only (same as the existing macOS/`archive` path). This stops
  R-CMD-check from failing when the servers are unreachable from CI runners.
  Adds `curl` to Suggests.

* `makeFireSenseLCC()`: remove commented-out dead code (an unused majority-NA
  block aggregation path) from the aggregation step.

# fireSenseUtils 0.1.5

* Fix namespace conflict warning: remove blanket `import(data.table)` in favour of
  explicit `importFrom` declarations; remove `importFrom(purrr, transpose)` which
  was silently overriding `data.table::transpose` (all call sites already use
  `purrr::transpose()` with explicit namespace).
* Add missing `importFrom(data.table, setorder)` used in `DEoptimIterative`.

# fireSenseUtils 0.0.5

* TODO
