# fireSenseUtils 0.2.3

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
