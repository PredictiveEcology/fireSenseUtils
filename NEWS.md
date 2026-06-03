# fireSenseUtils 0.2.3

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
