# fireSenseUtils 0.2.3

* `Firesense_LCC_flammability` vignette: probe the LCC (`ftp.maps.canada.ca`) and
  fire-polygon (`cwfis.cfs.nrcan.gc.ca`) data servers up front and skip the live
  download/analysis chunks when either is unreachable, degrading to
  documentation-only (same as the existing macOS/`archive` path). This stops
  R-CMD-check from failing when the servers are unreachable from CI runners.
  Adds `curl` to Suggests.

# fireSenseUtils 0.2.2

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
