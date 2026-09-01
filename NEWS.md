# fireSenseUtils 0.2.3

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
