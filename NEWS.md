# fireSenseUtils 0.2.2

* `makeTSD()` is now general purpose: the pixels to age from fire history and the
  flammable mask can be supplied directly via the new `pixToUpdate` and
  `flammablePixels` arguments, instead of only through `lcc` (the `landcoverDT`
  from `fireSense_dataPrepFit`). `lcc` is now optional and, when supplied, still
  derives both vectors (the `fireSense_dataPrepFit`-specific behaviour is
  preserved behind an `if`); explicit arguments take precedence. Fully backwards
  compatible (#18).

# fireSenseUtils 0.1.5

* Fix namespace conflict warning: remove blanket `import(data.table)` in favour of
  explicit `importFrom` declarations; remove `importFrom(purrr, transpose)` which
  was silently overriding `data.table::transpose` (all call sites already use
  `purrr::transpose()` with explicit namespace).
* Add missing `importFrom(data.table, setorder)` used in `DEoptimIterative`.

# fireSenseUtils 0.0.5

* TODO
