## reproducible::Cache digests only the called function's own code. A cached makeFireSenseLCC() call
## therefore keeps returning its old result when any function it CALLS changes, so callers must pass
## those functions in `.cacheExtra`. Only this package knows which they are -- the set depends on
## `lccSource`, which is a run-time option -- so a caller that hard-codes a list goes stale silently.
## That is exactly what happened: fireSense_dataPrepFit pinned prepInputs_NTEMS_LCC_FAO and kept it
## after the default moved to SCANFI, leaving SCANFI and CWIM changes invisible to the cache.

test_that("the dependency list holds the functions actually in the path, per lccSource", {
  ntems <- makeFireSenseLCCDeps("NTEMS")
  expect_true(is.list(ntems))
  expect_true(all(vapply(ntems, is.function, logical(1))))
  hasFn <- function(deps, fn) any(vapply(deps, function(d) identical(body(d), body(fn)), logical(1)))

  expect_true(hasFn(ntems, LandR::prepInputs_NTEMS_LCC_FAO))
  expect_false(hasFn(ntems, LandR::prepInputs_SCANFI_LCC_FAO))

  skip_if_not("prepInputs_CWIM" %in% getNamespaceExports("LandR"),
              "LandR without prepInputs_CWIM (PredictiveEcology/LandR#228)")
  scanfi <- makeFireSenseLCCDeps("SCANFI")
  expect_true(all(vapply(scanfi, is.function, logical(1))))
  expect_true(hasFn(scanfi, LandR::prepInputs_SCANFI_LCC_FAO))
  ## the wetland step is part of the SCANFI result, so it must be digested too
  expect_true(hasFn(scanfi, getExportedValue("LandR", "prepInputs_CWIM")))
  ## and the NTEMS function, which is NOT called, must not be
  expect_false(hasFn(scanfi, LandR::prepInputs_NTEMS_LCC_FAO))
})

test_that("switching lccSource changes the cache key it produces", {
  skip_if_not("prepInputs_CWIM" %in% getNamespaceExports("LandR"),
              "LandR without prepInputs_CWIM (PredictiveEcology/LandR#228)")
  ## the user-visible consequence: a cached call made under one source must not be reused under the
  ## other. With a hard-coded list both sources produced the same key.
  expect_false(identical(reproducible::.robustDigest(makeFireSenseLCCDeps("SCANFI")),
                         reproducible::.robustDigest(makeFireSenseLCCDeps("NTEMS"))))
})

test_that("lccSource is validated and case-insensitive, as in makeFireSenseLCC()", {
  expect_error(makeFireSenseLCCDeps("nope"), "SCANFI")
  expect_error(makeFireSenseLCCDeps(c("SCANFI", "NTEMS")), "SCANFI")
  expect_identical(makeFireSenseLCCDeps("ntems"), makeFireSenseLCCDeps("NTEMS"))
})
