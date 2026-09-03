## A run with no climate scenario carries NA -- and usually R's bare logical
## NA. Two consequences the plot helpers used to have:
##   * strsplit(NA, "_") errors with "non-character argument"
##   * pasting NA into a title or filename yields literal "NA" text, e.g.
##     burnSummary_someArea_NA.png

test_that(".scenarioParts never errors on a non-character scenario", {
  expect_identical(fireSenseUtils:::.scenarioParts(NA), character())            # logical NA
  expect_identical(fireSenseUtils:::.scenarioParts(NA_character_), character())
  expect_identical(fireSenseUtils:::.scenarioParts(character()), character())
  expect_error(strsplit(NA, "_"), "non-character")      # what the old code did
})

test_that(".scenarioParts splits a real scenario", {
  expect_identical(fireSenseUtils:::.scenarioParts("CNRM-ESM2-1_ssp370"),
                   c("CNRM-ESM2-1", "ssp370"))
  expect_identical(fireSenseUtils:::.scenarioParts("NRV"), "NRV")
})

test_that(".scenarioLabel calls an absent scenario NRV, and says so", {
  expect_message(expect_identical(fireSenseUtils:::.scenarioLabel(NA), "NRV"),
                 "declaring it to be NRV")
  expect_message(expect_identical(fireSenseUtils:::.scenarioLabel(NA_character_), "NRV"),
                 "declaring it to be NRV")
  expect_message(expect_identical(fireSenseUtils:::.scenarioLabel(""), "NRV"), "NRV")
  ## a filename built from it must never contain "NA"
  lbl <- suppressMessages(fireSenseUtils:::.scenarioLabel(NA))
  expect_false(grepl("_NA\\.png$", paste0("burnSummary_area_", lbl, ".png")))
})

test_that(".scenarioLabel passes a real scenario through silently", {
  expect_silent(
    expect_identical(fireSenseUtils:::.scenarioLabel("CNRM-ESM2-1_ssp370"),
                     "CNRM-ESM2-1_ssp370"))
})

test_that(".scenarioLabel errors on values that cannot be a scenario", {
  ## these mean a bug upstream, not an unset value; relabelling them NRV would
  ## hide it behind a plausible-looking filename
  expect_error(fireSenseUtils:::.scenarioLabel(NULL), "is missing")
  expect_error(fireSenseUtils:::.scenarioLabel(character()), "is missing")
  expect_error(fireSenseUtils:::.scenarioLabel(c("a_b", "c_d")), "must be length 1")
})
