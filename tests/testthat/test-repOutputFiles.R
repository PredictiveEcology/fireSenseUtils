## .repOutputFiles() locates one file per replicate. Callers can hand it the
## real paths (`simFiles`, e.g. the `file` column of SpaDES.core::outputs()),
## rather than having them rebuilt from a directory convention -- which does not
## hold in every project: some use `rep1`, not `rep01`.

test_that(".repOutputFiles selects by filename and orders replicates numerically", {
  sf <- c("/o/rep1/burn.csv", "/o/rep10/burn.csv", "/o/rep2/burn.csv", "/o/rep2/other.csv")
  r <- fireSenseUtils:::.repOutputFiles(sf, outputDir = NULL, reps = NULL, filename = "burn.csv")
  expect_identical(unname(r), c("/o/rep1/burn.csv", "/o/rep2/burn.csv", "/o/rep10/burn.csv"))
  expect_identical(names(r), c("1", "2", "10"))     # numeric order, not lexical
  expect_false(any(grepl("other.csv", r)))          # other filenames excluded
})

test_that(".repOutputFiles accepts both rep1 and rep01 styles", {
  sf <- c("/o/rep01/burn.csv", "/o/rep2/burn.csv")
  r <- fireSenseUtils:::.repOutputFiles(sf, NULL, NULL, "burn.csv")
  expect_identical(names(r), c("1", "2"))
})

test_that(".repOutputFiles errors rather than returning nothing", {
  ## silently returning character(0) would surface much later as a confusing
  ## failure inside fread()/rbindlist()
  expect_error(
    fireSenseUtils:::.repOutputFiles(c("/o/rep1/burn.csv"), NULL, NULL, "absent.csv"),
    "no file named"
  )
})

test_that(".repOutputFiles falls back to the repNN convention when simFiles is NULL", {
  r <- fireSenseUtils:::.repOutputFiles(NULL, outputDir = "/o", reps = 1:3, filename = "x.csv")
  expect_identical(names(r), c("1", "2", "3"))
  ## NOTE: the fallback zero-pads. Projects whose directories are `rep1` (not
  ## `rep01`) must pass `simFiles`; see the PR discussion.
  expect_identical(basename(dirname(r)), c("rep01", "rep02", "rep03"))
})

test_that(".stopOnMclapplyErrors reports the failures mclapply swallowed", {
  ok <- c(TRUE, FALSE, TRUE)
  res <- list(data.frame(a = 1), simpleError("boom"), data.frame(a = 2))
  expect_error(fireSenseUtils:::.stopOnMclapplyErrors(res, ok, "burn summaries"),
               "burn summaries")
  expect_silent(fireSenseUtils:::.stopOnMclapplyErrors(res, c(TRUE, TRUE, TRUE), "x"))
})
