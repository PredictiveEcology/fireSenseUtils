## .repOutputFiles() locates one file per replicate. Callers can hand it the
## real paths (`simFiles`, e.g. the `file` column of SpaDES.core::outputs()),
## rather than having them rebuilt from a directory convention -- which does not
## hold in every project: some use `rep1`, not `rep01`.

test_that(".repOutputFiles selects by filename and orders replicates numerically", {
  sf <- c("/o/rep1/burn.csv", "/o/rep10/burn.csv", "/o/rep2/burn.csv", "/o/rep2/other.csv")
  r <- fireSenseUtils:::.repOutputFiles(sf, filename = "burn.csv")
  expect_identical(unname(r), c("/o/rep1/burn.csv", "/o/rep2/burn.csv", "/o/rep10/burn.csv"))
  expect_identical(names(r), c("1", "2", "10"))     # numeric order, not lexical
  expect_false(any(grepl("other.csv", r)))          # other filenames excluded
})

test_that(".repOutputFiles accepts both rep1 and rep01 styles", {
  sf <- c("/o/rep01/burn.csv", "/o/rep2/burn.csv")
  r <- fireSenseUtils:::.repOutputFiles(sf, "burn.csv")
  expect_identical(names(r), c("1", "2"))
})

test_that(".repOutputFiles errors rather than returning nothing", {
  ## silently returning character(0) would surface much later as a confusing
  ## failure inside fread()/rbindlist()
  expect_error(
    fireSenseUtils:::.repOutputFiles(c("/o/rep1/burn.csv"), "absent.csv"),
    "no file named"
  )
})

test_that(".repOutputFiles requires simFiles rather than guessing paths", {
  ## reconstructing from a repNN convention is what this helper exists to
  ## avoid: it assumes zero-padded directories, which not every project uses,
  ## and yields paths that silently do not exist
  expect_error(fireSenseUtils:::.repOutputFiles(NULL, "x.csv"), "is required")
  expect_error(fireSenseUtils:::.repOutputFiles(character(), "x.csv"), "is required")
})

test_that(".stopOnMclapplyErrors reports the failures mclapply swallowed", {
  ok <- c(TRUE, FALSE, TRUE)
  res <- list(data.frame(a = 1), simpleError("boom"), data.frame(a = 2))
  expect_error(fireSenseUtils:::.stopOnMclapplyErrors(res, ok, "burn summaries"),
               "burn summaries")
  expect_silent(fireSenseUtils:::.stopOnMclapplyErrors(res, c(TRUE, TRUE, TRUE), "x"))
})
