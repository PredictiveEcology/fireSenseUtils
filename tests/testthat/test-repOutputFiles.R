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
  expect_identical(r, c("1" = "/o/rep01/burn.csv", "2" = "/o/rep2/burn.csv"))
})

test_that(".repOutputFiles takes each replicate's file from its own folder only", {
  ## A file of the same name in a subfolder of a replicate was taken as a second file of that
  ## replicate (two elements named "1"), and one in a folder named by a number that is not a
  ## replicate (`2023/`) as a replicate of its own.
  sf <- c("/o/rep1/burn.csv", "/o/rep1/figures/burn.csv", "/o/rep10/burn.csv",
          "/o/2023/burn.csv", "/o/burn.csv")
  r <- fireSenseUtils:::.repOutputFiles(sf, "burn.csv")
  expect_identical(r, c("1" = "/o/rep1/burn.csv", "10" = "/o/rep10/burn.csv"))
})

test_that(".repOutputFiles counts a path listed twice once", {
  sf <- c("/o/rep1/burn.csv", "/o/rep2/burn.csv", "/o/rep1/burn.csv")
  r <- fireSenseUtils:::.repOutputFiles(sf, "burn.csv")
  expect_identical(r, c("1" = "/o/rep1/burn.csv", "2" = "/o/rep2/burn.csv"))
})

test_that(".repOutputFiles stops on a replicate without the file, naming its folder", {
  ## rep2 has a burn map but no burn summary: it was left out of the figure, silently
  sf <- c("/o/rep1/burn.csv", "/o/rep1/map.tif", "/o/rep2/map.tif")
  expect_error(fireSenseUtils:::.repOutputFiles(sf, "burn.csv"), "burn.csv in /o/rep2")
  ## or has it only in a subfolder, which is not the replicate's own output
  expect_error(fireSenseUtils:::.repOutputFiles(c(sf, "/o/rep2/old/burn.csv"), "burn.csv"),
               "burn.csv in /o/rep2")
})

test_that(".repOutputFiles stops when no folder in simFiles is a replicate", {
  expect_error(fireSenseUtils:::.repOutputFiles(c("/o/burn.csv", "/o/sum/burn.csv"), "burn.csv"),
               "no replicate folder")
})

test_that(".repOutputFiles finds bare-number replicate folders, ordered numerically", {
  ## as SpaDES.project::pathBuild() writes them with the field label off
  sf <- c("/o/6.4/NRV/370/16/burn.csv", "/o/6.4/NRV/370/2/burn.csv", "/o/6.4/NRV/370/2/other.csv",
          "/o/6.4/NRV/370/10/burn.csv")
  r <- fireSenseUtils:::.repOutputFiles(sf, "burn.csv")
  expect_identical(names(r), c("2", "10", "16"))
  expect_identical(unname(r), sf[c(2, 4, 1)])
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
