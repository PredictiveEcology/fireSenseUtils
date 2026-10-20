## runELFs() drops ELFs that cannot be fitted before the queue is built from its names:
## the high Arctic (ecoprovinces 1.x and 2.x) and, data-driven, ELFs with no fire in the
## fitted years (sim$ELFsExcluded, from fireSense_ELFs). ELFs with no SCANFI tree species
## (3.2.1, 3.2.4) are NOT excluded any more: the dataPrep chain carries a zero-row sppEquiv
## and speciesLayers = NULL through, and fireSense fits nonForest fuel classes only.
##
## The module run is stubbed out: what is under test is which ELF
## names come back.

local_runELFs <- function(ids, env = parent.frame(), excluded = NULL) {
  byName <- stats::setNames(as.list(ids), ids)
  sim <- list(ELFs = list(rasWhole = byName, rasCentered = byName),
              spreadFitPreRun = stats::setNames(list(ids), polygonIDTxt),
              ELFsExcluded = excluded)
  local_mocked_bindings(simInitAndSpades2 = function(...) sim,
                        .package = "SpaDES.core", .env = env)
  list(modules = "fireSense_ELFs", paths = list(modulePath = withr::local_tempdir(.local_envir = env)),
       params = list(fireSense_ELFs = list()))
}

ids <- c("1.1", "2.3", "3.1.1", "3.2.1", "3.2.2", "3.2.4", "13.2.1")
kept <- c("3.1.1", "3.2.1", "3.2.2", "3.2.4", "13.2.1") # the no-tree ELFs 3.2.1 / 3.2.4 stay

test_that("runELFs() leaves out the Arctic but keeps the no-tree ELFs in the name lists", {
  prj <- local_runELFs(ids)
  expect_identical(runELFs(prj, whatOut = "allNames"), kept)
  expect_identical(runELFs(prj, whatOut = "fittedNamesOnly"), kept)
})

test_that("runELFs() keeps merged ELFs and leaves out the ELFs fireSense_ELFs could not fit", {
  ## After ELFmergePlan(): 3.1.1 and 3.1.2 are one ELF, 3.2.2 and 3.2.4 had too few fires even
  ## together, and an Arctic merge is still Arctic.
  merged <- c("1.1_2", "3.1.1_2", "3.2.2", "3.2.4", "3.2.25", "10.1")
  prj <- local_runELFs(merged, excluded = c("3.2.2", "3.2.4"))
  expect_identical(runELFs(prj, whatOut = "allNames"), c("3.1.1_2", "3.2.25", "10.1"))
})

test_that("runELFs() leaves out the Arctic but keeps the no-tree ELFs in the maps", {
  prj <- local_runELFs(ids)
  maps <- runELFs(prj, whatOut = "maps")
  expect_identical(names(maps$rasWhole), kept)
  expect_identical(names(maps$rasCentered), kept)
})

test_that("runELFs() leaves out the Prairie grassland ELFs 10.2.*, not 10.20 or 10.1.*", {
  prj <- local_runELFs(c("10.1.1", "10.2.1", "10.2.2", "10.20", "10.3", "3.1.1"))
  expect_identical(runELFs(prj, whatOut = "allNames"), c("10.1.1", "10.20", "10.3", "3.1.1"))
  expect_identical(names(runELFs(prj, whatOut = "maps")$rasWhole), c("10.1.1", "10.20", "10.3", "3.1.1"))
})
