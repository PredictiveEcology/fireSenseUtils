## runELFs() drops ELFs that cannot be fitted before the queue is built from its names:
## the high Arctic (ecoprovinces 1.x and 2.x) and, data-driven, ELFs with no fire in the
## fitted years (sim$ELFsExcluded, from fireSense_ELFs). ELFs with no SCANFI tree species
## (3.2.1, 3.2.4) are NOT excluded any more: the dataPrep chain carries a zero-row sppEquiv
## and speciesLayers = NULL through, and fireSense fits nonForest fuel classes only.
##
## The module run and the cloud calls are stubbed out: what is under test is which ELF
## names come back.

local_runELFs <- function(ids, env = parent.frame()) {
  byName <- stats::setNames(as.list(ids), ids)
  sim <- list(ELFs = list(rasWhole = byName, rasCentered = byName),
              spreadFitPreRun = stats::setNames(list(ids), polygonIDTxt))
  local_mocked_bindings(Cache = function(x, ...) x, .env = env)
  local_mocked_bindings(simInitAndSpades2 = function(...) sim,
                        .package = "SpaDES.core", .env = env)
  local_mocked_bindings(getRemoteMetadata = function(...) list(remoteHash = "hash"),
                        .package = "reproducible", .env = env)
  local_mocked_bindings(user = function(...) "notTheUploader",
                        .package = "SpaDES.project", .env = env)
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

test_that("runELFs() leaves out the Arctic but keeps the no-tree ELFs in the maps", {
  prj <- local_runELFs(ids)
  maps <- runELFs(prj, whatOut = "maps")
  expect_identical(names(maps$rasWhole), kept)
  expect_identical(names(maps$rasCentered), kept)
})
