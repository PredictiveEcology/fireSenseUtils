## runELFs() must find fireSense_ELFs whether the project lists it directly or lists a parent
## module (e.g. PredictiveEcology/fireSense@development) whose children include it.

makeToyModulePath <- function(env = parent.frame()) {
  mp <- withr::local_tempdir(.local_envir = env)
  suppressMessages({
    SpaDES.core::newModule("fam", mp, type = "parent",
                           children = c("fireSense_ELFs", "other"), open = FALSE)
    SpaDES.core::newModule("fireSense_ELFs", mp, open = FALSE)
    SpaDES.core::newModule("other", mp, open = FALSE)
    SpaDES.core::newModule("unrelated", mp, open = FALSE)
  })
  mp
}

test_that(".findModuleInProject sees through a listed parent module", {
  mp <- makeToyModulePath()
  expect_identical(
    .findModuleInProject("PredictiveEcology/fam@development", mp, pattern = "^fireSense_ELFs$"),
    "fireSense_ELFs")
  expect_identical(
    .findModuleInProject(c("PredictiveEcology/unrelated@development", "fam"), mp,
                         pattern = "^fireSense_ELFs$"),
    "fireSense_ELFs")
})

test_that(".findModuleInProject finds a module listed directly, spec or plain name", {
  mp <- makeToyModulePath()
  expect_identical(.findModuleInProject(c("unrelated", "fireSense_ELFs"), mp, "^fireSense_ELFs$"),
                   "fireSense_ELFs")
  expect_identical(.findModuleInProject("PredictiveEcology/fireSense_ELFs@development", mp,
                                        "^fireSense_ELFs$"),
                   "fireSense_ELFs")
})

test_that(".findModuleInProject stops, naming what was searched, when nothing matches", {
  mp <- makeToyModulePath()
  expect_error(.findModuleInProject(c("PredictiveEcology/unrelated@development", "other"), mp,
                                    "^fireSense_ELFs$"),
               "fireSense_ELFs.*unrelated.*other")
})

test_that(".findModuleInProject survives a parent that lists itself", {
  mp <- withr::local_tempdir()
  suppressMessages(SpaDES.core::newModule("loop", mp, type = "parent", children = "loop",
                                          open = FALSE))
  expect_error(.findModuleInProject("loop", mp, "^fireSense_ELFs$"), "fireSense_ELFs")
})

test_that("runELFs stops clearly when the project has no ELFs module", {
  mp <- makeToyModulePath()
  prj <- list(modules = "PredictiveEcology/unrelated@development", paths = list(modulePath = mp))
  expect_error(runELFs(prj), "fireSense_ELFs.*unrelated")
})

test_that(".elfOutputFile stops before upload when the ELF file is missing", {
  expect_error(.elfOutputFile(c("a/checkpoint.rds", "a/progress.png")), "No ELF output file")
  expect_identical(.elfOutputFile(c("a/checkpoint.rds", "a/ELFs.rds")), "a/ELFs.rds")
})
