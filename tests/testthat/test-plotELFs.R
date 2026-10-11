## plotELFs(which, fill) fills the cores of the named ELFs. makeELFs() and Cache() are
## stubbed with a two-ELF polygon layer; terra::plot() records what is drawn.

test_that("plotELFs fills the cores of the ELFs named in `which`", {
  sq <- function(x0) sprintf("POLYGON ((%1$s 0, %2$s 0, %2$s 1, %1$s 1, %1$s 0))", x0, x0 + 1)
  poly <- terra::vect(c(sq(0), sq(0), sq(2), sq(2)))
  poly$ID <- c("4.1", "4.1", "13.1", "13.1")
  poly$buffer <- c(1, 2, 1, 2)
  local_mocked_bindings(
    makeELFs = function(...) list(poly = poly),
    Cache = function(x, ...) x
  )
  filled <- list()
  local_mocked_bindings(
    plot = function(x, ..., col = NULL, add = FALSE) {
      if (isTRUE(add)) filled[[length(filled) + 1L]] <<- list(x = x, col = col)
    },
    text = function(...) NULL,
    .package = "terra"
  )

  noPlan <- data.table::data.table(action = character(0), ELF = character(0), members = list(), reason = character(0))
  local_mocked_bindings(legend = function(...) NULL, .package = "graphics")

  plotELFs(which = "13.1", fill = "red", plan = noPlan)
  expect_length(filled, 1L)
  expect_identical(filled[[1]]$x$ID, "13.1")
  expect_identical(filled[[1]]$x$buffer, 2)
  expect_identical(filled[[1]]$col, "red")

  filled <- list()
  expect_warning(plotELFs(which = c("4.1", "99.9"), plan = noPlan), "ELFs not found: 99.9")
  expect_identical(filled[[1]]$x$ID, "4.1")

  filled <- list()
  plotELFs(plan = noPlan)
  expect_length(filled, 0L)
})

## Three ELFs side by side, as makeELFs() returns them (rasters and `poly`). 11.1 and 11.2 merge, 11.3 is skipped.
post <- function() {
  whole <- bandELFs(list(`11.1` = 1:2, `11.2` = 3:4, `11.3` = 5:6))
  rasWhole <- stats::setNames(as.list(whole), names(whole))
  poly <- .ELFoutlines(rasWhole)
  list(rasWhole = rasWhole, rasCentered = rasWhole, poly = poly)
}
mergePlan <- function() {
  data.table::data.table(
    action = c("skip", "merge"), ELF = c(NA, "11.1_2"), members = list("11.3", c("11.1", "11.2")),
    naturalIgnitions = 1L, escapes = 1L, firePolygons = 1L, reason = c("too few fires", "smaller than 35000 km2"))
}
drawn <- function(...) {
  calls <- list(plots = list(), labels = character(0), legend = NULL)
  local_mocked_bindings(
    makeELFs = function(...) post(),
    Cache = function(x, ...) x
  )
  local_mocked_bindings(
    plot = function(x, ..., col = NULL, add = FALSE) calls$plots[[length(calls$plots) + 1L]] <<- list(x = x, col = col, add = add),
    text = function(x, labels, ...) calls$labels <<- c(calls$labels, labels),
    .package = "terra"
  )
  local_mocked_bindings(legend = function(x, legend, fill, ...) calls$legend <<- list(legend = legend, fill = fill),
                        .package = "graphics")
  plotELFs(...)
  calls
}

test_that("plotELFs draws the post-merge state: the merged ELF in place of its members", {
  out <- drawn(plan = mergePlan())
  expect_setequal(unique(out$plots[[1]]$x$ID), c("11.1_2", "11.3"))        # the outlines
  expect_setequal(out$labels, c("11.1_2", "11.3"))
  expect_false(any(c("11.1", "11.2") %in% out$labels))
})

test_that("the ELFs without enough fires are shaded: by default the plan's skips, or those given", {
  out <- drawn(plan = mergePlan(), fewFiresFill = "grey70")
  shaded <- Filter(function(p) identical(p$col, "grey70"), out$plots)
  expect_length(shaded, 1L)
  expect_identical(shaded[[1]]$x$ID, "11.3")                                  # the skipped ELF, core only
  expect_equal(shaded[[1]]$x$buffer, 2)
  expect_true("not enough fires" %in% out$legend$legend)
  expect_true("grey70" %in% out$legend$fill)

  out <- drawn(plan = mergePlan(), fewFires = "11.1_2", fewFiresFill = "grey70")  # explicit: overrides the skips
  shaded <- Filter(function(p) identical(p$col, "grey70"), out$plots)
  expect_length(shaded, 1L)
  expect_identical(shaded[[1]]$x$ID, "11.1_2")

  out <- drawn(plan = mergePlan(), fewFires = character(0), fewFiresFill = "grey70")  # none
  expect_length(Filter(function(p) identical(p$col, "grey70"), out$plots), 0L)
  expect_false("not enough fires" %in% out$legend$legend)
})

test_that("`which` is in the legend when it is used", {
  out <- drawn(plan = mergePlan(), which = "11.3", fill = "red", fewFires = character(0))
  expect_true("red" %in% out$legend$fill)
})

test_that("without a plan, plotELFs computes one with fewFireELFs()", {
  seen <- NULL
  local_mocked_bindings(
    fewFireELFs = function(ELFs, ...) { seen <<- list(...); list(ELFs = ELFs, plan = mergePlan()) },
    sizeRuleLandCover = function(...) "lcc.tif"
  )
  out <- drawn()
  expect_identical(seen$fireYears, defaultFireYears())
  expect_setequal(out$labels, c("11.1_2", "11.3"))
})

test_that("a given ELFs list is used, drawn from rasWhole when it has no outlines", {
  made <- FALSE
  local_mocked_bindings(makeELFs = function(...) { made <<- TRUE; post() }, Cache = function(x, ...) x)
  ELFs <- post()
  ELFs$poly <- NULL
  calls <- list(plots = list(), labels = character(0))
  local_mocked_bindings(
    plot = function(x, ..., col = NULL, add = FALSE) calls$plots[[length(calls$plots) + 1L]] <<- list(x = x, col = col),
    text = function(x, labels, ...) calls$labels <<- c(calls$labels, labels),
    .package = "terra")
  local_mocked_bindings(legend = function(...) NULL, .package = "graphics")
  plotELFs(plan = mergePlan(), ELFs = ELFs)
  expect_false(made)
  expect_setequal(calls$labels, c("11.1_2", "11.3"))
})
