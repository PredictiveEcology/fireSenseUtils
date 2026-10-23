## plotELFs(which, fill, labels, labelWhich, buffers, axes) fills the cores of the named ELFs. makeELFs() and Cache() are
## stubbed with a two-ELF polygon layer; terra::plot() records what is drawn.

## no merges and no skips
noPlan <- data.table::data.table(action = character(0), ELF = character(0), members = list(), reason = character(0))
## labels are placed with strwidth(), which needs an open device (terra::plot is mocked here)
local_plotDevice <- function(env = parent.frame()) {
  grDevices::pdf(NULL)
  withr::defer(grDevices::dev.off(), envir = env)
  graphics::plot.new()
  graphics::plot.window(c(0, 1e6), c(0, 1e6))
}

test_that("plotELFs fills the cores of the ELFs named in `which`", {
  sq <- function(x0) sprintf("POLYGON ((%1$s 0, %2$s 0, %2$s 1, %1$s 1, %1$s 0))", x0, x0 + 1)
  poly <- terra::vect(c(sq(0), sq(0), sq(2), sq(2)))
  poly$ID <- c("4.1", "4.1", "13.1", "13.1")
  poly$buffer <- c(1, 2, 1, 2)
  local_mocked_bindings(
    makeELFs = function(...) list(poly = poly),
    Cache = function(x, ...) x
  )
  ## terra::plot() is mocked, so give strwidth()/par() a (null) device with a plot region
  pdf(NULL)
  withr::defer(dev.off())
  graphics::plot.new()
  graphics::plot.window(c(0, 3), c(0, 1))
  filled <- list()
  local_mocked_bindings(
    plot = function(x, ..., col = NULL, add = FALSE) {
      if (isTRUE(add)) filled[[length(filled) + 1L]] <<- list(x = x, col = col)
    },
    text = function(...) NULL,
    .package = "terra"
  )

  local_mocked_bindings(legend = function(...) NULL, .package = "graphics")
  local_plotDevice()

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

## Four ELFs side by side, as makeELFs() returns them (rasters and `poly`). 11.1 and 11.2 merge, 11.3 is skipped,
## 1.1 is Arctic, 10.2.1 is Prairie grassland.
post <- function() {
  whole <- bandELFs(list(`11.1` = 1:2, `11.2` = 3:4, `11.3` = 5:6, `1.1` = 7:8, `10.2.1` = 9:10))
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
  local_plotDevice()
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
  expect_setequal(unique(out$plots[[1]]$x$ID), c("11.1_2", "11.3", "1.1", "10.2.1"))   # the outlines
  expect_equal(unique(out$plots[[1]]$x$buffer), 2)                           # of the cores only: no border twice
  expect_setequal(out$labels, c("11.1_2", "11.3", "1.1", "10.2.1"))
  expect_false(any(c("11.1", "11.2") %in% out$labels))
})

test_that("plotELFs(labels = 'name') names a merged ELF by its ecozone, keeping its merged code to tell it apart", {
  out <- drawn(plan = mergePlan(), labels = "name")
  expect_setequal(out$labels, c("Taiga Cordillera (11.1_2)", "Taiga Cordillera (11.3)", "Arctic Cordillera", "Prairie"))
  ## multi-level merged codes also resolve to their ecozone, never to the bare merged code
  expect_identical(elfLabels(c("3.1.1_2", "10.3.1_2"), "name"), c("Southern Arctic", "Prairie"))
})

test_that("Arctic ELFs have their own fill and a legend entry; merged ELFs are not filled, only labelled", {
  out <- drawn(plan = mergePlan(), arcticFill = "beige", fewFiresFill = "grey70")
  fillOf <- function(col) unique(unlist(lapply(Filter(function(p) identical(p$col, col), out$plots), function(p) p$x$ID)))
  expect_identical(fillOf("beige"), c("1.1", "10.2.1"))
  expect_identical(fillOf("grey70"), "11.3")
  filled <- unlist(lapply(Filter(function(p) isTRUE(p$add) && !is.null(p$col), out$plots), function(p) p$x$ID))
  expect_false("11.1_2" %in% filled)
  expect_true("11.1_2" %in% out$labels)
  expect_identical(out$legend$legend, c("not fitted (Arctic)", "not fitted (Prairie grassland)", "not enough fires"))
  expect_identical(out$legend$fill, c("beige", "beige", "grey70"))
})

test_that("plotELFs puts each label inside its ELF, even when the centroid falls outside it", {
  ## a U shape: its centroid lies in the notch, outside the ELF (as 10.2.1's does, in 10.3.1_2)
  u <- terra::vect("POLYGON ((0 0, 3 0, 3 3, 2 3, 2 1, 1 1, 1 3, 0 3, 0 0))")
  u$ID <- "10.2.1"
  u$buffer <- 2
  expect_false(terra::is.related(terra::centroids(u), u, "intersects"))
  local_mocked_bindings(makeELFs = function(...) list(poly = u), Cache = function(x, ...) x)
  pdf(NULL)
  withr::defer(dev.off())
  graphics::plot.new()
  graphics::plot.window(c(-1, 4), c(-1, 4))
  at <- NULL
  local_mocked_bindings(
    plot = function(...) NULL,
    text = function(x, labels, ...) at <<- x,
    .package = "terra"
  )
  plotELFs(plan = noPlan)
  expect_true(terra::is.related(at, u, "intersects"))
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
  expect_setequal(out$labels, c("11.1_2", "11.3", "1.1", "10.2.1"))
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
  local_plotDevice()
  plotELFs(plan = mergePlan(), ELFs = ELFs)
  expect_false(made)
  expect_setequal(calls$labels, c("11.1_2", "11.3", "1.1", "10.2.1"))
})

test_that("plotELFs rejects values outside the allowed options", {
  expect_error(plotELFs(labels = "ecoregion", plan = noPlan), "should be one of")
  expect_error(plotELFs(labelWhich = "some", plan = noPlan), "should be one of")
  expect_error(plotELFs(axes = "km", plan = noPlan), "should be one of")
})

test_that("elfLayer() drops the buffer rings when buffers = FALSE", {
  sq <- function(x0) sprintf("POLYGON ((%1$s 0, %2$s 0, %2$s 1, %1$s 1, %1$s 0))", x0, x0 + 1)
  poly <- terra::vect(c(sq(0), sq(0), sq(2), sq(2)))
  poly$ID <- c("4.1", "4.1", "13.1", "13.1")
  poly$buffer <- c(1, 2, 1, 2)
  expect_equal(NROW(elfLayer(poly, buffers = TRUE)), 4)
  cores <- elfLayer(poly, buffers = FALSE)
  expect_identical(cores$buffer, c(2, 2))
  expect_identical(cores$ID, c("4.1", "13.1"))
})

test_that("elfLabels() gives ecozone names for 'name', not codes", {
  ids <- c("13.1", "5.4", "14.3", "6.2.1")
  expect_identical(elfLabels(ids, "name"),
                   c("Pacific Maritime", "Taiga Shield", "Montane Cordillera", "Boreal Shield"))
  expect_identical(elfLabels(ids, "code"), ids)
  expect_identical(elfLabels(ids, "none"), rep("", 4L))
  ## two ELFs in one zone stay distinguishable
  expect_identical(elfLabels(c("6.2.1", "6.2.2", "5.4"), "name"),
                   c("Boreal Shield (6.2.1)", "Boreal Shield (6.2.2)", "Taiga Shield"))
  ## an unknown zone falls back to the code
  expect_identical(elfLabels("99.1", "name"), "99.1")
  expect_error(elfLabels("13.1", "zone"), "should be one of")
})

test_that("elfLabelXY() keeps separate labels put and moves overlapping ones apart", {
  ## far apart: untouched
  far <- elfLabelXY(c(0, 10), c(0, 0), w = c(2, 2), h = c(1, 1))
  expect_identical(far$x, c(0, 10)); expect_identical(far$y, c(0, 0))
  ## same spot: the second moves, and no two label boxes overlap afterwards
  near <- elfLabelXY(c(0, 0.5), c(0, 0), w = c(4, 4), h = c(1, 1))
  expect_identical(c(near$x[1], near$y[1]), c(0, 0))
  expect_true(abs(near$x[2] - near$x[1]) >= 4 || abs(near$y[2] - near$y[1]) >= 1)
  ## and a label never leaves the plot region
  edge <- elfLabelXY(c(0, 0.5), c(0, 0), w = c(4, 4), h = c(1, 1), usr = c(-2, 3, -3, 3))
  expect_true(all(edge$x - 2 >= -2 & edge$x + 2 <= 3 & edge$y - 0.5 >= -3 & edge$y + 0.5 <= 3))
})

test_that("plotELFs labels only the highlighted ELFs, with names, when asked", {
  sq <- function(x0) sprintf("POLYGON ((%1$s 0, %2$s 0, %2$s 1, %1$s 1, %1$s 0))", x0, x0 + 1)
  poly <- terra::vect(c(sq(0), sq(0), sq(2), sq(2)))
  poly$ID <- c("4.1", "4.1", "13.1", "13.1")
  poly$buffer <- c(1, 2, 1, 2)
  local_mocked_bindings(
    makeELFs = function(...) list(poly = poly),
    Cache = function(x, ...) x
  )
  pdf(NULL)
  withr::defer(dev.off())
  graphics::plot.new()
  graphics::plot.window(c(0, 3), c(0, 1))
  shown <- NULL
  drawn <- NULL
  local_mocked_bindings(
    plot = function(x, ..., add = FALSE) if (!isTRUE(add)) drawn <<- x,
    text = function(x, labels, ...) shown <<- labels,
    .package = "terra"
  )
  plotELFs(which = "13.1", labels = "name", labelWhich = "highlighted", buffers = FALSE, axes = "none", plan = noPlan)
  expect_identical(shown, "Pacific Maritime")
  expect_identical(drawn$buffer, c(2, 2))

  plotELFs(labels = "code", plan = noPlan)
  expect_identical(shown, c("4.1", "13.1"))
  expect_equal(NROW(drawn), 4)

  shown <- NULL
  plotELFs(which = "13.1", labels = "none", plan = noPlan)
  expect_null(shown)
})

test_that("plotELFs(axes = 'none' / 'longlat') switches off the map axes; only 'longlat' adds a graticule", {
  sq <- function(x0) sprintf("POLYGON ((%1$s 0, %2$s 0, %2$s 1, %1$s 1, %1$s 0))", x0, x0 + 1)
  poly <- terra::vect(c(sq(0), sq(2)))
  poly$ID <- c("4.1", "13.1")
  poly$buffer <- c(2, 2)
  gratFor <- NULL
  local_mocked_bindings(
    makeELFs = function(...) list(poly = poly),
    Cache = function(x, ...) x,
    elfGraticule = function(v, ...) gratFor <<- v
  )
  pdf(NULL)
  withr::defer(dev.off())
  graphics::plot.new()
  graphics::plot.window(c(0, 3), c(0, 1))
  plotArgs <- NULL
  local_mocked_bindings(
    plot = function(x, ..., add = FALSE) if (!isTRUE(add)) plotArgs <<- list(...),
    text = function(...) NULL,
    .package = "terra"
  )
  plotELFs(axes = "m", plan = noPlan)
  expect_null(plotArgs$axes)
  expect_null(gratFor)

  plotELFs(axes = "none", plan = noPlan)
  expect_false(plotArgs$axes)
  expect_false(plotArgs$box)
  expect_null(gratFor)

  plotELFs(axes = "longlat", plan = noPlan)
  expect_false(plotArgs$axes)
  expect_identical(plotArgs$mar, c(3, 3, 1, 1))
  expect_equal(NROW(gratFor), 2)
})

test_that("plotELFs widens the extent only when it labels the highlighted ELFs", {
  sq <- function(x0) sprintf("POLYGON ((%1$s 0, %2$s 0, %2$s 1, %1$s 1, %1$s 0))", x0, x0 + 1)
  poly <- terra::vect(c(sq(0), sq(2)))
  poly$ID <- c("4.1", "13.1")
  poly$buffer <- c(2, 2)
  local_mocked_bindings(makeELFs = function(...) list(poly = poly), Cache = function(x, ...) x)
  pdf(NULL)
  withr::defer(dev.off())
  graphics::plot.new()
  graphics::plot.window(c(0, 3), c(0, 1))
  plotArgs <- NULL
  local_mocked_bindings(
    plot = function(x, ..., add = FALSE) if (!isTRUE(add)) plotArgs <<- list(...),
    text = function(...) NULL,
    .package = "terra"
  )
  plotELFs(which = "4.1", labelWhich = "highlighted", plan = noPlan)
  expect_equal(unname(as.vector(plotArgs$ext)), c(-0.24, 3.24, -0.24, 1.24))
  plotELFs(which = "4.1", labelWhich = "highlighted", labels = "none", plan = noPlan)
  expect_null(plotArgs$ext)
  plotELFs(which = "4.1", labelWhich = "all", plan = noPlan)
  expect_null(plotArgs$ext)
})

test_that("plotELFs draws a leader line for a label moved off its ELF", {
  sq <- function(x0) sprintf("POLYGON ((%1$s 0, %2$s 0, %2$s 1, %1$s 1, %1$s 0))", x0, x0 + 0.1)
  poly <- terra::vect(c(sq(0), sq(0.05)))
  poly$ID <- c("4.1", "13.1")
  poly$buffer <- c(2, 2)
  local_mocked_bindings(makeELFs = function(...) list(poly = poly), Cache = function(x, ...) x)
  pdf(NULL)
  withr::defer(dev.off())
  graphics::plot.new()
  graphics::plot.window(c(-1, 2), c(-1, 2))
  segs <- NULL
  shown <- NULL
  local_mocked_bindings(segments = function(x0, y0, x1, y1, ...) segs <<- list(x0, y0, x1, y1),
                        .package = "graphics")
  local_mocked_bindings(
    plot = function(...) NULL,
    text = function(x, labels, ...) shown <<- list(xy = terra::crds(x), labels = labels),
    .package = "terra"
  )
  plotELFs(plan = noPlan)
  expect_identical(shown$labels, c("4.1", "13.1"))
  expect_length(segs[[1]], 1L) # one label moved, from its centroid ...
  expect_equal(unname(c(segs[[1]], segs[[2]])), c(0.1, 0.5))
  expect_false(isTRUE(all.equal(c(segs[[3]], segs[[4]]), c(segs[[1]], segs[[2]])))) # ... to a new spot
  expect_equal(unname(shown$xy[2, ]), c(segs[[3]], segs[[4]]))
})

test_that("elfGraticule() draws graticule lines and labels the edges in degrees", {
  ## a 4 x 4 degree-ish block in Canada Lambert (EPSG:3978), well inside 10-degree lines
  v <- terra::vect("POLYGON ((-2000000 -500000, 1500000 -500000, 1500000 2500000, -2000000 2500000, -2000000 -500000))",
                   crs = "EPSG:3978")
  pdf(NULL)
  withr::defer(dev.off())
  terra::plot(v, axes = FALSE)
  lines <- 0L
  axes <- list()
  local_mocked_bindings(
    plot = function(x, ..., add = FALSE) if (isTRUE(add)) lines <<- lines + 1L,
    .package = "terra"
  )
  local_mocked_bindings(
    axis = function(side, at, labels, ...) axes[[as.character(side)]] <<- list(at = at, labels = labels),
    .package = "graphics"
  )
  expect_null(elfGraticule(v))
  expect_gt(lines, 0L)
  ## meridians (side 1) and parallels (side 2) are labelled with N/S/E/W suffixes
  expect_true(all(grepl("^[0-9]+°[EW]$", axes[["1"]]$labels)))
  expect_true(all(grepl("^[0-9]+°[NS]$", axes[["2"]]$labels)))
  expect_length(axes[["1"]]$at, length(axes[["1"]]$labels))
  ## a label sits on the edge it names: parallels at the left edge x, inside the map's y range
  e <- terra::ext(v)
  expect_true(all(axes[["2"]]$at >= e$ymin & axes[["2"]]$at <= e$ymax))
  expect_true(all(axes[["1"]]$at >= e$xmin & axes[["1"]]$at <= e$xmax))
})

test_that("elfGraticule() copes with a map smaller than one graticule step and with lines that miss the map", {
  ## about 100 km square near the projection's central meridian: no multiple of 10 degrees lies inside it
  v <- terra::vect("POLYGON ((0 1000000, 100000 1000000, 100000 1100000, 0 1100000, 0 1000000))",
                   crs = "EPSG:3978")
  pdf(NULL)
  withr::defer(dev.off())
  terra::plot(v, axes = FALSE)
  axes <- 0L
  lines <- 0L
  local_mocked_bindings(axis = function(...) axes <<- axes + 1L, .package = "graphics")
  local_mocked_bindings(plot = function(x, ..., add = FALSE) lines <<- lines + 1L, .package = "terra")
  expect_null(elfGraticule(v))
  expect_identical(c(lines, axes), c(0L, 0L))

  ## a graticule line that does not reach the map edge gets no label
  big <- terra::vect("POLYGON ((-2000000 -500000, 1500000 -500000, 1500000 2500000, -2000000 2500000, -2000000 -500000))",
                     crs = "EPSG:3978")
  terra::plot(big, axes = FALSE)
  axes <- 0L
  local_mocked_bindings(crop = function(x, y, ...) terra::vect(), .package = "terra")
  expect_null(elfGraticule(big))
  expect_identical(axes, 0L)
})

test_that("plotELFs(metrics, metric) fills the cores that have a value, with their values", {
  m <- data.frame(ELF = c("11.1_2", "1.1"), areaBurnedHaPerYear = c(250, NA))
  out <- drawn(plan = mergePlan(), metrics = m, metric = "areaBurnedHaPerYear")
  valued <- Filter(function(p) "areaBurnedHaPerYear" %in% names(p$x), out$plots)
  expect_length(valued, 1L)
  expect_identical(valued[[1]]$x$ID, "11.1_2")                  # 1.1 has no value: it keeps its Arctic fill
  expect_equal(valued[[1]]$x$areaBurnedHaPerYear, 250)
  expect_true(valued[[1]]$add)
})

test_that("plotELFs(metrics) needs one row per ELF and a `metric` column, and warns of unknown ELFs", {
  expect_error(drawn(plan = mergePlan(), metrics = data.frame(ELF = c("11.3", "11.3"), v = 1:2), metric = "v"),
               "more than one row for ELF 11.3")
  expect_error(drawn(plan = mergePlan(), metrics = data.frame(ELF = "11.3", v = 1)), "`metric` must name")
  ## a pre-merge name is not on the map
  expect_warning(drawn(plan = mergePlan(), metrics = data.frame(ELF = c("11.1", "11.3"), v = 1:2), metric = "v"),
                 "ELFs in `metrics` not found: 11.1")
})

test_that("plotELFs(plotFires) fills the cores with ELFfireRates() on a green-to-red ramp", {
  seen <- NULL
  local_mocked_bindings(ELFfireRates = function(ELFs, plan, fireYears, pixelAreaHa, destinationPath) {
    seen <<- list(ELFs = names(ELFs$rasWhole), fireYears = fireYears, pixelAreaHa = pixelAreaHa)
    data.table::data.table(ELF = c("11.1_2", "11.3"), naturalFiresPerYear = c(2, 1), haBurnedPerYear = c(100, 10))
  })
  out <- drawn(plan = mergePlan(), plotFires = "area", fireYears = 2001:2003)
  expect_setequal(seen$ELFs, c("11.1", "11.2", "11.3", "1.1", "10.2.1"))   # before the merges
  expect_identical(seen$fireYears, 2001:2003)
  expect_identical(seen$pixelAreaHa, 5.76)
  valued <- Filter(function(p) "haBurnedPerYear" %in% names(p$x), out$plots)
  expect_length(valued, 1L)
  expect_equal(stats::setNames(valued[[1]]$x$haBurnedPerYear, valued[[1]]$x$ID)[c("11.1_2", "11.3")],
               c(`11.1_2` = 100, `11.3` = 10))
  ramp <- valued[[1]]$col
  expect_identical(ramp, grDevices::hcl.colors(100, "RdYlGn", rev = TRUE))
  rgbOf <- function(x) grDevices::col2rgb(x)[, 1]
  expect_gt(rgbOf(ramp[1])[["green"]], rgbOf(ramp[1])[["red"]])        # less: green
  expect_gt(rgbOf(ramp[100])[["red"]], rgbOf(ramp[100])[["green"]])    # more: red

  expect_false("not enough fires" %in% out$legend$legend)            # the fire values fill those ELFs

  out <- drawn(plan = mergePlan(), plotFires = "number")
  expect_length(Filter(function(p) "naturalFiresPerYear" %in% names(p$x), out$plots), 1L)
})

test_that("plotELFs(plotFires) takes 'number' or 'area', and not with `metrics`", {
  expect_error(drawn(plan = mergePlan(), plotFires = "count"), "should be one of")
  expect_error(drawn(plan = mergePlan(), plotFires = "area", metrics = data.frame(ELF = "11.3", v = 1), metric = "v"),
               "not both")
})
