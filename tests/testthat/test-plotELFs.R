## plotELFs(which, fill, labels, labelWhich, buffers, axes) fills the cores of the named ELFs. makeELFs() and Cache() are
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

  plotELFs(which = "13.1", fill = "red")
  expect_length(filled, 1L)
  expect_identical(filled[[1]]$x$ID, "13.1")
  expect_identical(filled[[1]]$x$buffer, 2)
  expect_identical(filled[[1]]$col, "red")

  filled <- list()
  expect_warning(plotELFs(which = c("4.1", "99.9")), "ELFs not found: 99.9")
  expect_identical(filled[[1]]$x$ID, "4.1")

  filled <- list()
  plotELFs()
  expect_length(filled, 0L)
})

test_that("plotELFs rejects values outside the allowed options", {
  expect_error(plotELFs(labels = "ecoregion"), "should be one of")
  expect_error(plotELFs(labelWhich = "some"), "should be one of")
  expect_error(plotELFs(axes = "km"), "should be one of")
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
  plotELFs(which = "13.1", labels = "name", labelWhich = "highlighted", buffers = FALSE, axes = "none")
  expect_identical(shown, "Pacific Maritime")
  expect_identical(drawn$buffer, c(2, 2))

  plotELFs(labels = "code")
  expect_identical(shown, c("4.1", "13.1"))
  expect_equal(NROW(drawn), 4)

  shown <- NULL
  plotELFs(which = "13.1", labels = "none")
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
  plotELFs(axes = "m")
  expect_null(plotArgs$axes)
  expect_null(gratFor)

  plotELFs(axes = "none")
  expect_false(plotArgs$axes)
  expect_false(plotArgs$box)
  expect_null(gratFor)

  plotELFs(axes = "longlat")
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
  plotELFs(which = "4.1", labelWhich = "highlighted")
  expect_equal(unname(as.vector(plotArgs$ext)), c(-0.24, 3.24, -0.24, 1.24))
  plotELFs(which = "4.1", labelWhich = "highlighted", labels = "none")
  expect_null(plotArgs$ext)
  plotELFs(which = "4.1", labelWhich = "all")
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
  plotELFs()
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
