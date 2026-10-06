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
