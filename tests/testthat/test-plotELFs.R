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
