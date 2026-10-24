## sizeLikT(): a Student-t density of sqrt(size), centred on the mean of the simulated sizes with their
## standard deviation as its scale. It was the objective's per-fire likelihood with `sizeLik = "t"`; the
## objective now uses kernelLogScore(), and sizeLikT() is kept for yearAreaNLL(), a diagnostic.

test_that("sizeLikT() is the t density of sqrt(size), and survives identical simulated sizes", {
  N <- c(16, 25, 36, 49)
  m <- mean(sqrt(N)); s <- sd(sqrt(N))
  expect_equal(sizeLikT(30, N, df = 5), dt((sqrt(30) - m) / s, 5) / s)
  expect_true(is.finite(sizeLikT(30, c(25, 25, 25), df = 5))) # sd is 0: the scale has a lower bound
})
