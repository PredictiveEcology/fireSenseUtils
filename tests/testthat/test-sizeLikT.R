## sizeLikT(): the Student-t likelihood of an observed size given simulated ones, used by the annual-area
## term with `sizeLik = "t"`.

test_that("sizeLikT() is the t density of sqrt(size), and survives identical simulated sizes", {
  N <- c(16, 25, 36, 49)
  m <- mean(sqrt(N)); s <- sd(sqrt(N))
  expect_equal(sizeLikT(30, N, df = 5), dt((sqrt(30) - m) / s, 5) / s)
  expect_true(is.finite(sizeLikT(30, c(25, 25, 25), df = 5))) # sd is 0: the scale has a lower bound
})
