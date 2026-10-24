## `adWeight` multiplies the Anderson-Darling statistic before it is added to the fire-size SNLL.
## The fixed 50 balanced the two terms only for the kernel-density likelihood with no size weight:
## measured across six ELFs, the AD term's share of the influence on the objective ranged from 0.16
## to 0.94 over the four sizeLik x weighted combinations. `adWeight = "auto"` scales with the number
## of fitted fires, which is what both terms grow with, and holds that share to 0.36-0.68.
## The objective no longer has the Anderson-Darling term (it depended on the buffers); adWeightAuto() stays
## for diagnostics.

test_that("adWeightAuto() is c * sqrt(nFires), with c from sizeLik and weighted", {
  expect_equal(adWeightAuto(100, "kde", FALSE), 1.897 * 10)
  expect_equal(adWeightAuto(100, "kde", "sqrt"), 5.098 * 10)
  expect_equal(adWeightAuto(100, "t", FALSE), 0.282 * 10)
  expect_equal(adWeightAuto(100, "t", "sqrt"), 0.804 * 10)
  ## TRUE means the log(size) weight, which was measured with the unweighted constants
  expect_equal(adWeightAuto(400, "kde", TRUE), adWeightAuto(400, "kde", FALSE))
  expect_gt(adWeightAuto(1624, "kde", FALSE), adWeightAuto(138, "kde", FALSE))
})
