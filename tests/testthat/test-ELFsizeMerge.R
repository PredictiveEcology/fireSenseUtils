## Regions smaller than minAreaKm2 merge with a similar neighbour of the same group, by an iterative rule:
## of the candidate merges, take the one that leaves the fewest small regions without a partner.

## stats table: area (km2), burn rate (%/yr) and land-cover shares (three classes) per region
sizeStats <- function(ELF, areaKm2, shares, burnRate = rep(1, length(ELF))) {
  data.table::data.table(ELF = ELF, areaKm2 = areaKm2, burnRate = burnRate,
                         landCover = lapply(shares, function(s) stats::setNames(s, c("a", "b", "c"))))
}
chain <- function(ids) neighboursOf(head(ids, -1), tail(ids, -1), rep(4000, length(ids) - 1))

test_that("the more similar pair is rejected when it orphans the third region (11.1 / 11.2 / 11.3)", {
  ## 11.2 and 11.3 are the most similar pair (0.15), but merging them leaves 11.1 (distance 0.4 to the
  ## merged region) with no partner. 11.1 + 11.2 (0.3) leaves nothing small: 11.3 is big.
  st <- sizeStats(c("11.1", "11.2", "11.3"), c(20000, 20000, 40000),
                  list(c(1, 0, 0), c(0.7, 0.3, 0), c(0.55, 0.45, 0)))
  plan <- ELFsizePlan(st, chain(st$ELF))
  expect_identical(plan$action, "merge")
  expect_identical(plan$ELF, "11.1_2")
  expect_identical(plan$members[[1]], c("11.1", "11.2"))
  expect_match(plan$reason, "^smaller than 35000 km2; merged with 11\\.[12] \\(land cover distance 0\\.30, burn ratio 1\\.0\\)$")
})

test_that("without the orphan problem the most similar pair wins", {
  ## all three small and similar: any merge leaves a small region that still has a partner
  st <- sizeStats(c("11.1", "11.2", "11.3"), c(20000, 20000, 20000),
                  list(c(1, 0, 0), c(0.7, 0.3, 0), c(0.55, 0.45, 0)))
  plan <- ELFsizePlan(st, chain(st$ELF), maxLandCoverDist = 0.5)
  ## 11.2 + 11.3 (0.15) merge first; the 40000 km2 result is not small, so 11.1 (0.4 to it) is an orphan
  ## only if it is small: it is, and 0.4 <= 0.5 gives it a partner, so it merges on the next step
  expect_identical(plan$ELF, "11.1_2_3")
  expect_identical(plan$members[[1]], c("11.1", "11.2", "11.3"))
})

test_that("no candidate: nothing changes", {
  st <- sizeStats(c("11.1", "11.2"), c(20000, 20000), list(c(1, 0, 0), c(0, 0, 1)))
  expect_identical(nrow(ELFsizePlan(st, chain(st$ELF))), 0L)                    # land cover too different
  st <- sizeStats(c("11.1", "11.2"), c(20000, 20000), list(c(1, 0, 0), c(1, 0, 0)), burnRate = c(1, 7))
  expect_identical(nrow(ELFsizePlan(st, chain(st$ELF))), 0L)                    # burn ratio 7 > 6
  expect_identical(nrow(ELFsizePlan(st, chain(st$ELF), maxBurnRatio = 8)), 1L)
  st <- sizeStats(c("11.1", "12.1"), c(20000, 20000), list(c(1, 0, 0), c(1, 0, 0)))
  expect_identical(nrow(ELFsizePlan(st, chain(st$ELF))), 0L)                    # different base
  st <- sizeStats(c("11.1", "11.2"), c(40000, 40000), list(c(1, 0, 0), c(1, 0, 0)))
  expect_identical(nrow(ELFsizePlan(st, chain(st$ELF))), 0L)                    # none small
  expect_identical(nrow(ELFsizePlan(st, chain(st$ELF), minAreaKm2 = NULL)), 0L)  # rule off
  expect_identical(nrow(ELFsizePlan(st[0, ], chain(st$ELF))), 0L)
})

test_that("a merged region still below the minimum can merge again, and independent pairs both merge", {
  st <- sizeStats(paste0("3.1.", 1:3), c(10000, 10000, 10000), rep(list(c(1, 0, 0)), 3))
  plan <- ELFsizePlan(st, chain(st$ELF))
  expect_identical(plan$ELF, "3.1.1_2_3")                    # 30000 km2 < 35000, but nothing left to merge
  expect_match(plan$reason, "; ")                            # both steps are in the reason
  ## groups are the ELF bases: 3.1 and 4.1 do not mix, though 3.1.2 and 4.1.1 touch
  st <- sizeStats(c("3.1.1", "3.1.2", "4.1.1", "4.1.2"), rep(20000, 4), rep(list(c(1, 0, 0)), 4))
  plan <- ELFsizePlan(st, chain(st$ELF))
  expect_identical(plan$ELF, c("3.1.1_2", "4.1.1_2"))
})

test_that("merged names from an earlier merge are accepted: members stay the original ids", {
  st <- sizeStats(c("3.1.1_2", "3.1.3"), c(20000, 10000), rep(list(c(1, 0, 0)), 2))
  st$members <- list(c("3.1.1", "3.1.2"), "3.1.3")
  plan <- ELFsizePlan(st, neighboursOf("3.1.1_2", "3.1.3", 4000))
  expect_identical(plan$ELF, "3.1.1_2_3")
  expect_identical(plan$members[[1]], c("3.1.1", "3.1.2", "3.1.3"))
})

test_that("ELFsizePlan output can be applied by mergeELFs()", {
  whole <- bandELFs(list(`11.1` = 1:2, `11.2` = 3:4, `11.3` = 5:6))
  st <- sizeStats(names(whole), c(20000, 20000, 40000),
                  list(c(1, 0, 0), c(0.7, 0.3, 0), c(0.55, 0.45, 0)))
  plan <- ELFsizePlan(st, ELFneighbours(whole))
  merged <- mergeELFs(list(rasWhole = stats::setNames(as.list(whole), names(whole)),
                              rasCentered = stats::setNames(as.list(whole), names(whole))), plan)
  expect_setequal(names(merged$rasWhole), c("11.1_2", "11.3"))
})

test_that("ELFregionStats gives area, land-cover shares and burn rate", {
  r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10000, ymin = 0, ymax = 10000, crs = "EPSG:3978")
  elfs <- c(terra::setValues(r, 2), terra::setValues(r, rep(c(2, 0), each = 50)))
  names(elfs) <- c("all", "half")
  ## land cover: coniferous (210) on the left half, water (20) on the right half
  lcc <- terra::setValues(r, rep(rep(c(210, 20), each = 5), 10))
  poly <- terra::vect("POLYGON ((0 0, 2000 0, 2000 1000, 0 1000, 0 0))", crs = "EPSG:3978")
  poly$YEAR <- 1985L
  st <- ELFregionStats(elfs, lcc, firePolys = poly, fireYears = 1985:1989)
  expect_identical(st$ELF, c("all", "half"))
  expect_equal(st$areaKm2, c(100, 50), tolerance = 1e-3)
  expect_equal(unname(st$landCover[[1]][c("20", "210")]), c(0.5, 0.5), tolerance = 0.02)
  expect_equal(sum(st$landCover[[1]]), 1)
  ## 2 km2 burned in 5 years, in the region "all" (100 km2): 0.4 %/yr
  expect_equal(st$burnRate[1], 100 * 2 / 100 / 5, tolerance = 1e-2)
  ## the polygon is in the bottom-left of "half" (cells 51-100 are the lower rows)? rows 1-5 are core
  expect_true(st$burnRate[2] %in% c(0, 100 * 2 / 50 / 5))
  ## deterministic
  expect_identical(ELFregionStats(elfs, lcc, firePolys = poly, fireYears = 1985:1989), st)
  ## without fire polygons the burn rate is NA
  expect_true(all(is.na(ELFregionStats(elfs, lcc)$burnRate)))
})
