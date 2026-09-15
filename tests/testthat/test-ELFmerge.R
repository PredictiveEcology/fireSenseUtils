## ELFs with too few fires are merged with a neighbour that shares their base, or left out
## (Eliot, 2026-09-14). These pin the decisions on small synthetic maps and tables.

## ELFs as vertical bands of core cells on one grid (1000 m cells), optionally with a buffer band.
bandELFs <- function(bands, nrow = 4L, buffers = list()) {
  ncol <- max(unlist(c(bands, buffers)))
  r <- terra::rast(nrows = nrow, ncols = ncol, xmin = 0, xmax = ncol * 1000,
                   ymin = 0, ymax = nrow * 1000, crs = "EPSG:3978")
  layers <- Map(cols = bands, nam = names(bands), function(cols, nam) {
    v <- matrix(0L, nrow, ncol)
    if (!is.null(buffers[[nam]])) v[, buffers[[nam]]] <- 1L
    v[, cols] <- 2L
    terra::setValues(r, as.vector(t(v)))
  })
  out <- terra::rast(unname(layers))
  names(out) <- names(bands)
  out
}

statusOf <- function(ELF, naturalIgnitions, firePolygons) {
  data.table::data.table(ELF = ELF, naturalIgnitions = naturalIgnitions, firePolygons = firePolygons,
                         status = ifelse(naturalIgnitions == 0 | firePolygons == 0, "zero",
                                         ifelse(naturalIgnitions < 50 | firePolygons < 50, "few", "ok")))
}

neighboursOf <- function(ELF1, ELF2, sharedLength) {
  data.table::data.table(ELF1 = ELF1, ELF2 = ELF2, sharedEdges = sharedLength / 1000,
                         sharedLength = sharedLength)
}

test_that("ELFneighbours measures the core border each pair of ELFs shares", {
  elfs <- bandELFs(list(`3.1.1` = 1:4, `3.1.2` = 5:8, `3.1.3` = 9:12))
  nb <- ELFneighbours(elfs)
  expect_setequal(paste(nb$ELF1, nb$ELF2), c("3.1.1 3.1.2", "3.1.2 3.1.3"))
  expect_true(all(nb$sharedEdges == 4L))   # 4 rows of adjacent cells
  expect_true(all(nb$sharedLength == 4000))
  ## a list of layers, as sim$ELFs$rasWhole holds them, gives the same answer
  expect_identical(ELFneighbours(as.list(elfs)), nb)
})

test_that("a buffer is not a border: only core cells count", {
  elfs <- bandELFs(list(`3.1.1` = 1:4, `3.1.2` = 7:10), buffers = list(`3.1.1` = 5:6))
  expect_identical(nrow(ELFneighbours(elfs)), 0L)
})

test_that("a thin piece merges with its sibling when together they have enough fires", {
  plan <- ELFmergePlan(statusOf(c("3.1.1", "3.1.2"), c(60, 10), c(60, 10)),
                       neighboursOf("3.1.1", "3.1.2", 4000))
  expect_identical(plan$action, "merge")
  expect_identical(plan$ELF, "3.1.1_2")
  expect_identical(plan$members[[1]], c("3.1.1", "3.1.2"))
  expect_identical(plan$naturalIgnitions, 70)
  expect_identical(ELFsSkipped(plan), character(0))
})

test_that("if the pair still has too few fires, neither is fitted", {
  plan <- ELFmergePlan(statusOf(c("3.1.1", "3.1.2"), c(20, 10), c(20, 10)),
                       neighboursOf("3.1.1", "3.1.2", 4000))
  expect_identical(plan$action, "skip")
  expect_true(is.na(plan$ELF))
  expect_setequal(ELFsSkipped(plan), c("3.1.1", "3.1.2"))
})

test_that("with several siblings, the partner is the one sharing the longest border", {
  status <- statusOf(c("3.2.1", "3.2.2", "3.2.4"), c(10, 100, 100), c(10, 100, 100))
  nb <- neighboursOf(c("3.2.1", "3.2.1"), c("3.2.2", "3.2.4"), c(2000, 5000))
  plan <- ELFmergePlan(status, nb)
  expect_identical(plan$ELF, "3.2.1_4")
})

test_that("a whole province merges only with another whole province of its ecozone", {
  status <- statusOf(c("12.1", "12.2", "12.3.1", "13.1"), c(10, 100, 100, 100), c(10, 100, 100, 100))
  nb <- neighboursOf(c("12.1", "12.1", "12.1"), c("12.2", "12.3.1", "13.1"), c(1000, 9000, 9000))
  plan <- ELFmergePlan(status, nb)
  expect_identical(plan$ELF, "12.1_2")
})

test_that("an ELF with no neighbour sharing its base is not fitted", {
  plan <- ELFmergePlan(statusOf(c("4.3", "5.1"), c(10, 100), c(10, 100)),
                       neighboursOf("4.3", "5.1", 7000))
  expect_identical(plan$action, "skip")
  expect_identical(ELFsSkipped(plan), "4.3")
})

test_that("an ELF whose core touches no other ELF is not fitted", {
  ## 4.3 has no row in neighbours at all, like an island or an ELF bordered only by the Arctic
  plan <- ELFmergePlan(statusOf(c("4.3", "5.1", "5.2"), c(10, 10, 100), c(10, 10, 100)),
                       neighboursOf("5.1", "5.2", 7000))
  expect_identical(plan$action, c("skip", "merge"))
  expect_identical(plan$reason[1], "too few fires; no neighbour shares its base")
  expect_identical(ELFsSkipped(plan), "4.3")
  ## and with no neighbours anywhere
  plan <- ELFmergePlan(statusOf("4.3", 10, 10), neighboursOf(character(0), character(0), numeric(0)))
  expect_identical(ELFsSkipped(plan), "4.3")
})

test_that("an ELF is part of at most one merge", {
  ## 3.1.1 and 3.1.3 are thin and both border only 3.1.2; 3.1.1 comes first and takes it
  status <- statusOf(c("3.1.1", "3.1.2", "3.1.3"), c(10, 100, 10), c(10, 100, 10))
  nb <- neighboursOf(c("3.1.1", "3.1.2"), c("3.1.2", "3.1.3"), c(4000, 4000))
  plan <- ELFmergePlan(status, nb)
  expect_identical(plan$action, c("merge", "skip"))
  expect_identical(plan$ELF[1], "3.1.1_2")
  expect_identical(ELFsSkipped(plan), "3.1.3")
})

test_that("ELFs with enough fires are left alone, and no thin ELF gives an empty plan", {
  plan <- ELFmergePlan(statusOf(c("6.1", "6.2"), c(100, 100), c(100, 100)),
                       neighboursOf("6.1", "6.2", 4000))
  expect_identical(nrow(plan), 0L)
  expect_identical(ELFsSkipped(plan), character(0))
})

test_that("ELFs of ecozones 1 and 2 are left out permanently: never merged, never partners", {
  expect_identical(ELFsArctic(c("1.1", "2.3.1", "10.1", "12.1", "3.1.1")), c("1.1", "2.3.1"))
  plan <- ELFmergePlan(statusOf(c("1.1", "1.2", "2.1", "2.2"), c(10, 100, 10, 10), c(10, 100, 10, 10)),
                       neighboursOf(c("1.1", "2.1"), c("1.2", "2.2"), c(4000, 4000)))
  expect_identical(nrow(plan), 0L)
})

test_that("merged names use the shared base and the members' last parts in numeric order", {
  expect_identical(ELFmergedName(c("3.2.4", "3.2.1")), "3.2.1_4")
  expect_identical(ELFmergedName(c("10.2", "10.1")), "10.1_2")
  expect_identical(ELFmergedName(c("3.2.10", "3.2.9")), "3.2.9_10")
  expect_error(ELFmergedName(c("3.1.1", "3.2.1")))
})

test_that("ELFrunName maps a merged member to the merged ELF and leaves others alone", {
  plan <- ELFmergePlan(statusOf(c("3.1.1", "3.1.2"), c(60, 10), c(60, 10)),
                       neighboursOf("3.1.1", "3.1.2", 4000))
  expect_identical(ELFrunName("3.1.2", plan), "3.1.1_2")
  expect_identical(ELFrunName("3.1.1", plan), "3.1.1_2")
  expect_identical(ELFrunName("4.3", plan), "4.3")
})

test_that("mergeELFs replaces the members by one ELF in every map", {
  whole <- bandELFs(list(`3.1.1` = 1:4, `3.1.2` = 6:9, `3.1.3` = 11:12),
                    buffers = list(`3.1.1` = 5, `3.1.2` = c(5, 10)))
  asPoly <- function(r, nam) {
    a <- r
    a[a[] == 0] <- NA
    vec <- terra::as.polygons(a)
    vec[, "ID"] <- nam
    vec[, "buffer"] <- vec[, nam]
    vec[, nam] <- NULL
    vec
  }
  ## named by ELF, as makeELFs() returns them
  layers <- stats::setNames(as.list(whole), names(whole))
  ELFs <- list(rasCentered = layers, rasWhole = layers,
               poly = Reduce(rbind, Map(asPoly, layers, names(layers))))
  plan <- ELFmergePlan(statusOf(c("3.1.1", "3.1.2", "3.1.3"), c(60, 10, 100), c(60, 10, 100)),
                       neighboursOf("3.1.1", "3.1.2", 4000))

  merged <- mergeELFs(ELFs, plan)
  expect_setequal(names(merged$rasWhole), c("3.1.3", "3.1.1_2"))
  expect_setequal(names(merged$rasCentered), c("3.1.3", "3.1.1_2"))
  expect_setequal(unique(merged$poly$ID), c("3.1.3", "3.1.1_2"))

  m <- matrix(terra::values(merged$rasWhole[["3.1.1_2"]], mat = FALSE), nrow = 4, byrow = TRUE)
  expect_true(all(m[, c(1:4, 6:9)] == 2))  # a cell that is core in either member is core
  expect_true(all(m[, c(5, 10)] == 1))     # buffer in the members, core in neither: still buffer
  expect_true(all(m[, 11:12] == 0))
  expect_identical(names(merged$rasCentered[["3.1.1_2"]]), "3.1.1_2")
  expect_true(2 %in% terra::values(merged$rasCentered[["3.1.1_2"]], mat = FALSE))
  ## an ELF not in the plan is untouched
  expect_identical(terra::values(merged$rasWhole[["3.1.3"]]), terra::values(ELFs$rasWhole[["3.1.3"]]))
})
