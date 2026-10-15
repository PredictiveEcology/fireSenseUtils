test_that("resolveDoFit: doFit NA, TRUE and FALSE pass through", {
  expect_identical(resolveDoFit(NA), NA)
  expect_identical(resolveDoFit(), NA)
  expect_identical(resolveDoFit(TRUE), TRUE)
  expect_identical(resolveDoFit(FALSE), FALSE)
})

test_that("resolveDoFit: the deprecated aliases map, and refitExisting wins", {
  q <- function(...) suppressMessages(resolveDoFit(...))
  expect_identical(q(refitExisting = TRUE), TRUE)
  expect_identical(q(stopIfNoPreRunFit = TRUE), FALSE)
  expect_identical(q(stopIfNoPreRunFit = FALSE), NA)       # only allowed a fit: what NA does
  expect_identical(q(refitExisting = FALSE), NA)
  expect_identical(q(refitExisting = TRUE, stopIfNoPreRunFit = TRUE), TRUE)
  expect_identical(q(refitExisting = FALSE, stopIfNoPreRunFit = TRUE), FALSE)
  expect_identical(q(doFit = FALSE, refitExisting = TRUE), FALSE)   # explicit doFit is kept
  expect_identical(q(doFit = TRUE, stopIfNoPreRunFit = TRUE), TRUE)
})

test_that("resolveDoFit: one message naming doFit and the module when an alias is set, none otherwise", {
  expect_message(resolveDoFit(refitExisting = TRUE, module = "mod"), "mod.*deprecated.*doFit")
  expect_message(resolveDoFit(stopIfNoPreRunFit = FALSE), "deprecated.*doFit")
  expect_no_message(resolveDoFit(doFit = NA))
  expect_no_message(resolveDoFit(refitExisting = TRUE, message = FALSE))
  msgs <- character()
  withCallingHandlers(resolveDoFit(refitExisting = TRUE, stopIfNoPreRunFit = TRUE),
                      message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  expect_length(msgs, 1L)
})

test_that("mustFit: doFit x ledger", {
  expect_true(mustFit(NA, FALSE)); expect_false(mustFit(NA, TRUE))
  expect_true(mustFit(TRUE, FALSE)); expect_true(mustFit(TRUE, TRUE))
  expect_false(mustFit(FALSE, TRUE))
  expect_error(mustFit(FALSE, FALSE, polygon = "9.9"), "polygon 9\\.9.*doFit")
})

test_that("mustFit does not consult the ledger when doFit is TRUE", {
  expect_true(mustFit(TRUE, stop("consulted")))
})
