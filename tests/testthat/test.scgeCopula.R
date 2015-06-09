library(scge)
context("scgeCopula")

defaultObject <- scgeCopula()

test_that("scgeCopula produces default fit", {
  expect_is(defaultObject, "scgeCopula")
})
