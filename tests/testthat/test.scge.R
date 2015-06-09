library(scge)
context("scge")

defaultObject <- scge()

test_that("scge produces default fit", {
  expect_is(defaultObject, "scge")
})
