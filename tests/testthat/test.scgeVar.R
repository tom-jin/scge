library(scge)
context("scgeVar")

defaultObject <- scgeVar()

test_that("scgeVar produces default fit", {
  expect_is(defaultObject, "scgeVar")
})
