library(scge)
context("scgeMean")

defaultObject <- scgeMean()

test_that("scgeMean produces default fit", {
  expect_is(defaultObject, "scgeMean")
  expect_true(is.na(defaultObject$data))
})

test_that("coef.scgeMean", {
  expect_equivalent(coef(defaultObject), c(5.5, 1.0))
})

