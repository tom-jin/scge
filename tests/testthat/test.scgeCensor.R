library(scge)
context("scgeCensor")

defaultObject <- scgeCensor()

test_that("scgeCensor produces default fit", {
  expect_is(defaultObject, "scgeCensor")
  expect_true(is.na(defaultObject$data))
  expect_true(is.na(defaultObject$geneMean))
  expect_true(is.na(defaultObject$geneCensor))
})

test_that("coef.scgeCensor", {
  expect_equivalent(coef(defaultObject), c(6, -1, 0.2))
})

test_that("predict.scgeCensor", {
  expect_true(predict(defaultObject, 42) - 0.9057088 < 1e-7)
})
