library(scge)
context("scgeCopula")
load(file = "testdata.RData")

defaultObject <- scgeCopula()
realdataObject <- scgeCopula(data = testdata)

test_that("scgeCopula produces default fit", {
  expect_is(defaultObject, "scgeCopula")
})

test_that("scgeCopula fits data", {
  expect_is(realdataObject, "scgeCopula")
  expect_equal(dim(realdataObject$chol), c(1000, 1000))
})

test_that("scgeCopula results match copula::normalCopula", {
  if (!require(copula)) {
    skip("Package copula not available")
  }

#   normal.cop <- normalCopula(dim = 100, dispstr = "un")
#   u <- pobs(testdata[, 1:100])
#   fit.ml <- fitCopula(normal.cop, u, method = "mpl")
#   getSigma(fit.ml$copula)
})