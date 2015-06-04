scgeCopula <- function(data) {
  cor_mat <- cor(data)
  fix_mat <- Matrix::nearPD(cor_mat, corr = TRUE)$mat
  cho_mat <- chol(fix_mat)
  object <- list(data = data, ncol = ncol(data), chol = cho_mat)
  class(object) <- "scgeCopula"
  return(object)
}

coef.scgeCopula <- function(object) {
  return(object$chol)
}

print.scgeCopula <- function(object) {
  message("An scgeCopulaobject from package scge.")
}

simulate.scgeCopula <- function(object, n = 1) {
  return(pnorm(matrix(rnorm(object$ncol * n), n, object$ncol) %*% object$chol))
}

summary.scgeCopula <- function(object) {
  message("Distribution: Gaussian copula")
}
