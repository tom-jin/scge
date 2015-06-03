scgeCopula <- function(data) {
  cor_mat <- cor(data)
  fix_mat <- nearPD(cor_mat, corr = TRUE)$mat
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

simulate.scgeCopula <- function(object, n) {
  return(object$chol %*% matrix(rnorm(object$ncol * n), object$ncol, n))
}

summary.scgeCopula <- function(object) {
  message("Distribution: Gaussian copula")
}
