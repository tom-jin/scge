#' Fitting Gene Correlation Models
#'
#' \code{scge} is used to fit gene correlation models.
#'
#' @param data A matrix of gene expression counts. Rows should represent cells
#' and columns represent genes.
#' @return \code{scgeCopula} returns an object of class "scgeCopula".
#' @seealso \code{\link{scge}}, \code{\link{scgeMean}}, \code{\link{scgeVar}}
#' and \code{\link{scgeCensor}}
#' @export
scgeCopula <- function(data = NULL) {
  if (is.null(data)) {
    object <- list(data = NA, ncol = 1000, chol = cho)
  } else {
    cor_mat <- cor(data)
    fix_mat <- Matrix::nearPD(cor_mat, corr = TRUE)$mat
    cho_mat <- chol(fix_mat)
    object <- list(data = data, ncol = ncol(data), chol = cho_mat)
  }
  class(object) <- "scgeCopula"
  return(object)
}

#' @export
coef.scgeCopula <- function(object, ...) {
  return(object$chol)
}

#' @export
print.scgeCopula <- function(x, ...) {
  message("An scgeCopulaobject from package scge.")
}

#' @export
simulate.scgeCopula <- function(object, nsim = 1, seed = NULL, ...) {
  return(pnorm(matrix(rnorm(object$ncol * nsim), nsim, object$ncol) %*% object$chol))
}

#' @export
summary.scgeCopula <- function(object, ...) {
  message("Gene Correlation Fit")
  message("Distribution: Gaussian copula")
  if (length(object$data) == 1) {
    message("Using default parameters.")
  }
}
