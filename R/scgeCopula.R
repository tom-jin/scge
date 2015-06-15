#' Fitting Gene Correlation Models
#'
#' \code{scge} is used to fit gene correlation models.
#'
#' @param data A matrix of gene expression counts. Rows should represent cells
#' and columns represent genes.
#' @return \code{scgeCopula} returns an object of class "scgeCopula".
#' @seealso \code{\link{scge}}, \code{\link{scgeMean}}, \code{\link{scgeVar}}
#' and \code{\link{scgeCensor}}
#' @references Pollen, A., Nowakowski, T., Shuga, J., Wang, X., Leyrat, A.,
#' Lui, J., Li, N., Szpankowski, L., Fowler, B., Chen, P., Ramalingam, N.,
#' Sun, G., Thu, M., Norris, M., Lebofsky, R., Toppani, D., Kemp, D., Wong, M.,
#' Clerkson, B., Jones, B., Wu, S., Knutsson, L., Alvarado, B., Wang, J.,
#' Weaver, L., May, A., Jones, R., Unger, M., Kriegstein, A. and West,
#' J. (2014). \emph{Low-coverage single-cell mRNA sequencing reveals cellular
#' heterogeneity and activated signaling pathways in developing cerebral cortex.}
#' Nat Biotechnol, [online] 32(10), pp.1053-1058. Available at:
#' \url{http://www.nature.com/nbt/journal/v32/n10/abs/nbt.2967.html}
#' @export
scgeCopula <- function(data = NULL) {
  if (is.null(data)) {
    object <- list(data = NA, ncol = 1000, chol = cho)
  } else {
    pdo_obs <- apply(data, 2, rank, na.last = NA)/(nrow(data) + 1)
    cor_mat <- cor(pdo_obs)
    fix_mat <- Matrix::nearPD(cor_mat, corr = TRUE)$mat
    cho_mat <- t(chol(fix_mat))
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
  return(pnorm(object$chol %*% matrix(rnorm(object$ncol * nsim), object$ncol, nsim)))
}

#' @export
summary.scgeCopula <- function(object, ...) {
  message("Gene Correlation Fit")
  message("Distribution: Gaussian copula")
  if (length(object$data) == 1) {
    message("Using default parameters.")
  }
}
