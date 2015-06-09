#' Fitting Gene Mean Expression Models
#'
#' \code{scge} is used to fit gene mean expression models.
#'
#' @param data A matrix of gene expression counts. Rows should represent cells
#' and columns represent genes.
#' @return \code{scgeMean} returns an object of class "scgeMean".
#' @seealso \code{\link{scge}}, \code{\link{scgeVar}}, \code{\link{scgeCensor}}
#' and \code{\link{scgeCopula}}
#' @export
scgeMean <- function(data) {
  geneLogMean <- apply(data, 2, function(x) {log(mean(x[x != 0]))})
  object <- list(data = data, geneLogMean = geneLogMean,
                 mean = mean(geneLogMean), sd = sd(geneLogMean))
  class(object) <- "scgeMean"
  return(object)
}

#' @export
coef.scgeMean <- function(object, ...) {
  return(c(mean = object$mean, sd = object$sd))
}

#' @export
plot.scgeMean <- function(x, ...) {
  object <- x
  hist(object$geneLogMean, freq = FALSE, main = "Log Mean Gene Expression Fit",
       xlab = "Log Mean Gene Expression", ...)
  support <- seq(min(object$geneLogMean), max(object$geneLogMean), length.out = 100)
  lines(support, dnorm(support, object$mean, object$sd), col = "blue")
  invisible()
}

#' @export
predict.scgeMean <- function(object, ...) {
  return(exp(object$mean))
}

#' @export
print.scgeMean <- function(x, ...) {
  message("An scgeMean object from package scge.")
}

#' @export
simulate.scgeMean <- function(object, nsim = 1, seed = NULL, ...) {
  exp(rnorm(nsim, object$mean, object$sd))
}

#' @export
summary.scgeMean <- function(object, ...) {
  message("Distribution: Log-normal")
  message("Location: ", object$mean)
  message("Scale: ", object$sd)
}
