#' Fitting Gene Mean-Variance Models
#'
#' \code{scge} is used to fit gene mean-variance models.
#'
#' @param data A matrix of gene expression counts. Rows should represent cells
#' and columns represent genes.
#' @return \code{scgeVar} returns an object of class "scgeVar".
#' @seealso \code{\link{scge}}, \code{\link{scgeMean}}, \code{\link{scgeCensor}}
#' and \code{\link{scgeCopula}}
#' @export
scgeVar <- function(data = NULL) {
  if (is.null(data)) {
    object <- list(data = NA, geneMean = NA, geneVar = NA, a = 2.8, b = 1.7,
                   noiseSD = 0.8)
  } else {
    geneMean <- apply(data, 2, function(x) {mean(x[x != 0])})
    geneVar <- apply(data, 2, function(x) {var(x[x != 0])})
    linearModel <- lm(I(log(geneVar)) ~ I(log(geneMean)))
    noiseSD <- sd(log(geneVar) - coef(linearModel)[2] * log(geneMean) -
                    coef(linearModel)[1], na.rm = TRUE)
    object <- list(data = data, geneMean = geneMean, geneVar = geneVar,
                   a = coef(linearModel)[1], b = coef(linearModel)[2],
                   noiseSD = noiseSD)
  }
  class(object) <- "scgeVar"
  return(object)
}

#' @export
coef.scgeVar <- function(object, ...) {
  return(c(a = object$a, b = object$b))
}

#' @export
plot.scgeVar <- function(x, ...) {
  object <- x
  if (length(object$data) == 1) {
    support <- exp(seq(log(5), log(50000),
                       length.out = 100))
    plot(support, exp(object$a)*support ^ object$b, type = "l", col = "blue",
         log = "xy", xlim = c(5, 50000), ylim = c(1, 1e10),
         xlab = "Gene Expression Mean", ylab = "Gene Expression Variance",
         main = "Gene Mean-Var of Non-Zero Data")
    abline(0, 1, col = "red", untf = TRUE)
    lines(support, exp(object$a + 2*object$noiseSD)*support ^ object$b, col = "green")
    lines(support, exp(object$a - 2*object$noiseSD)*support ^ object$b, col = "green")
  } else {
    plot(object$geneMean, object$geneVar, log = "xy", xlim = c(5, 50000),
         ylim = c(1, 1e10), xlab = "Gene Expression Mean",
         ylab = "Gene Expression Variance", main = "Gene Mean-Var of Non-Zero Data", ...)
    abline(0, 1, col = "red", untf = TRUE)
    support <- exp(seq(log(min(object$geneMean)), log(max(object$geneMean)),
                       length.out = 100))
    lines(support, exp(object$a)*support ^ object$b, col = "blue")
    lines(support, exp(object$a + 2*object$noiseSD)*support ^ object$b, col = "green")
    lines(support, exp(object$a - 2*object$noiseSD)*support ^ object$b, col = "green")
  }
  invisible()
}

#' @export
predict.scgeVar <- function(object, mean, ...) {
  return(exp(object$a) * mean ^ object$b)
}

#' @export
print.scgeVar <- function(x, ...) {
  message("An scgeVar object from package scge.")
}

#' @export
simulate.scgeVar <- function(object, nsim = length(mean), seed = NULL, mean, ...) {
  return(exp(object$a + rnorm(nsim, 0, object$noiseSD)) * mean ^ object$b)
}

#' @export
summary.scgeVar <- function(object, ...) {
  message("Gene Mean-Variance Fit")
  message("Distribution: Log-log linear")
  message("Intercept: ", object$a)
  message("Slope: ", object$b)
  if (length(object$data) == 1) {
    message("Using default parameters.")
  }
}
