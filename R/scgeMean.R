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
scgeMean <- function(data = NULL) {
  if (is.null(data)) {
    object <- list(data = NA, geneLogMean = NA, mean = 5.5, sd = 1.0)
  } else {
    geneLogMean <- apply(data, 2, function(x) {log(mean(x[x != 0]))})
    object <- list(data = data, geneLogMean = geneLogMean,
                   mean = mean(geneLogMean, na.rm = TRUE),
                   sd = sd(geneLogMean, na.rm = TRUE))
  }

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
  if (length(object$data) == 1) {
    support <- seq(0, 10, length.out = 100)
    plot(support, dnorm(support, object$mean, object$sd), type = "l",
         main = "Log Mean Gene Expression Fit", xlim = c(0, 10), ylim = c(0, 0.4),
         xlab = "Log Mean Gene Expression", ylab = "Density", col = "blue")
  } else {
    hist(object$geneLogMean, freq = FALSE, main = "Log Mean Gene Expression Fit",
         xlab = "Log Mean Gene Expression", xlim = c(0, 10), ylim = c(0, 0.4), ...)
    support <- seq(min(object$geneLogMean, na.rm = TRUE),
                   max(object$geneLogMean, na.rm = TRUE), length.out = 100)
    lines(support, dnorm(support, object$mean, object$sd), col = "blue")
  }

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
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if (is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }

  val <- exp(rnorm(nsim, object$mean, object$sd))
  attr(val, "seed") <- RNGstate
  return(val)
}

#' @export
summary.scgeMean <- function(object, ...) {
  message("Gene Mean Expression Fit")
  message("Distribution: Log-normal")
  message("Location: ", object$mean)
  message("Scale: ", object$sd)
  if (length(object$data) == 1) {
    message("Using default parameters.")
  }
}
