#' scge: A package for simulating single cell gene expression data.
#'
#' scge is a package to simulate realistic gene expression data, initially for
#' single cell genomics but generalisable to other sources of gene expression
#' data. Models are fitted to reproduce the mean expression, mean-variance,
#' mean-censorship and correlation relationships within the data. Parameters
#' from an example fit are provided should data be unavailable.
#'
#' @docType package
#' @name scge
NULL

#' Fitting Gene Expression Models
#'
#' \code{scge} is used to fit gene expression models. This can be performed
#' automatically by passing a matrix a counts or by passing components of the
#' model.
#'
#' @param data A matrix of gene expression counts. Rows should represent cells
#' and columns represent genes.
#' @param meanObject The option object from the \code{scgeMean} function.
#' @param varObject The option object from the \code{scgeVar} function.
#' @param censorObject The option object from the \code{scgeCensor} function.
#' @param copulaObject The option object from the \code{scgeCopula} function.
#' @return \code{scge} returns an object of class "scge".
#' @seealso \code{\link{scgeMean}}, \code{\link{scgeVar}},
#' \code{\link{scgeCensor}} and \code{\link{scgeCopula}}
#' @export
scge <- function(data = NULL, meanObject = NULL, varObject = NULL, censorObject = NULL,
                 copulaObject = NULL) {
  if (!is.null(data)) {
    data <- sanitise(data)
  }
  if (is.null(meanObject)) {
    meanObject <- scgeMean(data)
  }
  if (is.null(varObject)) {
    varObject <- scgeVar(data)
  }
  if (is.null(censorObject)) {
    censorObject <- scgeCensor(data)
  }
  if (is.null(copulaObject)) {
    copulaObject <- scgeCopula(data)
  }
  if (is.null(data)) {
    object <- list(data = data, mean = meanObject, var = varObject,
                   censor = censorObject, copula = copulaObject,
                   nGenes = 1000)
  } else {
    object <- list(data = data, mean = meanObject, var = varObject,
                   censor = censorObject, copula = copulaObject,
                   nGenes = ncol(data))
  }
  class(object) <- "scge"
  return(object)
}

#' @export
plot.scge <- function(x, ...) {
  plot(x$mean)
  null <- readline("Hit <Return> to see next plot: ")
  plot(x$var)
  null <- readline("Hit <Return> to see next plot: ")
  plot(x$censor)
}

#' @export
print.scge <- function(x, ...) {
  message("An scge object from package scge.")
}

#' @export
#' @importFrom stats simulate
simulate.scge <- function(object, nsim = 1, seed = NULL, ...) {
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

  if (is.na(object$mean$geneLogMean)) {
    geneMean <- simulate(object$mean, 1000)
  } else {
    geneMean <- object$var$geneMean
  }

  if (is.na(object$var$geneVar)) {
    geneVar <- simulate(object$var, mean = geneMean)
  } else {
    geneVar  <- object$var$geneVar
  }

  if (is.na(object$censor$geneCensor)) {
    geneCensor <- simulate.scgeCensor(object$censor, nsim , mean = geneMean)
  } else {
    geneCensor <- object$censor$geneCensor
  }

  geneQuant <- simulate.scgeCopula(object$copula, nsim)
  geneP <- param_p(geneMean, geneVar)
  geneN <- param_n(geneMean, geneP)

  adj_quantile <- (geneQuant - geneCensor) / (1 - geneCensor)

  if (nsim > 1) {
    adj_quantile <- apply(adj_quantile, 2, pmax, 0)
    adj_quantile <- apply(adj_quantile, 2, pmin, 1)
  } else {
    adj_quantile <- pmax(adj_quantile, 0)
    adj_quantile <- pmin(adj_quantile, 1)
  }

  output <- qnbinom(adj_quantile, geneN, geneP)
  output <- output * (geneQuant > geneCensor)
  attr(output, "seed") <- RNGstate
  return(output)
}

#' @export
summary.scge <- function(object, ...) {
  message("Gene Expression fit composed of:")
  summary(object$mean)
  summary(object$var)
  summary(object$censor)
  summary(object$copula)
}

param_p <- function(mean, var) {return(mean/var)}
param_n <- function(mean, p) {return(mean * p / (1 - p))}
