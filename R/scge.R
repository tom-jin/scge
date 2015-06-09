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
  geneMean <- simulate.scgeMean(object$mean, object$nGenes)
  geneVar  <- simulate.scgeVar(object$var, geneMean)
  geneCensor <- simulate.scgeCensor(object$censor, geneMean)
  geneQuant <- simulate.scgeCopula(object$copula, object$nGenes)
  geneP <- param_p(geneMean, geneVar)
  geneN <- param_n(geneMean, geneP)

  output <- rep(NA, object$nGenes)
  for (i in 1:object$nGenes) {
    if (geneQuant[i] < geneCensor[i]) {
      output[i] <- 0
    } else {
      adj_quantile <- pnbinom(0, geneN[i], geneP[i]) +
        ((1 - pnbinom(0, geneN[i], geneP[i])) * (geneQuant[i] - geneCensor[i]) /
           (1 - geneCensor[i]))
      output[i] <- qnbinom(adj_quantile, geneN[i], geneP[i])
    }
  }
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
