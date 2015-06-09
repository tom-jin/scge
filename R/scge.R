#' scge: A package for simulating single cell gene expression data.
#'
#' scge is a package to simulate realistic gene expression data, initially for
#' single cell genomics but generalisable to other sources of gene expression
#' data. Models are fitted to reproduce the mean expression, mean-variance,
#' mean-censorship and correlation relationships within the data. Parameters
#' from an example fit are provided should data be unavailable.
#'
#' @docType package
#' @name scge-package
NULL

scge <- function(data, meanObject = NULL, varObject = NULL, censorObject = NULL,
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

  object <- list(data = data, mean = meanObject, var = varObject,
                 censor = censorObject, copula = copulaObject,
                 nGenes = ncol(data))
  class(object) <- "scge"
  return(object)
}

plot.scge <- function(object) {
  plot(object$mean)
  null <- readline("Hit <Return> to see next plot: ")
  plot(object$var)
  null <- readline("Hit <Return> to see next plot: ")
  plot(object$censor)
}

print.scge <- function(object) {
  message("An scge object from package scge.")
}

simulate.scge <- function(object) {
  geneMean <- simulate(object$mean, object$nGenes)
  geneVar  <- simulate(object$var, geneMean)
  geneCensor <- simulate(object$censor, geneMean)
  geneQuant <- simulate(object$copula, object$nGenes)
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

summary.scge <- function(object) {
  message("An scge object from package scge.")
}

param_p <- function(mean, var) {return(mean/var)}
param_n <- function(mean, p) {return(mean * p / (1 - p))}
