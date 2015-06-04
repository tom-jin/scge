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

print.scge <- function(object) {
  message("An scge object from package scge.")
}

simulate.scge <- function(object) {
  param_p <- function(mean, var) {
    return(mean/var)
  }
  param_n <- function(mean, p) {
    return(mean * p / (1 - p))
  }

  geneMean <- simulate(object$meanObject, object$nGenes)
  geneVar  <- simulate(object$varObject, geneMean)
  geneCensor <- simulate(object$censorObject, geneMean)
  geneQuant <- simulate(object$copulaObject, object$nGenes)
  geneP <- param_p(geneMean, geneVar)
  geneN <- param_n(geneMean, geneP)

  output <- rep(NA, 15122)
  for (i in 1:15122) {
    if (geneQuant[i] < geneCensor[i]) {
      output[i] <- 0
    } else {
      adj_quantile <- pnbinom(0, geneN[i], geneP[i]) + ((1 - pnbinom(0, geneN[i], geneP[i])) * (geneQuant[i] - geneCensor[i]) / (1 - geneCensor[i]))
      output[i] <- qnbinom(adj_quantile, geneN[i], geneP[i])
    }
  }
  return(output)
}

summary.scge <- function(object) {
  message("An scge object from package scge.")
}
