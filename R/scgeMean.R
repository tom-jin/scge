scgeMean <- function(data) {
  geneLogMean <- apply(data, 2, function(x) {log(mean(x[x != 0]))})
  object <- list(data = data, geneLogMean = geneLogMean,
                 mean = mean(geneLogMean), sd = sd(geneLogMean))
  class(object) <- "scgeMean"
  return(object)
}

coef.scgeMean <- function(object) {
  return(c(mean = object$mean, sd = object$sd))
}

plot.scgeMean <- function(object) {
  hist(object$geneLogMean, freq = FALSE, main = "Log Mean Gene Expression Fit",
       xlab = "Log Mean Gene Expression")
  support <- seq(min(object$geneLogMean), max(object$geneLogMean), length.out = 100)
  lines(support, dnorm(support, object$mean, object$sd), col = "blue")
  invisible()
}

predict.scgeMean <- function(object) {
  return(exp(object$mean))
}

print.scgeMean <- function(object) {
  message("An scgeMean object from package scge.")
}

simulate.scgeMean <- function(object, n = 1) {
  exp(rnorm(n, object$mean, object$sd))
}

summary.scgeMean <- function(object) {
  message("Distribution: Log-normal")
  message("Location: ", object$mean)
  message("Scale: ", object$sd)
}
