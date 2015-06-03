scgeMean <- function(data) {
  geneLogMeans <- apply(data, 2, function(x) {log(mean(x[x != 0]))})
  object <- list(data = data, geneLogMeans = geneLogMeans,
                 mean = mean(geneLogMeans), sd = sd(geneLogMeans))
  class(object) <- "scgeMean"
  return(object)
}

coef.scgeMean <- function(object) {
  return(c(mean = object$mean, sd = object$sd))
}

plot.scgeMean <- function(object) {
  hist(object$geneLogMeans, freq = FALSE, main = "Log Mean Gene Expression Fit")
  support <- seq(min(object$geneLogMeans), max(object$geneLogMeans), length.out = 100)
  lines(support, dnorm(support, object$mean, object$sd), col = "blue")
  invisible()
}

pedict.scgeMean <- function(object) {
  return(exp(object$mean))
}

print.scgeMean <- function(object) {
  message("An scgeMean object from package scge.")
}

simulate.scgeMean <- function(object, n) {
  exp(rnorm(n, object$mean, object$sd))
}

summary.scgeMean <- function(object) {
  message("Distribution: Log-normal")
  message("Location: ", object$mean)
  message("Scale: ", object$sd)
}
