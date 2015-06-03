scgeVar <- function(data) {
  geneMean <- apply(data, 2, function(x) {mean(x[x != 0])})
  geneVar <- apply(data, 2, function(x) {var(x[x != 0])})
  linearModel <- lm(I(log(geneVar)) ~ I(log(geneMean)))
  object <- list(data = data, geneMean = geneMean, geneVar = geneVar,
                 a = coef(linearModel)[1], b = coef(linearModel)[2])
  class(object) <- "scgeVar"
  return(object)
}

coef.scgeVar <- function(object) {
  return(c(a = object$a, b = object$b))
}

plot.scgeVar <- function(object) {
  plot(object$geneMean, object$geneVar, log = "xy", xlab = "Gene Expression Mean",
       ylab = "Gene Expression Variance", main = "Gene Mean-Var of Non-Zero Data")
  abline(0, 1, col = "red", untf = TRUE)
  support <- exp(seq(log(min(object$geneMean)), log(max(object$geneMean)),
                     length.out = 100))
  lines(support, object$a*support ^ object$b, col = "blue")
}

pedict.scgeVar <- function(object, mean) {
  return(exp(object$a) * mean ^ object$b)
}

print.scgeVar <- function(object) {
  message("An scgeVar object from package scge.")
}

simulate.scgeVar <- function(object) {
  #TODO: Fit noise.
}

summary.scgeVar <- function(object) {
  message("Distribution: Log-log linear")
  message("Intercept: ", object$a)
  message("Slope: ", object$b)
}
