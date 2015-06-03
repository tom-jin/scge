scgeVar <- function(data) {
  geneMeans <- apply(data, 2, function(x) {mean(x[x != 0])})
  geneVars <- apply(data, 2, function(x) {var(x[x != 0])})
  linearModel <- lm(I(log(geneVars)) ~ I(log(geneMeans)))
  object <- list(data = data, geneMeans = geneMeans, geneVars = geneVars,
                 a = coef(linearModel)[1], b = coef(linearModel)[2])
  class(object) <- "scgeVar"
  return(object)
}

coef.scgeVar <- function(object) {
  return(c(a = object$a, b = object$b))
}

plot.scgeVar <- function(object) {
  plot(object$geneMeans, object$geneVars, log = "xy", xlab = "Gene Expression Mean",
       ylab = "Gene Expression Variance", main = "Gene Mean-Var of Non-Zero Data")
  abline(0, 1, col = "red", untf = TRUE)
  support <- exp(seq(log(min(object$geneMeans)), log(max(object$geneMeans)),
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
