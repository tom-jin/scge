scgeCensor <- function(data) {
  geneMean <- apply(data, 2, function(x) {mean(x[x != 0])})
  geneCensor <- apply(data, 2, function(x) {sum(x == 0)})
  #TODO: Fit sigmoid function.
  object <- list(data = data, geneMean = geneMean, geneCensor = geneCensor, sigmoidScale = -1.4,
                 position = 5, noiseScale = 0.5)
  class(object) <- "scgeCensor"
  return(object)
}

coef.scgeCensor <- function(object) {
  return(c(position = object$position, sigmoidScale = object$sigmoidScale,
           noiseScale = object$noiseScale))
}

plot.scgeCensor <- function(object) {
  plot(object$geneMean, object$geneCensor/nrow(object$data), log = "x",
       xlab = "Gene Expression Mean", ylab = "Gene Censorship",
       main = "Sigmoid Fit")
  support <- exp(seq(log(min(object$geneMean)), log(max(object$geneMean)),
                     length.out = 100))
  lines(support, sigmoid(support, object$sigmoidScale, object$position), col = "blue")
  lines(support, pmin(1, sigmoid(support, object$sigmoidScale, object$position) + secant(support, object$noiseScale, object$position)), col = "green")
  lines(support, pmax(0, sigmoid(support, object$sigmoidScale, object$position) - secant(support, object$noiseScale, object$position)), col = "green")
}

pedict.scgeCensor <- function(object, mean) {
  return(sigmoid(mean, object$sigmoidScale, object$position))
}

print.scgeCensor <- function(object) {
  message("An scgeCensor object from package scge.")
}

simulate.scgeCensor <- function(object, mean) {
  samples <- sapply(mean, function(m) {
    sample <- rnorm(1, mean = sigmoid(m, object$sigmoidScale, object$position),
                    sd = secant(m, object$noiseScale, object$position))
    sample <- pmax(0, sample)
    sample <- pmin(1, sample)
    sample
  })
  return(samples)
}

summary.scgeCensor <- function(object) {
  message("Distribution: Sigmoid with normal noise restricted to the [0,1] interval")
  message("Sigmoid Scale: ", object$sigmoidScale)
  message("Noise Scale: ", object$noiseScale)
  message("Position: ", object$position)
}

sigmoid <- function(x, scale, position) {1 / (1 + exp(scale * (position - log(x))))}
secant <- function(x, scale, position) {scale/(exp(log(x) - position) + exp(position - log(x)))}
