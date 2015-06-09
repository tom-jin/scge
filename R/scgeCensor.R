#' Fitting Gene Censorship Models
#'
#' \code{scge} is used to fit gene censorship models.
#'
#' @param data A matrix of gene expression counts. Rows should represent cells
#' and columns represent genes.
#' @return \code{scgeCensor} returns an object of class "scgeCensor".
#' @seealso \code{\link{scge}}, \code{\link{scgeMean}}, \code{\link{scgeVar}}
#' and \code{\link{scgeCopula}}
#' @export
scgeCensor <- function(data = NULL) {
  if (is.null(data)) {
    object <- list(data = NA, geneMean = NA, geneCensor = NA, scale = -1,
                   position = 6, sd = 0.2)
  } else {
    geneMean <- apply(data, 2, function(x) {mean(x[x != 0])})
    geneCensor <- apply(data, 2, function(x) {sum(x == 0)})

    d <- data.frame(mean = geneMean, censor = geneCensor/249,
                    logit.censor = logit(geneCensor/249))[geneCensor != 0,]
    q <- quantile(geneMean, probs = c(0.5, 1))
    f <- lm(logit.censor ~ I(log(mean)),d[d$mean > q[1] & d$mean < q[2], ])

    object <- list(data = data, geneMean = geneMean, geneCensor = geneCensor,
                   scale = coef(f)[2], position = coef(f)[1],
                   sd = sd(d$censor - sigmoid(coef(f)[1] + coef(f)[2]*log(d$mean))))
  }
  class(object) <- "scgeCensor"
  return(object)
}

#' @export
coef.scgeCensor <- function(object) {
  return(c(position = object$position, scale = object$scale, sd = object$sd))
}

#' @export
plot.scgeCensor <- function(object) {
  if (length(object$data) == 1) {
    plot(1, type = "n", log = "x",
         xlab = "Gene Expression Mean", ylab = "Gene Censorship",
         main = "Sigmoid Fit", xlim = c(5, 50000), ylim = c(0, 1))
    support <- exp(seq(log(5), log(50000), length.out = 100))
  } else {
    plot(object$geneMean, object$geneCensor/nrow(object$data), log = "x",
         xlab = "Gene Expression Mean", ylab = "Gene Censorship",
         main = "Sigmoid Fit", xlim = c(5, 50000), ylim = c(0, 1))
    support <- exp(seq(log(min(object$geneMean)), log(max(object$geneMean)),
                       length.out = 100))
  }
  lines(support, sigmoid(object$position + object$scale * log(support)),
        col = "blue")
  lines(support, pmin(1, sigmoid(object$position + object$scale *
                                   log(support)) + 2*object$sd), col = "green")
  lines(support, pmax(0, sigmoid(object$position + object$scale *
                                   log(support)) - 2*object$sd), col = "green")
}

#' @export
predict.scgeCensor <- function(object, mean) {
  return(sigmoid(object$position + object$scale * log(mean)))
}

#' @export
print.scgeCensor <- function(object) {
  message("An scgeCensor object from package scge.")
}

#' @export
simulate.scgeCensor <- function(object, mean) {
  samples <- sapply(mean, function(m) {
    sample <- rnorm(1, mean = sigmoid(object$position + object$scale * log(m)),
                    sd = object$sd)
    sample <- pmax(0, sample)
    sample <- pmin(1, sample)
    sample
  })
  return(samples)
}

#' @export
summary.scgeCensor <- function(object) {
  message("Distribution: Sigmoid with normal noise restricted to the [0,1] interval")
  message("Sigmoid Position: ", object$position)
  message("Sigmoid Scale: ", object$scale)
  message("Noise Scale: ", object$sd)
  if (length(object$data) == 1) {
    message("Using default parameters.")
  } else {
    message("Parameters fitted to ", nrow(object$data), " observations.")
  }
}

logit <- function(x) {-log((1/x) - 1)}
sigmoid <- function(x) {1 / (1 + exp(-x))}
