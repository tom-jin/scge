#' Fitting Gene Censorship Models
#'
#' \code{scge} is used to fit gene censorship models.
#'
#' @param data A matrix of gene expression counts. Rows should represent cells
#' and columns represent genes.
#' @return \code{scgeCensor} returns an object of class "scgeCensor".
#' @seealso \code{\link{scge}}, \code{\link{scgeMean}}, \code{\link{scgeVar}}
#' and \code{\link{scgeCopula}}
#' @examples
#' obj <- scgeCensor()
#' plot(obj)
#' predict(obj, 42)
#' simulate(obj, nsim = 10, mean = 42)
#' @export
scgeCensor <- function(data = NULL) {
  if (is.null(data)) {
    object <- list(data = NA, geneMean = NA, geneCensor = NA, scale = -1,
                   position = 6, sd = 0.2)
  } else {
    geneMean <- apply(data, 2, function(x) {mean(x[x != 0])})
    geneCensor <- apply(data, 2, function(x) {sum(x == 0)})
    rows <- nrow(data)
    d <- data.frame(mean = geneMean, censor = geneCensor/rows,
                    logit.censor = logit(geneCensor/rows))[geneCensor != 0,]
    d <- d[is.finite(d$mean),]
    d <- d[is.finite(d$logit.censor),]
    q <- quantile(geneMean, probs = c(0.5, 1), na.rm = TRUE)
    f <- lm(logit.censor ~ I(log(mean)),d[d$mean > q[1] & d$mean < q[2], ])

    object <- list(data = data, geneMean = geneMean, geneCensor = geneCensor,
                   scale = coef(f)[2], position = coef(f)[1],
                   sd = sd(d$censor - sigmoid(coef(f)[1] + coef(f)[2]*log(d$mean))))
  }
  class(object) <- "scgeCensor"
  return(object)
}

#' @export
coef.scgeCensor <- function(object, ...) {
  return(c(position = object$position, scale = object$scale, sd = object$sd))
}

#' @export
plot.scgeCensor <- function(x, ...) {
  object <- x
  if (length(object$data) == 1) {
    plot(1, type = "n", log = "x",
         xlab = "Gene Expression Mean", ylab = "Gene Censorship",
         main = "Sigmoid Fit", xlim = c(5, 50000), ylim = c(0, 1), ...)
    support <- exp(seq(log(5), log(50000), length.out = 100))
  } else {
    plot(object$geneMean, object$geneCensor/nrow(object$data), log = "x",
         xlab = "Gene Expression Mean", ylab = "Gene Censorship",
         main = "Sigmoid Fit", xlim = c(5, 50000), ylim = c(0, 1), ...)
    support <- exp(seq(log(5), log(max(object$geneMean)),
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
predict.scgeCensor <- function(object, mean, ...) {
  return(sigmoid(object$position + object$scale * log(mean)))
}

#' @export
print.scgeCensor <- function(x, ...) {
  message("An scgeCensor object from package scge.")
}

#' @export
simulate.scgeCensor <- function(object, nsim = 1, seed = NULL, mean, ...) {
  if (missing("mean"))
    stop("What mean gene expression do you want to simulate from?")
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

  sample <- sapply(mean, function(mean, n) {rnorm(n, mean = sigmoid(object$position + object$scale * log(mean)), sd = object$sd)}, n = nsim)

  if (nsim > 1) {
    sample <- apply(sample, 2, pmax, 0)
    sample <- apply(sample, 2, pmin, 1)
  } else {
    sample <- pmax(sample, 0)
    sample <- pmin(sample, 1)
  }

  attr(sample, "seed") <- RNGstate
  return(sample)
}

#' @export
summary.scgeCensor <- function(object, ...) {
  message("Gene Censorship Fit")
  message("Distribution: Sigmoid with normal noise restricted to the [0,1] interval")
  message("Sigmoid Position: ", object$position)
  message("Sigmoid Scale: ", object$scale)
  message("Noise Scale: ", object$sd)
  if (length(object$data) == 1) {
    message("Using default parameters.")
  }
}

logit <- function(x) {-log((1/x) - 1)}
sigmoid <- function(x) {1 / (1 + exp(-x))}
