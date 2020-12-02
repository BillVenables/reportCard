### rank discrete scores
###

#' @import graphics
NULL
#' @import stats
NULL

#' Aggregate scores
#'
#' Given a set of scores on [0,1] calculate an
#' aggregate score using Bayesian bootstrap methods
#'
#' @param scores numeric with values in [0,1]
#' @param bounds numeric grade boundaries
#' @param labels character grade labels
#' @param B numeric number of bootstrap simulations to use
#' @param x an object of class "agg_scores"
#' @param ... additional arguments
#'
#' @return A character string giving the aggregate grade, with
#'         several attributes containing details from the aggregation.
#' @export
#'
#' @examples
#' sc <- rbeta(6, 1.5, 3.5)
#' ag <- agg_scores(sc)
#' ag
agg_scores <- function(scores, bounds = (0:5)/5, labels = LETTERS[5:1], B = 100000) {
  stopifnot(is.numeric(scores), all(scores >= 0 & scores <= 1))
  X <- matrix(rexp(B*length(scores)), B)
  bootscores <- drop((X/rowSums(X)) %*% scores)
  fn <- function(a) mean(scores) - mean(bootscores^a)
  x1 <- x2 <- 1
  while(fn(x1) > 0) x1 <- x1/2
  while(fn(x2) < 0) x2 <- x2 + 1/2
  a <- uniroot(fn, c(x1, x2))$root
  bootscores <- bootscores^a
  boots <- factor(labels[findInterval(bootscores, bounds)], levels = labels)
  tab <- table(boots)/B*100
  combined <- names(tab)[which.max(tab)]
  structure(combined,
            bootscores = bootscores,
            grades = factor(labels[findInterval(scores, bounds, all.inside = TRUE)], levels = labels),
            errors = tab,
            class = "agg_scores")
}

#' @rdname agg_scores
#' @export
print.agg_scores <- function(x, ...) {
  cat("Original grades: ", paste(attr(x, "grades"), collapse = ", "), "\n")
  cat("Grade inclusion:\n")
  print(noquote(setNames(paste0(attr(x, "errors"), "%"), names(attr(x, "errors")))))
  cat("Aggregate: ", as.vector(x), "\n")
  invisible(x)
}

#' @rdname agg_scores
#' @export
plot.agg_scores <- function(x, ...) {
  freq <- attr(x, "errors")
  col <- rep("light blue", length(freq))
  col[which.max(freq)] <- "steel blue"
  b <- barplot(freq, col = col, border = "steel blue")
  grid()
  par(new = TRUE)
  barplot(freq, col = col, border = "steel blue")
  text(b, freq, labels = paste0(round(freq, 1), "%"), pos = 3, xpd = NA)
  abline(h = 0)
  invisible(x)
}

#' Extract grades
#'
#' Extract the individual grades from an aggregated scores object.
#'
#' @param ag  output of agg_scores()
#' @param ... other arguments, currently ignored
#'
#' @return A factor giving the grades
#' @export
#'
#' @examples
#' sc <- rbeta(6, 1.5, 3.5)
#' ag <- agg_scores(sc)
#' grades(ag)
grades <- function(ag, ...) {
  UseMethod("grades")
}

#' @rdname grades
#' @export
grades.agg_scores <- function(ag, ...) {
  attr(ag, "grades")
}


#' Distribution function and inverse
#'
#' Cumulative distribution function and its inverse
#' from the bootstrap simulations conducted as part
#' of the scores aggregation process.
#'
#' @param agg,x an object of class "agg_scores"
#' @param ... additional arguments
#' @param xlab,ylab,col graphics parameters
#'
#' @return A function giving either the cdf or its inverse
#' @export
#'
#' @examples
#' sc <- rbeta(6, 1.5, 3.5)
#' ag <- agg_scores(sc)
#' pag <- cdf(ag)
#' qag <- qdf(ag)
#' pp <- 1:4/5
#' plot(pag)
#' segments(qag(pp), 0, qag(pp), pp, lwd = 0.5)
#' segments(0, pp, qag(pp), pp, lwd = 0.5)
#' plot(qag)
#' segments(pp, 0, pp, qag(pp), lwd = 0.5)
#' segments(0, qag(pp), pp, qag(pp), lwd = 0.5)
cdf <- function(agg, ...) {
  UseMethod("cdf")
}

#' @rdname cdf
#' @export
cdf.agg_scores <- function(agg, ...) {
  boots <- attr(agg, "bootscores")
  y <- sort(boots)
  x <- (seq_along(y) - 0.5)/length(y)
  structure(approxfun(y, x, yleft = 0, yright = 1, rule = 2),
            class = "cdf")
}

#' @rdname cdf
#' @export
qdf <- function(agg, ...) {
  UseMethod("qdf")
}

#' @rdname cdf
#' @export
qdf.agg_scores <- function(agg, ...) {
  boots <- attr(agg, "bootscores")
  y <- sort(boots)
  x <- (seq_along(y) - 0.5)/length(y)
  structure(approxfun(x, y, yleft = 0, yright = 1, rule = 2),
    class = "qdf")
}

#' @rdname cdf
#' @export
plot.cdf <- function(x, ..., xlab = "score", ylab = "F(score)", col = "red") {
  fn <- x
  oldPar <- par(pty = "s", las = 1, xaxs = "i", yaxs = "i")
  on.exit(par(oldPar))
  x <- seq(0, 1, length.out = 1000)
  y <- fn(x)
  plot(x, y, type = "l", xlim = c(0,1), ylim = c(0,1), col = col,
       xlab = xlab, ylab = ylab, ..., panel.first = grid())
  invisible(fn)
}

#' @rdname cdf
#' @export
plot.qdf <- function(x, ..., xlab = "F(score)", ylab = "score", col = "blue") {
  fn <- x
  oldPar <- par(pty = "s", las = 1, xaxs = "i", yaxs = "i")
  on.exit(par(oldPar))
  x <- seq(0, 1, length.out = 1000)
  y <- fn(x)
  plot(x, y, type = "l", xlim = c(0, 1), ylim = c(0, 1), col = col,
       xlab = xlab, ylab = ylab, ..., panel.first = grid())
  invisible(fn)
}
