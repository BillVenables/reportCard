bc_old <- function(y, a, eps = 1.0e-4) {
  if(abs(a) > eps) {
    (y^a - 1)/a
  } else {
    ly <- log(y)
    aly <- a*ly
    ly*(1 + aly*(1 + aly*(1 + aly*(1 + aly*(1 + aly*(1 + aly/7)/6)/5)/4)/3)/2)
  }
}

bc <- function(y, a, eps = 1.0e-4) { ## vectorized wrt y and a
  n <- max(length(y), length(a))
  a <- rep_len(a, length.out = n)
  # y <- rep_len(y, length.out = n)
  ly <- log(y)
  aly <- a*ly
  ifelse(abs(a) > eps,  (exp(aly) - 1)/a,
         ly*(1 + aly*(1 + aly*(1 + aly*(1 + aly*(1 + aly*(1 + aly/7)/6)/5)/4)/3)/2))
}

bci <- function(z, a, eps = 1e-5) {
  n <- max(length(z), length(a))
  a <- rep_len(a, length.out = n)
  ifelse(abs(a) > eps, (1 + a*z)^(1/a),
    exp(z) * (1 - (z^2*a)/2 + ((3*z^4 + 8*z^3)*a^2)/24 -
                ((z^6 + 8*z^5 + 12*z^4)*a^3)/48 +
                ((15*z^8 + 240*z^7 + 1040*z^6 + 1152*z^5)*a^4)/5760 -
                ((3*z^10 + 80*z^9 + 680*z^8 + 2112*z^7 + 1920*z^6)*a^5)/11520))
}

find_as <- function(low, upp, ref = 1, s_low = 0.2, s_upp = 0.8,
                    qfn = qnorm) {
  stopifnot(low < ref, upp > ref, s_low < 0.5, s_upp > 0.5)
  y1 <- low/ref
  y2 <- upp/ref
  z1 <- qfn(s_low)
  z2 <- qfn(s_upp)
  a <- uniroot(function(a) z2*bc(y1, a) - z1*bc(y2, a), lower = -5, upper = 5)$root
  list(a = a, s = (bc(y1, a)/z1 + bc(y2, a)/z2)/2)
}

#' Find cross points
#'
#' Find the points in the value scale corresponding to grade boundaries
#' in the score scale, assuming the naive value to score transformation
#'
#'
#' @param guideline the guideline value, assumed corresponding to 0.5 score
#' @param score the score boundaries
#' @param increasing logical: is the score transformation increasing?
#'
#' @return A set of grade boundaries in the value scale
#' @export
#'
#' @examples
#' cross_points(10, c(0.5, (0:5)/5))
#' cross_points(10, c(0.5, (0:5)/5), TRUE)
cross_points <- function(guideline, score = c(0.2, 0.8), increasing = FALSE) {
  stopifnot(is.numeric(guideline), length(guideline) == 1)
  stopifnot(is.numeric(score), all(score >= 0), all(score <= 1))
  stopifnot(is.logical(increasing), length(increasing) == 1)
  score <- sort(score, decreasing = !increasing)
  mult <- if(increasing) 2^(2*score - 1) else 2^(1 - 2*score)
  structure(guideline*mult, names = format(score))
}

#' Scoring function maker
#'
#' Make a scoring function to enact a prescribed scoring scheme.
#'
#' @param A,E  cutoffs for the A- and E-grads in the value scale
#' @param guideline Guideline valuep
#' @param s_A,s_E cutoffs for the A- and E-grades in the score scale
#' @param pfn,qfn functions to enact the score transformation, and its inverse
#' @param inverse logical: should the function be the inverse mapping, score -> value?
#'
#' @return A function to enact the scoring scheme as determined by guideline and cutoffs
#' @export
#'
#' @examples
#' Ref <- 10
#' r <- cross_points(Ref)
#' fn <- make_score_function(r[1], r[2], Ref, pfn = plogis, qfn = qlogis)
#' oldPar <- par(yaxs = "i", mar = c(5,4,2,2)+0.1)
#' curve(fn, xlim = c(Ref/5, 2.5*Ref), ylim = 0:1,
#'       xlab = "value", ylab = "score")
#' abline(lty="dashed", h = (1:4)/5, v = r)
#' curve(1-(log(x/10, base = 2) + 1)/2, add = TRUE, col="blue", xlim = c(5,20))
#' text(x = rep(2, 5), y = (seq(1, 9, 2))/10, labels = LETTERS[5:1], font = 2)
#'
#' curve(fn, xlim = c(Ref/4,4*Ref), ylim = 0:1, log="x",
#'       xlab = "value (log scale)", ylab = "score")
#' abline(lty="dashed", h = (1:4)/5, v = r)
#' curve(1-(log(x/10, base = 2) + 1)/2, add = TRUE, col="red", xlim = c(5,20))
#' text(x = rep(2.5, 5), y = (seq(1, 9, 2))/10, labels = LETTERS[5:1], font = 2)
#' par(oldPar)
make_score_function <- function(A, E, guideline,
                                s_A = 0.8, # ifelse(increasing, 0.8, 0.2),
                                s_E = 0.2, # ifelse(increasing, 0.2, 0.8),
                                pfn = pnorm, qfn = qnorm, inverse = FALSE) {
  positive <- function(x) is.numeric(x) && length(x) == 1 && x > 0
  stopifnot(positive(A), positive(E), positive(guideline))
  increasing <- A > E
  if(increasing) {
    stopifnot(E < guideline && guideline < A)
  } else {
    stopifnot(A < guideline && guideline < E)
  }
  stopifnot(s_E < 0.5 && 0.5 < s_A)
  a_and_s <- if(increasing) {
    find_as(E, A, guideline, s_E, s_A, qfn = qfn)
  } else {
    find_as(A, E, guideline, 1-s_A, 1-s_E, qfn = qfn)
  }
  alpha <- a_and_s$a
  scale <- a_and_s$s
  force(pfn)
  if(inverse) {
    function(score) {
      qz <- qfn(score, lower.tail = increasing)*scale
      guideline * bci(qz, alpha)
    }
  } else {
    function(value) {
      z <- bc(value/guideline, alpha)/scale
      pfn(z, lower.tail = increasing)
    }
  }
}
