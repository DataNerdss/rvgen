#' Exponential random variables by inverse transform
#'
#' Generates draws from the Exponential distribution with rate \eqn{\lambda}
#' using the closed-form inverse CDF
#' \deqn{F(x) = 1 - e^{-\lambda x}, \qquad X = -\log(1 - U) / \lambda,}
#' where \eqn{U \sim \mathrm{Unif}(0,1)}.
#'
#' @param n Number of draws. A single non-negative whole number.
#' @param rate Rate parameter \eqn{\lambda > 0}. The mean is \eqn{1/\lambda}.
#'
#' @return A numeric vector of length `n` with non-negative values.
#'
#' @details The mean of the distribution is \eqn{1/\lambda} and the variance
#'   is \eqn{1/\lambda^2}.
#'
#' @seealso [stats::rexp()] for R's built-in generator, used as a reference.
#'
#' @examples
#' set.seed(1)
#' x <- r_exp_inv(10000, rate = 2)
#' mean(x)  # close to 0.5
#' var(x)   # close to 0.25
#' @export
r_exp_inv <- function(n, rate = 1) {
  n <- check_n(n)
  check_positive(rate, "rate")
  u <- stats::runif(n)
  -log(1 - u) / rate
}

#' Weibull random variables by inverse transform
#'
#' Generates draws from the Weibull distribution with shape \eqn{k} and
#' scale \eqn{\lambda} using the closed-form inverse CDF
#' \deqn{F(x) = 1 - e^{-(x/\lambda)^k}, \qquad
#'   X = \lambda\,[-\log(1 - U)]^{1/k},}
#' where \eqn{U \sim \mathrm{Unif}(0,1)}.
#'
#' @param n Number of draws. A single non-negative whole number.
#' @param shape Shape parameter \eqn{k > 0}.
#' @param scale Scale parameter \eqn{\lambda > 0}.
#'
#' @return A numeric vector of length `n` with non-negative values.
#'
#' @details The parameterization matches [stats::rweibull()]. The mean is
#'   \eqn{\lambda\,\Gamma(1 + 1/k)}. With `shape = 1` this is the
#'   Exponential distribution with rate `1 / scale`.
#'
#' @seealso [stats::rweibull()] for R's built-in generator, used as a
#'   reference.
#'
#' @examples
#' set.seed(1)
#' x <- r_weibull_inv(10000, shape = 2, scale = 3)
#' mean(x)                # close to 3 * gamma(1.5), about 2.66
#' @export
r_weibull_inv <- function(n, shape, scale = 1) {
  n <- check_n(n)
  check_positive(shape, "shape")
  check_positive(scale, "scale")
  u <- stats::runif(n)
  scale * (-log(1 - u))^(1 / shape)
}

#' Inverse-transform sampling with any quantile function
#'
#' Generates random draws from a continuous distribution by applying its
#' quantile function (inverse CDF) to uniform random numbers. If
#' \eqn{U \sim \mathrm{Unif}(0,1)} and \eqn{F} is a continuous CDF, then
#' \eqn{X = F^{-1}(U)} has CDF \eqn{F}.
#'
#' @param n Number of draws. A single non-negative whole number.
#' @param qfun A quantile function. Its first argument must be a vector of
#'   probabilities in (0, 1), and it must return one value per probability.
#' @param ... Further arguments passed to `qfun`, such as distribution
#'   parameters.
#'
#' @return A numeric vector of length `n`.
#'
#' @details Assumes the target distribution is continuous and that `qfun`
#'   is its exact inverse CDF. Only `stats::runif()` is used to supply
#'   randomness; the transformation is done by `qfun`.
#'
#' @seealso [r_exp_inv()] and [r_weibull_inv()] for closed-form versions.
#'
#' @examples
#' set.seed(1)
#' # Exponential(rate = 2) using R's quantile function
#' x <- r_inverse(1000, stats::qexp, rate = 2)
#' mean(x)  # close to 1/2
#'
#' # A user-written quantile function: Uniform(a, b)
#' qunif_ab <- function(p, a, b) a + (b - a) * p
#' y <- r_inverse(5, qunif_ab, a = 10, b = 20)
#' @export
r_inverse <- function(n, qfun, ...) {

  n <- check_n(n)
  if (!is.function(qfun)) {
    stop("`qfun` must be a function (a quantile function).", call. = FALSE)
  }
  u <- stats::runif(n)
  x <- qfun(u, ...)
  if (!is.numeric(x) || length(x) != n) {
    stop("`qfun` must return a numeric vector with one value per probability.",
         call. = FALSE)
  }
  if (any(is.na(x))) {
    stop("`qfun` returned NA or NaN values; check its parameters.", call. = FALSE)
  }
  x
}















