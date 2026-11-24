#' Acceptance-rejection sampling
#'
#' Generates draws from a target density \eqn{f} by proposing from an easier
#' density \eqn{g} and keeping only some proposals. A bound \eqn{M} with
#' \eqn{f(x) \le M g(x)} for all \eqn{x} is required. Each proposal
#' \eqn{Y \sim g} is accepted when
#' \deqn{U \le \frac{f(Y)}{M\,g(Y)}, \qquad U \sim \mathrm{Unif}(0,1).}
#' The accepted values have density exactly \eqn{f}, and the probability of
#' accepting a proposal is \eqn{1/M}.
#'
#' @param n Number of draws wanted. A single non-negative whole number.
#' @param f Target density function. Takes a numeric vector, returns
#'   non-negative values.
#' @param rg Proposal sampler: a function of one argument `k` that returns
#'   `k` draws from \eqn{g}.
#' @param g Proposal density function, evaluated at the proposals.
#' @param M Bound with \eqn{f(x) \le M g(x)} everywhere. Must be at least 1
#'   when both `f` and `g` are proper densities.
#' @param max_batches Iteration limit: the maximum number of proposal
#'   batches before giving up. Stops endless loops when almost nothing is
#'   accepted.
#'
#' @return A numeric vector of length `n` with these attributes:
#'   \describe{
#'     \item{`acceptance_rate`}{Accepted proposals divided by all proposals.
#'       Should be close to \eqn{1/M}.}
#'     \item{`n_proposed`}{Total number of proposals generated.}
#'     \item{`M`}{The bound used.}
#'     \item{`max_ratio`}{Largest observed \eqn{f(Y)/(M g(Y))}. A value
#'       above 1 means `M` was too small.}
#'   }
#'   Read them with, for example, `attr(x, "acceptance_rate")`.

#'
#' @details Proposals are generated in batches, sized from the expected
#'   acceptance rate, so the loop usually runs only once or twice. If any
#'   proposal has \eqn{f(Y) > M g(Y)}, the bound is invalid and the output
#'   does not follow \eqn{f}; a warning is given.
#'
#' @seealso [r_beta_ar()] for a ready-made Beta sampler.
#'
#' @examples
#' # Half-normal target, Exponential(1) proposal (sampled by inverse transform)
#' f <- function(x) 2 * dnorm(x)
#' g <- function(x) dexp(x)
#' M <- sqrt(2 * exp(1) / pi)           # about 1.315
#' set.seed(1)
#' x <- r_accept_reject(10000, f, rg = function(k) r_exp_inv(k), g = g, M = M)
#' mean(x)                              # close to sqrt(2/pi), about 0.798
#' attr(x, "acceptance_rate")           # close to 1/M, about 0.760
#' @export
r_accept_reject <- function(n, f, rg, g, M, max_batches = 1000) {
  # ---- 1. Check inputs ----
  n <- check_n(n)
  if (!is.function(f))  stop("`f` must be a function (the target density).", call. = FALSE)
  if (!is.function(rg)) stop("`rg` must be a function (the proposal sampler).", call. = FALSE)
  if (!is.function(g))  stop("`g` must be a function (the proposal density).", call. = FALSE)
  check_positive(M, "M")
  if (M < 1) {
    stop("`M` must be at least 1, because f <= M * g cannot hold for two ",
         "densities when M < 1.", call. = FALSE)
  }
  check_n(max_batches)
  if (max_batches < 1) stop("`max_batches` must be at least 1.", call. = FALSE)


  # ---- 2. Set up counters ----
  out <- numeric(0)      # accepted draws so far
  n_proposed <- 0        # total darts thrown
  n_accepted <- 0        # total darts kept
  max_ratio <- 0         # largest f/(M g) seen, to detect a bad M
  batches <- 0

  # ---- 3. Throw batches of darts until we have n ----
  while (length(out) < n) {
    if (batches >= max_batches) {
      stop("Stopped after ", max_batches, " batches (", n_proposed,
           " proposals) with only ", length(out), " of ", n,
           " draws accepted. Check that `f`, `g` and `M` are correct, ",
           "or increase `max_batches`.", call. = FALSE)
    }
    batches <- batches + 1
    k <- ceiling(1.1 * M * (n - length(out))) + 10   # expected need + margin

    y <- rg(k)                                       # x-positions of darts
    if (!is.numeric(y) || length(y) != k) {
      stop("`rg(k)` must return a numeric vector of length k.", call. = FALSE)
    }
    fy <- f(y)
    gy <- g(y)
    if (any(!is.finite(fy)) || any(fy < 0)) {
      stop("`f` returned negative, infinite or missing values.", call. = FALSE)
    }
    if (any(!is.finite(gy)) || any(gy <= 0)) {
      stop("`g` must be positive and finite wherever `rg` proposes values.",
           call. = FALSE)
    }


    ratio <- fy / (M * gy)                           # height of curve / height of box
    max_ratio <- max(max_ratio, ratio)
    u <- stats::runif(k)
    keep <- u <= ratio                               # under the curve?

    out <- c(out, y[keep])
    n_proposed <- n_proposed + k
    n_accepted <- n_accepted + sum(keep)
  }

  # ---- 4. Warn if the box did not cover the curve ----
  if (max_ratio > 1) {
    warning("f(x) > M * g(x) at some proposals (largest ratio ",
            signif(max_ratio, 4), "). `M` is too small, so the draws do not ",
            "follow the target density.", call. = FALSE)
  }

  # ---- 5. Return exactly n draws, with diagnostics attached ----
  x <- out[seq_len(n)]
  attr(x, "acceptance_rate") <- if (n_proposed > 0) n_accepted / n_proposed else NA_real_
  attr(x, "n_proposed") <- n_proposed
  attr(x, "M") <- M
  attr(x, "max_ratio") <- max_ratio
  x
}

#' Beta random variables by acceptance-rejection
#'
#' Generates Beta\eqn{(a, b)} draws by acceptance-rejection with a
#' Uniform(0, 1) proposal. The bound \eqn{M} is the largest value of the
#' Beta density, reached at the mode
#' \deqn{x^\star = \frac{a-1}{a+b-2}, \qquad M = f(x^\star),}
#' so the expected acceptance rate is \eqn{1/M}.
#'
#' @param n Number of draws. A single non-negative whole number.
#' @param a,b Shape parameters, both at least 1. For values below 1 the
#'   density is unbounded near 0 or 1, so no finite \eqn{M} exists for a
#'   uniform proposal.
#' @param max_batches Iteration limit passed to [r_accept_reject()].
#'
#' @return A numeric vector of length `n` in (0, 1), with the attributes
#'   described in [r_accept_reject()].
#'
#' @details The target density is evaluated with [stats::dbeta()]; all
#'   random numbers come from [stats::runif()]. The mean of
#'   Beta\eqn{(a,b)} is \eqn{a/(a+b)} and the variance is
#'   \eqn{ab / \{(a+b)^2(a+b+1)\}}.
#'
#' @seealso [stats::rbeta()] for R's built-in generator, used as a reference.
#'
#' @examples
#' set.seed(1)
#' x <- r_beta_ar(10000, a = 2, b = 2)
#' mean(x)                      # close to 0.5
#' attr(x, "M")                 # 1.5
#' attr(x, "acceptance_rate")   # close to 1/1.5 = 0.667
#'
#' y <- r_beta_ar(10000, a = 8, b = 8)

#' attr(y, "acceptance_rate")   # close to 1/3.14 = 0.318, about half as efficient
#' @export
r_beta_ar <- function(n, a, b, max_batches = 1000) {
  n <- check_n(n)
  check_positive(a, "a")
  check_positive(b, "b")
  if (a < 1 || b < 1) {
    stop("`a` and `b` must both be at least 1. For a < 1 or b < 1 the Beta ",
         "density is unbounded, so a uniform proposal cannot be used.",
         call. = FALSE)
  }
  mode <- if (a + b > 2) (a - 1) / (a + b - 2) else 0.5   # Beta(1,1) is flat
  M <- stats::dbeta(mode, a, b)

  r_accept_reject(
    n,
    f  = function(x) stats::dbeta(x, a, b),
    rg = function(k) stats::runif(k),
    g  = function(x) rep(1, length(x)),
    M  = M,
    max_batches = max_batches
  )
}
