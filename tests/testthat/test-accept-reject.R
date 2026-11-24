# Tests for the acceptance-rejection samplers.

# Exact variance of a Beta(a, b), used in several tests
beta_var <- function(a, b) a * b / ((a + b)^2 * (a + b + 1))

# ---- Known-answer / reference cases ---------------------------------------

test_that("r_beta_ar(2, 2) has the right mean, variance and distribution", {
  set.seed(2026)
  x <- r_beta_ar(1e5, 2, 2)
  expect_length(x, 1e5)
  expect_equal(mean(x), 0.5, tolerance = 0.01)
  expect_equal(var(x), beta_var(2, 2), tolerance = 0.03)
  expect_gt(suppressWarnings(ks.test(x[1:5000], "pbeta", 2, 2))$p.value, 0.01)
})

test_that("r_beta_ar works for a skewed Beta(2, 5)", {
  set.seed(2026)
  x <- r_beta_ar(1e5, 2, 5)
  expect_equal(mean(x), 2 / 7, tolerance = 0.01)
  expect_equal(var(x), beta_var(2, 5), tolerance = 0.03)
  expect_gt(suppressWarnings(ks.test(x[1:5000], "pbeta", 2, 5))$p.value, 0.01)
})

test_that("observed acceptance rate is close to 1/M", {
  set.seed(2026)
  x1 <- r_beta_ar(1e5, 2, 2)
  expect_equal(attr(x1, "M"), 1.5)
  expect_equal(attr(x1, "acceptance_rate"), 1 / 1.5, tolerance = 0.01)

  x2 <- r_beta_ar(1e5, 8, 8)
  M8 <- dbeta(0.5, 8, 8)

  expect_equal(attr(x2, "M"), M8)
  expect_equal(attr(x2, "acceptance_rate"), 1 / M8, tolerance = 0.02)
  # Peaked target is less efficient than the mild one
  expect_lt(attr(x2, "acceptance_rate"), attr(x1, "acceptance_rate"))
})

test_that("general sampler: half-normal from an Exponential(1) proposal", {
  f <- function(x) 2 * dnorm(x)
  g <- function(x) dexp(x)
  M <- sqrt(2 * exp(1) / pi)
  set.seed(2026)
  x <- r_accept_reject(5e4, f, rg = function(k) r_exp_inv(k), g = g, M = M)
  expect_equal(mean(x), sqrt(2 / pi), tolerance = 0.01)
  expect_equal(attr(x, "acceptance_rate"), 1 / M, tolerance = 0.02)
  expect_lte(attr(x, "max_ratio"), 1)
  half_normal_cdf <- function(q) 2 * pnorm(q) - 1
  expect_gt(suppressWarnings(ks.test(x[1:5000], half_normal_cdf))$p.value, 0.01)
})

test_that("results are reproducible with a seed", {
  set.seed(42); x <- r_beta_ar(100, 3, 4)
  set.seed(42); y <- r_beta_ar(100, 3, 4)
  expect_identical(x, y)
})

# ---- Edge cases -----------------------------------------------------------

test_that("Beta(1, 1) is flat: M = 1 and every proposal is accepted", {
  set.seed(1)
  x <- r_beta_ar(1000, 1, 1)
  expect_equal(attr(x, "M"), 1)
  expect_equal(attr(x, "acceptance_rate"), 1)

  expect_true(all(x > 0 & x < 1))
})

test_that("mode at a boundary works: Beta(1, 3) has M = 3 at x = 0", {
  set.seed(1)
  x <- r_beta_ar(5e4, 1, 3)
  expect_equal(attr(x, "M"), 3)
  expect_equal(mean(x), 1 / 4, tolerance = 0.02)
})

test_that("n = 0 returns an empty vector", {
  x <- r_beta_ar(0, 2, 2)
  expect_length(x, 0)
  expect_true(is.na(attr(x, "acceptance_rate")))
})

test_that("a bound M that is too small triggers a warning", {
  set.seed(1)
  expect_warning(
    r_accept_reject(1000, function(x) dbeta(x, 8, 8), runif,
                    function(x) rep(1, length(x)), M = 1.5),
    "too small"
  )
})

test_that("the iteration limit stops a sampler that never accepts", {
  never <- function(x) rep(0, length(x))
  expect_error(
    r_accept_reject(10, never, runif, function(x) rep(1, length(x)), M = 1,
                    max_batches = 3),
    "Stopped after 3 batches"
  )

})

# ---- Invalid inputs -------------------------------------------------------

test_that("r_beta_ar rejects invalid shape parameters", {
  expect_error(r_beta_ar(10, 0.5, 2), "at least 1")
  expect_error(r_beta_ar(10, 2, 0.9), "at least 1")
  expect_error(r_beta_ar(10, -1, 2), "`a` must be positive")
  expect_error(r_beta_ar(10, 2, "3"), "single finite number")
  expect_error(r_beta_ar(-5, 2, 2), "non-negative whole number")
})

test_that("r_accept_reject rejects invalid arguments", {
  f <- function(x) dbeta(x, 2, 2)
  g <- function(x) rep(1, length(x))
  expect_error(r_accept_reject(10, "f", runif, g, 1.5), "`f` must be a function")
  expect_error(r_accept_reject(10, f, 5, g, 1.5), "`rg` must be a function")
  expect_error(r_accept_reject(10, f, runif, NULL, 1.5), "`g` must be a function")
  expect_error(r_accept_reject(10, f, runif, g, 0), "`M` must be positive")
  expect_error(r_accept_reject(10, f, runif, g, 0.5), "at least 1")
  expect_error(r_accept_reject(10, f, function(k) runif(1), g, 1.5),
               "length k")
  expect_error(r_accept_reject(10, function(x) -x, runif, g, 1.5),
               "negative, infinite or missing")
  expect_error(r_accept_reject(10, f, runif, function(x) 0 * x, 1.5),
               "`g` must be positive")
})
