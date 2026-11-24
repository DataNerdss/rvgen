# Tests for the inverse-transform samplers.

# ---- Known-answer / reference cases ---------------------------------------

test_that("r_exp_inv applies the exact formula -log(1 - U) / rate", {
  set.seed(123)
  u <- runif(50)
  set.seed(123)
  x <- r_exp_inv(50, rate = 2)
  expect_equal(x, -log(1 - u) / 2)
  # Same draws as applying R's exponential quantile function to the same U
  expect_equal(x, qexp(u, rate = 2), tolerance = 1e-12)
})

test_that("r_exp_inv has the right mean and variance", {
  set.seed(2026)
  x <- r_exp_inv(1e5, rate = 2)
  expect_equal(mean(x), 0.5, tolerance = 0.02)
  expect_equal(var(x), 0.25, tolerance = 0.03)
})

test_that("r_exp_inv passes a Kolmogorov-Smirnov test against pexp", {
  set.seed(2026)
  x <- r_exp_inv(5000, rate = 3)
  expect_gt(suppressWarnings(ks.test(x, "pexp", rate = 3))$p.value, 0.01)
})

test_that("r_weibull_inv matches qweibull on the same uniforms", {
  set.seed(7)
  u <- runif(50)
  set.seed(7)
  x <- r_weibull_inv(50, shape = 2, scale = 3)
  expect_equal(x, qweibull(u, shape = 2, scale = 3), tolerance = 1e-12)
})

test_that("r_weibull_inv has the right mean and passes a KS test", {
  set.seed(2026)
  x <- r_weibull_inv(1e5, shape = 2, scale = 3)
  expect_equal(mean(x), 3 * gamma(1.5), tolerance = 0.02)
  expect_gt(
    suppressWarnings(ks.test(x[1:5000], "pweibull", shape = 2, scale = 3))$p.value,
    0.01
  )
})

test_that("r_inverse reproduces qexp and a user-written quantile function", {
  set.seed(11)
  u <- runif(20)
  set.seed(11)
  expect_equal(r_inverse(20, qexp, rate = 2), qexp(u, rate = 2))

  qunif_ab <- function(p, a, b) a + (b - a) * p
  set.seed(11)
  expect_equal(r_inverse(20, qunif_ab, a = 10, b = 20), 10 + 10 * u)
})

# ---- Edge cases -----------------------------------------------------------

test_that("n = 0 returns an empty numeric vector", {
  expect_identical(r_exp_inv(0), numeric(0))
  expect_identical(r_weibull_inv(0, shape = 2), numeric(0))
  expect_length(r_inverse(0, qexp), 0)
})

test_that("Weibull with shape = 1 equals Exponential with rate = 1/scale", {
  set.seed(99)
  w <- r_weibull_inv(100, shape = 1, scale = 4)
  set.seed(99)
  e <- r_exp_inv(100, rate = 1 / 4)
  expect_equal(w, e)
})

test_that("all exponential and Weibull draws are non-negative and finite", {
  set.seed(5)
  expect_true(all(is.finite(r_exp_inv(1e4, rate = 0.1))))
  expect_true(all(r_weibull_inv(1e4, shape = 0.5) >= 0))
})

# ---- Invalid inputs -------------------------------------------------------

test_that("invalid n is rejected with a clear message", {
  expect_error(r_exp_inv(-1), "non-negative whole number")
  expect_error(r_exp_inv(2.5), "non-negative whole number")
  expect_error(r_exp_inv(c(1, 2)), "single finite number")
  expect_error(r_exp_inv("10"), "single finite number")
  expect_error(r_exp_inv(NA), "single finite number")
  expect_error(r_exp_inv(Inf), "single finite number")
})

test_that("invalid parameters are rejected with a clear message", {
  expect_error(r_exp_inv(10, rate = 0), "`rate` must be positive")
  expect_error(r_exp_inv(10, rate = -2), "`rate` must be positive")
  expect_error(r_weibull_inv(10, shape = -1), "`shape` must be positive")
  expect_error(r_weibull_inv(10, shape = 2, scale = 0), "`scale` must be positive")
  expect_error(r_weibull_inv(10), "argument \"shape\" is missing")
})

test_that("r_inverse rejects a bad quantile function", {
  expect_error(r_inverse(10, "qexp"), "must be a function")
  expect_error(r_inverse(10, function(p) p[1]), "one value per probability")
  expect_error(suppressWarnings(r_inverse(10, qexp, rate = -1)), "NA or NaN")
})
