# Internal input-checking helpers. Not exported.

# Check that n is a single non-negative whole number.
check_n <- function(n) {
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || !is.finite(n)) {
    stop("`n` must be a single finite number.", call. = FALSE)
  }
  if (n < 0 || n != round(n)) {
    stop("`n` must be a non-negative whole number, not ", n, ".", call. = FALSE)
  }
  invisible(as.integer(n))
}

# Check that a parameter is a single positive finite number.
check_positive <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop("`", name, "` must be a single finite number.", call. = FALSE)
  }
  if (x <= 0) {
    stop("`", name, "` must be positive, not ", x, ".", call. = FALSE)
  }
  invisible(x)
}



