# =============================================================================
# rvgen: experiments for the STAT 8670 project report
#
# Run from the top folder of the package (the folder containing DESCRIPTION),
# after installing the package with devtools::install():
#
#     source("experiments/run_experiments.R")
#
# Outputs (overwritten on each run):
#   experiments/results/*.csv              tables used in the report
#   experiments/figures/*.pdf, *.png       figures used in the report
#   experiments/results/session_info.txt   software and hardware details
#
# Every experiment uses fixed seeds, so results are reproducible.
# Timing results depend on the computer.
# =============================================================================

library(rvgen)

# ---- Settings ---------------------------------------------------------------

out_res <- file.path("experiments", "results")
out_fig <- file.path("experiments", "figures")
dir.create(out_res, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)

R_REPS    <- 200            # Monte Carlo replications per condition
N_GRID    <- c(1e3, 1e5)    # sample sizes
BASE_SEED <- 8670

# ---- Helper functions -------------------------------------------------------


# Time a single call, in seconds.
time_it <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  force(expr)
  proc.time()[["elapsed"]] - t0
}

# Save the same base-R figure as PDF (for the report) and PNG (for slides).
save_fig <- function(name, width, height, plot_fun) {
  grDevices::pdf(file.path(out_fig, paste0(name, ".pdf")), width = width, height = height)
  plot_fun(); invisible(grDevices::dev.off())
  grDevices::png(file.path(out_fig, paste0(name, ".png")),
                 width = width, height = height, units = "in", res = 150)
  plot_fun(); invisible(grDevices::dev.off())
}

# One replication: generate a sample and record summary measures.
one_rep <- function(gen, cdf, true_mean, true_var) {
  tm <- time_it(x <- gen())
  ks <- suppressWarnings(stats::ks.test(x, cdf))
  acc <- attr(x, "acceptance_rate")
  data.frame(
    mean     = mean(x),
    var      = stats::var(x),
    mean_err = mean(x) - true_mean,
    var_err  = stats::var(x) - true_var,
    ks_stat  = unname(ks$statistic),
    ks_p     = ks$p.value,
    accept   = if (is.null(acc)) NA_real_ else acc,
    time_s   = tm
  )
}


# Run R_REPS replications of one method in one condition.
run_condition <- function(setting, method, n, gen_fun, cdf, true_mean, true_var) {
  rows <- vector("list", R_REPS)
  for (r in seq_len(R_REPS)) {
    set.seed(BASE_SEED + r)            # same seeds for every method
    rows[[r]] <- one_rep(function() gen_fun(n), cdf, true_mean, true_var)
  }
  cbind(setting = setting, method = method, n = n, rep = seq_len(R_REPS),
        do.call(rbind, rows))
}

# Summarise the replications: one row per setting x method x n.
summarise_runs <- function(raw) {
  key <- interaction(raw$setting, raw$method, raw$n, drop = TRUE, lex.order = TRUE)
  out <- lapply(split(raw, key), function(d) data.frame(
    setting        = d$setting[1],
    method         = d$method[1],
    n              = d$n[1],
    mean_bias      = mean(d$mean_err),
    rmse_mean      = sqrt(mean(d$mean_err^2)),
    mean_var       = mean(d$var),
    avg_ks_stat    = mean(d$ks_stat),
    ks_reject_5pct = mean(d$ks_p < 0.05),   # should be close to 0.05
    accept_rate    = mean(d$accept),
    median_time_ms = 1000 * stats::median(d$time_s)
  ))
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}


# ---- Experiment 1: Beta targets, r_beta_ar() vs rbeta() ---------------------

cat("Experiment 1: Beta targets by acceptance-rejection vs rbeta()\n")

beta_settings <- list(
  "Beta(2,2)" = c(a = 2, b = 2),   # Setting 1: mild target,   M = 1.5
  "Beta(8,8)" = c(a = 8, b = 8)    # Setting 2: peaked target, M = 3.14
)

raw_beta <- list()
for (s in names(beta_settings)) {
  a <- beta_settings[[s]][["a"]]; b <- beta_settings[[s]][["b"]]
  tm <- a / (a + b); tv <- a * b / ((a + b)^2 * (a + b + 1))
  cdf <- function(q) stats::pbeta(q, a, b)
  for (n in N_GRID) {
    cat("  ", s, " n =", format(n, scientific = FALSE), "\n")
    raw_beta[[length(raw_beta) + 1]] <- run_condition(
      s, "rvgen::r_beta_ar", n, function(k) r_beta_ar(k, a, b), cdf, tm, tv)
    raw_beta[[length(raw_beta) + 1]] <- run_condition(
      s, "stats::rbeta", n, function(k) stats::rbeta(k, a, b), cdf, tm, tv)
  }
}
raw_beta <- do.call(rbind, raw_beta)
sum_beta <- summarise_runs(raw_beta)
sum_beta$theory_accept <- ifelse(sum_beta$method == "rvgen::r_beta_ar",
                                 ifelse(sum_beta$setting == "Beta(2,2)", 1 / 1.5, 1 / stats::dbeta(0.5, 8, 8)), NA)

utils::write.csv(raw_beta, file.path(out_res, "beta_raw.csv"), row.names = FALSE)
utils::write.csv(sum_beta, file.path(out_res, "beta_summary.csv"), row.names = FALSE)

# ---- Experiment 2: inverse transform vs rexp() and rweibull() ---------------

cat("Experiment 2: inverse transform vs rexp() and rweibull()\n")

raw_inv <- list()
for (n in N_GRID) {
  cat("   n =", format(n, scientific = FALSE), "\n")
  # Exponential(rate = 2): mean 0.5, variance 0.25
  cdf_e <- function(q) stats::pexp(q, 2)
  raw_inv[[length(raw_inv) + 1]] <- run_condition(
    "Exp(2)", "rvgen::r_exp_inv", n, function(k) r_exp_inv(k, 2), cdf_e, 0.5, 0.25)
  raw_inv[[length(raw_inv) + 1]] <- run_condition(
    "Exp(2)", "stats::rexp", n, function(k) stats::rexp(k, 2), cdf_e, 0.5, 0.25)
  # Weibull(shape = 2, scale = 3)
  wm <- 3 * gamma(1.5); wv <- 9 * (gamma(2) - gamma(1.5)^2)
  cdf_w <- function(q) stats::pweibull(q, 2, 3)
  raw_inv[[length(raw_inv) + 1]] <- run_condition(
    "Weibull(2,3)", "rvgen::r_weibull_inv", n, function(k) r_weibull_inv(k, 2, 3), cdf_w, wm, wv)
  raw_inv[[length(raw_inv) + 1]] <- run_condition(
    "Weibull(2,3)", "stats::rweibull", n, function(k) stats::rweibull(k, 2, 3), cdf_w, wm, wv)
}
raw_inv <- do.call(rbind, raw_inv)
sum_inv <- summarise_runs(raw_inv)
utils::write.csv(raw_inv, file.path(out_res, "inverse_raw.csv"), row.names = FALSE)
utils::write.csv(sum_inv, file.path(out_res, "inverse_summary.csv"), row.names = FALSE)

# ---- Experiment 3: cost of acceptance-rejection as M grows ------------------

cat("Experiment 3: cost of acceptance-rejection as the target gets more peaked\n")

# Beta(a, a) for increasing a: the peak M grows, so cost should grow too.
a_grid <- c(1, 2, 4, 8, 16, 32, 64)

n_cost <- 1e5
cost <- do.call(rbind, lapply(a_grid, function(a) {
  M <- stats::dbeta(0.5, a, a)
  reps <- lapply(1:20, function(r) {
    set.seed(BASE_SEED + r)
    tm <- time_it(x <- r_beta_ar(n_cost, a, a))
    c(time = tm, acc = attr(x, "acceptance_rate"),
      props = attr(x, "n_proposed") / n_cost)
  })
  reps <- do.call(rbind, reps)
  data.frame(a = a, M = M, theory_accept = 1 / M,
             accept_rate = mean(reps[, "acc"]),
             proposals_per_draw = mean(reps[, "props"]),
             median_time_ms = 1000 * stats::median(reps[, "time"]))
}))
utils::write.csv(cost, file.path(out_res, "cost_vs_M.csv"), row.names = FALSE)


# ---- Experiment 4: batched sampler vs a one-proposal-at-a-time loop ---------

cat("Experiment 4: batched acceptance-rejection vs a one-at-a-time loop\n")

# Baseline: acceptance-rejection written as a loop with one proposal per
# iteration, as in typical teaching examples. Same target, proposal and M as
# r_beta_ar(), so the only difference is how the proposals are generated.
loop_beta_ar <- function(n, a, b) {
  mode <- if (a + b > 2) (a - 1) / (a + b - 2) else 0.5
  M <- stats::dbeta(mode, a, b)
  out <- numeric(n)
  k <- 0          # draws accepted so far
  tries <- 0      # proposals made so far
  while (k < n) {
    y <- stats::runif(1)
    u <- stats::runif(1)
    tries <- tries + 1
    if (u <= stats::dbeta(y, a, b) / M) {
      k <- k + 1
      out[k] <- y
    }
  }
  attr(out, "acceptance_rate") <- n / tries
  out
}

# Time per call (ms). Fast code is repeated until at least MIN_TIME seconds
# have passed and the total is divided by the number of calls, because the
# system clock (about 10 ms resolution on Windows) cannot time very short
# calls directly. The median over LOOP_REPS such measurements is returned,
# plus the mean and acceptance rate of one run as a correctness check.
LOOP_REPS <- 5
MIN_TIME  <- 0.25
time_method <- function(gen) {
  t <- numeric(LOOP_REPS)
  for (r in seq_len(LOOP_REPS)) {
    set.seed(BASE_SEED + r)
    calls <- 0
    t0 <- proc.time()[["elapsed"]]
    repeat {
      x <- gen()
      calls <- calls + 1
      elapsed <- proc.time()[["elapsed"]] - t0
      if (elapsed >= MIN_TIME) break
    }
    t[r] <- elapsed / calls
    if (r == 1) first <- x
  }
  acc <- attr(first, "acceptance_rate")
  c(time_ms = 1000 * stats::median(t), mean = mean(first),
    accept = if (is.null(acc)) NA_real_ else acc)
}
loop_res <- list()
for (s in names(beta_settings)) {
  a <- beta_settings[[s]][["a"]]; b <- beta_settings[[s]][["b"]]
  for (n in c(1e4, 1e5)) {
    cat("  ", s, " n =", format(n, scientific = FALSE), "\n")
    lp <- time_method(function() loop_beta_ar(n, a, b))
    bt <- time_method(function() r_beta_ar(n, a, b))
    rb <- time_method(function() stats::rbeta(n, a, b))
    loop_res[[length(loop_res) + 1]] <- data.frame(
      setting         = s,
      n               = n,
      loop_ms         = lp[["time_ms"]],
      batched_ms      = bt[["time_ms"]],
      rbeta_ms        = rb[["time_ms"]],
      speedup_vs_loop = lp[["time_ms"]] / bt[["time_ms"]],
      loop_mean       = lp[["mean"]],
      batched_mean    = bt[["mean"]],
      loop_accept     = lp[["accept"]],
      batched_accept  = bt[["accept"]]
    )
  }
}
loop_res <- do.call(rbind, loop_res)
utils::write.csv(loop_res, file.path(out_res, "loop_vs_batched.csv"), row.names = FALSE)


# ---- Figures ----------------------------------------------------------------

cat("Figures\n")

# Figure 1: histograms with the true density, one large sample per setting.
save_fig("fig1_beta_histograms", 8, 3.6, function() {
  op <- graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1))
  for (s in names(beta_settings)) {
    a <- beta_settings[[s]][["a"]]; b <- beta_settings[[s]][["b"]]
    set.seed(BASE_SEED)
    x <- r_beta_ar(1e5, a, b)
    graphics::hist(x, breaks = 60, freq = FALSE, col = "grey85", border = "white",
                   main = sprintf("%s  (M = %.2f, accepted %.1f%%)", s,
                                  attr(x, "M"), 100 * attr(x, "acceptance_rate")),
                   xlab = "x", cex.main = 0.95,
                   ylim = c(0, 1.3 * attr(x, "M")))
    graphics::curve(stats::dbeta(x, a, b), add = TRUE, lwd = 2, col = "firebrick")
    graphics::legend("top", c("r_beta_ar draws", "true density"),
                     pch = c(15, NA), pt.cex = 1.6, lty = c(NA, 1), lwd = c(NA, 2),
                     col = c("grey75", "firebrick"), horiz = TRUE,
                     bty = "n", cex = 0.8)
  }
  graphics::par(op)
})

# Figure 2: QQ plots against the exact Beta quantiles.
save_fig("fig2_beta_qq", 8, 3.6, function() {
  op <- graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1))
  for (s in names(beta_settings)) {
    a <- beta_settings[[s]][["a"]]; b <- beta_settings[[s]][["b"]]
    set.seed(BASE_SEED)
    x <- r_beta_ar(1e4, a, b)

    p <- stats::ppoints(length(x))
    graphics::plot(stats::qbeta(p, a, b), sort(x), pch = 16, cex = 0.3,
                   col = grDevices::adjustcolor("navy", 0.4),
                   xlab = "Theoretical Beta quantile", ylab = "Sample quantile",
                   main = paste("QQ plot,", s), cex.main = 0.95)
    graphics::abline(0, 1, col = "firebrick", lwd = 2)
  }
  graphics::par(op)
})

# Figure 3: KS p-values across replications. If the sampler is correct they are
# Uniform(0,1), exactly as for rbeta().
save_fig("fig3_ks_pvalues", 8, 3.6, function() {
  op <- graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1))
  for (s in names(beta_settings)) {
    d <- raw_beta[raw_beta$setting == s & raw_beta$n == 1e5, ]
    br <- seq(0, 1, by = 0.1)
    h1 <- graphics::hist(d$ks_p[d$method == "rvgen::r_beta_ar"], breaks = br, plot = FALSE)
    h2 <- graphics::hist(d$ks_p[d$method == "stats::rbeta"], breaks = br, plot = FALSE)
    graphics::barplot(rbind(h1$counts, h2$counts), beside = TRUE,
                      names.arg = sprintf("%.1f", br[-1]), col = c("navy", "grey70"),
                      border = NA, xlab = "KS p-value (upper bin edge)",
                      ylab = "Replications",
                      ylim = c(0, 1.45 * max(h1$counts, h2$counts)),
                      main = paste0(s, ", n = 100,000, ", R_REPS, " reps"), cex.main = 0.95)
    # dashed line = count expected in each bin if p-values are Uniform(0,1)
    graphics::abline(h = R_REPS / 10, lty = 2)
    graphics::legend("top", c("r_beta_ar", "rbeta", "expected"),
                     pch = c(15, 15, NA), pt.cex = 1.4, lty = c(NA, NA, 2),
                     col = c("navy", "grey70", "black"), ncol = 3,
                     x.intersp = 0.6, bty = "n", cex = 0.7)
  }

  graphics::par(op)
})

# Figure 4: cost grows with M.
save_fig("fig4_cost_vs_M", 8, 3.6, function() {
  op <- graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1))
  graphics::plot(cost$M, cost$accept_rate, pch = 16, col = "navy", log = "x",
                 xlab = "Bound M (log scale)", ylab = "Acceptance rate",
                 main = "Acceptance rate vs M, Beta(a, a)", cex.main = 0.95)
  graphics::curve(1 / x, add = TRUE, col = "firebrick", lwd = 2)
  graphics::text(cost$M, cost$accept_rate, paste0("a=", cost$a), pos = 4, cex = 0.7)
  graphics::legend("topright", c("observed", "theory 1/M"), pch = c(16, NA),
                   lwd = c(NA, 2), col = c("navy", "firebrick"), bty = "n", cex = 0.8)
  graphics::plot(cost$M, cost$median_time_ms, pch = 16, col = "navy", log = "x",
                 xlab = "Bound M (log scale)", ylab = "Median time (ms)",
                 main = "Time for 100,000 draws vs M", cex.main = 0.95)
  graphics::lines(cost$M, cost$median_time_ms, col = "navy")
  graphics::par(op)
})

# ---- Software and hardware details (for the report) -------------------------

sink(file.path(out_res, "session_info.txt"))
cat("Run date:", format(Sys.time()), "\n")
cat("Replications:", R_REPS, " Base seed:", BASE_SEED, "\n")
cat("CPU cores:", parallel::detectCores(), "\n")
cat("Machine:", paste(Sys.info()[c("sysname", "release", "machine")], collapse = " "), "\n\n")
print(utils::sessionInfo())
sink()

cat("\nDone. Tables in", out_res, " Figures in", out_fig, "\n")
