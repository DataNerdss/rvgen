# rvgen

**Random variable generation by inverse transform and acceptance-rejection.**

`rvgen` is an R package written for the STAT 8670 (Computational Methods in
Statistics, Fall 2026) course project. It shows how to turn uniform random
numbers into draws from other distributions using two classical methods, and
evaluates their accuracy and computing cost against R's built-in generators.

**Group:** Group Nine
**Members:** Moses Kwame Hovor, Juliet Seadey

## What the package does

| Function | Method | Target |
|---|---|---|
| `r_inverse(n, qfun, ...)` | Inverse transform, any quantile function | Any continuous distribution |
| `r_exp_inv(n, rate)` | Inverse transform, closed form | Exponential |
| `r_weibull_inv(n, shape, scale)` | Inverse transform, closed form | Weibull |
| `r_accept_reject(n, f, rg, g, M)` | Acceptance-rejection, user-supplied proposal | Any density with a known bound |
| `r_beta_ar(n, a, b)` | Acceptance-rejection, uniform proposal | Beta (a, b ≥ 1) |

All randomness comes from `stats::runif()` (plus a user-supplied proposal
sampler in `r_accept_reject()`); the transformation and accept/reject steps
are written in this package. Acceptance-rejection results carry diagnostics:
the bound `M`, the observed acceptance rate (theory: `1/M`), the number of
proposals, and a check that `f <= M g` held.

## Installation

From GitHub:


```r
# install.packages("remotes")
remotes::install_github("DataNerdss/rvgen", build_vignettes = TRUE)
```

From a local copy of this repository, in R, with the working directory set
to the package folder:

```r
devtools::install(build_vignettes = TRUE)
```

## Minimal example

```r
library(rvgen)
set.seed(8670)

x <- r_exp_inv(10000, rate = 2)
mean(x)                       # about 0.5 = 1/rate

b <- r_beta_ar(10000, a = 8, b = 8)
mean(b)                       # about 0.5
attr(b, "M")                  # 3.14: height of the Beta(8,8) density
attr(b, "acceptance_rate")    # about 0.318 = 1/M
```

The full worked example is in the vignette: `vignette("rvgen")`.

## Running the tests and checks

From the package folder:


```r
devtools::test()     # 69 tests: known-answer, edge-case and invalid-input
devtools::check()    # 0 errors, 0 warnings, 0 notes
```

The check log for the submitted version is in `check-log.txt`.

## Reproducing the report results

From the package folder, after installing the package and restarting R:

```r
source("experiments/run_experiments.R")
```

This writes all tables to `experiments/results/` and all figures to
`experiments/figures/`, using fixed seeds (`set.seed(8670 + r)` for
replication `r`, the same for every method). It takes about 4 minutes.

| Experiment | Question | Output |
|---|---|---|
| 1 | Does `r_beta_ar()` reproduce Beta(2,2) and Beta(8,8) as well as `rbeta()`? (n = 1,000 and 100,000; 200 replications) | `beta_summary.csv`, `beta_raw.csv` |
| 2 | Do `r_exp_inv()` and `r_weibull_inv()` match `rexp()` and `rweibull()`? | `inverse_summary.csv`, `inverse_raw.csv` |
| 3 | How do the acceptance rate and time change as M grows? (Beta(a,a), a = 1 to 64) | `cost_vs_M.csv` |
| 4 | How much faster is batching than a one-proposal-at-a-time loop? (n = 10,000 and 100,000; median of 5 runs) | `loop_vs_batched.csv` |

Figures: `fig1_beta_histograms`, `fig2_beta_qq`, `fig3_ks_pvalues` and
`fig4_cost_vs_M`, each saved as PDF and PNG. Software and hardware details
are saved in `results/session_info.txt`. Timings below about 10 ms are under
the timer resolution and may show as 0.


## References

- Devroye, L. (1986). *Non-Uniform Random Variate Generation*. Springer.
- Rizzo, M. L. (2019). *Statistical Computing with R*, 2nd ed. Chapman & Hall/CRC.
- Ross, S. M. (2013). *Simulation*, 5th ed. Academic Press.

- von Neumann, J. (1951). Various techniques used in connection with random
  digits. *National Bureau of Standards Applied Mathematics Series*, 12, 36–38.

The Beta target density is evaluated with `stats::dbeta()`; comparisons use
`stats::rbeta()`, `stats::rexp()` and `stats::rweibull()` as reference
generators. No external code was copied.

## License

MIT. See `LICENSE.md`.
