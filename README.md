# ivcrtest

# Summary of the CR test
We develop a Correlation Restriction (CR) test for assessing the validity of instrumental variables when the direction of endogeneity is known. Consistent with Gunsilius (2020), the test demonstrates that additional structural assumptions can render the instrument validity testable, even in models with continuous endogenous variables. Building on DiTraglia and Garcia-Jimeno (2021), our approach exploits the joint correlation structure among the instrument, regressor, and outcome to construct confidence intervals for partially identified parameters using a modified union-bound procedure. Monte Carlo simulations and empirical applications---including analyses of returns to education, recidivism, and development---demonstrate that the CR test reliably detects invalid instruments and characterises the range of exogeneity, thereby enabling formal empirical evaluation of instrument validity.

# The package
Testing the validity of IVs using the Correlation Restriction
The package computes the union bound CIs (simple and Bei (2024))  membership probability

Pr(r_zu=0 is a in C_n | r_xu in D), 

To install the package:

```r
remotes::install_git("https://github.com/Ainura-tur/ivcrtest.git")
library(ivcrtest)
```

## Reproducing Table 4 (Card 1993)

```r
card <- wooldridge::card
H <- c("exper", "expersq", "black", "south", "smsa", "smsa66",
       paste0("reg66", 1:8))

# Ability bias implies corr(educ, u) > 0, so k = +1 and D lies in [0, 1).
res <- iv_cr_test(card, X = "educ", Y = "lwage", H = H,
                  Z = c("fatheduc", "motheduc", "nearc2", "nearc4"),
                  k = 1, alpha = 0.05, rxu_range = c(0, 0.8), K = 80)
print(res)
```

Each instrument is evaluated on its own complete-case sample, giving
n = 2320, 2657, 3010 and 3010 respectively. Do **not** pre-filter with
`complete.cases(card)`: that drops to n = 1600, because `card` carries
missing values in columns the specification never uses.

`k` sets the sign of the maintained endogeneity restriction and must agree
with `rxu_range`. With `rxu_range = NULL` the domain defaults to `c(0, 0.8)`
for `k = 1` and `c(-0.8, 0)` for `k = -1`.

## Diagnostics

`iv_cr_test()` returns the sign condition of Lemma 2.3, the endpoint-gap
statistic, the grid discretization measures and the eigenvalues of the
reparameterised covariance, all under `res$diagnostics`. `cr_valid_range()`
computes the valid and breakdown ranges for the endogeneity magnitude.

## Method note

Inference uses the exact two-dimensional linear reparameterisation of the
finite-grid map, so the Gaussian simulation involves no estimated coefficient
matrix. Set `route = "primitive"` to use the older gradient-based scheme; the
two differ by o_p(1) and give identical grid-level standard errors.
