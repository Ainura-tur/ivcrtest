# ivcrtest

Testing instrumental variable validity by restricting the sign and magnitude of
endogeneity. R implementation of the Correlation Restriction (CR) test.

The CR test asks whether the orthogonality restriction \(\rho_{zu}=0\) remains
compatible with the data once the researcher is willing to restrict the
direction and magnitude of endogeneity to a domain \(\mathcal D\). Inference is
by the modified conditional union-bound interval of Bei (2024), with
Imbens–Manski/Stoye and naive selected-endpoint intervals reported alongside as
diagnostics.

## Install

```r
remotes::install_git("https://github.com/Ainura-tur/ivcrtest.git")
library(ivcrtest)
```

Depends on `MASS`, `Matrix` and `dplyr`; `wooldridge` is needed only for the
example below.

## Quick start

Card (1993), with four candidate instruments for schooling:

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

```
        Z        plug_in         CI_MCUB Zero_in_CI_MCUB p_adm
 fatheduc [0.020, 0.239] [-0.013, 0.274]               ✓ 0.177
 motheduc [0.038, 0.240]  [0.003, 0.272]               × 0.033
   nearc2 [0.033, 0.043]  [0.001, 0.076]               × 0.050
   nearc4 [0.020, 0.065] [-0.010, 0.093]               ✓ 0.152
```

A ✓ means zero lies in the interval, so the test does not reject validity at
the 5% level under the maintained domain.

Two things that trip people up. Each instrument is evaluated on its own
complete-case sample, giving n = 2320, 2657, 3010 and 3010; do **not**
pre-filter with `complete.cases(card)`, which drops to n = 1600 because `card`
carries missing values in columns the specification never uses. And `k` sets the
sign of the maintained endogeneity restriction and must agree with
`rxu_range` — with `rxu_range = NULL` the domain defaults to `c(0, 0.8)` for
`k = 1` and `c(-0.8, 0)` for `k = -1`.

## Sensitivity

How strong an endogeneity assumption can the instrument survive?

```r
# largest rho* such that assuming rho_xu >= rho* still fails to reject
cr_switch_point(df, rho_bar = 0.8, k = 1, direction = "below")

# cheaper screen: single assumed values compatible with validity
cr_pointwise_range(df, rxu_range = c(0, 0.8), k = 1)
```

`df` holds `x`, `y` and `z` already residualised on any controls.
`cr_switch_point()` bisects over MCUB fits and is defined through the test
itself; `cr_pointwise_range()` needs no MCUB fit and agrees closely.

## Diagnostics

`iv_cr_test()` returns, under `res$diagnostics`, the sign condition of
Lemma 2.3, the endpoint-gap statistic, the grid discretization measures and the
eigenvalues of the reparameterised covariance.

## Method note

Inference uses the exact two-dimensional linear reparameterisation of the
finite-grid map, so the Gaussian simulation involves no estimated coefficient
matrix. Set `route = "primitive"` for the older gradient-based scheme; the two
differ by o_p(1) and give identical grid-level standard errors.

## Background

We develop a Correlation Restriction test for assessing the validity of
instrumental variables when the direction of endogeneity is known. Consistent
with Gunsilius (2020), the test demonstrates that additional structural
assumptions can render instrument validity testable, even in models with
continuous endogenous variables. Building on DiTraglia and García-Jimeno
(2021), the approach exploits the joint correlation structure among the
instrument, regressor and outcome to construct confidence intervals for
partially identified parameters using a modified union-bound procedure. Monte
Carlo simulations and empirical applications, including returns to education,
recidivism and development, show that the CR test reliably detects invalid
instruments and characterises the range of exogeneity, enabling formal
empirical evaluation of instrument validity.

## References

Bei, X. (2024). Local linearization based subvector inference in moment
inequality models.

DiTraglia, F. J. and C. García-Jimeno (2021). A framework for eliciting,
incorporating, and disciplining identification beliefs in linear models.

Gunsilius, F. (2020). A path-sampling method to partially identify causal
effects in instrumental variable models.

Imbens, G. W. and C. F. Manski (2004). Confidence intervals for partially
identified parameters. Stoye, J. (2009). More on confidence intervals for
partially identified parameters.

## Citation

Dzhumashev, R. and A. Tursunalieva. A test for instrumental variable validity
using correlation restrictions.
