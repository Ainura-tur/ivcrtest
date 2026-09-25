# ivcrtest

Testing instrumental variable validity by restricting the sign and magnitude of
endogeneity. R implementation of the Correlation Restriction (CR) test.

The CR test asks whether the orthogonality restriction $\rho_{zu}=0$ stays
compatible with the data once you are willing to restrict the direction and
magnitude of endogeneity to a domain $\mathcal{D}$. Inference uses the modified
conditional union-bound interval of Bei (2024), with Imbens–Manski/Stoye and
naive selected-endpoint intervals reported alongside as diagnostics.

## Install

```r
remotes::install_git("https://github.com/Ainura-tur/ivcrtest.git", ref = "v0.2.0")
library(ivcrtest)
```

Version 0.2.0 is the version used in the paper; drop `ref` for the latest
version. Imports `MASS` and `Matrix`. `wooldridge` is needed only for the
examples below.

## Quick start

Card (1993), four candidate instruments for schooling. About 80 seconds.

```r
card <- wooldridge::card
H <- c("exper", "expersq", "black", "south", "smsa", "smsa66",
       paste0("reg66", 1:8))

# Ability bias implies corr(educ, u) > 0, so k = +1 and D lies in [0, 1).
res <- iv_cr_test(card, X = "educ", Y = "lwage", H = H,
                  Z = c("fatheduc", "motheduc", "nearc2", "nearc4"),
                  k = 1, alpha = 0.05, rxu_range = c(0, 0.8), K = 80,
                  compute_pvalue = FALSE)
print(res)
```

```
Instrumental Variable Correlation Restriction test
  n = varies by instrument, see n, D = [0.00, 0.80], K = 80, alpha = 0.050, route = reparam

        Z    n        plug_in         CI_MCUB           CI_IM        CI_naive
 fatheduc 2320 [0.020, 0.239] [-0.013, 0.274] [-0.013, 0.273] [-0.020, 0.280]
 motheduc 2657 [0.038, 0.240]  [0.005, 0.273]  [0.006, 0.272] [-0.001, 0.278]
   nearc2 3010 [0.033, 0.043]  [0.001, 0.075]  [0.003, 0.074] [-0.003, 0.078]
   nearc4 3010 [0.020, 0.065] [-0.012, 0.095] [-0.011, 0.095] [-0.017, 0.100]
 Zero_in_CI_MCUB Zero_in_CI_IM p_adm_label      G_hat status
               ✓             ✓        <NA> 7.92273343     ok
               ×             ×        <NA> 6.90677141     ok
               ×             ×        <NA> 0.01382591     ok
               ✓             ✓        <NA> 1.32814879     ok

MCUB is the primary procedure; IM/Stoye and the naive endpoint interval
are reported as diagnostics. p_adm is censored to the admissible range.
```

A ✓ means zero lies in the interval, so the test does not reject validity at the
5% level under the maintained domain. Here parental education survives for
fathers but not mothers, and proximity to a four-year college survives while
proximity to a two-year college does not.

The call also prints two warnings per instrument. One notes that the sign
condition of Lemma 2.3 fails, and points to the discretization diagnostics
described below. The other reports the share of grid pairs that the tolerance
`tol_r` screens out of the truncation bounds; the paper uses the same default,
`tol_r = 0.001`.

Each instrument uses its own complete-case sample, which is why `n` differs
across rows. Pass `common_sample = TRUE` to force the intersection instead,
which makes the columns directly comparable at the cost of observations — here
it would drop all four to n = 2220.

`k` sets the sign of the maintained endogeneity restriction and must agree with
`rxu_range`. With `rxu_range = NULL` the domain defaults to `c(0, 0.8)` when
`k = 1` and `c(-0.8, 0)` when `k = -1`; a mismatch is an error rather than a
silent override.

### Inverted-CI p-values

`compute_pvalue = TRUE`, the default, adds `p_adm`, the level at which the
interval first excludes zero, censored to the theorem-backed range
$(0.01, 0.5)$. It costs roughly 100 seconds per instrument, so it is worth
turning off while exploring.

The bootstrap seed depends on an instrument's position in `Z`, so results can
differ slightly with how the instruments are grouped. The Quick start grouping
reproduces the intervals in the paper's Table 2, whose p-values come from one
call per instrument:

```
one call per instrument   fatheduc 0.177   motheduc 0.033   nearc2 0.050   nearc4 0.152
all four in one call      fatheduc 0.177   motheduc 0.031   nearc2 0.049   nearc4 0.158
```

## Sensitivity

How strong an endogeneity assumption can an instrument survive? Both functions
take a data frame of `x`, `y` and `z` already residualised on the controls:

```r
v  <- c("lwage", "educ", "fatheduc", H)
d  <- card[complete.cases(card[, v]), v]
rr <- function(nm) resid(lm(reformulate(H, nm), data = d))
df <- data.frame(x = rr("educ"), y = rr("lwage"), z = rr("fatheduc"))

# largest rho* for which assuming rho_xu >= rho* still fails to reject
cr_switch_point(df, rho_bar = 0.8, k = 1, direction = "below")
#> tau 0.049, status "interior"

# cheaper screen: the single assumed values compatible with validity
cr_pointwise_range(df, rxu_range = c(0, 0.8), k = 1)
#> range [0.000, 0.070], contiguous
```

`cr_switch_point()` bisects over MCUB fits and is defined through the test
itself; it brackets automatically and reports `accepts throughout` or
`rejects throughout` rather than returning a number when no switch exists.
`cr_pointwise_range()` needs no MCUB fit and runs in seconds.

## Diagnostics

`res$diagnostics` holds, per instrument, the sign condition of Lemma 2.3 in the
paper, the endpoint-gap statistic, the grid discretization measures and the
eigenvalues of the reparameterised covariance:

```r
d <- res$diagnostics[["fatheduc"]]
d$sign_condition$holds        # FALSE: h_r = -0.019 with D in (0, 1)
d$endpoint_gap$G              # 7.92
d$discretization$hausdorff    # 0
```

A small `G_hat` flags a near-tie between the selected endpoint and its
neighbour, which happens when the extremum is interior and the map is locally
flat. `nearc2` above is the example, at 0.014.

## Method note

Inference uses the exact two-dimensional linear reparameterisation of the
finite-grid map, so the Gaussian simulation involves no estimated coefficient
matrix. Set `route = "primitive"` for the older gradient-based scheme; the two
differ by $o_p(1)$ and give identical grid-level standard errors.

## Relation to the paper

The Quick start output reproduces the plug-in, MCUB and IM/Stoye columns of
Table 2 in Dzhumashev and Tursunalieva (*Econometric Reviews*, forthcoming),
and the one-call-per-instrument p-values above reproduce its p-value row. The
naive column can differ in the third decimal, because the paper's replication
scripts estimate that interval's covariance from a separate bootstrap draw. The
replication package that accompanies the paper is the authoritative source of
the published numbers; it includes a copy of version 0.2.0 of this package.

## Background

The CR test assesses instrument validity when the direction of endogeneity is
known. Consistent with Gunsilius (2021), it shows that additional structural
assumptions can render validity testable even with continuous endogenous
variables. Building on DiTraglia and García-Jimeno (2021), it exploits the joint
correlation structure among instrument, regressor and outcome to construct
confidence intervals for partially identified parameters through a modified
union-bound procedure. Simulations and applications to returns to education,
recidivism and development show that the test detects invalid instruments and
characterises the range of exogeneity over which validity survives.

## References

Bei, X. (2024). Inference on union bounds. SSRN Working Paper 5023605.
https://doi.org/10.2139/ssrn.5023605

Card, D. (1993). Using geographic variation in college proximity to estimate
the return to schooling. NBER Working Paper 4483.

DiTraglia, F. J. and C. García-Jimeno (2021). A framework for eliciting,
incorporating, and disciplining identification beliefs in linear models.
*Journal of Business & Economic Statistics* 39(4), 1038–1053.

Gunsilius, F. F. (2021). Nontestability of instrument validity under
continuous treatments. *Biometrika* 108(4), 989–995.

Imbens, G. W. and C. F. Manski (2004). Confidence intervals for partially
identified parameters. *Econometrica* 72(6), 1845–1857.

Stoye, J. (2009). More on confidence intervals for partially identified
parameters. *Econometrica* 77(4), 1299–1315.

## Citation

Dzhumashev, R. and A. Tursunalieva. A test for instrumental variable validity
using a partially identified correlation restriction. *Econometric Reviews*,
forthcoming.

## License

MIT. See `LICENSE`.
