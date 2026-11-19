# ivcrtest for Python (in development)

This folder holds a Python port, in progress, of the R package ivcrtest. It
will implement the correlation-restriction (CR) test of instrument validity
from Dzhumashev and Tursunalieva (2026), *Econometric Reviews*,
https://doi.org/10.1080/07474938.2026.2741548.

Nothing is usable yet. So far the folder has the package layout, reference
results recorded from the R package, and tests that check those results.

## Layout

```
python/
  pyproject.toml           package metadata (built with hatchling)
  src/ivcrtest/            the package
  tests/                   pytest tests
  tests/fixtures/          reference results from R, as JSON
  tools/make_fixtures.R    the R script that writes the fixtures
```

The R package ignores this folder: `.Rbuildignore` excludes `python/`.

## Running the tests

```
cd python
python -m venv .venv
source .venv/bin/activate        # on Windows: .venv\Scripts\activate
pip install -e ".[test]"
pytest
```

## Reference results from R

The port is tested against version 0.2.2 of the R package. The script
`tools/make_fixtures.R` runs the R functions and saves their inputs and
results. Each fixture is one of three kinds:

* **exact**: deterministic results. The port must match them to about 1e-10.
* **replay**: results that depend on random draws. The script records the
  draws R used, so the port can be fed the same draws. Bootstrap indices come
  from `sample.int()` and are 1-based, one row per replication. Normal draws
  come from `MASS::mvrnorm()`, in call order.
* **tolerance**: complete runs. The port, using its own random numbers, must
  agree within Monte Carlo error and reach the same decisions.

| File | Kind | Covers |
|---|---|---|
| `algebra.json` | exact | `cr_h`, `cr_delta`, `cr_jacobian`, `cr_Amat`, `cr_omega_delta`, `cr_grid_map`, `g_fun`, `g_grad`, `g_xu_safe`, `local_compute_gradient_safe`, the four diagnostics, `solve_c_im` |
| `nearpd.json` | exact | `Matrix::nearPD()` as the bootstrap uses it, with its settings |
| `data.json` | exact | simulated data with missing values; residualisation on controls; sample correlations; `cr_pointwise_range` for a given `Sigma_rho` |
| `bootstrap.json` | replay | `estimate_cov_corr_boot` and `im_stoye_inference`, with their bootstrap indices; the delta-method IM interval (exact) |
| `union_bounds.json` | exact | `CIcon_TNbounds`, `CIcon` and `CIproj_p` for given inputs and draws, on both routes |
| `cihybrid.json` | replay | `CIhybrid` with all its normal draws, on both routes, plus R's optimiser path |
| `end_to_end.json` | tolerance | `iv_cr_test` on simulated data (both routes, with p-values), the sensitivity functions, and the Card (1993) application |

To regenerate the fixtures after changing the R package, install the new
version and run this from the root of the repository:

```
Rscript python/tools/make_fixtures.R
```

The script needs the R packages ivcrtest and jsonlite. It records the Card
application only if wooldridge is installed. It stops if the installed
ivcrtest does not match the version in `DESCRIPTION`. It takes about two
minutes.

## Plan

Port in this order, testing each step against the fixtures:

1. The exact algebra: the reparameterisation, the grid maps, the diagnostics,
   `solve_c_im`, residualisation and correlations, and the pointwise range.
2. The bootstrap covariance, including `nearPD`, and the IM/Stoye interval.
3. The union-bound routines: `CIcon_TNbounds`, `CIcon` and `CIproj_p`, then
   `CIhybrid`.
4. The user-facing functions: `iv_cr_test`, `check_compatibility`, the
   sensitivity functions and a result object.

## Design notes

* **Random numbers.** Functions take a numpy `Generator` or a seed, and never
  set a global seed. All draws go through one small interface, so that tests
  can replay the draws R recorded.
* **The optimiser.** In `CIhybrid`, R calls `optim(method = "L-BFGS-B")`
  without a gradient. R then uses central differences with step 1e-3 and keeps
  5 corrections. SciPy's defaults are forward differences with step 1e-8 and
  10 corrections. The objective is an average of indicators over fixed draws,
  so it changes in jumps, and SciPy's tiny step would usually see a zero
  gradient. The port therefore passes R's gradient rule as `jac` and sets
  `maxcor=5`. The two L-BFGS-B implementations still differ, so `CIhybrid`
  results are compared with a tolerance. `cihybrid.json` records R's path
  through the optimiser, to show where a difference starts.
* **nearPD.** SciPy has no equivalent of `Matrix::nearPD()`. The port will
  implement the same algorithm (Higham 2002) with the same settings.
* **A known difference from R 0.2.2.** In `CIhybrid` (R/CIhybrid.R, line 80),
  R sets `c_bd` after dividing the draws by their standard deviations with
  `/`, which recycles the vector down the rows. Bei's original code and
  ivcrtest's own `CIproj_p` divide column by column. The port follows Bei.
  Once the R package is fixed, regenerate the fixtures.
* **Quantiles.** R's default quantile (type 7) matches numpy's default
  (`method="linear"`).
* **Residualisation.** R's `lm()` handles factor controls and drops collinear
  ones. The port will start with numeric controls.

## Licence

MIT; see `LICENSE`. The union-bound routines derive from Xinyue Bei's
UnionBounds R package, also MIT-licensed; see `COPYRIGHTS`.
