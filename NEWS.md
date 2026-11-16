# ivcrtest 0.2.3

* `CIhybrid()` now computes the bound `c_bd` as the paper's Online Appendix
  describes. It divides each column of the simulated draws by that column's
  standard deviation, then takes the 1 - `eta` quantile of each draw's
  largest absolute value. Earlier versions divided the draws by a vector
  that R recycled down the rows, so many entries were scaled by another
  column's standard deviation. The fix matches the standardisation that
  `CIproj_p()` already uses when it checks draws against `c_bd`. `c_bd` sets
  the optimiser's search box and the draws kept in the auxiliary set, so the
  projection interval, and with it the MCUB interval, can move slightly. The
  conditional interval does not use `c_bd` and is unchanged. In the paper's
  applications, no test decision changes.

# ivcrtest 0.2.2

* `iv_cr_test()`, `estimate_cov_corr_boot()` and `im_stoye_inference()` now
  restore the session's random-number state when they exit, as `CIhybrid()`
  already did. Results are unchanged. With `seed = NULL`,
  `estimate_cov_corr_boot()` now draws from the session's current stream
  instead of reseeding at random.
* Every exported function now has an example, and the help pages of the
  lower-level functions describe what they compute. The `iv_cr_test()`
  example no longer fails when wooldridge is not installed.
* Xinyue Bei is now listed as a contributor and copyright holder, because
  `CIhybrid()`, `CIcon()`, `CIcon_TNbounds()` and `CIproj_p()` are adapted
  from her UnionBounds package. See `inst/COPYRIGHTS`.
* The package now cites the published article, Dzhumashev and Tursunalieva
  (2026, *Econometric Reviews*), and Bei (2024), each with its DOI.
  `citation("ivcrtest")` returns the article.

# ivcrtest 0.2.1

* `cr_valid_range()` now returns the valid range of Section 3.5. Version
  0.2.0 took each endpoint from the wrong family of domains and returned
  the outer cut points whenever any domain accepted. The accept/reject
  table was correct. Results in the paper are unaffected.
* `cr_valid_range()` gains `Sigma_rho`, `eps`, explicit MCUB tuning
  arguments, and the fields `full_domain_accepts`, `lower_bracket`,
  `upper_bracket` and `status`. `...` and `one_sided` are removed, and the
  default `seed` is now 123.
* `cr_valid_range()` and `cr_switch_point()` share one MCUB fit routine, so
  with the same inputs the bisection always lies inside the grid bracket.
* `cr_switch_point()` now stops when an MCUB fit fails, instead of
  counting the failure as a rejection.
* `check_compatibility()` and `iv_cr_test()` gain `boot_seed`, which seeds
  the covariance and IM/Stoye bootstraps. The default reproduces 0.2.0.
* Zero-inclusion marks are now UTF-8 in every locale, results gain logical
  columns `zero_in_mcub`, `zero_in_im` and `zero_in_naive`, and `print()`
  shows "yes" and "no" in sessions that are not UTF-8.
* Corrected the help pages of `cr_valid_range()`, `cr_pointwise_range()`,
  `check_compatibility()` and `iv_cr_test()`.