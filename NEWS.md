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