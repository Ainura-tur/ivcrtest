# Pushing this to github.com/Ainura-tur/ivcrtest

Version 0.2.0. This is a new repository. It supersedes
`github.com/ratbekd/ivcrtest` at commit `4d05337` (2026-07-28), which the
manuscript no longer cites.

The repository name must stay `ivcrtest`, matching the package name, or the
`remotes::install_git` line in Appendix B.1 and in the README will need
changing to match.

## Before the first commit

`man/` in this archive still holds the upstream `.Rd` files and has no entries
for the thirteen new functions. Regenerate it, which also rewrites `NAMESPACE`
from the `@export` tags:

```r
setwd("path/to/ivcrtest")
roxygen2::roxygenise(".", clean = TRUE)
devtools::install(".", upgrade = FALSE)
devtools::test(".")          # 6 blocks, 23 assertions, 0 failures expected
```

Expect 23 `.Rd` files and 31 exports plus `S3method(print, iv_cr_test)`.

## Push

```bash
cd path/to/ivcrtest
git init                      # new repository
git add -A
git commit -m "v0.2.0: exact two-dimensional reparameterisation; correctness fixes

- iv_cr_test(): drop the projection of z onto (x, y); residualise on H only
- iv_cr_test(): forward alpha, rxu_range, k, seed instead of discarding them
- check_compatibility(): K = 80; single shared grid and bootstrap; route switch
- CIcon_TNbounds(): fix the polyhedral truncation, which had collapsed to Tl
- ci_naive_endpoint(): z_{1-alpha/2} in place of the union-bound critical value
- pvalue_mcub_zero(): remove the alpha/2 halving; censor to (0.01, 0.5)
- CIhybrid(): drop library() calls; save and restore the caller's RNG state
- new: reparam_delta, cr_diagnostics, cr_valid_range, cr_pointwise_range
- package now installs on a clean machine; regression tests pinned to Table 4"
git remote add origin https://github.com/Ainura-tur/ivcrtest.git
git branch -M main
git push -u origin main
```

## What changed and why it matters

Six defects in the previous release, five of which are not in the scripts that
produced the paper, which is why the published tables reproduce.

1. `iv_cr_test()` replaced the residualised instrument with its projection onto
   (x, y), twice, the second call idempotent. This put z in span(x, y), making
   the correlation matrix exactly singular and the bootstrap covariance
   numerically singular.
2. `alpha`, `rxu_range`, `k` and `seed` were accepted, documented and then
   ignored; the call to `check_compatibility()` hard-coded them.
3. The grid was K = 50 against K = 80 in the paper.
4. The naive benchmark used `qnorm(1 - alpha/(4*K))`, about 3.60, where the
   paper specifies `z_{1-alpha/2}` = 1.96.
5. The p-value inversion halved alpha and searched outside the admissible range
   without censoring.
6. `check_compatibility()` opened with five `library()` calls to packages in
   Suggests that nothing used, so the package did not run on a clean install.

The sixth defect, in `CIcon_TNbounds()`, is shared with the paper scripts: the
polyhedral truncation used the row minimum rather than the j-th statistic, so
the expression collapsed to `Tl` identically and the conditional component of
the hybrid never bound.

## Verification on a clean install

```
version 0.2.0                             PASS
dplyr verbs not exported                  PASS
31 exports + print.iv_cr_test             PASS
Lemma 3.1(i): A delta == g(rho; r)        PASS   (0 error)
Lemma 3.1(ii): rank(A) == 2               PASS
eq (12): A Omega A' == grad Sigma grad'   PASS   (5.2e-18)
eq (22): naive critical value == 1.96     PASS
Card fatheduc n=2320 [0.020, 0.239]       PASS
Card motheduc n=2657 [0.038, 0.240]       PASS
Card nearc2   n=3010 [0.033, 0.043]       PASS
Card nearc4   n=3010 [0.020, 0.065]       PASS
```

## Note for the co-author

The repository now lives under the Ainura-tur account; `DESCRIPTION` still lists
Ratbek Dzhumashev as maintainer and both authors as `aut`, since the GitHub
account and the package maintainer field are independent. If Ratbek would rather
host it, fork rather than re-upload, so the two copies keep a shared history and
the manuscript URL can be repointed with a single edit.
`dplyr` moved from Suggests to Imports, since the namespace-protection block at
the top of each file in `R/` references it. `ivreg`, `lmtest`, `sandwich` and
`AER` were dropped from Suggests, as nothing in the package uses them.

Consider archiving a release on Zenodo and citing the DOI alongside the URL in
Appendix B.1, so the reference survives any later rename or transfer.
