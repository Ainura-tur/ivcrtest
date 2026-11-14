## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

* Possibly misspelled words in DESCRIPTION: Bei, Dzhumashev, Imbens,
  Manski, Stoye and Tursunalieva are surnames; endogeneity and
  reparameterisation are technical terms.

## Test environments

* local: Windows 11, R 4.5.2
* win-builder: R-devel
* macOS builder: macOS 26.6 (Apple M1), R 4.6.1

## Other comments

* The \donttest example of iv_cr_test() uses the suggested package
  wooldridge. It is skipped when wooldridge is not installed, and it
  takes about 15 seconds when it runs.
