#' ivcrtest: Correlation-Restriction Test for Instrument Validity
#'
#' Tests whether the orthogonality restriction on an instrument remains
#' compatible with the data once the sign and magnitude of endogeneity are
#' restricted.
#'
#' @references Dzhumashev, R. and A. Tursunalieva (2026). A test for
#'   instrumental variable validity using a partially identified correlation
#'   restriction. \emph{Econometric Reviews}.
#'   \doi{10.1080/07474938.2026.2741548}
#' @keywords internal
#' @importFrom stats as.formula coef complete.cases cor cov lm median optim
#' @importFrom stats pnorm qnorm quantile resid rnorm sd var
#' @importFrom MASS mvrnorm
"_PACKAGE"
