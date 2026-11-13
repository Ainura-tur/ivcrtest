
#' Bootstrap covariance of the three sample correlations
#'
#' Estimates \eqn{\hat\Sigma_\rho = \widehat{Var}(\hat\rho)}, the sampling
#' covariance of \eqn{\hat\rho = (\hat\rho_{xy}, \hat\rho_{xz},
#' \hat\rho_{yz})'}, by the nonparametric bootstrap. Rows with a missing value
#' in \code{x}, \code{y} or \code{z} are dropped first. The estimate is on the
#' estimator scale, so it is not divided by n. It is symmetrised and, if needed,
#' moved to the nearest positive semi-definite matrix with
#' \code{\link[Matrix]{nearPD}}.
#'
#' @param x Numeric vector, the residualised endogenous regressor.
#' @param y Numeric vector, the residualised outcome.
#' @param z Numeric vector, the residualised instrument.
#' @param B Integer number of bootstrap replications.
#' @param seed Integer seed for the bootstrap draws, or NULL to use the
#'   current random-number stream. The caller's random-number state is
#'   restored on exit.
#' @return A 3x3 covariance matrix, ordered (rho_xy, rho_xz, rho_yz).
#' @export
#' @examples
#' set.seed(1)
#' n <- 300
#' x <- rnorm(n); y <- 0.5 * x + rnorm(n); z <- 0.4 * x + rnorm(n)
#' estimate_cov_corr_boot(x, y, z, B = 200)

estimate_cov_corr_boot <- function(x, y, z, B = 800, seed = 123) {
  .cr_local_seed(seed)
  dat <- cbind(x, y, z)
  dat <- dat[complete.cases(dat), , drop = FALSE]
  n <- nrow(dat)

  Rb <- matrix(NA_real_, nrow = B, ncol = 3)
  for (b in 1:B) {
    idx <- sample.int(n, n, replace = TRUE)
    xb <- dat[idx, 1]; yb <- dat[idx, 2]; zb <- dat[idx, 3]
    Rb[b, ] <- c(cor(xb, yb), cor(xb, zb), cor(yb, zb))
  }

  S <- cov(Rb, use = "complete.obs")
  S <- (S + t(S)) / 2  # symmetrize

  # optional: make PSD in case bootstrap noise makes it slightly indefinite
  if (requireNamespace("Matrix", quietly = TRUE)) {
    S <- as.matrix(Matrix::nearPD(S, corr = FALSE)$mat)
  }
  S
}
