
#' Estimates bootstrap covariance of correlation vector
#' @param x Numeric vector.
#' @param y Numeric vector.
#' @param z Numeric vector.
#' @param B Integer number of bootstrap replications.
#' @param seed Integer seed for reproducibility.
#' @return A 3x3 covariance matrix.
#' @export

estimate_cov_corr_boot <- function(x, y, z, B = 800, seed = 123) {
  set.seed(seed)
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
