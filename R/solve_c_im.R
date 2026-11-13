
#' Imbens-Manski critical value
#'
#' Solves \eqn{\Phi(c + t) - \Phi(-c) = 1 - \alpha} for \eqn{c \ge 0}, where
#' \eqn{t} is the width of the estimated identified set divided by the larger
#' endpoint standard error. With \eqn{t = 0} the solution is
#' \eqn{z_{1-\alpha/2}}; as \eqn{t} grows it falls towards
#' \eqn{z_{1-\alpha}}. \code{\link{im_stoye_inference}} uses it to form the
#' Imbens-Manski/Stoye interval.
#'
#' @param alpha Numeric significance level.
#' @param t Non-negative numeric scalar: the set width divided by the larger
#'   endpoint standard error.
#' @return Numeric scalar, the critical value.
#' @references Imbens, G. W. and C. F. Manski (2004). Confidence intervals for
#'   partially identified parameters. \emph{Econometrica} 72(6), 1845--1857.
#' @export
#' @examples
#' solve_c_im(0.05, 0)    # equals qnorm(0.975)
#' solve_c_im(0.05, 2)
#' solve_c_im(0.05, 10)   # close to qnorm(0.95)
solve_c_im <- function(alpha, t) {
  # Solve: Phi(c+t) - Phi(-c) = 1 - alpha for c >= 0
  # t >= 0 (typically t = Delta_hat / se_max)
  if (!is.finite(t) || t < 0) {
    stop("t must be finite and nonnegative.")
  }

  # For t=0, solution is z_{1-alpha/2}
  if (abs(t) < 1e-12) {
    return(stats::qnorm(1 - alpha / 2))
  }

  # Define objective function
  f <- function(c) stats::pnorm(c + t) - stats::pnorm(-c) - (1 - alpha)

  # Bracket the root
  lo <- 0
  hi <- max(2, stats::qnorm(1 - alpha / 2) + t + 2)

  while (f(hi) < 0 && hi < 1e6) {
    hi <- hi * 2
  }

  if (f(hi) < 0) {
    stop("Failed to bracket the IM critical value root (check inputs).")
  }

  stats::uniroot(f, c(lo, hi), tol = 1e-10)$root
}
