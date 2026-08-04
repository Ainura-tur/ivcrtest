
#' Inverted-CI membership result for the null value zero
#'
#' Computes \eqn{p_0^{adm} = \inf\{\alpha \in (0.01, 0.5) : 0 \notin
#' CI_{MCUB,1-\alpha}\}} by inverting the MCUB interval over the theorem-backed
#' range only, under the maintained tuning rule \eqn{\alpha_c = 0.8\alpha} with
#' \eqn{\eta} held fixed.
#'
#' Following Remark 3.1, the result is censored rather than extrapolated. If
#' zero is already excluded at the lower boundary the result is reported as
#' \code{"<=lower"}; if zero remains included throughout the admissible range it
#' is reported as \code{">upper"}. These are censored inversion results, not
#' exact p-values of zero or one.
#'
#' The admissible range follows from the tuning conditions of Proposition 3.2:
#' \eqn{\alpha < 1/2} and \eqn{\eta < (\alpha - \alpha_c)/2 = \alpha/10}, so
#' with \eqn{\eta = 0.001} the lower boundary is \eqn{\alpha = 0.01}.
#'
#' @param deltahat,deltaSigma,Al,Au,g Arguments passed to \code{\link{CIhybrid}},
#'   in whichever parameterisation the caller is using.
#' @param eta Lower-truncation confidence level. The admissible lower boundary
#'   is \code{10 * eta} unless \code{alpha_lo} is larger.
#' @param B,Blarge Numbers of Gaussian draws.
#' @param tol,tol_r Tolerances passed to \code{\link{CIhybrid}}.
#' @param alpha_lo,alpha_hi Endpoints of the admissible inversion range.
#' @param alpha_grid Optional bracketing grid; must lie inside the admissible
#'   range. Defaults to a twelve-point grid spanning it.
#' @param refine_steps Integer number of bisection refinement steps.
#' @param seed Integer seed passed to \code{\link{CIhybrid}}.
#' @return A list with \code{p_adm}, a character \code{flag} taking the values
#'   \code{"<=lower"}, \code{"interior"} or \code{">upper"}, a formatted
#'   \code{label}, and the admissible range used.
#' @export
#' @importFrom stats qnorm

pvalue_mcub_zero <- function(deltahat, deltaSigma, Al, Au, g,
                             eta          = 0.001,
                             B            = 300L,
                             Blarge       = 3000L,
                             tol          = 1e-3,
                             tol_r        = 1e-3,
                             alpha_lo     = 0.01,
                             alpha_hi     = 0.5,
                             alpha_grid   = NULL,
                             refine_steps = 8L,
                             seed         = 123L) {

  # Tuning admissibility: eta < (alpha - alphac)/2 with alphac = 0.8 * alpha.
  alpha_lo <- max(alpha_lo, 10 * eta)
  if (!(alpha_lo < alpha_hi)) {
    stop("pvalue_mcub_zero: empty admissible range; check eta, alpha_lo and alpha_hi.")
  }
  if (alpha_hi >= 0.5 + 1e-12) {
    stop("pvalue_mcub_zero: Bei (2024) requires alpha < 1/2.")
  }

  if (is.null(alpha_grid)) {
    # Twelve points strictly inside the admissible range, denser at the bottom
    # where the decision usually switches.
    u <- c(0.002, 0.01, 0.02, 0.04, 0.08, 0.14, 0.22, 0.33, 0.46, 0.62, 0.80, 0.998)
    alpha_grid <- alpha_lo + u * (alpha_hi - alpha_lo)
  }
  alpha_grid <- sort(unique(alpha_grid))
  if (any(alpha_grid <= alpha_lo - 1e-12) || any(alpha_grid >= alpha_hi + 1e-12)) {
    stop("pvalue_mcub_zero: alpha_grid must lie inside [alpha_lo, alpha_hi]; ",
         "evaluations outside the admissible range are computational ",
         "extrapolations with no coverage interpretation.")
  }

  inside0 <- function(alpha) {
    res <- CIhybrid(deltahat, deltaSigma, Al, Au,
                    alpha = alpha, alphac = 0.8 * alpha,
                    eta = eta, B = B, Blarge = Blarge,
                    tol = tol, tol_r = tol_r, index = NULL, g = g,
                    seed = seed)
    CI <- res$CI_h
    isTRUE(CI[1] <= 0 && 0 <= CI[2])
  }

  ins <- vapply(alpha_grid, inside0, logical(1))

  # Zero already excluded at the lower boundary of the admissible range.
  if (!ins[1]) {
    return(list(p_adm = alpha_lo, flag = "<=lower",
                label = sprintf("<=%.2f", alpha_lo),
                alpha_lo = alpha_lo, alpha_hi = alpha_hi,
                alpha_grid = alpha_grid, inside = ins))
  }

  # Zero still included at the top of the admissible range.
  if (ins[length(ins)]) {
    return(list(p_adm = alpha_hi, flag = ">upper",
                label = sprintf(">%.2f", alpha_hi),
                alpha_lo = alpha_lo, alpha_hi = alpha_hi,
                alpha_grid = alpha_grid, inside = ins))
  }

  j  <- which(!ins)[1]
  lo <- alpha_grid[j - 1]   # zero included
  hi <- alpha_grid[j]       # zero excluded

  for (k in seq_len(refine_steps)) {
    mid <- 0.5 * (lo + hi)
    if (inside0(mid)) lo <- mid else hi <- mid
  }

  list(p_adm = hi, flag = "interior",
       label = sprintf("%.3f", hi),
       alpha_lo = alpha_lo, alpha_hi = alpha_hi,
       alpha_grid = alpha_grid, inside = ins)
}

#' Deprecated alias for the inverted-CI membership result
#'
#' The former \code{pvalue_mcub_zero_fast} halved the nominal level before
#' inversion and searched a grid extending outside the theorem-backed range
#' (0.01, 0.5) without censoring. This alias forwards to
#' \code{\link{pvalue_mcub_zero}} and returns the scalar \code{p_adm}.
#'
#' @inheritParams pvalue_mcub_zero
#' @param B_fast,Blarge_fast Deprecated; use \code{B} and \code{Blarge}.
#' @return Numeric scalar \code{p_adm}.
#' @export
pvalue_mcub_zero_fast <- function(deltahat, deltaSigma, Al, Au, g,
                                  eta = 0.001,
                                  B_fast = 300L, Blarge_fast = 3000L,
                                  tol = 1e-3, tol_r = 1e-3,
                                  alpha_grid = NULL,
                                  refine_steps = 8L,
                                  seed = 123L) {
  warning("pvalue_mcub_zero_fast() is deprecated and now forwards to ",
          "pvalue_mcub_zero(), which no longer halves alpha and censors the ",
          "result to the admissible range (0.01, 0.5). Results will differ.",
          call. = FALSE)
  pvalue_mcub_zero(deltahat, deltaSigma, Al, Au, g,
                   eta = eta, B = B_fast, Blarge = Blarge_fast,
                   tol = tol, tol_r = tol_r, alpha_grid = alpha_grid,
                   refine_steps = refine_steps, seed = seed)$p_adm
}
