# Namespace protection (avoid MASS/stats masking dplyr verbs)
select    <- dplyr::select;      filter    <- dplyr::filter
mutate    <- dplyr::mutate;      slice     <- dplyr::slice
recode    <- dplyr::recode;      rename    <- dplyr::rename
summarise <- dplyr::summarise;   summarize <- dplyr::summarize
arrange   <- dplyr::arrange;     count     <- dplyr::count;   lag <- dplyr::lag

#' Solves for Imbens-Manski critical value
#' @param alpha Numeric significance level.
#' @param t Numeric ratio of identified set width to max SE.
#' @return Numeric critical value.
#' @export
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
