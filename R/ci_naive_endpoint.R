# Namespace protection (avoid MASS/stats masking dplyr verbs)
select    <- dplyr::select;      filter    <- dplyr::filter
mutate    <- dplyr::mutate;      slice     <- dplyr::slice
recode    <- dplyr::recode;      rename    <- dplyr::rename
summarise <- dplyr::summarise;   summarize <- dplyr::summarize
arrange   <- dplyr::arrange;     count     <- dplyr::count;   lag <- dplyr::lag

#' Naive selected-endpoint benchmark interval
#'
#' Implements equation (22): the plug-in finite-grid endpoints expanded by
#' pointwise normal critical values \eqn{z_{1-\alpha/2}} using endpoint-specific
#' delta-method grid-level standard errors. This interval treats the selected
#' endpoint indices as fixed. It does not use the full finite-grid union-bound
#' structure and does not apply the Imbens-Manski/Stoye critical-value
#' correction; it is reported only as a benchmark.
#'
#' Grid-level standard errors are computed from the exact representation
#' \eqn{\hat\Sigma_\lambda = A\hat\Omega_\delta A'}. By equation (12) these are
#' numerically identical to the delta-method expression in \eqn{\rho}.
#'
#' @param x Numeric vector, the residualised endogenous regressor.
#' @param y Numeric vector, the residualised outcome.
#' @param z Numeric vector, the residualised instrument.
#' @param rxu_range Numeric vector of length 2 giving the admissible domain D.
#' @param K Integer number of equally spaced grid points. Default 80.
#' @param alpha Numeric significance level.
#' @param Sigma_rho Optional 3x3 sampling covariance of rho_hat on the estimator
#'   scale. Supply the same object used by the MCUB routine so that the two
#'   intervals rest on one bootstrap.
#' @param B_boot Integer bootstrap replications used only when Sigma_rho is NULL.
#' @param seed Integer seed used only when Sigma_rho is NULL.
#' @param tol Numeric clamping tolerance.
#' @return A list with the interval, plug-in endpoints, selected indices, grid
#'   standard errors and the critical value used.
#' @export
#' @importFrom stats complete.cases cor qnorm

ci_naive_endpoint <- function(x, y, z,
                              rxu_range = c(0, 0.8),
                              K         = 80L,
                              alpha     = 0.05,
                              Sigma_rho = NULL,
                              B_boot    = 800L,
                              seed      = 123L,
                              tol       = 1e-10) {

  dat <- data.frame(x = x, y = y, z = z)
  dat <- dat[complete.cases(dat), , drop = FALSE]
  n <- nrow(dat)
  if (n < 30L) stop("ci_naive_endpoint: need at least 30 complete observations.")

  rho_hat <- c(rho_xy = cor(dat$x, dat$y),
               rho_xz = cor(dat$x, dat$z),
               rho_yz = cor(dat$y, dat$z))
  if (any(!is.finite(rho_hat))) {
    stop("ci_naive_endpoint: non-finite sample correlations (check for constant variables).")
  }

  if (is.null(Sigma_rho)) {
    Sigma_rho <- estimate_cov_corr_boot(dat$x, dat$y, dat$z, B = B_boot, seed = seed)
  }
  Sigma_rho <- as.matrix(Sigma_rho)
  if (any(dim(Sigma_rho) != c(3L, 3L))) {
    stop("ci_naive_endpoint: Sigma_rho must be 3x3 for (rho_xy, rho_xz, rho_yz).")
  }

  r_grid  <- seq(rxu_range[1], rxu_range[2], length.out = K)
  Amat    <- cr_Amat(r_grid, tol)
  delta   <- cr_delta(rho_hat, tol)
  Omega   <- cr_omega_delta(rho_hat, Sigma_rho, tol)

  lambda_hat  <- as.numeric(Amat %*% delta)
  Sigma_lam   <- Amat %*% Omega %*% t(Amat)
  sigma_b     <- sqrt(pmax(diag(Sigma_lam), 0))

  b_l <- which.min(lambda_hat)
  b_u <- which.max(lambda_hat)
  plug_in <- c(lambda_hat[b_l], lambda_hat[b_u])

  # Equation (22): pointwise normal critical value at the selected endpoints.
  crit <- qnorm(1 - alpha / 2)

  CI <- c(plug_in[1] - crit * sigma_b[b_l],
          plug_in[2] + crit * sigma_b[b_u])

  list(CI         = CI,
       plug_in    = plug_in,
       r_star     = c(r_grid[b_l], r_grid[b_u]),
       b_l        = b_l,
       b_u        = b_u,
       lambda_hat = lambda_hat,
       sigma_b    = sigma_b,
       rho_hat    = rho_hat,
       delta_hat  = delta,
       Omega_delta = Omega,
       Sigma_rho  = Sigma_rho,
       r_grid     = r_grid,
       crit       = crit)
}

#' Deprecated alias for the naive selected-endpoint benchmark
#'
#' The former \code{ci_simple_union} applied a union-bound critical value
#' \code{qnorm(1 - alpha / (4 * K))}, which is not the construction in equation
#' (22). This alias forwards to \code{\link{ci_naive_endpoint}} and therefore
#' returns a different, narrower interval than the old function did.
#'
#' @inheritParams ci_naive_endpoint
#' @param grid_length Deprecated; use \code{K}.
#' @param cov_method Deprecated and ignored; the bootstrap is always used.
#' @param make_psd Deprecated and ignored.
#' @return See \code{\link{ci_naive_endpoint}}.
#' @export
ci_simple_union <- function(x, y, z,
                            rxu_range   = c(0, 0.8),
                            grid_length = 80L,
                            alpha       = 0.05,
                            Sigma_rho   = NULL,
                            cov_method  = NULL,
                            B_boot      = 800L,
                            seed        = 123L,
                            tol         = 1e-10,
                            make_psd    = TRUE) {
  warning("ci_simple_union() is deprecated and now forwards to ci_naive_endpoint(), ",
          "which uses the z_{1-alpha/2} critical value of equation (22) rather than ",
          "the former union-bound critical value. Results will differ.",
          call. = FALSE)
  ci_naive_endpoint(x, y, z,
                    rxu_range = rxu_range, K = grid_length, alpha = alpha,
                    Sigma_rho = Sigma_rho, B_boot = B_boot, seed = seed, tol = tol)
}
