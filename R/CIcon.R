# Adapted from CIcon() in the UnionBounds package by Xinyue Bei
# (https://github.com/xinyuebei-econ/UnionBounds), MIT licence,
# Copyright (c) 2024 Xinyue Bei. See inst/COPYRIGHTS.

#' Conditional confidence interval for a union bound
#'
#' Computes the conditional confidence interval of Bei (2024) for a parameter
#' whose identified set is the union bound \eqn{[\min_b \lambda_{l,b},
#' \max_b \lambda_{u,b}]}. Each endpoint starts at the plug-in bound widened
#' by \code{c_LF} standard errors and moves inward in steps of \code{tol}
#' until the conditional test stops rejecting or the endpoint reaches the
#' plug-in bound. The test uses the truncation points from
#' \code{\link{CIcon_TNbounds}}. \code{\link{CIhybrid}} combines this interval
#' with the calibrated projection interval.
#'
#' @param deltahat Numeric vector, the point estimate \eqn{\hat\delta}.
#' @param deltaSigma Sampling covariance matrix of \code{deltahat}.
#' @param Al,Au Matrices with one row per grid point that map \eqn{\delta} to
#'   the lower and upper bounds; equal in the degenerate-bound case.
#' @param c_LF Critical value used to widen the plug-in bounds and to cap the
#'   truncation points; \code{\link{CIhybrid}} uses \code{qnorm(1 - eta)}.
#' @param alphac Level of the conditional test.
#' @param tol Step by which each endpoint moves inward.
#' @param tol_r Tolerance for the pairwise-correlation screen in
#'   \code{\link{CIcon_TNbounds}}.
#' @param g Function mapping a matrix of \eqn{\delta} values, one per row, to
#'   the grid-level bounds, one column per grid point; for example
#'   \code{cr_grid_map(A)}.
#' @return Numeric vector of length 2, the lower and upper limits of the
#'   conditional interval.
#' @source Adapted from \code{CIcon()} in the UnionBounds package by Xinyue
#'   Bei, \url{https://github.com/xinyuebei-econ/UnionBounds}, MIT licence.
#'   See the file \file{COPYRIGHTS} installed with ivcrtest.
#' @references Bei, X. (2024). Inference on union bounds. SSRN Working Paper
#'   5023605. Version of October 6, 2025. Available at
#'   \url{https://ssrn.com/abstract=5023605} or \doi{10.2139/ssrn.5023605}.
#'   Implemented in the UnionBounds R package.
#' @export
#' @examples
#' set.seed(1)
#' n <- 500
#' z <- rnorm(n); v <- rnorm(n); u <- 0.5 * v + rnorm(n)
#' x <- 0.6 * z + v; y <- 1 + 0.5 * x + u
#' rho <- c(cor(x, y), cor(x, z), cor(y, z))
#' Sigma_rho <- estimate_cov_corr_boot(x, y, z, B = 200)
#'
#' A <- cr_Amat(seq(0, 0.8, length.out = 20))
#' CIcon(cr_delta(rho), cr_omega_delta(rho, Sigma_rho), A, A,
#'       c_LF = qnorm(1 - 0.001), alphac = 0.04, tol = 0.001, tol_r = 0.001,
#'       g = cr_grid_map(A))

CIcon <- function(deltahat, deltaSigma, Al, Au,
                  c_LF, alphac, tol, tol_r, g) {
  # -------------------------------
  # Equivalent to MATLAB CIcon.m
  # -------------------------------

  kk <- nrow(Al)

  # ----- 1. Lambda covariance and correlations -----
  Lambda <- rbind(Al, Au)
  lambdaSigma <- Lambda %*% deltaSigma %*% t(Lambda)
  lambda_sigma <- sqrt(diag(lambdaSigma))

  corr_all <- diag(1 / lambda_sigma) %*% lambdaSigma %*% diag(1 / lambda_sigma)
  corr_m <- corr_all[1:kk, (kk + 1):(2 * kk)]
  corr_l <- corr_all[1:kk, 1:kk]
  corr_u <- corr_all[(kk + 1):(2 * kk), (kk + 1):(2 * kk)]

  sigma_l <- lambda_sigma[1:kk]
  sigma_u <- lambda_sigma[(kk + 1):(2 * kk)]

  # ----- 2. Compute lambda bounds -----
  gval <- as.numeric(g(matrix(deltahat, nrow = 1)))

  lb <- gval - sigma_l * c_LF
  ub <- gval + sigma_u * c_LF
  lb <- min(lb)
  ub <- max(ub)

  mid <- (lb + ub) / 2
  lb_pt <- min(gval)
  ub_pt <- max(gval)

  # ----- 3. CI lower bound -----
  rej <- 1
  theta <- lb

  while (rej == 1 && theta <= min(mid, lb_pt)) {
    Tcl <- min((gval - theta) / sigma_l)
    Tcu <- min((theta - gval) / sigma_u)
    Tc <- max(Tcl, Tcu)

    th_bounds <- CIcon_TNbounds(theta, matrix(deltahat, nrow = 1),
                                Al, Au, sigma_l, sigma_u,
                                corr_m, corr_l, corr_u,
                                c_LF, c_LF, tol_r, g)
    th_1 <- th_bounds[[1]]
    th_2 <- th_bounds[[2]]

    t <- qnorm((1 - alphac) * pnorm(th_2) + alphac * pnorm(th_1))
    rej <- as.numeric(Tc > t)
    theta <- theta + tol
  }
  CI_lower <- theta

  # ----- 4. CI upper bound -----
  rej <- 1
  theta <- ub

  while (rej == 1 && theta >= max(mid, ub_pt)) {
    Tcl <- min((gval - theta) / sigma_l)
    Tcu <- min((theta - gval) / sigma_u)
    Tc <- max(Tcl, Tcu)

    th_bounds <- CIcon_TNbounds(theta, matrix(deltahat, nrow = 1),
                                Al, Au, sigma_l, sigma_u,
                                corr_m, corr_l, corr_u,
                                c_LF, c_LF, tol_r, g)
    th_1 <- th_bounds[[1]]
    th_2 <- th_bounds[[2]]

    t <- qnorm((1 - alphac) * pnorm(th_2) + alphac * pnorm(th_1))
    rej <- as.numeric(Tc > t)
    theta <- theta - tol
  }
  CI_upper <- theta

  # ----- 5. Return final CI -----
  CIcon <- c(CI_lower, CI_upper)
  return(CIcon)
}
