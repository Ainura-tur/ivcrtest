# Adapted from CIproj_p() in the UnionBounds package by Xinyue Bei
# (https://github.com/xinyuebei-econ/UnionBounds), MIT licence,
# Copyright (c) 2024 Xinyue Bei. See inst/COPYRIGHTS.

#' Simulated rejection probability for calibrating the projection critical value
#'
#' Estimates, by simulation, the rejection probability of the modified
#' conditional test of Bei (2024) at the point \code{delta} for a candidate
#' projection critical value \code{c}. The draws are \code{delta} plus the rows
#' of \code{deltastar_demean}. Only draws in the auxiliary set count: those
#' whose largest standardised deviation is at most \code{c_bd}. At the lower
#' endpoint, the midpoint and the upper endpoint of the plug-in bound, a draw
#' rejects when both the projection test with critical value \code{c} and the
#' conditional test at level \code{alphac} reject. The result is the larger
#' of two rejection rates, one for the lower endpoint or the midpoint and one
#' for the upper endpoint or the midpoint, plus \code{eta}.
#' \code{\link{CIhybrid}} searches for the \code{c} at which this probability
#' equals \code{alpha} at the least favourable \code{delta}.
#'
#' @param c Candidate critical value of the projection test.
#' @param delta Numeric vector, the point at which the probability is
#'   evaluated.
#' @param alphac Level of the conditional test.
#' @param c_bd Radius of the auxiliary set in standardised units.
#' @param deltastar_demean Matrix of mean-zero normal draws with covariance
#'   \code{deltaSigma}, one per row.
#' @param Al,Au Matrices with one row per grid point that map \eqn{\delta} to
#'   the lower and upper bounds.
#' @param sigma_l,sigma_u Standard deviations of the lower and upper
#'   grid-level bounds.
#' @param deltaSigma Sampling covariance matrix of \eqn{\hat\delta}, used to
#'   standardise the draws.
#' @param corr_m Correlation matrix between the lower bounds (rows) and the
#'   upper bounds (columns).
#' @param corr_l,corr_u Correlation matrices of the lower bounds and of the
#'   upper bounds.
#' @param eta Lower-truncation level, added to the rejection rate.
#' @param tol_r Tolerance for the pairwise-correlation screen in
#'   \code{\link{CIcon_TNbounds}}.
#' @param g Function mapping a matrix of \eqn{\delta} values, one per row, to
#'   the grid-level bounds, one column per grid point; for example
#'   \code{cr_grid_map(A)}.
#' @return Numeric scalar, the estimated rejection probability, or NA when it
#'   is not finite.
#' @source Adapted from \code{CIproj_p()} in the UnionBounds package by Xinyue
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
#' Omega <- cr_omega_delta(rho, estimate_cov_corr_boot(x, y, z, B = 200))
#'
#' A <- cr_Amat(seq(0, 0.8, length.out = 20))
#' V <- A %*% Omega %*% t(A)
#' s <- sqrt(diag(V))
#' C <- V / outer(s, s)
#'
#' draws <- MASS::mvrnorm(200, mu = c(0, 0), Sigma = Omega)
#' dev <- apply(sweep(abs(draws), 2, sqrt(diag(Omega)), "/"), 1, max)
#' c_bd <- quantile(dev, 1 - 0.001)
#'
#' # Rejection probability at the point estimate for c = 1.5
#' CIproj_p(c = 1.5, delta = cr_delta(rho), alphac = 0.04, c_bd = c_bd,
#'          deltastar_demean = draws, Al = A, Au = A,
#'          sigma_l = s, sigma_u = s, deltaSigma = Omega,
#'          corr_m = C, corr_l = C, corr_u = C,
#'          eta = 0.001, tol_r = 0.001, g = cr_grid_map(A))

CIproj_p <- function(c, delta, alphac, c_bd,
                     deltastar_demean,
                     Al, Au,
                     sigma_l, sigma_u,
                     deltaSigma,
                     corr_m, corr_l, corr_u,
                     eta, tol_r, g) {
  # -------------------------------
  # Equivalent to MATLAB CIproj_p.m
  # -------------------------------
  # g : function(delta_matrix) -> lambda_matrix
  # delta, deltastar_demean: row or column vectors/matrices
  # Al, Au: constraint matrices (lower/upper)

  # ----- Lambda at point delta -----
  lambda_l <- g(matrix(delta, nrow = 1))
  lambda_u <- g(matrix(delta, nrow = 1))
  lb <- min(lambda_l)
  ub <- max(lambda_u)
  mb <- (lb + ub) / 2

  # ----- Randomized samples -----
  deltastar <- sweep(deltastar_demean, 2, delta, "+")
  lambdastar_l <- g(deltastar)
  lambdastar_u <- g(deltastar)

  # ----- Condition Δ -----
  delta_sigma <- sqrt(diag(deltaSigma))
  maxdev <- apply(abs(deltastar_demean / matrix(delta_sigma, nrow = nrow(deltastar_demean), ncol = length(delta_sigma), byrow = TRUE)), 1, max)
  ind_Delta <- as.numeric(maxdev <= c_bd)

  # ----- Projection statistics -----
  Tstar_l <- pmax(
    apply((lambdastar_l - lb) / matrix(sigma_l, nrow = nrow(lambdastar_l), ncol = length(sigma_l), byrow = TRUE), 1, min),
    apply((lb - lambdastar_u) / matrix(sigma_u, nrow = nrow(lambdastar_u), ncol = length(sigma_u), byrow = TRUE), 1, min)
  )
  Tstar_m <- pmax(
    apply((lambdastar_l - mb) / matrix(sigma_l, nrow = nrow(lambdastar_l), ncol = length(sigma_l), byrow = TRUE), 1, min),
    apply((mb - lambdastar_u) / matrix(sigma_u, nrow = nrow(lambdastar_u), ncol = length(sigma_u), byrow = TRUE), 1, min)
  )
  Tstar_u <- pmax(
    apply((lambdastar_l - ub) / matrix(sigma_l, nrow = nrow(lambdastar_l), ncol = length(sigma_l), byrow = TRUE), 1, min),
    apply((ub - lambdastar_u) / matrix(sigma_u, nrow = nrow(lambdastar_u), ncol = length(sigma_u), byrow = TRUE), 1, min)
  )

  ind_proj_l <- as.numeric(Tstar_l <= c)
  ind_proj_m <- as.numeric(Tstar_m <= c)
  ind_proj_u <- as.numeric(Tstar_u <= c)

  # ----- Condition cLF -----
  cLF <- qnorm(1 - eta)

  # === Helper function for conditional parts ===
  cond_block <- function(bound, select, lambdastar_l, lambdastar_u, deltastar, sigma_l, sigma_u, label) {
    if (sum(select) == 0) {
      return(rep(0, nrow(deltastar)))
    }
    sub_deltastar <- deltastar[select == 1, , drop = FALSE]
    sub_l_l <- lambdastar_l[select == 1, , drop = FALSE]
    sub_l_u <- lambdastar_u[select == 1, , drop = FALSE]

    Tc1 <- apply((sub_l_l - bound) / matrix(sigma_l, nrow = nrow(sub_l_l), ncol = length(sigma_l), byrow = TRUE), 1, min)
    Tc2 <- apply((bound - sub_l_u) / matrix(sigma_u, nrow = nrow(sub_l_u), ncol = length(sigma_u), byrow = TRUE), 1, min)
    Tc <- pmax(Tc1, Tc2)

    # Call to the truncated-normal bounds function
    th_bounds <- CIcon_TNbounds(bound, sub_deltastar, Al, Au, sigma_l, sigma_u,
                                corr_m, corr_l, corr_u, cLF, cLF, tol_r, g)
    th1 <- th_bounds[[1]]; th2 <- th_bounds[[2]]

    phi <- (pnorm(Tc) - pnorm(th1)) / (pnorm(th2) - pnorm(th1))
    sub_ind <- as.numeric(phi <= 1 - alphac)

    ind_c <- rep(0, nrow(deltastar))
    ind_c[select == 1] <- sub_ind
    ind_c
  }

  # ----- l / m / u conditions -----
  selectl <- (1 - ind_proj_l) * ind_Delta
  ind_c_l <- cond_block(lb, selectl, lambdastar_l, lambdastar_u, deltastar, sigma_l, sigma_u, "l")

  selectm <- (1 - ind_proj_m) * ind_Delta
  ind_c_m <- cond_block(mb, selectm, lambdastar_l, lambdastar_u, deltastar, sigma_l, sigma_u, "m")

  selectu <- (1 - ind_proj_u) * ind_Delta
  ind_c_u <- cond_block(ub, selectu, lambdastar_l, lambdastar_u, deltastar, sigma_l, sigma_u, "u")

  # ----- Combine indicators -----
  ind_l <- 1 - (1 - ind_c_l) * (1 - ind_proj_l)
  ind_m <- 1 - (1 - ind_c_m) * (1 - ind_proj_m)
  ind_u <- 1 - (1 - ind_c_u) * (1 - ind_proj_u)

  # p1 <- mean((1 - ind_l * ind_m) * ind_Delta)
  # p2 <- mean((1 - ind_u * ind_m) * ind_Delta)
  # p <- max(p1, p2) + eta
  p1  <- mean((1 - ind_l * ind_m) * ind_Delta, na.rm = TRUE)
  p2  <- mean((1 - ind_u * ind_m) * ind_Delta, na.rm = TRUE)
  p   <- max(p1, p2, na.rm = TRUE) + eta

  if (!is.finite(p)) p <- NA_real_
  return(p)

  return(p)
}
