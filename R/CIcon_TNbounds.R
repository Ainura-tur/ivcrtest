# Adapted from CIcon_TNbounds() in the UnionBounds package by Xinyue Bei
# (https://github.com/xinyuebei-econ/UnionBounds), MIT licence,
# Copyright (c) 2024 Xinyue Bei. See inst/COPYRIGHTS.

#' Truncation points for the conditional test
#'
#' For a hypothesised value \code{theta} and each row of \code{deltahat},
#' computes the truncation points \eqn{\theta_1} and \eqn{\theta_2} of the
#' truncated-normal distribution used by the conditional test of Bei (2024).
#' The test statistic is the larger of the smallest standardised distances
#' from \code{theta} to the lower and upper grid-level bounds, and the
#' truncation conditions on the grid points that attain those minima. Pairs
#' of grid points whose correlation lies within \code{tol_r} of -1 or 1 are
#' screened out; \code{\link{cr_corr_screen}} reports how many pairs this
#' removes.
#'
#' @param theta Numeric scalar, the hypothesised value of the parameter.
#' @param deltahat Matrix of \eqn{\delta} values, one per row: a single row
#'   for the point estimate, or one row per simulated draw.
#' @param Al,Au Matrices with one row per grid point that map \eqn{\delta} to
#'   the lower and upper bounds.
#' @param sigma_l,sigma_u Standard deviations of the lower and upper
#'   grid-level bounds.
#' @param corr_m Correlation matrix between the lower bounds (rows) and the
#'   upper bounds (columns).
#' @param corr_l,corr_u Correlation matrices of the lower bounds and of the
#'   upper bounds.
#' @param cLFl,cLFu Caps on the upper truncation point when the lower or the
#'   upper side binds; \code{\link{CIcon}} uses \code{qnorm(1 - eta)} for both.
#' @param tol_r Tolerance for the pairwise-correlation screen.
#' @param g Function mapping a matrix of \eqn{\delta} values, one per row, to
#'   the grid-level bounds, one column per grid point; for example
#'   \code{cr_grid_map(A)}.
#' @return A list with \code{th_1} and \code{th_2}, the lower and upper
#'   truncation points, each with one element per row of \code{deltahat}.
#' @source Adapted from \code{CIcon_TNbounds()} in the UnionBounds package by
#'   Xinyue Bei, \url{https://github.com/xinyuebei-econ/UnionBounds}, MIT
#'   licence. See the file \file{COPYRIGHTS} installed with ivcrtest.
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
#' # With Al = Au = A, every correlation block is the same matrix.
#' A <- cr_Amat(seq(0, 0.8, length.out = 20))
#' V <- A %*% Omega %*% t(A)
#' s <- sqrt(diag(V))
#' C <- V / outer(s, s)
#' CIcon_TNbounds(theta = 0, deltahat = matrix(cr_delta(rho), nrow = 1),
#'                Al = A, Au = A, sigma_l = s, sigma_u = s,
#'                corr_m = C, corr_l = C, corr_u = C,
#'                cLFl = qnorm(0.999), cLFu = qnorm(0.999), tol_r = 0.001,
#'                g = cr_grid_map(A))
CIcon_TNbounds <- function(theta, deltahat,
                           Al, Au,
                           sigma_l, sigma_u,
                           corr_m, corr_l, corr_u,
                           cLFl, cLFu, tol_r, g) {
  B   <- nrow(deltahat)
  k_l <- nrow(Al)
  k_u <- nrow(Au)

  lambdahat_l <- g(deltahat)
  lambdahat_u <- lambdahat_l

  Tlb <- sweep(lambdahat_l - theta, 2, sigma_l, "/")   # B x k_l
  Tub <- sweep(theta - lambdahat_u, 2, sigma_u, "/")   # B x k_u

  Tl <- apply(Tlb, 1, min)         # length B
  bl <- apply(Tlb, 1, which.min)   # argmin per row (1..k_l)
  Tu <- apply(Tub, 1, min)
  bu <- apply(Tub, 1, which.min)   # argmin per row (1..k_u)

  # corr_m rows/cols at those argmins
  corr_m_blb <- corr_m[bl, , drop = FALSE]           # B x k_u
  corr_m_bub <- t(corr_m[, bu, drop = FALSE])        # B x k_l

  ## ----- Case Tl >= Tu : “lower” side -----
  # TS part
  tTS_l  <- matrix(1e10, nrow = B, ncol = k_u)
  # Broadcast Tl and Tub row-wise
  tTS_ll <- (1 + corr_m_blb)^(-1) * (Tub + corr_m_blb * Tl)
  maskTS_l <- (1 + corr_m_blb) > tol_r
  tTS_l[maskTS_l] <- tTS_ll[maskTS_l]

  th_1_1 <- apply(tTS_l, 1, min) * (Tl >= Tu)
  th_1_1[th_1_1 == 1e10] <- -1e10
  th_1_1 <- th_1_1 * (Tl >= Tu)

  # LF part
  tLFl <- rep(cLFl, B)

  # B part: polyhedral upper truncation for the selected minimum.
  # Constraint T_{l,j} - T_l >= 0 decomposed along the direction of T_l gives
  #   taux_1[i,j] = ( Tlb[i,j] - corr_l[bl[i],j] * Tl[i] ) / ( 1 - corr_l[bl[i],j] ).
  # The numerator must use the j-th statistic Tlb[i,j], not the row minimum.
  # Using the minimum collapses the expression to Tl[i] for every j.
  CL_bl  <- corr_l[bl, , drop = FALSE]                                  # B x k_l
  Tl_rep <- matrix(Tl, nrow = B, ncol = k_l)
  taux_1 <- (1 - CL_bl)^(-1) * (Tlb - CL_bl * Tl_rep)                   # B x k_l

  tB_l2   <- matrix(1e10, nrow = B, ncol = k_l)
  maskB_l <- (1 > CL_bl + tol_r)
  tB_l2[maskB_l] <- taux_1[maskB_l]

  th_2_1 <- apply(cbind(tLFl, tB_l2), 1, min) * (Tl >= Tu)

  ## ----- Case Tl < Tu : “upper” side -----
  # TS part
  tTS_u  <- matrix(1e10, nrow = B, ncol = k_l)
  tTS_u1 <- (1 + corr_m_bub)^(-1) * (Tlb + corr_m_bub * Tu)
  maskTS_u <- (1 + corr_m_bub) > tol_r
  tTS_u[maskTS_u] <- tTS_u1[maskTS_u]

  th_1_2 <- apply(tTS_u, 1, min) * (Tl < Tu)
  th_1_2[th_1_2 == 1e10] <- -1e10
  th_1_2 <- th_1_2 * (Tl < Tu)

  # LF part
  tLFu <- rep(cLFu, B)

  # B part: same construction on the upper side.
  CU_bu  <- corr_u[bu, , drop = FALSE]                                  # B x k_u
  Tu_rep <- matrix(Tu, nrow = B, ncol = k_u)
  taux_2 <- (1 - CU_bu)^(-1) * (Tub - CU_bu * Tu_rep)                   # B x k_u

  tB_u2   <- matrix(1e10, nrow = B, ncol = k_u)
  maskB_u <- (1 > CU_bu + tol_r)
  tB_u2[maskB_u] <- taux_2[maskB_u]

  th_2_2 <- apply(cbind(tLFu, tB_u2), 1, min) * (Tl < Tu)

  th_1 <- th_1_1 + th_1_2
  th_2 <- th_2_1 + th_2_2
  list(th_1 = th_1, th_2 = th_2)
}
