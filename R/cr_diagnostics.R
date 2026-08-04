
#' Diagnostics for the correlation-restriction test
#'
#' Four diagnostics reported in the paper but not previously implemented:
#' the sign condition of Lemmas 2.3 and 2.4, the endpoint-gap statistic of
#' Section 3.4 and Table 7, the grid discretization measures of Remark A.2 and
#' Table 6, and a screen reporting how many grid pairs the \code{tol_r}
#' tolerance removes from the truncated-normal bounds.
#'
#' @param rho Numeric vector of length 3, ordered (rho_xy, rho_xz, rho_yz).
#' @param rxu_range Numeric vector of length 2 giving the admissible domain D.
#' @param K Integer number of equally spaced grid points.
#' @param Sigma_rho 3x3 sampling covariance of rho_hat on the estimator scale.
#' @param tol Numeric clamping tolerance.
#' @name cr_diagnostics
NULL

#' Sign condition and curvature of the grid map
#'
#' Evaluates \eqn{\hat h_r = \hat\rho_{xy}\hat\rho_{xz} - \hat\rho_{yz}} and
#' checks the condition of Lemma 2.3, namely \eqn{h_r > 0} with
#' \eqn{D \subset (0,1)} or \eqn{h_r < 0} with \eqn{D \subset (-1,0)}. When it
#' holds, Lemma 2.4(ii) gives that \eqn{g(\rho;\cdot)} is strictly increasing on
#' D, the extrema are attained at the grid corners, and \eqn{C_K = C} exactly
#' for every \eqn{K \ge 2}. The critical point \eqn{r^*} and the monotonicity
#' margin \eqn{\min_D (\partial g/\partial r)/\rho_{xz}} are also returned.
#'
#' @rdname cr_diagnostics
#' @return \code{cr_sign_condition} returns a list with \code{h_r},
#'   \code{r_star}, \code{holds}, \code{monotonicity_margin} and \code{note}.
#' @export
cr_sign_condition <- function(rho, rxu_range, tol = 1e-10) {
  rho <- as.numeric(rho)
  rho_xy <- min(max(rho[1], -1 + tol), 1 - tol)
  rho_xz <- rho[2]
  rho_yz <- rho[3]

  h_r <- rho_xy * rho_xz - rho_yz
  h   <- h_r / sqrt(1 - rho_xy^2)

  # Critical point r* = c / sqrt(1 + c^2) with c = -rho_xz sqrt(1-rho_xy^2)/h_r.
  r_star <- if (abs(h_r) < .Machine$double.eps) {
    NA_real_
  } else {
    cc <- -rho_xz * sqrt(1 - rho_xy^2) / h_r
    cc / sqrt(1 + cc^2)
  }

  D_pos <- all(rxu_range >= 0)
  D_neg <- all(rxu_range <= 0)
  holds <- (D_pos && h_r > 0) || (D_neg && h_r < 0)

  # dg/dr = rho_xz + h * r / sqrt(1 - r^2), evaluated over D.
  r  <- seq(rxu_range[1], rxu_range[2], length.out = 2001L)
  r  <- pmin(pmax(r, -1 + tol), 1 - tol)
  dg <- rho_xz + h * r / sqrt(1 - r^2)
  margin <- if (abs(rho_xz) < .Machine$double.eps) NA_real_ else min(dg) / rho_xz

  note <- if (holds) {
    "Sign condition holds: g strictly increasing on D, extrema at the grid corners, C_K = C exactly."
  } else if (!D_pos && !D_neg) {
    "D straddles zero; Lemma 2.3 is stated for D contained in (0,1) or in (-1,0)."
  } else {
    "Sign condition fails: Lemma 2.4(iii) applies, one extremum may be interior at r_star."
  }

  list(h_r = h_r, h = h, r_star = r_star, holds = holds,
       r_star_in_D = isTRUE(r_star >= min(rxu_range) && r_star <= max(rxu_range)),
       monotonicity_margin = margin, note = note)
}

#' Standardized endpoint-gap diagnostic
#'
#' Computes \eqn{\hat G_\ell}, \eqn{\hat G_u} and \eqn{\hat G = \min\{\hat
#' G_\ell, \hat G_u\}} of Section 3.4, the standardized separation between each
#' selected endpoint index and its nearest competing grid point. Standard errors
#' of the differences use \eqn{(A_b - A_{\hat b})\hat\Omega_\delta(A_b -
#' A_{\hat b})'}, which equals the gradient-difference expression by equation
#' (12). Large values indicate that Assumption 6 is empirically plausible; the
#' benchmark 1.96 is a familiar numerical reference, not a size-calibrated
#' critical value.
#'
#' @rdname cr_diagnostics
#' @return \code{cr_endpoint_gap} returns a list with \code{G_l}, \code{G_u},
#'   \code{G}, the raw gaps and the selected indices.
#' @export
cr_endpoint_gap <- function(rho, Sigma_rho, rxu_range, K = 80L, tol = 1e-10) {
  r_grid <- seq(rxu_range[1], rxu_range[2], length.out = K)
  A      <- cr_Amat(r_grid, tol)
  Omega  <- cr_omega_delta(rho, Sigma_rho, tol)
  lambda <- as.numeric(A %*% cr_delta(rho, tol))

  b_l <- which.min(lambda)
  b_u <- which.max(lambda)

  se_diff <- function(b, bref) {
    d <- A[b, ] - A[bref, ]
    sqrt(max(as.numeric(t(d) %*% Omega %*% d), 0))
  }

  idx_l <- setdiff(seq_len(K), b_l)
  idx_u <- setdiff(seq_len(K), b_u)

  raw_l <- lambda[idx_l] - lambda[b_l]
  raw_u <- lambda[b_u] - lambda[idx_u]

  s_l <- vapply(idx_l, se_diff, numeric(1), bref = b_l)
  s_u <- vapply(idx_u, se_diff, numeric(1), bref = b_u)

  ok_l <- s_l > 0
  ok_u <- s_u > 0

  G_l <- if (any(ok_l)) min(raw_l[ok_l] / s_l[ok_l]) else NA_real_
  G_u <- if (any(ok_u)) min(raw_u[ok_u] / s_u[ok_u]) else NA_real_

  list(G_l = G_l, G_u = G_u, G = min(G_l, G_u, na.rm = FALSE),
       gap_l = min(raw_l), gap_u = min(raw_u),
       b_l = b_l, b_u = b_u,
       corner_solution = (b_l == 1L && b_u == K) || (b_l == K && b_u == 1L))
}

#' Grid discretization diagnostics
#'
#' Returns the grid spacing \eqn{\delta_K = (r_u - r_\ell)/(K-1)}, the Lipschitz
#' constant \eqn{\hat L = \max_D |\partial g/\partial r|}, the conservative bound
#' \eqn{\hat L \delta_K} of Remark A.2, and the exact Hausdorff distance between
#' the continuous identified set and its finite-grid approximation. Under the
#' sign condition the exact distance is zero to machine precision.
#'
#' @rdname cr_diagnostics
#' @param n_fine Integer number of points used to approximate the continuous
#'   optimum when computing the exact Hausdorff distance.
#' @return \code{cr_discretization} returns a one-row data frame.
#' @export
cr_discretization <- function(rho, rxu_range, K = 80L, n_fine = 200001L, tol = 1e-10) {
  rho <- as.numeric(rho)
  rho_xy <- min(max(rho[1], -1 + tol), 1 - tol)
  h <- (rho_xy * rho[2] - rho[3]) / sqrt(1 - rho_xy^2)

  delta_K <- diff(range(rxu_range)) / (K - 1)

  r_fine <- seq(rxu_range[1], rxu_range[2], length.out = n_fine)
  r_fine <- pmin(pmax(r_fine, -1 + tol), 1 - tol)
  L_hat  <- max(abs(rho[2] + h * r_fine / sqrt(1 - r_fine^2)))

  g_of <- function(r) as.numeric(cr_Amat(r, tol) %*% cr_delta(rho, tol))
  g_K    <- g_of(seq(rxu_range[1], rxu_range[2], length.out = K))
  g_fine <- g_of(r_fine)

  dH <- max(abs(min(g_fine) - min(g_K)), abs(max(g_fine) - max(g_K)))

  data.frame(K = K, delta_K = delta_K, L_hat = L_hat,
             L_delta_K = L_hat * delta_K, hausdorff = dH,
             stringsAsFactors = FALSE)
}

#' Pairwise-correlation screen in the truncated-normal bounds
#'
#' \code{CIcon_TNbounds} discards a grid pair whenever its implied correlation
#' exceeds \code{1 - tol_r}. Lemma 3.1(iii) shows the boundary cases never arise
#' in exact arithmetic, but on a fine grid adjacent points have one minus
#' correlation of order \eqn{\delta_K^2}, so a fixed \code{tol_r} does remove
#' pairs, and the share removed varies with both K and D. This function reports
#' that share so the choice of \code{tol_r} can be made deliberately.
#'
#' @rdname cr_diagnostics
#' @param tol_r Tolerance used by \code{\link{CIcon_TNbounds}}.
#' @return \code{cr_corr_screen} returns a list with the share of off-diagonal
#'   pairs retained and the smallest one-minus-correlation on the grid.
#' @export
cr_corr_screen <- function(rho, Sigma_rho, rxu_range, K = 80L, tol_r = 1e-3, tol = 1e-10) {
  r_grid <- seq(rxu_range[1], rxu_range[2], length.out = K)
  A      <- cr_Amat(r_grid, tol)
  Omega  <- cr_omega_delta(rho, Sigma_rho, tol)
  V      <- A %*% Omega %*% t(A)
  s      <- sqrt(pmax(diag(V), 0))
  C      <- V / outer(s, s)

  off <- !diag(TRUE, K)
  list(share_retained = mean(C[off] < 1 - tol_r),
       share_screened = mean(C[off] >= 1 - tol_r),
       min_one_minus_corr = min(1 - C[off]),
       tol_r = tol_r, K = K)
}
