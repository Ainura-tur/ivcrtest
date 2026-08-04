# Namespace protection (avoid MASS/stats masking dplyr verbs)
select    <- dplyr::select;      filter    <- dplyr::filter
mutate    <- dplyr::mutate;      slice     <- dplyr::slice
recode    <- dplyr::recode;      rename    <- dplyr::rename
summarise <- dplyr::summarise;   summarize <- dplyr::summarize
arrange   <- dplyr::arrange;     count     <- dplyr::count;   lag <- dplyr::lag

#' Computes the CIs and probabilities
#' @param r_xu Numeric scalar, value of correlation coefficient rho_xu.
#' @param delta Numeric vector for observed correlations coefficient: r_xy, r_xz, r_yz.
#' @param tol tolerance level
#' @return values of function g(rho, r_xu).
#' @export
g_xu_safe <- function(r_xu, delta, tol = 1e-10) {
  # ---------------------------------------------------------
  # Safe version of g_xu that accepts both vector and matrix delta
  # ---------------------------------------------------------
  if (is.null(dim(delta))) {
    # plain numeric vector of length 3
    delta <- matrix(delta, nrow = 1)
  } else if (ncol(delta) != 3 && nrow(delta) == 3) {
    # 3×1 column -> convert to 1×3 row
    delta <- t(delta)
  } else if (ncol(delta) != 3) {
    stop("delta must have 3 elements: [rho_xy, rho_xz, rho_yz]")
  }

  rho_xy <- delta[, 1]
  rho_xz <- delta[, 2]
  rho_yz <- delta[, 3]

  B <- nrow(delta)
  K1 <- length(r_xu)

  rho_xy_mat <- matrix(rho_xy, nrow = B, ncol = K1)
  rho_xz_mat <- matrix(rho_xz, nrow = B, ncol = K1)
  rho_yz_mat <- matrix(rho_yz, nrow = B, ncol = K1)
  r_xu_mat   <- matrix(r_xu,   nrow = B, ncol = K1, byrow = TRUE)

  # Clip safely
  rho_xy_mat <- pmin(pmax(rho_xy_mat, -1 + tol), 1 - tol)
  r_xu_mat   <- pmin(pmax(r_xu_mat,   -1 + tol), 1 - tol)

  denom <- sqrt(pmax(1 - rho_xy_mat^2, tol))
  sqrt_term <- sqrt(pmax(1 - r_xu_mat^2, 0))

  gval <- rho_xz_mat * r_xu_mat -
    (rho_xy_mat * rho_xz_mat - rho_yz_mat) / denom * sqrt_term

  gval[!is.finite(gval)] <- NA_real_
  return(gval)
}
