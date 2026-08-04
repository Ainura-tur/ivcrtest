# Namespace protection (avoid MASS/stats masking dplyr verbs)
select    <- dplyr::select;      filter    <- dplyr::filter
mutate    <- dplyr::mutate;      slice     <- dplyr::slice
recode    <- dplyr::recode;      rename    <- dplyr::rename
summarise <- dplyr::summarise;   summarize <- dplyr::summarize
arrange   <- dplyr::arrange;     count     <- dplyr::count;   lag <- dplyr::lag

#' Computes gradient of g_xu(r_xu, rho_xy, rho_xz, rho_yz)
#
# with NA/Inf safety.
#
# Arguments:
#   r_xu, rho_xy, rho_xz, rho_yz : numeric scalars or vectors
#   tol : tolerance for numerical stability
#
# Returns:
#   grad : 3 x length(r_xu) matrix (each column is a gradient vector)
# ---------------------------------------------------------
#' @param rho_xy observed correlation values.
#' @param rho_xz observed correlation values.
#' @param rho_yz observed correlation values.
#' @param r_xu assumed correlation  values.
#' @param tol  Numeric scalar for error tolerance.
#' @return 3 x length(r_xu) matrix (each column is a gradient vector).
#' @export

local_compute_gradient_safe <- function(r_xu, rho_xy, rho_xz, rho_yz, tol = 1e-10) {


  # Clamp correlations safely inside (-1, 1)
  rho_xy <- pmin(pmax(rho_xy, -1 + tol), 1 - tol)
  r_xu   <- pmin(pmax(r_xu,   -1 + tol), 1 - tol)

  # Common sqrt term
  sqrt_term <- sqrt(pmax((1 - r_xu^2) / (1 - rho_xy^2), 0))

  # Partial derivatives
  dg_drho_xy <- -rho_xz * sqrt_term -
    (rho_xy * rho_xz - rho_yz) * sqrt_term * rho_xy / (1 - rho_xy^2)

  dg_drho_xz <- r_xu - rho_xy * sqrt_term
  dg_drho_yz <- sqrt_term

  grad <- rbind(dg_drho_xy, dg_drho_xz, dg_drho_yz)
  grad[!is.finite(grad)] <- NA_real_

  return(grad)
}
