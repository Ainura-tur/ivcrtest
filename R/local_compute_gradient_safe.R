
#' Gradient of the grid map on a grid, with clamping
#'
#' Vectorised version of \code{\link{g_grad}}. Returns the gradient of
#' \eqn{g(\rho; r)} with respect to (rho_xy, rho_xz, rho_yz) at each value in
#' \code{r_xu}, clamping \eqn{\rho_{xy}} and \eqn{r} inside \eqn{(-1, 1)}.
#' Under the primitive route, \code{\link{check_compatibility}} uses these
#' gradients as the constraint matrices passed to \code{\link{CIhybrid}}.
#'
#' @param r_xu Numeric vector of grid points for \eqn{\rho_{xu}}.
#' @param rho_xy,rho_xz,rho_yz Observed correlations.
#' @param tol Clamping tolerance.
#' @return A matrix with three rows, the derivatives with respect to rho_xy,
#'   rho_xz and rho_yz, and one column per grid point. Non-finite values are
#'   returned as NA.
#' @export
#' @examples
#' local_compute_gradient_safe(c(0, 0.4, 0.8), rho_xy = 0.5, rho_xz = 0.4,
#'                             rho_yz = 0.3)

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
