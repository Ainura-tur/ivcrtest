
#' Gradient of the grid map with respect to the observed correlations
#'
#' Returns the gradient of \code{\link{g_fun}} with respect to
#' (rho_xy, rho_xz, rho_yz) at a single value of \eqn{\rho_{xu}}.
#' \code{\link{local_compute_gradient_safe}} is a vectorised version that
#' clamps its inputs.
#'
#' @param r Numeric scalar, the assumed value of \eqn{\rho_{xu}}.
#' @param rho_xy,rho_xz,rho_yz Observed correlations.
#' @return Numeric vector of length 3: the derivatives with respect to
#'   rho_xy, rho_xz and rho_yz.
#' @export
#' @examples
#' g_grad(0.4, rho_xy = 0.5, rho_xz = 0.4, rho_yz = 0.3)

g_grad <- function(r, rho_xy, rho_xz, rho_yz) {
  S <- sqrt((1 - r^2) / (1 - rho_xy^2))
  c(
    -rho_xz * S - (rho_xy * rho_xz - rho_yz) * S * rho_xy / (1 - rho_xy^2), # ∂/∂rho_xy
    r - rho_xy * S,                                                         # ∂/∂rho_xz
    S                                                                       # ∂/∂rho_yz
  )
}
