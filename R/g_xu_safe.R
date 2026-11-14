
#' Grid map on a grid of endogeneity values, with clamping
#'
#' Vectorised version of \code{\link{g_fun}}. Evaluates \eqn{g(\rho; r)} for
#' every row of \code{delta} and every value in \code{r_xu}, clamping
#' \eqn{\rho_{xy}} and \eqn{r} inside \eqn{(-1, 1)}. Under the primitive
#' route, \code{\link{check_compatibility}} passes it to \code{\link{CIhybrid}}
#' as the map \code{g}.
#'
#' @param r_xu Numeric vector of grid points for \eqn{\rho_{xu}}.
#' @param delta Correlations ordered (rho_xy, rho_xz, rho_yz): a vector of
#'   length 3, or a matrix with one such row per draw.
#' @param tol Clamping tolerance.
#' @return A matrix with one row per row of \code{delta} and one column per
#'   grid point. Non-finite values are returned as NA.
#' @export
#' @examples
#' g_xu_safe(seq(0, 0.8, by = 0.2), c(0.5, 0.4, 0.3))
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
