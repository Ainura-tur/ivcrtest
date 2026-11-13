
#' Evaluate a grid map at every grid point
#'
#' Calls \code{g_fun(rho_vec, grid_r_xu)}. If \code{g_fun} is not vectorised
#' over the grid, it is evaluated one grid point at a time instead.
#'
#' @param g_fun Function \code{g_fun(rho, r)} returning the grid map at the
#'   correlation vector \code{rho} and the grid point or points \code{r}.
#' @param rho_vec Numeric vector of correlations, ordered (rho_xy, rho_xz,
#'   rho_yz).
#' @param grid_r_xu Numeric vector of grid points.
#' @return Numeric vector with one value per grid point.
#' @export
#' @examples
#' g <- function(rho, r) as.numeric(cr_Amat(r) %*% cr_delta(rho))
#' eval_g_on_grid(g, c(0.5, 0.4, 0.3), seq(0, 0.8, by = 0.2))
eval_g_on_grid <- function(g_fun, rho_vec, grid_r_xu) {
  # Evaluate g function on grid, handling both vectorized and non-vectorized cases
  out <- g_fun(rho_vec, grid_r_xu)

  # If g_fun is not vectorized, evaluate point-by-point
  if (length(out) == 1L && length(grid_r_xu) > 1L) {
    out <- vapply(grid_r_xu, function(r) g_fun(rho_vec, r), numeric(1))
  } else {
    out <- as.numeric(out)
  }

  # Validate output length
  if (length(out) != length(grid_r_xu)) {
    stop("g_fun(rho, grid_r_xu) must return length(grid_r_xu) values.")
  }

  out
}
