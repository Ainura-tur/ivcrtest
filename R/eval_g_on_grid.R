
#' Evaluates g function on a grid
#' @param g_fun Function mapping rho and r_xu to lambda.
#' @param rho_vec Named numeric vector of correlations.
#' @param grid_r_xu Numeric vector of grid points.
#' @return Numeric vector of g values.
#' @export
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
