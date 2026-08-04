
#' Estimates rho vector from data
#' @param data Data frame.
#' @param x_col Name of x column.
#' @param y_col Name of y column.
#' @param z_col Name of z column.
#' @return Named numeric vector of correlations.
#' @export

estimate_rho <- function(data, x_col = "x", y_col = "y", z_col = "z") {
  # Validate column presence
  if (!all(c(x_col, y_col, z_col) %in% names(data))) {
    stop("Data must contain columns: ", paste(c(x_col, y_col, z_col), collapse = ", "))
  }

  # Use complete cases only
  d <- data[complete.cases(data[, c(x_col, y_col, z_col)]), c(x_col, y_col, z_col)]
  if (nrow(d) < 5) {
    stop("Not enough complete observations to estimate correlations.")
  }

  c(
    rho_xy = suppressWarnings(stats::cor(d[[x_col]], d[[y_col]])),
    rho_xz = suppressWarnings(stats::cor(d[[x_col]], d[[z_col]])),
    rho_yz = suppressWarnings(stats::cor(d[[y_col]], d[[z_col]]))
  )
}
