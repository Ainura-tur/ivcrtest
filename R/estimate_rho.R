
#' Sample correlations of the regressor, outcome and instrument
#'
#' Computes \eqn{\hat\rho = (\hat\rho_{xy}, \hat\rho_{xz}, \hat\rho_{yz})'} on
#' the rows with no missing value in the three columns. At least five such
#' rows are required.
#'
#' @param data Data frame.
#' @param x_col,y_col,z_col Names of the regressor, outcome and instrument
#'   columns.
#' @return Named numeric vector \code{c(rho_xy, rho_xz, rho_yz)}.
#' @export
#' @examples
#' set.seed(1)
#' x <- rnorm(100)
#' d <- data.frame(x = x, y = x + rnorm(100), z = 0.5 * x + rnorm(100))
#' estimate_rho(d)

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
