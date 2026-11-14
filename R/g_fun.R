
#' Grid map from the assumed endogeneity correlation to rho_zu
#'
#' Evaluates \eqn{g(\rho; r) = \rho_{xz} r - (\rho_{xy}\rho_{xz} - \rho_{yz})
#' \sqrt{(1 - r^2)/(1 - \rho_{xy}^2)}}, the value of \eqn{\rho_{zu}} implied
#' by the observed correlations when \eqn{\rho_{xu} = r}. The identified set
#' for \eqn{\rho_{zu}} is the range of this map over the domain D.
#' \code{\link{g_xu_safe}} is a vectorised version that clamps its inputs.
#'
#' @param r Numeric vector of assumed values of \eqn{\rho_{xu}}.
#' @param rho_xy,rho_xz,rho_yz Observed correlations.
#' @return Numeric vector of implied \eqn{\rho_{zu}} values, one per element
#'   of \code{r}.
#' @export
#' @examples
#' r <- seq(0, 0.8, by = 0.2)
#' g_fun(r, rho_xy = 0.5, rho_xz = 0.4, rho_yz = 0.3)

g_fun <- function(r, rho_xy, rho_xz, rho_yz) {
  rho_xz * r - (rho_xy * rho_xz - rho_yz) * sqrt((1 - r^2) / (1 - rho_xy^2))
}
