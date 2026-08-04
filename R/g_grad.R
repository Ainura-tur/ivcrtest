# Namespace protection (avoid MASS/stats masking dplyr verbs)
select    <- dplyr::select;      filter    <- dplyr::filter
mutate    <- dplyr::mutate;      slice     <- dplyr::slice
recode    <- dplyr::recode;      rename    <- dplyr::rename
summarise <- dplyr::summarise;   summarize <- dplyr::summarize
arrange   <- dplyr::arrange;     count     <- dplyr::count;   lag <- dplyr::lag

#' Return value of gradient of g function
#' @param r r_xu values
#' @param rho_xy observed correlation values
#' @param rho_xz observed correlation values
#' @param rho_yz observed correlation values

#' @return Return value of gradient of g function
#' @export
#'

g_grad <- function(r, rho_xy, rho_xz, rho_yz) {
  S <- sqrt((1 - r^2) / (1 - rho_xy^2))
  c(
    -rho_xz * S - (rho_xy * rho_xz - rho_yz) * S * rho_xy / (1 - rho_xy^2), # ∂/∂rho_xy
    r - rho_xy * S,                                                         # ∂/∂rho_xz
    S                                                                       # ∂/∂rho_yz
  )
}
