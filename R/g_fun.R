# Namespace protection (avoid MASS/stats masking dplyr verbs)
select    <- dplyr::select;      filter    <- dplyr::filter
mutate    <- dplyr::mutate;      slice     <- dplyr::slice
recode    <- dplyr::recode;      rename    <- dplyr::rename
summarise <- dplyr::summarise;   summarize <- dplyr::summarize
arrange   <- dplyr::arrange;     count     <- dplyr::count;   lag <- dplyr::lag

#'Return value of r_zu
#' @param r r_xu values
#' @param rho_xy observed correlation values
#' @param rho_xz observed correlation values
#' @param rho_yz observed correlation values
#' @return Return value of r_zu
#' @export

g_fun <- function(r, rho_xy, rho_xz, rho_yz) {
  rho_xz * r - (rho_xy * rho_xz - rho_yz) * sqrt((1 - r^2) / (1 - rho_xy^2))
}
