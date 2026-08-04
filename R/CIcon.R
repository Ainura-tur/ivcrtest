# Namespace protection (avoid MASS/stats masking dplyr verbs)
select    <- dplyr::select;      filter    <- dplyr::filter
mutate    <- dplyr::mutate;      slice     <- dplyr::slice
recode    <- dplyr::recode;      rename    <- dplyr::rename
summarise <- dplyr::summarise;   summarize <- dplyr::summarize
arrange   <- dplyr::arrange;     count     <- dplyr::count;   lag <- dplyr::lag

#' Returns final CI
#' @param deltahat (B x k) matrix
#' @param Al constraint matrices
#' @param Au constraint matrices
#' @param deltaSigma (B x k) matrix
#' @param alphac adjusted alpha level
#' @param c_LF critical value
#' @param tol tolerance for correlation threshold
#' @param tol_r tolerance for correlation threshold
#' @param g mapping function returning lambda estimates
#' @return Return final CI
#' @export

CIcon <- function(deltahat, deltaSigma, Al, Au,
                  c_LF, alphac, tol, tol_r, g) {
  # -------------------------------
  # Equivalent to MATLAB CIcon.m
  # -------------------------------

  kk <- nrow(Al)

  # ----- 1. Lambda covariance and correlations -----
  Lambda <- rbind(Al, Au)
  lambdaSigma <- Lambda %*% deltaSigma %*% t(Lambda)
  lambda_sigma <- sqrt(diag(lambdaSigma))

  corr_all <- diag(1 / lambda_sigma) %*% lambdaSigma %*% diag(1 / lambda_sigma)
  corr_m <- corr_all[1:kk, (kk + 1):(2 * kk)]
  corr_l <- corr_all[1:kk, 1:kk]
  corr_u <- corr_all[(kk + 1):(2 * kk), (kk + 1):(2 * kk)]

  sigma_l <- lambda_sigma[1:kk]
  sigma_u <- lambda_sigma[(kk + 1):(2 * kk)]

  # ----- 2. Compute lambda bounds -----
  gval <- as.numeric(g(matrix(deltahat, nrow = 1)))

  lb <- gval - sigma_l * c_LF
  ub <- gval + sigma_u * c_LF
  lb <- min(lb)
  ub <- max(ub)

  mid <- (lb + ub) / 2
  lb_pt <- min(gval)
  ub_pt <- max(gval)

  # ----- 3. CI lower bound -----
  rej <- 1
  theta <- lb

  while (rej == 1 && theta <= min(mid, lb_pt)) {
    Tcl <- min((gval - theta) / sigma_l)
    Tcu <- min((theta - gval) / sigma_u)
    Tc <- max(Tcl, Tcu)

    th_bounds <- CIcon_TNbounds(theta, matrix(deltahat, nrow = 1),
                                Al, Au, sigma_l, sigma_u,
                                corr_m, corr_l, corr_u,
                                c_LF, c_LF, tol_r, g)
    th_1 <- th_bounds[[1]]
    th_2 <- th_bounds[[2]]

    t <- qnorm((1 - alphac) * pnorm(th_2) + alphac * pnorm(th_1))
    rej <- as.numeric(Tc > t)
    theta <- theta + tol
  }
  CI_lower <- theta

  # ----- 4. CI upper bound -----
  rej <- 1
  theta <- ub

  while (rej == 1 && theta >= max(mid, ub_pt)) {
    Tcl <- min((gval - theta) / sigma_l)
    Tcu <- min((theta - gval) / sigma_u)
    Tc <- max(Tcl, Tcu)

    th_bounds <- CIcon_TNbounds(theta, matrix(deltahat, nrow = 1),
                                Al, Au, sigma_l, sigma_u,
                                corr_m, corr_l, corr_u,
                                c_LF, c_LF, tol_r, g)
    th_1 <- th_bounds[[1]]
    th_2 <- th_bounds[[2]]

    t <- qnorm((1 - alphac) * pnorm(th_2) + alphac * pnorm(th_1))
    rej <- as.numeric(Tc > t)
    theta <- theta - tol
  }
  CI_upper <- theta

  # ----- 5. Return final CI -----
  CIcon <- c(CI_lower, CI_upper)
  return(CIcon)
}
