
#' # Return final projection indicator
#' @param delta, d row or column vectors/matrices.
#' @param deltastar_demean row or column vectors/matrices.
#' @param Al constraint matrices
#' @param Au constraint matrices
#' @param sigma_l  standard deviations of lambda_hat_l / u
#' @param sigma_u standard deviations of lambda_hat_l / u
#' @param deltaSigma (B x k) matrix
#' @param alphac adjusted critical level
#' @param eta critical area level
#' @param c modifier
#' @param c_bd critical value
#' @param corr_m correlation blocks
#' @param corr_l correlation blocks
#' @param corr_u correlation blocks
#' @param tol_r tolerance for correlation threshold
#' @param g mapping function returning lambda estimates
#' @return Return final projection indicator
#' @export

CIproj_p <- function(c, delta, alphac, c_bd,
                     deltastar_demean,
                     Al, Au,
                     sigma_l, sigma_u,
                     deltaSigma,
                     corr_m, corr_l, corr_u,
                     eta, tol_r, g) {
  # -------------------------------
  # Equivalent to MATLAB CIproj_p.m
  # -------------------------------
  # g : function(delta_matrix) -> lambda_matrix
  # delta, deltastar_demean: row or column vectors/matrices
  # Al, Au: constraint matrices (lower/upper)

  # ----- Lambda at point delta -----
  lambda_l <- g(matrix(delta, nrow = 1))
  lambda_u <- g(matrix(delta, nrow = 1))
  lb <- min(lambda_l)
  ub <- max(lambda_u)
  mb <- (lb + ub) / 2

  # ----- Randomized samples -----
  deltastar <- sweep(deltastar_demean, 2, delta, "+")
  lambdastar_l <- g(deltastar)
  lambdastar_u <- g(deltastar)

  # ----- Condition Δ -----
  delta_sigma <- sqrt(diag(deltaSigma))
  maxdev <- apply(abs(deltastar_demean / matrix(delta_sigma, nrow = nrow(deltastar_demean), ncol = length(delta_sigma), byrow = TRUE)), 1, max)
  ind_Delta <- as.numeric(maxdev <= c_bd)

  # ----- Projection statistics -----
  Tstar_l <- pmax(
    apply((lambdastar_l - lb) / matrix(sigma_l, nrow = nrow(lambdastar_l), ncol = length(sigma_l), byrow = TRUE), 1, min),
    apply((lb - lambdastar_u) / matrix(sigma_u, nrow = nrow(lambdastar_u), ncol = length(sigma_u), byrow = TRUE), 1, min)
  )
  Tstar_m <- pmax(
    apply((lambdastar_l - mb) / matrix(sigma_l, nrow = nrow(lambdastar_l), ncol = length(sigma_l), byrow = TRUE), 1, min),
    apply((mb - lambdastar_u) / matrix(sigma_u, nrow = nrow(lambdastar_u), ncol = length(sigma_u), byrow = TRUE), 1, min)
  )
  Tstar_u <- pmax(
    apply((lambdastar_l - ub) / matrix(sigma_l, nrow = nrow(lambdastar_l), ncol = length(sigma_l), byrow = TRUE), 1, min),
    apply((ub - lambdastar_u) / matrix(sigma_u, nrow = nrow(lambdastar_u), ncol = length(sigma_u), byrow = TRUE), 1, min)
  )

  ind_proj_l <- as.numeric(Tstar_l <= c)
  ind_proj_m <- as.numeric(Tstar_m <= c)
  ind_proj_u <- as.numeric(Tstar_u <= c)

  # ----- Condition cLF -----
  cLF <- qnorm(1 - eta)

  # === Helper function for conditional parts ===
  cond_block <- function(bound, select, lambdastar_l, lambdastar_u, deltastar, sigma_l, sigma_u, label) {
    if (sum(select) == 0) {
      return(rep(0, nrow(deltastar)))
    }
    sub_deltastar <- deltastar[select == 1, , drop = FALSE]
    sub_l_l <- lambdastar_l[select == 1, , drop = FALSE]
    sub_l_u <- lambdastar_u[select == 1, , drop = FALSE]

    Tc1 <- apply((sub_l_l - bound) / matrix(sigma_l, nrow = nrow(sub_l_l), ncol = length(sigma_l), byrow = TRUE), 1, min)
    Tc2 <- apply((bound - sub_l_u) / matrix(sigma_u, nrow = nrow(sub_l_u), ncol = length(sigma_u), byrow = TRUE), 1, min)
    Tc <- pmax(Tc1, Tc2)

    # Call to the truncated-normal bounds function
    th_bounds <- CIcon_TNbounds(bound, sub_deltastar, Al, Au, sigma_l, sigma_u,
                                corr_m, corr_l, corr_u, cLF, cLF, tol_r, g)
    th1 <- th_bounds[[1]]; th2 <- th_bounds[[2]]

    phi <- (pnorm(Tc) - pnorm(th1)) / (pnorm(th2) - pnorm(th1))
    sub_ind <- as.numeric(phi <= 1 - alphac)

    ind_c <- rep(0, nrow(deltastar))
    ind_c[select == 1] <- sub_ind
    ind_c
  }

  # ----- l / m / u conditions -----
  selectl <- (1 - ind_proj_l) * ind_Delta
  ind_c_l <- cond_block(lb, selectl, lambdastar_l, lambdastar_u, deltastar, sigma_l, sigma_u, "l")

  selectm <- (1 - ind_proj_m) * ind_Delta
  ind_c_m <- cond_block(mb, selectm, lambdastar_l, lambdastar_u, deltastar, sigma_l, sigma_u, "m")

  selectu <- (1 - ind_proj_u) * ind_Delta
  ind_c_u <- cond_block(ub, selectu, lambdastar_l, lambdastar_u, deltastar, sigma_l, sigma_u, "u")

  # ----- Combine indicators -----
  ind_l <- 1 - (1 - ind_c_l) * (1 - ind_proj_l)
  ind_m <- 1 - (1 - ind_c_m) * (1 - ind_proj_m)
  ind_u <- 1 - (1 - ind_c_u) * (1 - ind_proj_u)

  # p1 <- mean((1 - ind_l * ind_m) * ind_Delta)
  # p2 <- mean((1 - ind_u * ind_m) * ind_Delta)
  # p <- max(p1, p2) + eta
  p1  <- mean((1 - ind_l * ind_m) * ind_Delta, na.rm = TRUE)
  p2  <- mean((1 - ind_u * ind_m) * ind_Delta, na.rm = TRUE)
  p   <- max(p1, p2, na.rm = TRUE) + eta

  if (!is.finite(p)) p <- NA_real_
  return(p)

  return(p)
}
