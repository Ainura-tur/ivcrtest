# Namespace protection (avoid MASS/stats masking dplyr verbs)
select    <- dplyr::select;      filter    <- dplyr::filter
mutate    <- dplyr::mutate;      slice     <- dplyr::slice
recode    <- dplyr::recode;      rename    <- dplyr::rename
summarise <- dplyr::summarise;   summarize <- dplyr::summarize
arrange   <- dplyr::arrange;     count     <- dplyr::count;   lag <- dplyr::lag

#' Computes list(th_1, th_2): lower / upper bounds vectors
#' @param theta scalar
#' @param deltahat (B x k) matrix
#' @param Al constraint matrice
#' @param Au constraint matrices
#' @param sigma_l standard deviations of lambda_hat_l / u
#' @param sigma_u standard deviations of lambda_hat_l / u
#' @param corr_m correlation blocks
#' @param corr_l correlation blocks
#' @param corr_u correlation blocks
#' @param cLFl critical values
#' @param cLFu critical values
#' @param tol_r tolerance for correlation threshold
#' @param g mapping function returning lambda estimates
#' @return list(th_1, th_2): lower / upper bounds vectors
#' @export
CIcon_TNbounds <- function(theta, deltahat,
                           Al, Au,
                           sigma_l, sigma_u,
                           corr_m, corr_l, corr_u,
                           cLFl, cLFu, tol_r, g) {
  B   <- nrow(deltahat)
  k_l <- nrow(Al)
  k_u <- nrow(Au)

  lambdahat_l <- g(deltahat)
  lambdahat_u <- lambdahat_l

  Tlb <- sweep(lambdahat_l - theta, 2, sigma_l, "/")   # B x k_l
  Tub <- sweep(theta - lambdahat_u, 2, sigma_u, "/")   # B x k_u

  Tl <- apply(Tlb, 1, min)         # length B
  bl <- apply(Tlb, 1, which.min)   # argmin per row (1..k_l)
  Tu <- apply(Tub, 1, min)
  bu <- apply(Tub, 1, which.min)   # argmin per row (1..k_u)

  # corr_m rows/cols at those argmins
  corr_m_blb <- corr_m[bl, , drop = FALSE]           # B x k_u
  corr_m_bub <- t(corr_m[, bu, drop = FALSE])        # B x k_l

  ## ----- Case Tl >= Tu : “lower” side -----
  # TS part
  tTS_l  <- matrix(1e10, nrow = B, ncol = k_u)
  # Broadcast Tl and Tub row-wise
  tTS_ll <- (1 + corr_m_blb)^(-1) * (Tub + corr_m_blb * Tl)
  maskTS_l <- (1 + corr_m_blb) > tol_r
  tTS_l[maskTS_l] <- tTS_ll[maskTS_l]

  th_1_1 <- apply(tTS_l, 1, min) * (Tl >= Tu)
  th_1_1[th_1_1 == 1e10] <- -1e10
  th_1_1 <- th_1_1 * (Tl >= Tu)

  # LF part
  tLFl <- rep(cLFl, B)

  # B part: polyhedral upper truncation for the selected minimum.
  # Constraint T_{l,j} - T_l >= 0 decomposed along the direction of T_l gives
  #   taux_1[i,j] = ( Tlb[i,j] - corr_l[bl[i],j] * Tl[i] ) / ( 1 - corr_l[bl[i],j] ).
  # The numerator must use the j-th statistic Tlb[i,j], not the row minimum.
  # Using the minimum collapses the expression to Tl[i] for every j.
  CL_bl  <- corr_l[bl, , drop = FALSE]                                  # B x k_l
  Tl_rep <- matrix(Tl, nrow = B, ncol = k_l)
  taux_1 <- (1 - CL_bl)^(-1) * (Tlb - CL_bl * Tl_rep)                   # B x k_l

  tB_l2   <- matrix(1e10, nrow = B, ncol = k_l)
  maskB_l <- (1 > CL_bl + tol_r)
  tB_l2[maskB_l] <- taux_1[maskB_l]

  th_2_1 <- apply(cbind(tLFl, tB_l2), 1, min) * (Tl >= Tu)

  ## ----- Case Tl < Tu : “upper” side -----
  # TS part
  tTS_u  <- matrix(1e10, nrow = B, ncol = k_l)
  tTS_u1 <- (1 + corr_m_bub)^(-1) * (Tlb + corr_m_bub * Tu)
  maskTS_u <- (1 + corr_m_bub) > tol_r
  tTS_u[maskTS_u] <- tTS_u1[maskTS_u]

  th_1_2 <- apply(tTS_u, 1, min) * (Tl < Tu)
  th_1_2[th_1_2 == 1e10] <- -1e10
  th_1_2 <- th_1_2 * (Tl < Tu)

  # LF part
  tLFu <- rep(cLFu, B)

  # B part: same construction on the upper side.
  CU_bu  <- corr_u[bu, , drop = FALSE]                                  # B x k_u
  Tu_rep <- matrix(Tu, nrow = B, ncol = k_u)
  taux_2 <- (1 - CU_bu)^(-1) * (Tub - CU_bu * Tu_rep)                   # B x k_u

  tB_u2   <- matrix(1e10, nrow = B, ncol = k_u)
  maskB_u <- (1 > CU_bu + tol_r)
  tB_u2[maskB_u] <- taux_2[maskB_u]

  th_2_2 <- apply(cbind(tLFu, tB_u2), 1, min) * (Tl < Tu)

  th_1 <- th_1_1 + th_1_2
  th_2 <- th_2_1 + th_2_2
  list(th_1 = th_1, th_2 = th_2)
}
