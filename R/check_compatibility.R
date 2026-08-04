# Namespace protection (avoid MASS/stats masking dplyr verbs)
select    <- dplyr::select;      filter    <- dplyr::filter
mutate    <- dplyr::mutate;      slice     <- dplyr::slice
recode    <- dplyr::recode;      rename    <- dplyr::rename
summarise <- dplyr::summarise;   summarize <- dplyr::summarize
arrange   <- dplyr::arrange;     count     <- dplyr::count;   lag <- dplyr::lag

#' Compute the CR test intervals for one instrument
#'
#' Constructs the MCUB confidence interval for the scalar partially identified
#' parameter \eqn{\theta_0 = \rho_{zu}}, together with the Imbens-Manski/Stoye
#' diagnostic interval and the naive selected-endpoint benchmark, on a fixed
#' finite grid over the admissible domain D.
#'
#' Two routes are available. \code{route = "reparam"} is the default and
#' implements Remark A.1: the Gaussian simulation is carried out entirely in the
#' exact two-dimensional reparameterised experiment, with the known K x 2 matrix
#' \eqn{A} and \eqn{\hat\Omega_\delta = J(\hat\rho)\hat\Sigma_\rho J(\hat\rho)'}.
#' \code{route = "primitive"} draws in the three-dimensional space of
#' \eqn{\hat\rho} and evaluates the nonlinear map, with the estimated gradient
#' array in place of \eqn{A}. The two schemes differ by \eqn{o_p(1)} and give
#' numerically identical grid-level standard errors by equation (12), but only
#' the first matches the theory exactly, since \eqn{A} involves no estimated
#' quantity.
#'
#' All covariance quantities are on the estimator scale, \eqn{\hat\Sigma_{\rho,n}
#' = \widehat{Var}(\hat\rho)}, so no further division by n is applied anywhere.
#'
#' @param df Data frame containing columns \code{x}, \code{y}, \code{z}.
#' @param i Integer label for the instrument, used in messages and seeds.
#' @param n Ignored; determined internally from \code{df}.
#' @param k Sign of the maintained endogeneity restriction, +1 or -1. Used only
#'   to validate \code{rxu_range}.
#' @param alpha Numeric significance level.
#' @param rxu_range Numeric vector of length 2 giving the admissible domain D.
#' @param K Integer number of equally spaced grid points. Default 80, matching
#'   Remark 3.3 and Table 6.
#' @param route Either \code{"reparam"} or \code{"primitive"}.
#' @param B_boot Integer bootstrap replications for the covariance of rho_hat.
#' @param B,Blarge Numbers of Gaussian draws in the MCUB calculation.
#' @param eta Lower-truncation confidence level.
#' @param tol,tol_r Tolerances passed to the MCUB routines.
#' @param B_im Integer bootstrap replications for the IM/Stoye endpoint standard
#'   errors.
#' @param seed Integer seed.
#' @param compute_pvalue Logical; whether to invert the interval over the
#'   admissible range of alpha. This is the expensive step.
#' @param verbose Logical; print progress.
#' @return A one-row data frame of results, with diagnostics attached as the
#'   \code{"diagnostics"} attribute.
#' @export
#' @importFrom stats complete.cases cor qnorm

check_compatibility <- function(df, i = 1L, n = NULL, k = 1,
                                alpha          = 0.05,
                                rxu_range      = c(0, 0.8),
                                K              = 80L,
                                route          = c("reparam", "primitive"),
                                B_boot         = 800L,
                                B              = 500L,
                                Blarge         = 5000L,
                                eta            = 0.001,
                                tol            = 1e-3,
                                tol_r          = 1e-3,
                                B_im           = 999L,
                                seed           = 1L,
                                compute_pvalue = TRUE,
                                verbose        = TRUE) {

  route <- match.arg(route)

  if (!all(c("x", "y", "z") %in% names(df))) {
    stop("check_compatibility: df must contain columns 'x', 'y', 'z'.")
  }
  df <- df[complete.cases(df[, c("x", "y", "z")]), , drop = FALSE]
  n  <- nrow(df)
  if (n < 30L) {
    warning("check_compatibility: small sample size (", n, " observations).", call. = FALSE)
    return(NULL)
  }

  if (length(rxu_range) != 2L || !all(is.finite(rxu_range))) {
    stop("check_compatibility: rxu_range must be a finite numeric vector of length 2.")
  }
  rxu_range <- sort(rxu_range)
  if (any(abs(rxu_range) >= 1)) {
    stop("check_compatibility: rxu_range must lie strictly inside (-1, 1); ",
         "Assumption 2 rules out perfect correlation.")
  }
  if (k > 0 && rxu_range[1] < 0) {
    stop("check_compatibility: k = ", k, " asserts corr(x, u) > 0, so D must lie in [0, 1).")
  }
  if (k < 0 && rxu_range[2] > 0) {
    stop("check_compatibility: k = ", k, " asserts corr(x, u) < 0, so D must lie in (-1, 0].")
  }
  if (K < 2L) stop("check_compatibility: K must be at least 2.")
  if (!(alpha > 0 && alpha < 0.5)) {
    stop("check_compatibility: Bei (2024) requires alpha in (0, 1/2).")
  }
  alphac <- 0.8 * alpha
  if (!(eta < (alpha - alphac) / 2)) {
    stop("check_compatibility: tuning parameters violate Proposition 3.2; ",
         "eta must be below (alpha - alphac)/2 = ", (alpha - alphac) / 2, ".")
  }

  say <- function(...) if (isTRUE(verbose)) cat(sprintf(...))

  # ---- sample correlations and one shared bootstrap covariance ----
  rho_hat <- c(rho_xy = cor(df$x, df$y),
               rho_xz = cor(df$x, df$z),
               rho_yz = cor(df$y, df$z))
  Sigma_rho <- estimate_cov_corr_boot(df$x, df$y, df$z, B = B_boot, seed = 1000L + i)

  r_grid <- seq(rxu_range[1], rxu_range[2], length.out = K)

  # ---- diagnostics ----
  sign_diag  <- cr_sign_condition(rho_hat, rxu_range)
  gap_diag   <- cr_endpoint_gap(rho_hat, Sigma_rho, rxu_range, K = K)
  disc_diag  <- cr_discretization(rho_hat, rxu_range, K = K)
  screen_diag <- cr_corr_screen(rho_hat, Sigma_rho, rxu_range, K = K, tol_r = tol_r)

  if (!sign_diag$holds) {
    warning("check_compatibility: the sign condition of Lemma 2.3 fails for instrument ", i,
            " (h_r = ", signif(sign_diag$h_r, 3), "). C_K need not equal C; ",
            "see the discretization diagnostics.", call. = FALSE)
  }
  if (screen_diag$share_screened > 0.05) {
    warning("check_compatibility: tol_r = ", tol_r, " screens ",
            round(100 * screen_diag$share_screened, 1),
            "% of grid pairs out of the truncated-normal bounds at K = ", K,
            ". Consider a smaller tol_r.", call. = FALSE)
  }

  # ---- MCUB inputs, by route ----
  Amat  <- cr_Amat(r_grid)
  Omega <- cr_omega_delta(rho_hat, Sigma_rho)

  if (route == "reparam") {
    dhat  <- cr_delta(rho_hat)
    dSig  <- Omega
    Al    <- Amat
    g_map <- cr_grid_map(Amat)
  } else {
    dhat  <- as.numeric(rho_hat)
    dSig  <- Sigma_rho
    Al    <- t(vapply(r_grid, function(r) {
      as.numeric(local_compute_gradient_safe(r, rho_hat[1], rho_hat[2], rho_hat[3]))
    }, numeric(3)))
    g_map <- function(delta) g_xu_safe(r_grid, delta)
  }
  Au <- Al

  # Degenerate lower and upper bounds: lambda_{l,b} = lambda_{u,b} = lambda_b.
  lambda_hat <- as.numeric(Amat %*% cr_delta(rho_hat))
  sigma_b    <- sqrt(pmax(diag(Amat %*% Omega %*% t(Amat)), 0))
  plug_in    <- c(min(lambda_hat), max(lambda_hat))

  status <- character(0)

  # ---- naive selected-endpoint benchmark, equation (22) ----
  CI_naive <- c(NA_real_, NA_real_)
  res_naive <- tryCatch({
    say("  Instrument %d: naive endpoint benchmark...\n", i)
    ci_naive_endpoint(df$x, df$y, df$z, rxu_range = rxu_range, K = K,
                      alpha = alpha, Sigma_rho = Sigma_rho)
  }, error = function(e) {
    status <<- c(status, paste0("naive:", conditionMessage(e)))
    NULL
  })
  if (!is.null(res_naive)) CI_naive <- res_naive$CI

  # ---- MCUB interval ----
  CI_mcub <- c(NA_real_, NA_real_)
  res_mcub <- tryCatch({
    say("  Instrument %d: MCUB interval (%s route)...\n", i, route)
    CIhybrid(dhat, dSig, Al, Au,
             alpha = alpha, alphac = alphac, eta = eta,
             B = B, Blarge = Blarge, tol = tol, tol_r = tol_r,
             index = NULL, g = g_map, seed = seed)
  }, error = function(e) {
    status <<- c(status, paste0("mcub:", conditionMessage(e)))
    NULL
  })
  if (!is.null(res_mcub)) CI_mcub <- res_mcub$CI_h

  # ---- inverted-CI membership result ----
  p_adm <- NA_real_
  p_lab <- NA_character_
  if (isTRUE(compute_pvalue) && !is.null(res_mcub)) {
    res_p <- tryCatch({
      say("  Instrument %d: inverting over the admissible range...\n", i)
      pvalue_mcub_zero(dhat, dSig, Al, Au, g_map,
                       eta = eta, B = max(1L, B %/% 2L), Blarge = max(1L, Blarge %/% 2L),
                       tol = tol, tol_r = tol_r, seed = seed)
    }, error = function(e) {
      status <<- c(status, paste0("pvalue:", conditionMessage(e)))
      NULL
    })
    if (!is.null(res_p)) {
      p_adm <- res_p$p_adm
      p_lab <- res_p$label
    }
  }

  # ---- Imbens-Manski/Stoye diagnostic interval ----
  CI_IM <- c(NA_real_, NA_real_)
  reject_IM <- NA
  res_IM <- tryCatch({
    say("  Instrument %d: IM/Stoye diagnostic...\n", i)
    im_stoye_inference(data = df, grid_r_xu = r_grid,
                       g_fun = function(rho_vec, r_xu) {
                         as.numeric(cr_Amat(r_xu) %*% cr_delta(rho_vec))
                       },
                       x_col = "x", y_col = "y", z_col = "z",
                       se_method = "bootstrap", B = B_im, alpha = alpha,
                       seed = 1000L + i)
  }, error = function(e) {
    status <<- c(status, paste0("im:", conditionMessage(e)))
    NULL
  })
  if (!is.null(res_IM)) {
    CI_IM     <- as.numeric(res_IM$CI_IM)
    reject_IM <- res_IM$reject_H0_0_in_C
  }

  zero_in <- function(ci) {
    if (any(!is.finite(ci))) return(NA)
    ci[1] <= 0 && 0 <= ci[2]
  }
  mark <- function(v) if (is.na(v)) NA_character_ else if (isTRUE(v)) "\u2713" else "\u00d7"

  fmt <- function(ci) {
    if (any(!is.finite(ci))) NA_character_ else sprintf("[%.3f, %.3f]", ci[1], ci[2])
  }

  out <- data.frame(
    Z               = paste0("Z", i),
    n               = n,
    route           = route,
    K               = K,
    D               = sprintf("[%.2f, %.2f]", rxu_range[1], rxu_range[2]),
    plug_in         = sprintf("[%.3f, %.3f]", plug_in[1], plug_in[2]),
    CI_MCUB         = fmt(CI_mcub),
    CI_IM           = fmt(CI_IM),
    CI_naive        = fmt(CI_naive),
    Zero_in_CI_MCUB = mark(zero_in(CI_mcub)),
    Zero_in_CI_IM   = mark(if (is.na(reject_IM)) NA else !reject_IM),
    Zero_in_naive   = mark(zero_in(CI_naive)),
    p_adm           = p_adm,
    p_adm_label     = p_lab,
    G_hat           = gap_diag$G,
    sign_condition  = sign_diag$holds,
    status          = if (length(status) == 0L) "ok" else paste(status, collapse = "; "),
    stringsAsFactors = FALSE
  )

  attr(out, "diagnostics") <- list(
    rho_hat        = rho_hat,
    Sigma_rho      = Sigma_rho,
    delta_hat      = cr_delta(rho_hat),
    Omega_delta    = Omega,
    eig_Omega      = eigen(Omega, symmetric = TRUE, only.values = TRUE)$values,
    r_grid         = r_grid,
    lambda_hat     = lambda_hat,
    sigma_b        = sigma_b,
    sign_condition = sign_diag,
    endpoint_gap   = gap_diag,
    discretization = disc_diag,
    corr_screen    = screen_diag,
    mcub           = res_mcub
  )
  out
}
