# Namespace protection (avoid MASS/stats masking dplyr verbs)
select    <- dplyr::select;      filter    <- dplyr::filter
mutate    <- dplyr::mutate;      slice     <- dplyr::slice
recode    <- dplyr::recode;      rename    <- dplyr::rename
summarise <- dplyr::summarise;   summarize <- dplyr::summarize
arrange   <- dplyr::arrange;     count     <- dplyr::count;   lag <- dplyr::lag

#' Modified conditional union-bound confidence interval
#'
#' Returns the modified conditional interval of Bei (2024), obtained as the
#' union of the calibrated projection interval and the conditional interval,
#' corresponding to the critical value \eqn{\hat c_m(\theta;\alpha) =
#' \max\{\hat c_c(\theta;\alpha_c), \hat c_t\}}.
#'
#' The routine is parameterisation-agnostic. Under the reparameterised route of
#' Remark A.1 it is called with \code{deltahat} the two-dimensional
#' \eqn{\hat\delta}, \code{deltaSigma} the 2x2 \eqn{\hat\Omega_\delta}, and
#' \code{Al = Au = A} the known K x 2 grid matrix, with \code{g} the linear map
#' \eqn{\delta \mapsto A\delta}. Under the primitive route it is called with the
#' three-dimensional \eqn{\hat\rho}, the 3x3 \eqn{\hat\Sigma_\rho}, the estimated
#' gradient array, and the nonlinear \code{g}.
#'
#' @param deltahat Numeric vector, the point estimate being perturbed.
#' @param deltaSigma Sampling covariance matrix of \code{deltahat}.
#' @param Al,Au Constraint matrices; equal in the degenerate-bound case.
#' @param alpha Nominal level.
#' @param alphac Conditional level; admissible when in (alpha/2, alpha).
#' @param eta Lower-truncation confidence level for the auxiliary set.
#' @param B,Blarge Numbers of Gaussian draws.
#' @param tol Bisection tolerance.
#' @param tol_r Tolerance for the pairwise-correlation screen in the truncation
#'   bounds. See \code{\link{cr_corr_screen}} for a diagnostic of how many grid
#'   pairs this removes at a given grid size.
#' @param index Optional indices of \code{deltahat} to hold at zero; may be NULL.
#' @param g Mapping returning the grid-level lambda values.
#' @param seed Integer seed. The RNG state of the calling session is saved on
#'   entry and restored on exit.
#' @return A list with the hybrid, conditional and projection intervals, the
#'   least-favourable point reached, and the calibrated critical value.
#' @export
#' @importFrom stats qnorm quantile optim
#' @importFrom MASS mvrnorm

CIhybrid <- function(deltahat, deltaSigma, Al, Au,
                     alpha, alphac, eta,
                     B, Blarge,
                     tol, tol_r, index, g,
                     seed = 0L) {

  # ----- 0. Save and restore the caller's RNG stream -----
  has_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (has_seed) get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        rm(".Random.seed", envir = globalenv())
      }
    } else {
      assign(".Random.seed", old_seed, envir = globalenv())
    }
  }, add = TRUE)

  set.seed(seed)

  # ----- 1. Random draws and basic quantities -----
  k_dim <- length(deltahat)
  deltastar_demean_large <- MASS::mvrnorm(Blarge, mu = rep(0, k_dim), Sigma = deltaSigma)
  deviation <- apply(abs(deltastar_demean_large / sqrt(diag(deltaSigma))), 1, max)
  c_bd <- quantile(deviation, 1 - eta)

  sigma_delta <- sqrt(diag(deltaSigma))

  # ----- 2. Construct lambdaSigma and related components -----
  Lambda <- rbind(Al, Au)
  lambdaSigma <- Lambda %*% deltaSigma %*% t(Lambda)
  lambdasigma <- sqrt(diag(lambdaSigma))

  kk <- nrow(Al)
  sigma_l <- lambdasigma[1:kk]
  sigma_u <- lambdasigma[(kk + 1):(2 * kk)]

  corr_all <- diag(1 / lambdasigma) %*% lambdaSigma %*% diag(1 / lambdasigma)
  corr_l <- corr_all[1:kk, 1:kk]
  corr_u <- corr_all[(kk + 1):(2 * kk), (kk + 1):(2 * kk)]
  corr_m <- corr_all[1:kk, (kk + 1):(2 * kk)]

  # ----- 3. Feasible bounds -----
  lb <- deltahat - sigma_delta * c_bd
  ub <- deltahat + sigma_delta * c_bd
  delta1 <- deltahat

  if (!is.null(index) && length(index) > 0L) {
    lb[index] <- 0
    ub[index] <- 0
    delta1[index] <- 0
  }

  cl <- 0
  cu <- qnorm(1 - alpha / 2)
  c <- (cl + cu) / 2

  delta_fea <- list()
  c_fea <- numeric()

  # ----- 4. Objective wrapper -----
  obj_large <- function(delta, c_check) {
    (alpha - CIproj_p(c_check, delta, alphac, c_bd,
                      deltastar_demean_large,
                      Al, Au, sigma_l, sigma_u,
                      deltaSigma, corr_m, corr_l, corr_u,
                      eta, tol_r, g)) * 100
  }

  # ----- 5. Iterative bisection loop -----
  K1 <- 1
  while ((cu - cl) > tol) {
    set.seed(seed + K1)
    K1 <- K1 + 1

    deltastar_demean <- MASS::mvrnorm(B, mu = rep(0, k_dim), Sigma = deltaSigma)

    obj <- function(delta) {
      (alpha - CIproj_p(c, delta, alphac, c_bd, deltastar_demean,
                        Al, Au, sigma_l, sigma_u, deltaSigma,
                        corr_m, corr_l, corr_u, eta, tol_r, g)) * 100
    }

    p1 <- obj_large(delta1, c)

    if (!is.finite(p1)) {
      warning("CIhybrid: non-finite objective at delta1 = ",
              paste(round(delta1, 3), collapse = ", "),
              ", c = ", round(c, 3), call. = FALSE)
    }

    if (p1 >= 0) {
      f_optim <- tryCatch({
        optim(par = delta1, fn = obj, method = "L-BFGS-B",
              lower = lb, upper = ub, control = list(maxit = 1000))
      }, error = function(e) NULL)

      if (!is.null(f_optim)) {
        delta1 <- f_optim$par
      }
      p1 <- obj_large(delta1, c)

      if (p1 >= 0) {
        cu <- c
      } else {
        cl <- c
        delta_fea[[length(delta_fea) + 1]] <- delta1
        c_fea <- c(c, c_fea)
      }
    } else {
      cl <- c
      delta_fea[[length(delta_fea) + 1]] <- delta1
      c_fea <- c(c, c_fea)
    }

    c <- (cl + cu) / 2
  }

  # ----- 6. Final refinement -----
  cl <- c
  cu <- qnorm(1 - alpha / 2)

  while ((cu - cl) > tol) {
    c <- (cl + cu) / 2
    p <- obj_large(delta1, c)
    if (p >= 0) {
      cu <- c
    } else {
      cl <- c
    }
  }

  # ----- 7. Compute confidence intervals -----
  lambdahat_l <- g(deltahat)
  lambdahat_u <- lambdahat_l

  c_LF <- qnorm(1 - eta)
  CI_p <- c(min(lambdahat_l - c * sigma_l), max(lambdahat_u + c * sigma_u))
  CI_c <- CIcon(deltahat, deltaSigma, Al, Au, c_LF, alphac, tol, tol_r, g)

  CI_h <- c(min(CI_p[1], CI_c[1]), max(CI_p[2], CI_c[2]))

  list(CI_h = CI_h,
       CI_c = CI_c,
       CI_p = CI_p,
       delta1 = delta1,
       c = c)
}
