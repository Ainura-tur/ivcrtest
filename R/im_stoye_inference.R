
#' Imbens-Manski/Stoye interval for the identified set
#'
#' Computes the Imbens and Manski (2004) and Stoye (2009) confidence interval
#' for a parameter whose identified set is \eqn{[\theta_l, \theta_u]}, where
#' \eqn{\theta_l} and \eqn{\theta_u} are the minimum and maximum of
#' \code{g_fun(rho_hat, r)} over the grid. Endpoint standard errors come from
#' the bootstrap or the delta method. The critical value solves
#' \eqn{\Phi(c + \hat\Delta/\hat s) - \Phi(-c) = 1 - \alpha}, where
#' \eqn{\hat\Delta} is the width of the estimated set and \eqn{\hat s} the
#' larger endpoint standard error; see \code{\link{solve_c_im}}. The function
#' also reports whether zero lies outside the interval.
#' \code{\link{check_compatibility}} reports this interval as a diagnostic
#' alongside the MCUB interval.
#'
#' @param data Data frame containing the three variables.
#' @param grid_r_xu Numeric vector of grid points for \eqn{\rho_{xu}}.
#' @param g_fun Function \code{g_fun(rho, r)} returning the grid map at the
#'   correlation vector \code{rho}, ordered (rho_xy, rho_xz, rho_yz), and the
#'   grid point or points \code{r}.
#' @param alpha Numeric significance level.
#' @param se_method Either \code{"bootstrap"} or \code{"delta"}.
#' @param B Integer number of bootstrap replications.
#' @param grad_fun Used when \code{se_method = "delta"}: a function
#'   \code{grad_fun(rho, r)} returning the gradient of the grid map with
#'   respect to \code{rho}.
#' @param cov_rho Used when \code{se_method = "delta"}: an optional 3x3
#'   covariance of \eqn{\hat\rho} on the estimator scale. It is bootstrapped
#'   when NULL.
#' @param seed Integer seed for the bootstrap draws, or NULL to use the
#'   current random-number stream. The caller's random-number state is
#'   restored on exit.
#' @param x_col,y_col,z_col Names of the regressor, outcome and instrument
#'   columns in \code{data}.
#' @param return_boot Logical; whether to return the bootstrap endpoint draws.
#' @return A list with the interval \code{CI_IM}; the logical
#'   \code{reject_H0_0_in_C}, TRUE when zero lies outside the interval; the
#'   plug-in endpoints \code{theta_l_hat} and \code{theta_u_hat}; their
#'   standard errors \code{se_l} and \code{se_u}; the critical value
#'   \code{c_alpha}; and the inputs used.
#' @references Imbens, G. W. and C. F. Manski (2004). Confidence intervals for
#'   partially identified parameters. \emph{Econometrica} 72(6), 1845--1857.
#'
#'   Stoye, J. (2009). More on confidence intervals for partially identified
#'   parameters. \emph{Econometrica} 77(4), 1299--1315.
#' @export
#' @examples
#' set.seed(1)
#' n <- 300
#' z <- rnorm(n); v <- rnorm(n); u <- 0.5 * v + rnorm(n)
#' x <- 0.6 * z + v; y <- 1 + 0.5 * x + u
#' df <- data.frame(x = x, y = y, z = z)
#' g <- function(rho, r) as.numeric(cr_Amat(r) %*% cr_delta(rho))
#' res <- im_stoye_inference(df, seq(0, 0.8, length.out = 20), g_fun = g,
#'                           B = 199)
#' res$CI_IM
#' res$reject_H0_0_in_C
im_stoye_inference <- function(data,
                               grid_r_xu,
                               g_fun,                 # function(rho_vec, r_xu) -> lambda
                               alpha = 0.05,
                               se_method = c("bootstrap", "delta"),
                               B = 999L,
                               grad_fun = NULL,       # function(rho_vec, r_xu) -> gradient wrt rho (for delta)
                               cov_rho = NULL,        # Var(rho_hat) (not scaled by n); if NULL and delta, bootstrap it
                               seed = 1L,
                               x_col = "x",
                               y_col = "y",
                               z_col = "z",
                               return_boot = FALSE) {

  se_method <- match.arg(se_method)

  # Validate and prepare data
  if (!all(c(x_col, y_col, z_col) %in% names(data))) {
    stop("Data must contain columns: ", paste(c(x_col, y_col, z_col), collapse = ", "))
  }

  d <- data[complete.cases(data[, c(x_col, y_col, z_col)]), ]
  n <- nrow(d)

  if (n < 30) {
    warning(sprintf("Small sample size (n=%d). Normal/IM approximations may be rough.", n))
  }

  # Step 1: Estimate rho
  rho_hat <- estimate_rho(d, x_col, y_col, z_col)

  # Step 2: Evaluate lambda on grid and find endpoints
  lambda_hat <- eval_g_on_grid(g_fun, rho_hat, grid_r_xu)

  theta_l_hat <- min(lambda_hat[is.finite(lambda_hat)])
  theta_u_hat <- max(lambda_hat[is.finite(lambda_hat)])

  if (!is.finite(theta_l_hat) || !is.finite(theta_u_hat)) {
    stop("Non-finite endpoint(s). Check g_fun() and grid_r_xu for invalid points.")
  }

  b_l_hat <- which.min(lambda_hat)
  b_u_hat <- which.max(lambda_hat)
  Delta_hat <- theta_u_hat - theta_l_hat

  # Step 3: Compute standard errors for endpoints
  boot_endpoints <- NULL

  .cr_local_seed(seed)

  if (se_method == "bootstrap") {
    # Bootstrap standard errors
    boot_endpoints <- matrix(NA_real_, nrow = B, ncol = 2,
                             dimnames = list(NULL, c("theta_l", "theta_u")))

    for (bb in seq_len(B)) {
      idx <- sample.int(n, n, replace = TRUE)
      rho_b <- estimate_rho(d[idx, ], x_col, y_col, z_col)
      lam_b <- eval_g_on_grid(g_fun, rho_b, grid_r_xu)
      boot_endpoints[bb, ] <- c(min(lam_b, na.rm = TRUE), max(lam_b, na.rm = TRUE))
    }

    se_l <- stats::sd(boot_endpoints[, "theta_l"], na.rm = TRUE)
    se_u <- stats::sd(boot_endpoints[, "theta_u"], na.rm = TRUE)

  } else {
    # Delta method standard errors
    if (is.null(grad_fun)) {
      stop("For se_method='delta', you must supply grad_fun(rho, r_xu).")
    }

    # If cov_rho not provided, estimate Var(rho_hat) by bootstrap
    if (is.null(cov_rho)) {
      rho_boot <- matrix(NA_real_, nrow = B, ncol = length(rho_hat))
      colnames(rho_boot) <- names(rho_hat)

      for (bb in seq_len(B)) {
        idx <- sample.int(n, n, replace = TRUE)
        rho_boot[bb, ] <- estimate_rho(d[idx, ], x_col, y_col, z_col)
      }

      cov_rho <- stats::cov(rho_boot, use = "pairwise.complete.obs")
    }

    # Compute gradients at endpoints
    grad_l <- as.numeric(grad_fun(rho_hat, grid_r_xu[b_l_hat]))
    grad_u <- as.numeric(grad_fun(rho_hat, grid_r_xu[b_u_hat]))

    # Validate dimensions
    if (!is.matrix(cov_rho) || ncol(cov_rho) != length(grad_l) || ncol(cov_rho) != length(grad_u)) {
      stop("Dimension mismatch: cov_rho must be a square matrix matching length(grad_fun()).")
    }

    # Delta method SEs
    se_l <- sqrt(drop(t(grad_l) %*% cov_rho %*% grad_l))
    se_u <- sqrt(drop(t(grad_u) %*% cov_rho %*% grad_u))
  }

  # Step 4: Compute IM critical value and confidence interval
  se_max <- max(se_l, se_u)

  if (!is.finite(se_max) || se_max <= 0) {
    # Fallback to standard normal quantile
    c_alpha <- stats::qnorm(1 - alpha / 2)
    ci_lower <- theta_l_hat
    ci_upper <- theta_u_hat
  } else {
    t_val <- Delta_hat / se_max
    c_alpha <- solve_c_im(alpha, t_val)
    ci_lower <- theta_l_hat - c_alpha * se_l
    ci_upper <- theta_u_hat + c_alpha * se_u
  }

  # Step 5: Membership test (reject H0: 0 in C iff 0 not in CI_IM)
  reject_H0 <- !(ci_lower <= 0 && 0 <= ci_upper)

  # Return results
  out <- list(
    n = n,
    rho_hat = rho_hat,
    grid_r_xu = grid_r_xu,
    lambda_hat = lambda_hat,
    b_l_hat = b_l_hat,
    b_u_hat = b_u_hat,
    theta_l_hat = theta_l_hat,
    theta_u_hat = theta_u_hat,
    Delta_hat = Delta_hat,
    se_l = se_l,
    se_u = se_u,
    c_alpha = c_alpha,
    CI_IM = c(lower = ci_lower, upper = ci_upper),
    reject_H0_0_in_C = reject_H0,
    alpha = alpha,
    se_method = se_method
  )

  if (return_boot) out$boot_endpoints <- boot_endpoints

  out
}
