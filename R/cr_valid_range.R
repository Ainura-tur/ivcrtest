
#' Valid range and breakdown range for the endogeneity magnitude
#'
#' Implements the sensitivity analysis of Section 3.5. For a sequence of cut
#' points \eqn{\rho_k}, the MCUB interval is recomputed on the restricted domains
#' \eqn{D_k = (\rho_k, \bar\rho)} and \eqn{D'_k = (0, \rho_k)} and zero-inclusion
#' is recorded. The reported valid range is
#' \eqn{[\rho^*_{xu}, \rho^{*}_{xu}]}, where the lower endpoint is the smallest
#' \eqn{\rho_k} at which restricting from below yields acceptance and the upper
#' endpoint is the largest \eqn{\rho_k} at which restricting from above yields
#' acceptance. The breakdown range is the complement within \eqn{(0, \bar\rho)}.
#'
#' A narrow valid range means acceptance is sensitive to the assumed magnitude
#' of endogeneity; a wide one means validity is robust across plausible
#' magnitudes. A breakdown range reaching very small values means the instrument
#' is rejected even under mild endogeneity assumptions.
#'
#' The negative branch is handled by mirroring: with \code{k = -1} the domains
#' are \eqn{(-\bar\rho, \rho_k)} and \eqn{(\rho_k, 0)} for negative cut points.
#'
#' @param df Data frame with columns \code{x}, \code{y}, \code{z}, already
#'   residualised on any exogenous controls.
#' @param rho_bar Outer bound of the feasible domain, positive. Default 0.8.
#' @param cut_grid Numeric vector of interior cut points, strictly between 0 and
#'   \code{rho_bar}. Defaults to \code{seq(0.1, 0.7, by = 0.1)}.
#' @param k Sign of the maintained endogeneity restriction, +1 or -1.
#' @param alpha Numeric significance level.
#' @param K Integer number of grid points per domain.
#' @param route Either \code{"reparam"} or \code{"primitive"}.
#' @param B_boot Integer bootstrap replications for the covariance of rho_hat.
#' @param seed Integer seed.
#' @param verbose Logical; print progress.
#' @param ... Further arguments passed to \code{\link{check_compatibility}}.
#' @return A list with the per-cut-point \code{table}, the \code{valid_range},
#'   the \code{breakdown_range} and the settings used.
#' @export
#' @importFrom stats complete.cases

cr_valid_range <- function(df,
                           rho_bar  = 0.8,
                           cut_grid = seq(0.1, 0.7, by = 0.1),
                           k        = 1,
                           alpha    = 0.05,
                           K        = 80L,
                           route    = c("reparam", "primitive"),
                           B_boot   = 800L,
                           seed     = 1L,
                           verbose  = TRUE,
                           ...) {

  route <- match.arg(route)
  if (!all(c("x", "y", "z") %in% names(df))) {
    stop("cr_valid_range: df must contain columns 'x', 'y', 'z'.")
  }
  if (!(rho_bar > 0 && rho_bar < 1)) {
    stop("cr_valid_range: rho_bar must lie strictly between 0 and 1.")
  }
  cut_grid <- sort(unique(as.numeric(cut_grid)))
  if (any(cut_grid <= 0) || any(cut_grid >= rho_bar)) {
    stop("cr_valid_range: cut_grid must lie strictly between 0 and rho_bar.")
  }

  sgn <- sign(k)

  accepts <- function(lo, hi, tag) {
    D <- sgn * sort(c(lo, hi))
    D <- sort(D)
    res <- tryCatch(
      check_compatibility(df, i = 1L, k = k, alpha = alpha, rxu_range = D,
                          K = K, route = route, B_boot = B_boot, seed = seed,
                          compute_pvalue = FALSE, verbose = FALSE, ...),
      error = function(e) NULL
    )
    if (is.null(res)) return(NA)
    if (isTRUE(verbose)) {
      cat(sprintf("  %-14s D = [%6.3f, %6.3f]  CI_MCUB = %-18s zero: %s\n",
                  tag, D[1], D[2], res$CI_MCUB, res$Zero_in_CI_MCUB))
    }
    identical(res$Zero_in_CI_MCUB, "\u2713")
  }

  if (isTRUE(verbose)) cat("Restricting from below, D = (rho_k, rho_bar):\n")
  acc_from_below <- vapply(cut_grid, function(rk) {
    accepts(rk, rho_bar, sprintf("rho_k = %.2f", rk))
  }, logical(1))

  if (isTRUE(verbose)) cat("Restricting from above, D = (0, rho_k):\n")
  acc_from_above <- vapply(cut_grid, function(rk) {
    accepts(0, rk, sprintf("rho_k = %.2f", rk))
  }, logical(1))

  lower <- if (any(acc_from_below, na.rm = TRUE)) {
    min(cut_grid[which(acc_from_below)])
  } else NA_real_

  upper <- if (any(acc_from_above, na.rm = TRUE)) {
    max(cut_grid[which(acc_from_above)])
  } else NA_real_

  tab <- data.frame(rho_k             = sgn * cut_grid,
                    accept_from_below = acc_from_below,
                    accept_from_above = acc_from_above,
                    stringsAsFactors  = FALSE)

  # sort AFTER applying the sign: with k = -1, sgn * sort(c(lower, upper))
  # returns the endpoints in descending order.
  # The paper notes the valid range may be one-sided. If no domain of the form
  # (rho*, rho_bar) accepts, the lower endpoint is not identified from that
  # direction but the range is still informative from above, and conversely.
  # Collapsing either case to NA would discard that.
  one_sided <- "no"
  if (is.na(lower) && !is.na(upper)) {
    lower <- 0
    one_sided <- "upper only: no domain of the form (rho*, rho_bar) accepts"
  } else if (!is.na(lower) && is.na(upper)) {
    upper <- rho_bar
    one_sided <- "lower only: no domain of the form (0, rho*) accepts"
  }

  valid <- if (any(is.na(c(lower, upper)))) c(NA_real_, NA_real_) else
             sort(sgn * c(lower, upper))

  breakdown <- if (any(is.na(valid))) {
    "entire domain (no cut point yields acceptance in either direction)"
  } else {
    sprintf("(%.2f, %.2f) union (%.2f, %.2f)",
            min(0, sgn * rho_bar), min(valid), max(valid), max(0, sgn * rho_bar))
  }

  list(table           = tab,
       valid_range     = valid,
       one_sided       = one_sided,
       breakdown_range = breakdown,
       rho_bar         = sgn * rho_bar,
       k               = k,
       alpha           = alpha,
       K               = K,
       route           = route)
}

#' Pointwise compatibility range for the endogeneity magnitude
#'
#' Reports the set of single assumed values \eqn{\rho_{xu} = t} at which the
#' orthogonality restriction remains compatible with the data, that is
#' \eqn{\{t : 0 \in [g(\rho;t) \pm z_{1-\alpha/2}\,\sigma(t)]\}} with
#' \eqn{\sigma(t)^2 = A(t)\hat\Omega_\delta A(t)'} and
#' \eqn{A(t) = (t, -\sqrt{1-t^2})}.
#'
#' This is the quantity the Table 5 note describes, the subinterval of assumed
#' \eqn{\rho_{xu}} values for which zero lies in the interval. It differs from
#' \code{\link{cr_valid_range}}, which implements the two-sided infimum and
#' supremum construction over nested domains. Because identified sets are nested
#' in the domain, that construction returns the endpoints of the search grid
#' whenever either direction admits any acceptance, so it cannot produce the
#' narrow ranges the sensitivity discussion refers to. The pointwise version
#' discriminates, needs no MCUB fit, and is therefore fast.
#'
#' At a single \eqn{t} the identified set is the point \eqn{g(\rho;t)}, so the
#' relevant interval is the ordinary normal one and no union bound arises.
#'
#' @param df Data frame with columns \code{x}, \code{y}, \code{z}, already
#'   residualised on any exogenous controls.
#' @param rxu_range Numeric length-2 search domain. Defaults to \code{c(0, 0.8)}
#'   when \code{k > 0} and \code{c(-0.8, 0)} when \code{k < 0}.
#' @param k Sign of the maintained endogeneity restriction, +1 or -1.
#' @param alpha Numeric significance level.
#' @param n_grid Integer number of points at which to evaluate.
#' @param Sigma_rho Optional 3x3 sampling covariance of rho_hat.
#' @param B_boot Integer bootstrap replications, used only when Sigma_rho is NULL.
#' @param seed Integer seed for the bootstrap.
#' @param tol Numeric clamping tolerance.
#' @return A list with \code{range}, a logical \code{contiguous}, the
#'   \code{segments} making up the set, \code{breakdown_range}, and the
#'   evaluated \code{grid}, \code{g}, \code{se} and \code{included} vectors.
#' @export
cr_pointwise_range <- function(df,
                               rxu_range = NULL,
                               k         = 1,
                               alpha     = 0.05,
                               n_grid    = 2001L,
                               Sigma_rho = NULL,
                               B_boot    = 800L,
                               seed      = 123L,
                               tol       = 1e-10) {

  if (!all(c("x", "y", "z") %in% names(df))) {
    stop("cr_pointwise_range: df must contain columns 'x', 'y', 'z'.")
  }
  d <- df[complete.cases(df[, c("x", "y", "z")]), , drop = FALSE]

  if (is.null(rxu_range)) rxu_range <- if (k > 0) c(0, 0.8) else c(-0.8, 0)
  rxu_range <- sort(as.numeric(rxu_range))
  if (any(abs(rxu_range) >= 1)) {
    stop("cr_pointwise_range: rxu_range must lie strictly inside (-1, 1).")
  }

  rho_hat <- c(rho_xy = cor(d$x, d$y),
               rho_xz = cor(d$x, d$z),
               rho_yz = cor(d$y, d$z))
  if (is.null(Sigma_rho)) {
    Sigma_rho <- estimate_cov_corr_boot(d$x, d$y, d$z, B = B_boot, seed = seed)
  }
  Omega <- cr_omega_delta(rho_hat, Sigma_rho, tol)
  delta <- cr_delta(rho_hat, tol)

  # Keep the endpoints strictly inside the domain; A(t) is degenerate at |t| = 1.
  eps  <- diff(rxu_range) / (2 * n_grid)
  grid <- seq(rxu_range[1] + eps, rxu_range[2] - eps, length.out = n_grid)

  A  <- cr_Amat(grid, tol)
  g  <- as.numeric(A %*% delta)
  se <- sqrt(pmax(rowSums((A %*% Omega) * A), 0))   # avoids an n_grid^2 matrix

  crit <- qnorm(1 - alpha / 2)
  inc  <- abs(g) <= crit * se

  if (!any(inc)) {
    return(list(range = c(NA_real_, NA_real_), contiguous = NA,
                segments = data.frame(from = numeric(0), to = numeric(0)),
                breakdown_range = sprintf("entire domain [%.3f, %.3f]",
                                          rxu_range[1], rxu_range[2]),
                grid = grid, g = g, se = se, included = inc,
                rho_hat = rho_hat, Omega_delta = Omega, crit = crit))
  }

  # contiguous runs of TRUE
  rl   <- rle(inc)
  ends <- cumsum(rl$lengths)
  begs <- ends - rl$lengths + 1L
  keep <- rl$values
  segs <- data.frame(from = grid[begs[keep]], to = grid[ends[keep]])

  rng <- c(min(segs$from), max(segs$to))
  brk <- paste(c(if (rng[1] > rxu_range[1] + 2 * eps)
                   sprintf("[%.3f, %.3f)", rxu_range[1], rng[1]),
                 if (rng[2] < rxu_range[2] - 2 * eps)
                   sprintf("(%.3f, %.3f]", rng[2], rxu_range[2])),
               collapse = " union ")
  if (!nzchar(brk)) brk <- "empty"

  list(range = rng, contiguous = nrow(segs) == 1L, segments = segs,
       breakdown_range = brk, grid = grid, g = g, se = se, included = inc,
       rho_hat = rho_hat, Omega_delta = Omega, crit = crit)
}

#' One-sided switching points for the maintained domain
#'
#' Locates the point at which a one-sided restriction on the endogeneity
#' magnitude flips the CR test between acceptance and rejection.
#'
#' Two families are considered. Restricting from below, \eqn{\mathcal D =
#' (\rho^*, \bar\rho)}, identified sets shrink as \eqn{\rho^*} rises, so
#' acceptance is downward closed and the switch is
#' \eqn{\tau_{\mathrm{below}} = \sup\{\rho^* : 0 \in \mathrm{CI}_{\mathrm{MCUB}}\}}:
#' assuming \eqn{\rho_{xu} \ge \rho^*} rejects the instrument once \eqn{\rho^*}
#' passes it. Restricting from above, \eqn{\mathcal D = (0, \rho^*)}, acceptance
#' is upward closed and the switch is the corresponding infimum.
#'
#' Unlike \code{\link{cr_pointwise_range}} this is defined through the test
#' itself rather than through a point approximation to it, at the cost of a
#' bisection over MCUB fits. The two agree closely in practice, with the
#' switching point the more conservative of the two.
#'
#' \eqn{\hat\Omega_\delta} does not depend on the domain, so the bootstrap runs
#' once and every bisection step reuses it.
#'
#' @param df Data frame with columns \code{x}, \code{y}, \code{z}, residualised
#'   on any exogenous controls.
#' @param rho_bar Outer bound of the feasible domain, positive.
#' @param k Sign of the maintained endogeneity restriction, +1 or -1.
#' @param direction \code{"below"} for the \eqn{(\rho^*, \bar\rho)} family,
#'   \code{"above"} for \eqn{(0, \rho^*)}.
#' @param bracket Optional length-2 bracket in absolute value. When NULL a
#'   coarse scan on \code{scan_grid} locates one.
#' @param scan_grid Grid used to bracket the switch when \code{bracket} is NULL.
#' @param steps Integer bisection steps. Resolution is
#'   \code{diff(bracket) / 2^steps}.
#' @param alpha,K,route,B,Blarge,eta,tol,tol_r Passed through to the MCUB fit.
#' @param Sigma_rho Optional 3x3 covariance of rho_hat; bootstrapped when NULL.
#' @param B_boot,seed Bootstrap size and seed.
#' @param verbose Logical; report each evaluation.
#' @return A list with \code{tau} (signed), \code{direction}, the
#'   \code{bracket} used, \code{n_fits}, and a \code{status} of
#'   \code{"interior"}, \code{"accepts throughout"} or \code{"rejects throughout"}.
#' @export
cr_switch_point <- function(df,
                            rho_bar   = 0.8,
                            k         = 1,
                            direction = c("below", "above"),
                            bracket   = NULL,
                            scan_grid = seq(0.05, 0.75, by = 0.10),
                            steps     = 6L,
                            alpha     = 0.05,
                            K         = 80L,
                            route     = c("reparam", "primitive"),
                            Sigma_rho = NULL,
                            B_boot    = 800L,
                            seed      = 123L,
                            B         = 500L,
                            Blarge    = 5000L,
                            eta       = 0.001,
                            tol       = 1e-3,
                            tol_r     = 1e-3,
                            verbose   = FALSE) {

  direction <- match.arg(direction)
  route     <- match.arg(route)
  if (!all(c("x", "y", "z") %in% names(df))) {
    stop("cr_switch_point: df must contain columns 'x', 'y', 'z'.")
  }
  d   <- df[complete.cases(df[, c("x", "y", "z")]), , drop = FALSE]
  sgn <- sign(k)

  rho_hat <- c(cor(d$x, d$y), cor(d$x, d$z), cor(d$y, d$z))
  if (is.null(Sigma_rho)) {
    Sigma_rho <- estimate_cov_corr_boot(d$x, d$y, d$z, B = B_boot, seed = seed)
  }
  delta <- cr_delta(rho_hat)
  Omega <- cr_omega_delta(rho_hat, Sigma_rho)   # free of the domain

  n_fits <- 0L
  accepts <- function(rs) {
    D <- sort(sgn * if (direction == "below") c(rs, rho_bar) else c(0, rs))
    A <- cr_Amat(seq(D[1], D[2], length.out = K))
    if (route == "primitive") {
      stop("cr_switch_point: only route = \"reparam\" is implemented.")
    }
    r <- try(CIhybrid(delta, Omega, A, A, alpha = alpha, alphac = 0.8 * alpha,
                      eta = eta, B = B, Blarge = Blarge, tol = tol,
                      tol_r = tol_r, index = NULL, g = cr_grid_map(A),
                      seed = 0L), silent = TRUE)
    n_fits <<- n_fits + 1L
    if (inherits(r, "try-error")) return(NA)
    out <- isTRUE(r$CI_h[1] <= 0 && 0 <= r$CI_h[2])
    if (isTRUE(verbose)) {
      cat(sprintf("    rho* = %.4f  D = [%.3f, %.3f]  accept: %s\n",
                  rs, D[1], D[2], out))
    }
    out
  }

  if (is.null(bracket)) {
    sc  <- sort(scan_grid)
    acc <- vapply(sc, function(v) isTRUE(accepts(v)), logical(1))
    if (all(acc)) {
      return(list(tau = NA_real_, direction = direction, bracket = range(sc),
                  n_fits = n_fits, status = "accepts throughout"))
    }
    if (!any(acc)) {
      return(list(tau = NA_real_, direction = direction, bracket = range(sc),
                  n_fits = n_fits, status = "rejects throughout"))
    }
    # below: acceptance is downward closed; above: upward closed
    j <- if (direction == "below") max(which(acc)) else min(which(acc))
    bracket <- if (direction == "below") c(sc[j], sc[j + 1L]) else c(sc[j - 1L], sc[j])
  }

  lo <- bracket[1]; hi <- bracket[2]

  # Validate the bracket. Bisection silently converges to an endpoint if the
  # two sides do not straddle the switch, which would return the bracket floor
  # dressed up as a switching point.
  a_lo <- isTRUE(accepts(lo)); a_hi <- isTRUE(accepts(hi))
  want <- if (direction == "below") c(TRUE, FALSE) else c(FALSE, TRUE)
  if (identical(c(a_lo, a_hi), c(TRUE, TRUE))) {
    return(list(tau = NA_real_, direction = direction, bracket = sgn * sort(bracket),
                n_fits = n_fits, status = "accepts throughout"))
  }
  if (identical(c(a_lo, a_hi), c(FALSE, FALSE))) {
    return(list(tau = NA_real_, direction = direction, bracket = sgn * sort(bracket),
                n_fits = n_fits, status = "rejects throughout"))
  }
  if (!identical(c(a_lo, a_hi), want)) {
    stop("cr_switch_point: acceptance runs the wrong way across the bracket; ",
         "check `direction`.")
  }

  for (i in seq_len(steps)) {
    m <- (lo + hi) / 2
    a <- isTRUE(accepts(m))
    if (direction == "below") { if (a) lo <- m else hi <- m }
    else                      { if (a) hi <- m else lo <- m }
  }

  list(tau = sgn * (lo + hi) / 2, direction = direction,
       bracket = sgn * sort(bracket), n_fits = n_fits, status = "interior")
}
