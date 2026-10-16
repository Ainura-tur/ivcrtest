#' Valid range and breakdown range for the endogeneity magnitude
#'
#' Implements the sensitivity analysis of Section 3.5 on a grid of cut points.
#' For each cut point \eqn{\rho_k} the MCUB interval is recomputed on two
#' restricted domains, and zero-inclusion is recorded:
#' \itemize{
#'   \item restricting from below, \eqn{\mathcal D_k = (\rho_k, \bar\rho)}.
#'     Identified sets shrink as \eqn{\rho_k} rises, so acceptance is downward
#'     closed in \eqn{\rho_k};
#'   \item restricting from above, \eqn{\mathcal D'_k = (0, \rho_k)}.
#'     Acceptance is upward closed in \eqn{\rho_k}.
#' }
#' Following Section 3.5, the upper endpoint of the valid range is
#' \eqn{\bar\rho^*_{xu} = \sup\{\rho^* : 0 \in \mathrm{CI}_{\mathrm{MCUB}}
#' \textrm{ when } \mathcal D = (\rho^*, \bar\rho)\}}, taken from the
#' "from below" family, and the lower endpoint is
#' \eqn{\underline\rho^*_{xu} = \inf\{\rho^* : 0 \in \mathrm{CI}_{\mathrm{MCUB}}
#' \textrm{ when } \mathcal D = (0, \rho^*)\}}, taken from the "from above"
#' family. Assuming \eqn{\rho_{xu} \ge \rho^*} rejects the instrument once
#' \eqn{\rho^*} exceeds the upper endpoint, and assuming
#' \eqn{\rho_{xu} \le \rho^*} rejects it once \eqn{\rho^*} falls below the
#' lower endpoint.
#'
#' On a grid each endpoint is known only up to the bracket between two
#' adjacent cut points. The reported endpoint is the grid value that still
#' accepts, and the bracket is returned alongside it. Use
#' \code{\link{cr_switch_point}} to refine an endpoint inside its bracket; with
#' the same \code{Sigma_rho} and \code{seed} the two functions use identical
#' MCUB fits and therefore agree.
#'
#' The negative branch is handled by mirroring: with \code{k = -1} the domains
#' are \eqn{(-\bar\rho, -\rho_k)} and \eqn{(-\rho_k, 0)}, and the valid range is
#' reported with its sign.
#'
#' @param df Data frame with columns \code{x}, \code{y}, \code{z}, already
#'   residualised on any exogenous controls.
#' @param rho_bar Outer bound of the feasible domain, positive. Default 0.8.
#' @param cut_grid Numeric vector of interior cut points, strictly between 0 and
#'   \code{rho_bar}. Defaults to \code{seq(0.1, 0.7, by = 0.1)}.
#' @param k Sign of the maintained endogeneity restriction, +1 or -1.
#' @param alpha Numeric significance level.
#' @param K Integer number of grid points per domain.
#' @param route Only \code{"reparam"} is implemented.
#' @param Sigma_rho Optional 3x3 sampling covariance of rho_hat. When NULL it is
#'   bootstrapped once and reused for every domain, since it does not depend on
#'   the domain.
#' @param B_boot Integer bootstrap replications, used only when
#'   \code{Sigma_rho} is NULL.
#' @param seed Integer seed for the bootstrap of \code{Sigma_rho}.
#' @param B,Blarge,eta,tol,tol_r Passed through to the MCUB fit
#'   (\code{\link{CIhybrid}}).
#' @param eps Width of the edge domains \eqn{(0, \epsilon)} and
#'   \eqn{(\bar\rho - \epsilon, \bar\rho)}, fitted only when the smallest cut
#'   already accepts from above or the largest cut still accepts from below.
#'   They decide whether an endpoint reaches 0 or \eqn{\bar\rho}, as in the
#'   bisection routine used for Table 3 of the paper.
#' @param verbose Logical; print each evaluation.
#' @return A list with
#'   \item{table}{one row per cut point: \code{rho_k} (signed),
#'     \code{accept_from_below}, \code{accept_from_above};}
#'   \item{full_domain_accepts}{zero-inclusion on the whole domain
#'     \eqn{(0, \bar\rho)}, mirrored when \code{k = -1};}
#'   \item{valid_range}{the signed valid range, in increasing order, or
#'     \code{c(NA, NA)} when zero is excluded on the whole domain;}
#'   \item{lower_bracket, upper_bracket}{signed brackets locating each
#'     endpoint on the grid;}
#'   \item{status}{\code{"range"}; \code{"rejects throughout"} when zero is
#'     excluded on the full domain; or \code{"below grid resolution"} when both
#'     endpoints fall between the same pair of cut points, in which case
#'     \code{valid_range} is \code{NA} and the brackets locate it;}
#'   \item{breakdown_range}{the complement of the valid range within the
#'     signed domain, as text;}
#'   and the settings used.
#' @seealso \code{\link{cr_switch_point}} to refine an endpoint by bisection.
#' @export
#' @importFrom stats complete.cases
cr_valid_range <- function(df,
                           rho_bar   = 0.8,
                           cut_grid  = seq(0.1, 0.7, by = 0.1),
                           k         = 1,
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
                           eps       = 1e-3,
                           verbose   = TRUE) {

  route <- match.arg(route)
  if (route == "primitive") {
    stop("cr_valid_range: only route = \"reparam\" is implemented.")
  }
  if (!(eps > 0 && eps < min(cut_grid, rho_bar - max(cut_grid)))) {
    stop("cr_valid_range: eps must be positive and smaller than the gap ",
         "between the cut grid and 0 or rho_bar.")
  }
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

  d   <- df[complete.cases(df[, c("x", "y", "z")]), , drop = FALSE]
  sgn <- sign(k)

  # One covariance for every domain: Omega_delta does not depend on D.
  rho_hat <- c(cor(d$x, d$y), cor(d$x, d$z), cor(d$y, d$z))
  if (is.null(Sigma_rho)) {
    Sigma_rho <- estimate_cov_corr_boot(d$x, d$y, d$z, B = B_boot, seed = seed)
  }
  delta <- cr_delta(rho_hat)
  Omega <- cr_omega_delta(rho_hat, Sigma_rho)

  accepts <- function(lo, hi, tag) {
    D   <- sort(sgn * c(lo, hi))
    out <- .cr_mcub_accepts(delta, Omega, D, K = K, alpha = alpha, eta = eta,
                            B = B, Blarge = Blarge, tol = tol, tol_r = tol_r)
    if (isTRUE(verbose)) {
      cat(sprintf("  %-14s D = [%6.3f, %6.3f]  zero in CI_MCUB: %s\n",
                  tag, D[1], D[2], out))
    }
    out
  }

  full <- accepts(0, rho_bar, "full domain")

  if (isTRUE(verbose)) cat("Restricting from below, D = (rho_k, rho_bar):\n")
  acc_from_below <- vapply(cut_grid, function(rk) {
    accepts(rk, rho_bar, sprintf("rho_k = %.2f", rk))
  }, logical(1))

  if (isTRUE(verbose)) cat("Restricting from above, D = (0, rho_k):\n")
  acc_from_above <- vapply(cut_grid, function(rk) {
    accepts(0, rk, sprintf("rho_k = %.2f", rk))
  }, logical(1))

  # Edge fits, as in zero_inclusion_range() of the replication package: when
  # the smallest cut already accepts from above, check (0, eps) to tell a
  # lower endpoint of 0 from one inside (0, cut_1]; likewise at rho_bar.
  edge_lo <- if (isTRUE(full) && isTRUE(acc_from_above[1]))
    accepts(0, eps, sprintf("edge (0, %.3f)", eps)) else NA
  edge_hi <- if (isTRUE(full) && isTRUE(acc_from_below[length(cut_grid)]))
    accepts(rho_bar - eps, rho_bar, sprintf("edge (%.3f, rho_bar)", rho_bar - eps)) else NA

  ends <- .cr_range_from_table(cut_grid, acc_from_below, acc_from_above,
                               full, rho_bar, eps, edge_lo, edge_hi)

  tab <- data.frame(rho_k             = sgn * cut_grid,
                    accept_from_below = acc_from_below,
                    accept_from_above = acc_from_above,
                    stringsAsFactors  = FALSE)

  signed <- function(v) sort(sgn * v)
  if (ends$status == "range") {
    valid     <- signed(c(ends$lower, ends$upper))
    breakdown <- .cr_breakdown_text(ends$lower, ends$upper, rho_bar, sgn)
  } else if (ends$status == "below grid resolution") {
    valid     <- c(NA_real_, NA_real_)
    breakdown <- "not resolved on this grid: see lower_bracket and upper_bracket"
  } else {
    valid     <- c(NA_real_, NA_real_)
    breakdown <- "entire domain: zero is excluded even on the full domain"
  }

  list(table               = tab,
       full_domain_accepts = full,
       valid_range         = valid,
       lower_bracket       = signed(ends$lower_bracket),
       upper_bracket       = signed(ends$upper_bracket),
       status              = ends$status,
       breakdown_range     = breakdown,
       rho_bar             = sgn * rho_bar,
       k                   = k,
       alpha               = alpha,
       K                   = K,
       route               = route)
}

# One MCUB fit on domain D: TRUE when zero lies in CI_MCUB, NA when the fit
# fails. Shared by cr_valid_range() and cr_switch_point() so that the two use
# identical fits for the same inputs.
.cr_mcub_accepts <- function(delta, Omega, D, K, alpha, eta, B, Blarge,
                             tol, tol_r, seed = 0L) {
  A <- cr_Amat(seq(D[1], D[2], length.out = K))
  r <- try(CIhybrid(delta, Omega, A, A, alpha = alpha, alphac = 0.8 * alpha,
                    eta = eta, B = B, Blarge = Blarge, tol = tol,
                    tol_r = tol_r, index = NULL, g = cr_grid_map(A),
                    seed = seed), silent = TRUE)
  if (inherits(r, "try-error") || any(!is.finite(r$CI_h))) return(NA)
  r$CI_h[1] <= 0 && 0 <= r$CI_h[2]
}

# Endpoints of the valid range from the accept/reject table, on the |rho| scale.
# from_below: acceptance on (rho_k, rho_bar), downward closed in rho_k.
# from_above: acceptance on (0, rho_k), upward closed in rho_k.
# full:       acceptance on (0, rho_bar).
# The upper endpoint comes from from_below and the lower from from_above
# (Section 3.5). Each endpoint is the accepting grid value nearest the switch,
# with the bracket in which the switch lies.
.cr_range_from_table <- function(cut_grid, from_below, from_above, full,
                                 rho_bar, eps = 1e-3,
                                 edge_lo = NA, edge_hi = NA) {
  if (anyNA(c(from_below, from_above, full))) {
    stop("cr_valid_range: an MCUB fit failed; the table contains NA.")
  }
  na_out <- list(lower = NA_real_, upper = NA_real_,
                 lower_bracket = c(NA_real_, NA_real_),
                 upper_bracket = c(NA_real_, NA_real_))
  if (!isTRUE(full)) return(c(list(status = "rejects throughout"), na_out))

  # Nesting implies a single switch in each family. Monte Carlo noise in the
  # critical values can break this near the switch; say so rather than guess.
  monotone_down <- function(a) !any(diff(as.integer(a)) > 0)
  if (!monotone_down(from_below) || !monotone_down(rev(from_above))) {
    warning("cr_valid_range: acceptance is not monotone across cut points; ",
            "endpoints are taken from the outermost switch. Refine with ",
            "cr_switch_point().", call. = FALSE)
  }

  n <- length(cut_grid)

  # Upper endpoint (Section 3.5): sup of rho* with (rho*, rho_bar) accepting.
  # The reported value is the largest accepting cut, so it never overstates.
  if (!any(from_below)) {
    upper <- 0;               upper_bracket <- c(0, cut_grid[1])
  } else {
    j <- max(which(from_below))
    if (j < n) {
      upper <- cut_grid[j];   upper_bracket <- c(cut_grid[j], cut_grid[j + 1L])
    } else if (isTRUE(edge_hi)) {
      upper <- rho_bar;       upper_bracket <- c(rho_bar, rho_bar)
    } else {
      upper <- cut_grid[n];   upper_bracket <- c(cut_grid[n], rho_bar - eps)
    }
  }

  # Lower endpoint: inf of rho* with (0, rho*) accepting.
  # The reported value is the smallest accepting cut, so it never understates.
  if (!any(from_above)) {
    lower <- rho_bar;         lower_bracket <- c(cut_grid[n], rho_bar)
  } else {
    i <- min(which(from_above))
    if (i > 1L) {
      lower <- cut_grid[i];   lower_bracket <- c(cut_grid[i - 1L], cut_grid[i])
    } else if (isTRUE(edge_lo)) {
      lower <- 0;             lower_bracket <- c(0, 0)
    } else {
      lower <- cut_grid[1];   lower_bracket <- c(eps, cut_grid[1])
    }
  }

  status <- if (lower <= upper) "range" else "below grid resolution"
  list(status = status, lower = lower, upper = upper,
       lower_bracket = lower_bracket, upper_bracket = upper_bracket)
}

.cr_breakdown_text <- function(lower, upper, rho_bar, sgn) {
  fmt  <- function(a, b, open_a, open_b) {
    v <- sort(sgn * c(a, b))
    sprintf("%s%.3f, %.3f%s", if (open_a) "(" else "[", v[1], v[2],
            if (open_b) ")" else "]")
  }
  parts <- c(if (lower > 0)       fmt(0, lower, TRUE, TRUE),
             if (upper < rho_bar) fmt(upper, rho_bar, TRUE, TRUE))
  if (!length(parts)) "empty: the instrument survives the whole domain" else
    paste(parts, collapse = " union ")
}

#' Pointwise compatibility range for the endogeneity magnitude
#'
#' Reports the set of single assumed values \eqn{\rho_{xu} = t} at which the
#' orthogonality restriction remains compatible with the data, that is
#' \eqn{\{t : 0 \in [g(\rho;t) \pm z_{1-\alpha/2}\,\sigma(t)]\}} with
#' \eqn{\sigma(t)^2 = A(t)\hat\Omega_\delta A(t)'} and
#' \eqn{A(t) = (t, -\sqrt{1-t^2})}.
#'
#' It is a fast screen, not the valid range of Section 3.5. That range is
#' defined through the MCUB test on nested domains: \code{\link{cr_valid_range}}
#' computes it on a grid of cut points and \code{\link{cr_switch_point}} refines
#' an endpoint by bisection. The pointwise set needs no MCUB fit and usually
#' lies close to the valid range, which is the more conservative of the two.
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
#'   Spaced finely near zero: a strong instrument with small
#'   \eqn{\hat\Omega_\delta} can switch below 0.05, and a grid starting there
#'   would report \code{"rejects throughout"} rather than locating it.
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
                            scan_grid = c(0.005, 0.01, 0.02, 0.05, 0.10,
                                          0.20, 0.35, 0.50, 0.65, 0.75),
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

  if (route == "primitive") {
    stop("cr_switch_point: only route = \"reparam\" is implemented.")
  }

  n_fits <- 0L
  accepts <- function(rs) {
    D   <- sort(sgn * if (direction == "below") c(rs, rho_bar) else c(0, rs))
    out <- .cr_mcub_accepts(delta, Omega, D, K = K, alpha = alpha, eta = eta,
                            B = B, Blarge = Blarge, tol = tol, tol_r = tol_r)
    n_fits <<- n_fits + 1L
    # A failed fit used to count silently as a rejection, which could move the
    # bracket. Stop instead.
    if (is.na(out)) {
      stop(sprintf("cr_switch_point: the MCUB fit failed on D = [%.3f, %.3f].",
                   D[1], D[2]))
    }
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
