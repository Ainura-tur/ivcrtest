# Namespace protection (avoid MASS/stats masking dplyr verbs)
select    <- dplyr::select;      filter    <- dplyr::filter
mutate    <- dplyr::mutate;      slice     <- dplyr::slice
recode    <- dplyr::recode;      rename    <- dplyr::rename
summarise <- dplyr::summarise;   summarize <- dplyr::summarize
arrange   <- dplyr::arrange;     count     <- dplyr::count;   lag <- dplyr::lag

#' Perform the Instrumental Variable Correlation Restriction test
#'
#' Tests whether the orthogonality value \eqn{\rho_{zu} = 0} lies in a
#' parameter-coverage confidence interval for the scalar partially identified
#' parameter \eqn{\theta_0 = \rho_{zu}}, given a researcher-specified admissible
#' domain D for the unobserved endogeneity correlation \eqn{\rho_{xu}}.
#'
#' Exogenous controls are removed by residualisation only, as in Appendix A.1:
#' \eqn{\tilde y = M_H y}, \eqn{\tilde x = M_H x} and \eqn{\tilde z = M_H z} with
#' \eqn{M_H = I - H(H'H)^{-1}H'}. No further transformation is applied to the
#' instrument.
#'
#' If the residualised instrument has \eqn{\hat\rho_{xz} < 0} it is replaced by
#' its negative, which is without loss of generality under Lemma 2.3 since
#' \eqn{E[z'u] = 0} if and only if \eqn{E[(-z)'u] = 0}. This negates the
#' identified set and leaves the zero-inclusion decision unchanged.
#'
#' @param data A data frame containing all variables.
#' @param X Character. Name of the endogenous regressor.
#' @param Y Character. Name of the dependent variable.
#' @param H Character vector of exogenous control names; may be NULL or empty,
#'   in which case the variables are only demeaned.
#' @param Z Character vector of one or more candidate instrument names. Each is
#'   tested separately.
#' @param n Ignored; determined internally from \code{data}.
#' @param k Sign of the maintained endogeneity restriction: \code{+1} when
#'   corr(x, u) is believed positive, \code{-1} when negative. Note that the
#'   default changed to \code{+1}, since the previous default of \code{-1} was
#'   inconsistent with the default positive \code{rxu_range}.
#' @param alpha Numeric significance level. Forwarded, not ignored.
#' @param seed Integer seed. Forwarded, not ignored.
#' @param rxu_range Numeric vector of length 2 giving D. When NULL it defaults to
#'   \code{c(0, 0.8)} if \code{k > 0} and \code{c(-0.8, 0)} if \code{k < 0}.
#' @param K Integer number of equally spaced grid points. Default 80.
#' @param route Either \code{"reparam"} (default, Remark A.1) or
#'   \code{"primitive"}.
#' @param B_boot Integer bootstrap replications for the covariance of rho_hat.
#' @param compute_pvalue Logical; whether to invert over the admissible range.
#' @param common_sample Logical. When FALSE, the default, each instrument is
#'   evaluated on its own complete-case sample, so a missing value in one
#'   candidate instrument does not shrink the sample used for the others. This
#'   is what the published applications do. Set TRUE to force the intersection
#'   sample across all instruments in \code{Z}, which makes the columns
#'   directly comparable but discards observations.
#' @param verbose Logical; print progress.
#' @param ... Further arguments passed to \code{\link{check_compatibility}}.
#' @return An object of class \code{"iv_cr_test"}: a list with \code{results}
#'   (a data frame, one row per instrument), \code{diagnostics}, \code{domain},
#'   \code{route} and the call.
#' @export
#' @importFrom stats complete.cases as.formula lm resid cor
#' @examples
#' \donttest{
#' data <- wooldridge::card
#' H <- c("exper", "expersq", "black", "south", "smsa", "smsa66",
#'        "reg661", "reg662", "reg663", "reg664",
#'        "reg665", "reg666", "reg667", "reg668")
#' # Card (1993): ability bias implies corr(educ, u) > 0, hence k = +1.
#' res <- iv_cr_test(data, X = "educ", Y = "lwage", H = H,
#'                   Z = c("fatheduc", "motheduc", "nearc2", "nearc4"),
#'                   k = 1, alpha = 0.05, rxu_range = c(0, 0.8))
#' print(res)
#' }

iv_cr_test <- function(data, X, Y, H, Z,
                       n              = NULL,
                       k              = 1,
                       alpha          = 0.05,
                       seed           = 123L,
                       rxu_range      = NULL,
                       K              = 80L,
                       route          = c("reparam", "primitive"),
                       B_boot         = 800L,
                       compute_pvalue = TRUE,
                       common_sample  = FALSE,
                       verbose        = TRUE,
                       ...) {

  route <- match.arg(route)
  cl    <- match.call()

  if (!is.data.frame(data)) stop("iv_cr_test: 'data' must be a data frame.")
  H <- if (is.null(H)) character(0) else as.character(H)
  Z <- as.character(Z)
  if (length(Z) < 1L) stop("iv_cr_test: supply at least one instrument name in Z.")
  if (!(length(k) == 1L && k != 0)) stop("iv_cr_test: k must be a single nonzero value.")

  vars_needed <- unique(c(Y, X, H, Z))
  missing_vars <- setdiff(vars_needed, names(data))
  if (length(missing_vars) > 0L) {
    stop("iv_cr_test: variables not found in data: ", paste(missing_vars, collapse = ", "))
  }

  # ---- admissible domain, tied to the sign restriction ----
  if (is.null(rxu_range)) {
    rxu_range <- if (k > 0) c(0, 0.8) else c(-0.8, 0)
  }
  rxu_range <- sort(as.numeric(rxu_range))
  if (k > 0 && rxu_range[1] < 0) {
    stop("iv_cr_test: k = ", k, " asserts corr(x, u) > 0, so rxu_range must lie in [0, 1).")
  }
  if (k < 0 && rxu_range[2] > 0) {
    stop("iv_cr_test: k = ", k, " asserts corr(x, u) < 0, so rxu_range must lie in (-1, 0].")
  }

  set.seed(seed)

  # ---- residualisation on H only, per Appendix A.1 ----
  rhs <- if (length(H) > 0L) paste(H, collapse = " + ") else "1"
  resid_on_H <- function(v, dat) {
    resid(lm(as.formula(paste0("`", v, "` ~ ", rhs)), data = dat))
  }

  # Sample selection. By default each instrument keeps its own complete cases;
  # pooling them would let a missing value in one candidate shrink the sample
  # used for every other, which is not what the applications do.
  common_rows <- if (isTRUE(common_sample)) {
    complete.cases(data[, vars_needed, drop = FALSE])
  } else NULL

  results <- list()
  diags   <- list()
  n_by_Z  <- integer(0)

  for (j in seq_along(Z)) {
    vj  <- unique(c(Y, X, H, Z[j]))
    sel <- if (is.null(common_rows)) complete.cases(data[, vj, drop = FALSE]) else common_rows
    dat <- data[sel, vj, drop = FALSE]
    n   <- nrow(dat)
    if (n < 30L) {
      warning("iv_cr_test: fewer than 30 complete observations for ", Z[j],
              "; skipping.", call. = FALSE)
      next
    }
    n_by_Z[Z[j]] <- n

    x  <- resid_on_H(X, dat)
    y  <- resid_on_H(Y, dat)
    zj <- resid_on_H(Z[j], dat)

    # Lemma 2.3 normalisation: work with the orientation giving rho_xz > 0.
    flipped <- FALSE
    if (cor(x, zj) < 0) {
      zj <- -zj
      flipped <- TRUE
    }

    df <- data.frame(x = x, y = y, z = zj)

    if (isTRUE(verbose)) {
      cat(sprintf("\n=== CR test: instrument %s (n = %d, D = [%.2f, %.2f], K = %d, %s) ===\n",
                  Z[j], n, rxu_range[1], rxu_range[2], K, route))
    }

    res <- check_compatibility(df, i = j, k = k,
                               alpha = alpha, rxu_range = rxu_range,
                               K = K, route = route, B_boot = B_boot,
                               seed = seed, compute_pvalue = compute_pvalue,
                               verbose = verbose, ...)
    if (is.null(res)) next

    res$Z <- Z[j]
    res$sign_flipped <- flipped
    diags[[Z[j]]] <- attr(res, "diagnostics")
    attr(res, "diagnostics") <- NULL
    results[[j]] <- res
  }

  if (!length(results)) stop("iv_cr_test: no instrument had enough complete observations.")

  out <- list(results     = do.call(rbind, results),
              diagnostics = diags,
              n_by_Z      = n_by_Z,
              domain      = rxu_range,
              k           = k,
              alpha       = alpha,
              K           = K,
              route       = route,
              n           = if (length(unique(n_by_Z)) == 1L) unname(n_by_Z[1]) else n_by_Z,
              call        = cl)
  class(out) <- "iv_cr_test"
  out
}

#' Print method for CR test results
#'
#' @param x An object of class \code{"iv_cr_test"}.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @export
print.iv_cr_test <- function(x, ...) {
  cat("Instrumental Variable Correlation Restriction test\n")
  nlab <- if (length(x$n) == 1L) format(x$n) else "varies by instrument, see n"
  cat(sprintf("  n = %s, D = [%.2f, %.2f], K = %d, alpha = %.3f, route = %s\n\n",
              nlab, x$domain[1], x$domain[2], x$K, x$alpha, x$route))
  cols <- c("Z", "n", "plug_in", "CI_MCUB", "CI_IM", "CI_naive",
            "Zero_in_CI_MCUB", "Zero_in_CI_IM", "p_adm_label", "G_hat", "status")
  cols <- intersect(cols, names(x$results))
  print(x$results[, cols, drop = FALSE], row.names = FALSE)
  cat("\nMCUB is the primary procedure; IM/Stoye and the naive endpoint interval",
      "\nare reported as diagnostics. p_adm is censored to the admissible range.\n")
  invisible(x)
}
