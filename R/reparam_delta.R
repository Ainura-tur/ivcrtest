
#' Exact linear reparameterisation of the finite-grid map
#'
#' Implements Lemma 3.1 and Remark A.1 of the CR test paper: the grid map
#' \eqn{g(\rho; r)} is exactly linear in the two-dimensional reparameterisation
#' \eqn{\delta = (\rho_{xz}, h(\rho))'} with a coefficient matrix \eqn{A} that
#' depends on the grid alone and contains no estimated quantity.
#'
#' \code{cr_h} evaluates \eqn{h(\rho) = (\rho_{xy}\rho_{xz} - \rho_{yz}) /
#' \sqrt{1 - \rho_{xy}^2}}; \code{cr_delta} stacks it with \eqn{\rho_{xz}};
#' \code{cr_jacobian} returns \eqn{J(\rho) = \partial\delta/\partial\rho'} of
#' equation (11); \code{cr_Amat} builds the known \eqn{K \times 2} matrix with
#' rows \eqn{(r_b, -\sqrt{1 - r_b^2})}; \code{cr_omega_delta} forms
#' \eqn{\Omega_\delta = J \Sigma_\rho J'}; and \code{cr_grid_map} returns the
#' linear map \eqn{\delta \mapsto A\delta} in the interface expected by
#' \code{CIhybrid}.
#'
#' @param rho Numeric vector of length 3, ordered (rho_xy, rho_xz, rho_yz).
#' @param r_grid Numeric vector of grid points in (-1, 1).
#' @param Sigma_rho 3x3 sampling covariance matrix of rho_hat, on the estimator
#'   scale (see the scale convention note in Appendix B.1).
#' @param Amat A K x 2 matrix as returned by \code{cr_Amat}.
#' @param tol Numeric tolerance used to clamp correlations inside (-1, 1).
#' @return \code{cr_h} a scalar; \code{cr_delta} a named length-2 vector;
#'   \code{cr_jacobian} a 2x3 matrix; \code{cr_Amat} a K x 2 matrix;
#'   \code{cr_omega_delta} a 2x2 matrix; \code{cr_grid_map} a function.
#' @name reparam_delta
NULL

#' @rdname reparam_delta
#' @export
cr_h <- function(rho, tol = 1e-10) {
  rho <- as.numeric(rho)
  if (length(rho) != 3L) stop("cr_h: rho must have 3 elements (rho_xy, rho_xz, rho_yz).")
  rho_xy <- min(max(rho[1], -1 + tol), 1 - tol)
  (rho_xy * rho[2] - rho[3]) / sqrt(1 - rho_xy^2)
}

#' @rdname reparam_delta
#' @export
cr_delta <- function(rho, tol = 1e-10) {
  rho <- as.numeric(rho)
  c(rho_xz = rho[2], h = cr_h(rho, tol))
}

#' @rdname reparam_delta
#' @export
cr_jacobian <- function(rho, tol = 1e-10) {
  rho <- as.numeric(rho)
  if (length(rho) != 3L) stop("cr_jacobian: rho must have 3 elements.")
  rho_xy <- min(max(rho[1], -1 + tol), 1 - tol)
  rho_xz <- rho[2]
  rho_yz <- rho[3]
  s <- sqrt(1 - rho_xy^2)

  J <- matrix(0, nrow = 2L, ncol = 3L,
              dimnames = list(c("rho_xz", "h"),
                              c("rho_xy", "rho_xz", "rho_yz")))
  J[1L, ] <- c(0, 1, 0)
  J[2L, ] <- c(rho_xz / s + (rho_xy * rho_xz - rho_yz) * rho_xy / s^3,
               rho_xy / s,
               -1 / s)
  J
}

#' @rdname reparam_delta
#' @export
cr_Amat <- function(r_grid, tol = 1e-10) {
  r <- pmin(pmax(as.numeric(r_grid), -1 + tol), 1 - tol)
  A <- cbind(r, -sqrt(pmax(1 - r^2, 0)))
  dimnames(A) <- list(NULL, c("rho_xz", "h"))
  A
}

#' @rdname reparam_delta
#' @export
cr_omega_delta <- function(rho, Sigma_rho, tol = 1e-10) {
  Sigma_rho <- as.matrix(Sigma_rho)
  if (any(dim(Sigma_rho) != c(3L, 3L))) {
    stop("cr_omega_delta: Sigma_rho must be a 3x3 matrix for (rho_xy, rho_xz, rho_yz).")
  }
  J  <- cr_jacobian(rho, tol)
  Om <- J %*% Sigma_rho %*% t(J)
  (Om + t(Om)) / 2
}

#' @rdname reparam_delta
#' @export
cr_grid_map <- function(Amat) {
  Amat <- as.matrix(Amat)
  force(Amat)
  function(delta) {
    D <- if (is.null(dim(delta))) matrix(as.numeric(delta), nrow = 1L) else as.matrix(delta)
    if (ncol(D) != ncol(Amat)) D <- t(D)
    if (ncol(D) != ncol(Amat)) {
      stop("cr_grid_map: delta must have ", ncol(Amat), " columns.")
    }
    D %*% t(Amat)
  }
}
