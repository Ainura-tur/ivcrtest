# Regression tests for the CR test.
#
# The previous file ran one call and asserted only that the result inherited
# from data.frame, which gave no protection against numerical regressions.
# These pin the algebraic identities of Lemma 3.1 and the published Table 4.

test_that("Lemma 3.1: the grid map is exactly linear in delta", {
  rho <- c(0.94, 0.578, 0.568)
  r   <- seq(0, 0.8, length.out = 80)
  A   <- cr_Amat(r)
  expect_equal(as.numeric(A %*% cr_delta(rho)),
               as.numeric(g_xu_safe(r, rho)), tolerance = 1e-12)
  expect_equal(qr(A)$rank, 2L)
})

test_that("equation (12): both routes give identical grid standard errors", {
  rho   <- c(0.94, 0.578, 0.568)
  r     <- seq(0, 0.8, length.out = 80)
  Sigma <- matrix(c(4, 1, 1, 1, 9, 2, 1, 2, 6), 3, 3) / 1e4
  A     <- cr_Amat(r)
  grad  <- t(vapply(r, function(rr)
    as.numeric(local_compute_gradient_safe(rr, rho[1], rho[2], rho[3])), numeric(3)))
  expect_equal(grad, unname(A %*% cr_jacobian(rho)), tolerance = 1e-12)
  expect_equal(A %*% cr_omega_delta(rho, Sigma) %*% t(A),
               grad %*% Sigma %*% t(grad), tolerance = 1e-12)
})

test_that("the naive benchmark uses the equation (22) critical value", {
  set.seed(1); n <- 400
  x <- rnorm(n); u <- 0.5*x + rnorm(n); y <- 2*x + u; z <- 0.5*x + rnorm(n)
  res <- ci_naive_endpoint(x, y, z, rxu_range = c(0, 0.8), K = 80,
                           alpha = 0.05, B_boot = 100)
  expect_equal(res$crit, qnorm(0.975), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(res$crit, qnorm(1 - 0.05/(4*80)))))
})

test_that("the admissible inversion range is enforced", {
  A <- cr_Amat(seq(0, 0.8, length.out = 5))
  expect_error(
    pvalue_mcub_zero(c(0, 0), diag(2), A, A, cr_grid_map(A),
                     alpha_grid = c(0.001, 0.2)), "admissible range")
})

test_that("arguments are forwarded rather than discarded", {
  set.seed(2); n <- 300
  x <- rnorm(n); u <- 0.5*x + rnorm(n); y <- 2*x + u; z <- 0.5*x + rnorm(n)
  d <- data.frame(y = y, x = x, z = z)
  r1 <- iv_cr_test(d, X="x", Y="y", H=NULL, Z="z", k=1, rxu_range=c(0.4,0.6),
                   K=80, B_boot=100, compute_pvalue=FALSE, verbose=FALSE)
  expect_equal(r1$domain, c(0.4, 0.6)); expect_equal(r1$K, 80)
  r2 <- iv_cr_test(d, X="x", Y="y", H=NULL, Z="z", k=-1, K=80, B_boot=100,
                   compute_pvalue=FALSE, verbose=FALSE)
  expect_equal(r2$domain, c(-0.8, 0))
  expect_error(iv_cr_test(d, X="x", Y="y", H=NULL, Z="z", k=-1,
                          rxu_range=c(0,0.8), verbose=FALSE), "must lie in")
})

test_that("Card (1993) reproduces the published plug-in sets", {
  skip_on_cran(); skip_if_not_installed("wooldridge")
  card <- wooldridge::card
  H <- c("exper","expersq","black","south","smsa","smsa66", paste0("reg66",1:8))
  # Per-instrument samples. complete.cases() on the whole frame gives 1600.
  exp_n  <- c(fatheduc=2320, motheduc=2657, nearc2=3010, nearc4=3010)
  exp_lo <- c(fatheduc=0.020, motheduc=0.038, nearc2=0.033, nearc4=0.020)
  exp_hi <- c(fatheduc=0.239, motheduc=0.240, nearc2=0.043, nearc4=0.065)
  for (z in names(exp_n)) {
    res <- iv_cr_test(card, X="educ", Y="lwage", H=H, Z=z, k=1, alpha=0.05,
                      rxu_range=c(0,0.8), K=80, B_boot=200,
                      compute_pvalue=FALSE, verbose=FALSE)
    expect_equal(res$n, unname(exp_n[z]), info = z)
    lam <- res$diagnostics[[z]]$lambda_hat
    expect_equal(round(min(lam),3), unname(exp_lo[z]), info = z)
    expect_equal(round(max(lam),3), unname(exp_hi[z]), info = z)
  }
})
