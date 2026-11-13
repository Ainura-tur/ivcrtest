# Functions with a `seed` argument must not change the caller's random-number
# stream (CRAN policy: packages must not modify the global environment).

sim_xyz <- function(n = 200) {
  x <- rnorm(n); u <- 0.5 * x + rnorm(n); y <- 2 * x + u; z <- 0.5 * x + rnorm(n)
  data.frame(x = x, y = y, z = z)
}
g_lin <- function(rho, r) as.numeric(cr_Amat(r) %*% cr_delta(rho))

test_that("seeded functions leave the caller's random-number stream alone", {
  set.seed(1); d <- sim_xyz()
  set.seed(42); ref <- runif(3)

  set.seed(42)
  invisible(estimate_cov_corr_boot(d$x, d$y, d$z, B = 20))
  expect_identical(runif(3), ref)

  set.seed(42)
  invisible(im_stoye_inference(d, c(0, 0.4, 0.8), g_fun = g_lin, B = 20))
  expect_identical(runif(3), ref)

  set.seed(42)
  invisible(suppressWarnings(
    iv_cr_test(d, X = "x", Y = "y", H = NULL, Z = "z", K = 5, B_boot = 20,
               compute_pvalue = FALSE, verbose = FALSE,
               B = 20, Blarge = 50, B_im = 19)))
  expect_identical(runif(3), ref)
})

test_that("a seed still makes the draws reproducible", {
  set.seed(1); d <- sim_xyz()
  set.seed(7);  a <- estimate_cov_corr_boot(d$x, d$y, d$z, B = 20, seed = 5)
  set.seed(99); b <- estimate_cov_corr_boot(d$x, d$y, d$z, B = 20, seed = 5)
  expect_identical(a, b)
})

test_that("seed = NULL uses and advances the session's stream", {
  set.seed(1); d <- sim_xyz()
  set.seed(3); a <- estimate_cov_corr_boot(d$x, d$y, d$z, B = 20, seed = NULL)
  after_a <- runif(1)
  set.seed(3); b <- estimate_cov_corr_boot(d$x, d$y, d$z, B = 20, seed = NULL)
  expect_identical(a, b)
  set.seed(3); expect_false(identical(runif(1), after_a))
})

test_that("no random-number state is left behind when none existed", {
  set.seed(1); d <- sim_xyz()
  genv <- globalenv()
  old <- get(".Random.seed", envir = genv)
  on.exit(assign(".Random.seed", old, envir = genv), add = TRUE)
  rm(".Random.seed", envir = genv)
  invisible(estimate_cov_corr_boot(d$x, d$y, d$z, B = 20))
  expect_false(exists(".Random.seed", envir = genv, inherits = FALSE))
})
