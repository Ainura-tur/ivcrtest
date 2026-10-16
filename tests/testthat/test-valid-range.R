# Regression tests for the valid range of Section 3.5.
#
# ivcrtest 0.2.0 took the lower endpoint from the "restrict from below" family
# and the upper endpoint from the "restrict from above" family. Because
# acceptance is nested, that returned the outer cut points whenever anything
# accepted. The tables below are the accept/reject tables observed on real
# data; the endpoint rule is checked on them without any MCUB fit.

cuts <- seq(0.1, 0.7, by = 0.1)

test_that("Thornton (2008), any incentive: endpoints follow Section 3.5", {
  # Observed: (rho_k, 0.8) accepts only at 0.1; (0, rho_k) accepts throughout.
  below <- c(TRUE, rep(FALSE, 6))
  above <- rep(TRUE, 7)
  r <- ivcrtest:::.cr_range_from_table(cuts, below, above, full = TRUE,
                                       rho_bar = 0.8, edge_lo = TRUE)
  expect_equal(r$status, "range")
  expect_equal(r$upper, 0.1)
  expect_equal(r$upper_bracket, c(0.1, 0.2))  # bisection gives 0.182
  expect_equal(r$lower, 0)
  # The 0.2.0 rule returned c(0.1, 0.7) here.
  expect_false(isTRUE(all.equal(c(r$lower, r$upper), c(0.1, 0.7))))
})

test_that("Thornton (2008), incentive amount: switch between 0.2 and 0.3", {
  below <- c(TRUE, TRUE, rep(FALSE, 5))
  r <- ivcrtest:::.cr_range_from_table(cuts, below, rep(TRUE, 7), TRUE, 0.8,
                                       edge_lo = TRUE)
  expect_equal(r$upper, 0.2)
  expect_equal(r$upper_bracket, c(0.2, 0.3))  # bisection gives 0.218
})

test_that("an instrument that survives everywhere spans the whole domain", {
  # NHEFS price instruments: every fit accepts.
  r <- ivcrtest:::.cr_range_from_table(cuts, rep(TRUE, 7), rep(TRUE, 7), TRUE,
                                       0.8, edge_lo = TRUE, edge_hi = TRUE)
  expect_equal(c(r$lower, r$upper), c(0, 0.8))
  # Without the edge fits the endpoints are only bracketed, never overstated.
  r2 <- ivcrtest:::.cr_range_from_table(cuts, rep(TRUE, 7), rep(TRUE, 7), TRUE,
                                        0.8, edge_lo = FALSE, edge_hi = FALSE)
  expect_equal(c(r2$lower, r2$upper), c(0.1, 0.7))
  expect_equal(r2$lower_bracket, c(0.001, 0.1))
  expect_equal(r2$upper_bracket, c(0.7, 0.799))
})

test_that("a range away from zero is reported from both sides", {
  # Accepts only when the domain covers roughly 0.3 to 0.5.
  below <- cuts <= 0.5
  above <- cuts >= 0.3
  r <- ivcrtest:::.cr_range_from_table(cuts, below, above, TRUE, 0.8)
  expect_equal(c(r$lower, r$upper), c(0.3, 0.5))
  expect_equal(r$lower_bracket, c(0.2, 0.3))
  expect_equal(r$upper_bracket, c(0.5, 0.6))
})

test_that("rejection on the full domain is reported, not hidden", {
  r <- ivcrtest:::.cr_range_from_table(cuts, rep(FALSE, 7), rep(FALSE, 7),
                                       full = FALSE, rho_bar = 0.8)
  expect_equal(r$status, "rejects throughout")
  expect_true(all(is.na(c(r$lower, r$upper))))
})

test_that("a range narrower than one grid cell is flagged", {
  # Full domain accepts, but no cut accepts from below and the first accepting
  # cut from above is 0.1: both endpoints lie in (0, 0.1).
  r <- ivcrtest:::.cr_range_from_table(cuts, rep(FALSE, 7), rep(TRUE, 7), TRUE,
                                       0.8, edge_lo = FALSE)
  expect_equal(r$status, "below grid resolution")
})

test_that("a failed fit stops rather than counting as rejection", {
  expect_error(ivcrtest:::.cr_range_from_table(cuts, c(NA, rep(TRUE, 6)),
                                               rep(TRUE, 7), TRUE, 0.8),
               "MCUB fit failed")
})

test_that("non-monotone tables warn", {
  expect_warning(ivcrtest:::.cr_range_from_table(
    cuts, c(TRUE, FALSE, TRUE, rep(FALSE, 4)), rep(TRUE, 7), TRUE, 0.8,
    edge_lo = TRUE), "not monotone")
})

test_that("grid endpoints agree with cr_switch_point on the same inputs", {
  skip_on_cran()
  set.seed(7); n <- 1000
  e <- rnorm(n); z <- rnorm(n)
  x <- 0.5 * z + e
  y <- 0.5 * x + 0.3 * e + rnorm(n)            # valid instrument, rho_xu > 0
  H <- function(v) v - mean(v)
  df <- data.frame(x = H(x), y = H(y), z = H(z))
  S  <- estimate_cov_corr_boot(df$x, df$y, df$z, B = 300, seed = 11)
  vr <- cr_valid_range(df, cut_grid = c(0.1, 0.3, 0.5, 0.7), k = 1,
                       Sigma_rho = S, verbose = FALSE)
  expect_equal(vr$status, "range")
  if (vr$upper_bracket[1] < vr$upper_bracket[2] &&
      vr$upper_bracket[2] < 0.8 - 1e-3) {
    sp <- cr_switch_point(df, k = 1, direction = "below",
                          bracket = vr$upper_bracket, Sigma_rho = S, steps = 3L)
    expect_equal(sp$status, "interior")
    expect_true(sp$tau >= vr$upper_bracket[1] && sp$tau <= vr$upper_bracket[2])
  }
  # Mirrored sign gives a mirrored domain.
  vr_neg <- cr_valid_range(data.frame(x = df$x, y = df$y, z = df$z),
                           cut_grid = c(0.3), k = -1, Sigma_rho = S,
                           verbose = FALSE)
  expect_true(all(vr_neg$table$rho_k < 0))
})

# ---- boot_seed: explicit control of the bootstrap draws ----------------------

sim_xyz <- function(n = 300, s = 3) {
  set.seed(s)
  x <- rnorm(n); u <- 0.5 * x + rnorm(n); y <- 2 * x + u; z <- 0.5 * x + rnorm(n)
  data.frame(x = x, y = y, z = z)
}
fit_cc <- function(df, ...) suppressWarnings(
  check_compatibility(df, k = 1, rxu_range = c(0, 0.8), B_boot = 100,
                      B = 100, Blarge = 1000, B_im = 99,
                      compute_pvalue = FALSE, verbose = FALSE, ...))

test_that("the default boot_seed reproduces version 0.2.0 (1000 + i)", {
  df <- sim_xyz()
  r  <- fit_cc(df, i = 2L)
  S0 <- estimate_cov_corr_boot(df$x, df$y, df$z, B = 100, seed = 1002L)
  expect_equal(attr(r, "diagnostics")$Sigma_rho, S0)
})

test_that("a supplied boot_seed sets the bootstrap, and seed does not", {
  df <- sim_xyz()
  r5 <- fit_cc(df, i = 1L, boot_seed = 5L)
  expect_equal(attr(r5, "diagnostics")$Sigma_rho,
               estimate_cov_corr_boot(df$x, df$y, df$z, B = 100, seed = 5L))
  ra <- fit_cc(df, i = 1L, seed = 1L)
  rb <- fit_cc(df, i = 1L, seed = 99L)
  expect_equal(attr(ra, "diagnostics")$Sigma_rho, attr(rb, "diagnostics")$Sigma_rho)
  expect_error(fit_cc(df, boot_seed = c(1, 2)), "boot_seed")
})

test_that("iv_cr_test gives instrument j the bootstrap seed boot_seed + j", {
  df <- sim_xyz()
  df$z2 <- df$z + rnorm(nrow(df))
  r <- suppressWarnings(iv_cr_test(df, X = "x", Y = "y", H = NULL,
                                   Z = c("z", "z2"), k = 1, B_boot = 100,
                                   boot_seed = 40L, compute_pvalue = FALSE,
                                   verbose = FALSE, B = 100, Blarge = 1000,
                                   B_im = 99))
  d2 <- data.frame(x = df$x - mean(df$x), y = df$y - mean(df$y),
                   z = df$z2 - mean(df$z2))
  if (r$results$sign_flipped[2]) d2$z <- -d2$z
  expect_equal(r$diagnostics[["z2"]]$Sigma_rho,
               estimate_cov_corr_boot(d2$x, d2$y, d2$z, B = 100, seed = 42L),
               tolerance = 1e-10)
})

# ---- check marks: logical columns and a plain-text display -------------------

test_that("logical zero-inclusion columns match the marks", {
  r <- fit_cc(sim_xyz(), i = 1L)
  expect_type(r$zero_in_mcub, "logical")
  expect_identical(r$Zero_in_CI_MCUB == "✓", r$zero_in_mcub)
  expect_identical(r$Zero_in_naive == "✓", r$zero_in_naive)
})

test_that("marks print as words when the locale is not UTF-8", {
  tab <- data.frame(Z = c("a", "b", "c"),
                    Zero_in_CI_MCUB = c("✓", "×", NA),
                    stringsAsFactors = FALSE)
  plain <- ivcrtest:::.cr_marks_for_display(tab, utf8 = FALSE)
  expect_identical(plain$Zero_in_CI_MCUB, c("yes", "no", NA))
  expect_identical(ivcrtest:::.cr_marks_for_display(tab, utf8 = TRUE), tab)
})

test_that("stored marks are UTF-8 whatever the install locale", {
  r <- fit_cc(sim_xyz(), i = 1L)
  expect_identical(Encoding(r$Zero_in_CI_MCUB), "UTF-8")
  expect_identical(utf8ToInt(r$Zero_in_CI_MCUB),
                   if (isTRUE(r$zero_in_mcub)) 0x2713L else 0x00D7L)
})
