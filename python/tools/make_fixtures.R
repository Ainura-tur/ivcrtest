#!/usr/bin/env Rscript
# Record reference results from the R package ivcrtest, so that the Python
# port can be tested against R.
#
# Usage, from the root of the repository:
#
#   Rscript python/tools/make_fixtures.R
#
# The script writes JSON files to python/tests/fixtures/. It needs the R
# package ivcrtest installed at the version in DESCRIPTION, plus jsonlite.
# If wooldridge is installed, it also records the Card (1993) application.
#
# Each file is one of three kinds:
#
#   exact      Deterministic results. Python must match them to about 1e-10.
#   replay     Results that depend on random draws. The script records the
#              draws that R used, so that Python can replay them: bootstrap
#              indices from sample.int() (1-based, one row per replication)
#              and normal draws from MASS::mvrnorm(), in call order.
#   tolerance  End-to-end results. Python, using its own random numbers, must
#              agree within Monte Carlo error and reach the same decisions.
#
# Numbers are written with 17 significant digits, so they read back exactly.
# Matrices are written as arrays of rows. Missing values are written as null.

suppressPackageStartupMessages({
  library(ivcrtest)
  library(jsonlite)
})

out_dir <- file.path("python", "tests", "fixtures")
if (!file.exists("DESCRIPTION") || !dir.exists(out_dir)) {
  stop("Run this script from the root of the ivcrtest repository.")
}
version <- read.dcf("DESCRIPTION", fields = "Version")[1, 1]
if (as.character(packageVersion("ivcrtest")) != version) {
  stop("The installed ivcrtest is version ", packageVersion("ivcrtest"),
       ", but DESCRIPTION says ", version, ". Install this version first.")
}
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

meta <- list(ivcrtest_version = version,
             r_version = paste(R.version$major, R.version$minor, sep = "."),
             generator = "python/tools/make_fixtures.R")

write_fixture <- function(name, kind, ...) {
  obj  <- list(meta = c(meta, list(kind = kind)), ...)
  path <- file.path(out_dir, paste0(name, ".json"))
  write_json(obj, path, digits = I(17), auto_unbox = TRUE, matrix = "rowmajor",
             dataframe = "columns", na = "null", null = "null")
  cat(sprintf("wrote %-42s %7.1f KB\n", path, file.size(path) / 1024))
}

# Keep numeric vectors as JSON arrays, even when they have length 1.
vec <- function(x) I(as.numeric(x))

# Data frames as lists of columns, each kept as an array.
cols <- function(d) lapply(as.list(d), I)

# ---------------------------------------------------------------------------
# Recording random draws
# ---------------------------------------------------------------------------
# record() evaluates `expr` and returns its value together with, in call
# order, every result of MASS::mvrnorm() and sample.int() and every call to
# optim() made from inside ivcrtest. The optim() record is for debugging: it
# shows where a port's optimiser first takes a different path.

.rec <- new.env()
.rec_push <- function(slot, value) {
  .rec[[slot]][[length(.rec[[slot]]) + 1L]] <- value
  invisible(NULL)
}

swap_binding <- function(env, name, value) {
  unlockBinding(name, env)
  assign(name, value, envir = env)
  lockBinding(name, env)
}

record <- function(expr) {
  .rec$mvrnorm <- list()
  .rec$sample_int <- list()
  .rec$optim <- list()

  mass_ns  <- asNamespace("MASS")
  imports  <- parent.env(asNamespace("ivcrtest"))
  mvrnorm0 <- get("mvrnorm", envir = mass_ns)
  optim0   <- get("optim", envir = imports)

  swap_binding(mass_ns, "mvrnorm", function(n = 1, mu, Sigma, ...) {
    draws <- mvrnorm0(n = n, mu = mu, Sigma = Sigma, ...)
    .rec_push("mvrnorm", list(n = n, draws = as.matrix(draws)))
    draws
  })
  swap_binding(imports, "optim", function(par, fn, ...) {
    args <- list(...)
    fit  <- optim0(par = par, fn = fn, ...)
    .rec_push("optim", list(start = vec(par), lower = vec(args$lower),
                            upper = vec(args$upper), par = vec(fit$par),
                            value = fit$value, counts = vec(fit$counts),
                            convergence = fit$convergence))
    fit
  })
  suppressMessages(trace("sample.int", where = baseenv(), print = FALSE,
                         exit = quote(.rec_push("sample_int", returnValue()))))
  on.exit({
    swap_binding(mass_ns, "mvrnorm", mvrnorm0)
    swap_binding(imports, "optim", optim0)
    suppressMessages(untrace("sample.int", where = baseenv()))
  }, add = TRUE)

  value <- force(expr)
  list(value = value, mvrnorm = .rec$mvrnorm, sample_int = .rec$sample_int,
       optim = .rec$optim)
}

# Bootstrap indices as a matrix: one row per replication.
index_matrix <- function(rec) do.call(rbind, rec$sample_int)

# ---------------------------------------------------------------------------
# 1. algebra.json (exact): reparameterisation, grid maps, diagnostics
# ---------------------------------------------------------------------------

Sigma_fixed <- matrix(c(4, 1, 1, 1, 9, 2, 1, 2, 6), 3, 3) / 1e4
K_alg <- 20L

algebra_case <- function(name, rho, D) {
  r     <- seq(D[1], D[2], length.out = K_alg)
  A     <- cr_Amat(r)
  delta <- cr_delta(rho)
  delta2 <- rbind(delta, delta + c(0.01, -0.02))
  rho2   <- rbind(rho, rho + c(0.01, -0.01, 0.02))
  list(name = name, rho = vec(rho), D = vec(D), K = K_alg, r = vec(r),
       Sigma_rho = Sigma_fixed,
       h = cr_h(rho), delta = vec(delta), jacobian = cr_jacobian(rho),
       A = A, Omega = cr_omega_delta(rho, Sigma_fixed),
       lambda = vec(A %*% delta),
       grid_map_input = delta2, grid_map = cr_grid_map(A)(delta2),
       g_fun = vec(g_fun(r, rho[1], rho[2], rho[3])),
       g_grad_r = r[5], g_grad = vec(g_grad(r[5], rho[1], rho[2], rho[3])),
       g_xu_safe_input = rho2, g_xu_safe = g_xu_safe(r, rho2),
       gradient_safe = local_compute_gradient_safe(r, rho[1], rho[2], rho[3]),
       sign_condition = cr_sign_condition(rho, D),
       endpoint_gap = cr_endpoint_gap(rho, Sigma_fixed, D, K = K_alg),
       discretization = as.list(cr_discretization(rho, D, K = K_alg)),
       corr_screen = cr_corr_screen(rho, Sigma_fixed, D, K = K_alg, tol_r = 1e-3))
}

c_im <- expand.grid(alpha = c(0.01, 0.05, 0.10), t = c(0, 0.25, 1, 2, 5, 10))
c_im$c <- mapply(solve_c_im, c_im$alpha, c_im$t)

write_fixture(
  "algebra", "exact",
  cases = list(
    algebra_case("sign_condition_fails", c(0.5, 0.4, 0.3), c(0, 0.8)),
    algebra_case("sign_condition_holds", c(0.5, 0.4, 0.1), c(0, 0.8)),
    algebra_case("negative_domain", c(-0.3, 0.45, 0.2), c(-0.8, 0))
  ),
  solve_c_im = cols(c_im)
)

# ---------------------------------------------------------------------------
# 2. nearpd.json (exact): Matrix::nearPD() as estimate_cov_corr_boot() uses it
# ---------------------------------------------------------------------------

nearpd_case <- function(name, S) {
  fit <- Matrix::nearPD(S, corr = FALSE)
  list(name = name, input = S, output = as.matrix(fit$mat),
       iterations = fit$iterations, converged = fit$converged)
}
v <- c(1, 2, 3) / 100

write_fixture(
  "nearpd", "exact",
  settings = list(corr = FALSE, keepDiag = FALSE, do2eigen = TRUE,
                  doDykstra = TRUE, eig_tol = 1e-6, conv_tol = 1e-7,
                  posd_tol = 1e-8, maxit = 100L, conv_norm_type = "I"),
  cases = list(
    nearpd_case("positive_definite", Sigma_fixed),
    nearpd_case("indefinite",
                matrix(c(1, 0.9, 0.3, 0.9, 1, 0.9, 0.3, 0.9, 1), 3, 3) / 1e4),
    nearpd_case("near_singular", outer(v, v) + diag(1e-12, 3))
  )
)

# ---------------------------------------------------------------------------
# 3. data.json (exact): simulated data, residualisation on the controls,
#    sample correlations, and the pointwise range for a given Sigma_rho
# ---------------------------------------------------------------------------

set.seed(20261114)
n_sim <- 300L
w1 <- rnorm(n_sim); w2 <- rbinom(n_sim, 1, 0.4)
z1 <- rnorm(n_sim); v_ <- rnorm(n_sim); u <- 0.5 * v_ + rnorm(n_sim)
x  <- 0.6 * z1 + 0.3 * w1 - 0.2 * w2 + v_
y  <- 1 + 0.5 * x + 0.2 * w1 + 0.1 * w2 + u
z2 <- 0.4 * z1 + 0.3 * u + rnorm(n_sim)   # correlated with u: invalid
z3 <- -0.5 * z1 + rnorm(n_sim)            # negatively related to x: flipped
z2[c(5L, 17L, 101L)] <- NA                # exercise complete cases
w1[42L] <- NA
sim <- data.frame(y = y, x = x, z1 = z1, z2 = z2, z3 = z3, w1 = w1, w2 = w2)
H_sim <- c("w1", "w2")

# Mirrors the residualisation inside iv_cr_test().
residualise <- function(zname) {
  vj  <- unique(c("y", "x", H_sim, zname))
  sel <- complete.cases(sim[, vj])
  dat <- sim[sel, vj]
  rhs <- paste(H_sim, collapse = " + ")
  res <- function(v) unname(resid(lm(as.formula(paste0("`", v, "` ~ ", rhs)), data = dat)))
  xr <- res("x"); yr <- res("y"); zr <- res(zname)
  flipped <- cor(xr, zr) < 0
  if (flipped) zr <- -zr
  list(instrument = zname, rows = I(which(sel)), n = sum(sel),
       x = vec(xr), y = vec(yr), z = vec(zr), flipped = flipped,
       rho = vec(c(cor(xr, yr), cor(xr, zr), cor(yr, zr))))
}
resid_cases <- lapply(c("z1", "z2", "z3"), residualise)
names(resid_cases) <- c("z1", "z2", "z3")

# Check against iv_cr_test() itself, so that this script cannot drift from it.
chk <- suppressWarnings(iv_cr_test(sim, X = "x", Y = "y", H = H_sim,
                                   Z = c("z1", "z2", "z3"), K = 10L, B_boot = 50L,
                                   compute_pvalue = FALSE, verbose = FALSE,
                                   B = 50L, Blarge = 200L, B_im = 49L))
for (z in names(resid_cases)) {
  stopifnot(isTRUE(all.equal(unname(chk$diagnostics[[z]]$rho_hat),
                             as.numeric(resid_cases[[z]]$rho), tolerance = 1e-14)))
}

d1  <- resid_cases$z1
df1 <- data.frame(x = as.numeric(d1$x), y = as.numeric(d1$y), z = as.numeric(d1$z))

write_fixture(
  "data", "exact",
  data = cols(sim),
  controls = I(H_sim),
  residualised = resid_cases,
  estimate_rho = vec(estimate_rho(df1)),
  pointwise_range = list(
    instrument = "z1", rxu_range = vec(c(0, 0.8)), k = 1, alpha = 0.05,
    n_grid = 401L, Sigma_rho = Sigma_fixed,
    result = cr_pointwise_range(df1, rxu_range = c(0, 0.8), k = 1, alpha = 0.05,
                                n_grid = 401L, Sigma_rho = Sigma_fixed)
  )
)

# ---------------------------------------------------------------------------
# 4. bootstrap.json (replay): the covariance and IM/Stoye bootstraps
# ---------------------------------------------------------------------------

B_small <- 40L
rec_cov <- record(estimate_cov_corr_boot(df1$x, df1$y, df1$z, B = B_small, seed = 7L))
idx_cov <- index_matrix(rec_cov)
stopifnot(nrow(idx_cov) == B_small, ncol(idx_cov) == nrow(df1))
Rb <- t(apply(idx_cov, 1, function(i) {
  c(cor(df1$x[i], df1$y[i]), cor(df1$x[i], df1$z[i]), cor(df1$y[i], df1$z[i]))
}))
S_raw <- cov(Rb)
S_raw <- (S_raw + t(S_raw)) / 2
stopifnot(isTRUE(all.equal(as.matrix(Matrix::nearPD(S_raw, corr = FALSE)$mat),
                           rec_cov$value, tolerance = 1e-14)))

S_default <- estimate_cov_corr_boot(df1$x, df1$y, df1$z)

grid20 <- seq(0, 0.8, length.out = 20L)
g_lin  <- function(rho, r) as.numeric(cr_Amat(r) %*% cr_delta(rho))
rec_im <- record(im_stoye_inference(df1, grid20, g_fun = g_lin, B = 30L,
                                    seed = 5L, return_boot = TRUE))
im_delta <- im_stoye_inference(df1, grid20, g_fun = g_lin, se_method = "delta",
                               grad_fun = function(rho, r) {
                                 as.numeric(cr_Amat(r) %*% cr_jacobian(rho))
                               },
                               cov_rho = Sigma_fixed)

write_fixture(
  "bootstrap", "replay",
  note = paste("Data are data.json residualised$z1. Indices are 1-based.",
               "cov_default and im_delta need no replay."),
  cov_small = list(B = B_small, seed = 7L, indices = idx_cov, Rb = Rb,
                   S_raw = S_raw, S = rec_cov$value),
  cov_default = list(B = 800L, seed = 123L, S = S_default,
                     kind = "tolerance"),
  im_bootstrap = list(grid = vec(grid20), alpha = 0.05, B = 30L, seed = 5L,
                      indices = index_matrix(rec_im), result = rec_im$value),
  im_delta = list(grid = vec(grid20), alpha = 0.05, cov_rho = Sigma_fixed,
                  result = im_delta, kind = "exact")
)

# ---------------------------------------------------------------------------
# 5. union_bounds.json (exact): CIcon_TNbounds(), CIcon() and CIproj_p() for
#    given inputs, on the reparameterised and the primitive route
# ---------------------------------------------------------------------------

rho1 <- as.numeric(d1$rho)

# Mirrors steps 2 and 3 of CIhybrid().
union_setup <- function(route, K = 20L) {
  r_grid <- seq(0, 0.8, length.out = K)
  if (route == "reparam") {
    A    <- cr_Amat(r_grid)
    dhat <- as.numeric(cr_delta(rho1))
    dSig <- cr_omega_delta(rho1, S_default)
    Al   <- A
    g    <- cr_grid_map(A)
  } else {
    dhat <- rho1
    dSig <- S_default
    Al   <- t(vapply(r_grid, function(r) {
      as.numeric(local_compute_gradient_safe(r, rho1[1], rho1[2], rho1[3]))
    }, numeric(3)))
    g    <- function(delta) g_xu_safe(r_grid, delta)
  }
  Lambda   <- rbind(Al, Al)
  LS       <- Lambda %*% dSig %*% t(Lambda)
  ls       <- sqrt(diag(LS))
  corr_all <- diag(1 / ls) %*% LS %*% diag(1 / ls)
  list(route = route, r_grid = r_grid, deltahat = dhat, deltaSigma = dSig,
       Al = Al, Au = Al, g = g,
       sigma_l = ls[1:K], sigma_u = ls[(K + 1):(2 * K)],
       corr_l = corr_all[1:K, 1:K], corr_u = corr_all[(K + 1):(2 * K), (K + 1):(2 * K)],
       corr_m = corr_all[1:K, (K + 1):(2 * K)])
}

union_case <- function(route) {
  s   <- union_setup(route)
  k   <- length(s$deltahat)
  eta <- 0.001
  set.seed(11L)
  draws <- MASS::mvrnorm(200L, rep(0, k), s$deltaSigma)
  dev   <- apply(sweep(abs(draws), 2, sqrt(diag(s$deltaSigma)), "/"), 1, max)
  c_bd  <- unname(quantile(dev, 1 - eta))
  lam   <- as.numeric(s$g(s$deltahat))
  thetas <- c(min(lam) - 0.05, min(lam), 0, mean(range(lam)), max(lam))
  rows   <- rbind(s$deltahat, sweep(draws[1:5, , drop = FALSE], 2, s$deltahat, "+"))

  tn <- lapply(thetas, function(th) {
    b <- CIcon_TNbounds(th, rows, s$Al, s$Au, s$sigma_l, s$sigma_u,
                        s$corr_m, s$corr_l, s$corr_u,
                        qnorm(1 - eta), qnorm(1 - eta), 1e-3, s$g)
    list(theta = th, th_1 = vec(b$th_1), th_2 = vec(b$th_2))
  })
  offsets <- rbind(rep(0, k), s$deltahat * 0 + 0.01, s$deltahat * 0 - 0.01)
  pp <- list()
  for (cc in c(0.5, 1, 1.5, 2)) for (j in seq_len(nrow(offsets))) {
    dl <- s$deltahat + offsets[j, ]
    pp[[length(pp) + 1L]] <- list(
      c = cc, delta = vec(dl),
      p = CIproj_p(cc, dl, 0.04, c_bd, draws, s$Al, s$Au, s$sigma_l, s$sigma_u,
                   s$deltaSigma, s$corr_m, s$corr_l, s$corr_u, eta, 1e-3, s$g))
  }
  list(route = route, r_grid = vec(s$r_grid), deltahat = vec(s$deltahat),
       deltaSigma = s$deltaSigma, Al = s$Al, Au = s$Au,
       sigma_l = vec(s$sigma_l), sigma_u = vec(s$sigma_u),
       corr_l = s$corr_l, corr_u = s$corr_u, corr_m = s$corr_m,
       alphac = 0.04, eta = eta, tol_r = 1e-3,
       draws = draws, c_bd = c_bd,
       tnbounds_rows = rows, tnbounds = tn,
       cicon = list(c_LF = qnorm(1 - eta), tol = 1e-3,
                    result = vec(CIcon(s$deltahat, s$deltaSigma, s$Al, s$Au,
                                       qnorm(1 - eta), 0.04, 1e-3, 1e-3, s$g))),
       ciproj_p = pp)
}

write_fixture(
  "union_bounds", "exact",
  note = paste("Inputs come from data.json residualised$z1 and bootstrap.json",
               "cov_default$S. The map g is cr_grid_map(Al) on the reparam",
               "route and g_xu_safe(r_grid, delta) on the primitive route.",
               "c_bd is computed column by column, as in Bei's code."),
  cases = list(union_case("reparam"), union_case("primitive"))
)

# ---------------------------------------------------------------------------
# 6. cihybrid.json (replay): CIhybrid() with its normal draws recorded
# ---------------------------------------------------------------------------

cihybrid_case <- function(route) {
  s <- union_setup(route)
  settings <- list(alpha = 0.05, alphac = 0.04, eta = 0.001, B = 100L,
                   Blarge = 1000L, tol = 1e-3, tol_r = 1e-3, seed = 1L)
  rec <- record(CIhybrid(s$deltahat, s$deltaSigma, s$Al, s$Au,
                         alpha = settings$alpha, alphac = settings$alphac,
                         eta = settings$eta, B = settings$B, Blarge = settings$Blarge,
                         tol = settings$tol, tol_r = settings$tol_r, index = NULL,
                         g = s$g, seed = settings$seed))
  list(route = route, r_grid = vec(s$r_grid), deltahat = vec(s$deltahat),
       deltaSigma = s$deltaSigma, Al = s$Al, Au = s$Au, settings = settings,
       draws = lapply(rec$mvrnorm, function(m) m$draws),
       optim_calls = rec$optim,
       result = list(CI_h = vec(rec$value$CI_h), CI_c = vec(rec$value$CI_c),
                     CI_p = vec(rec$value$CI_p), delta1 = vec(rec$value$delta1),
                     c = rec$value$c))
}

write_fixture(
  "cihybrid", "replay",
  note = paste("draws[[1]] is the Blarge draw made once at the start; the rest",
               "are the B draws made at each bisection step, in order.",
               "R 0.2.2 standardises the large draws row-wise when it sets",
               "c_bd (CIhybrid.R line 80); Bei's code does it column-wise.",
               "Results may also differ slightly because R's L-BFGS-B differs",
               "from SciPy's; optim_calls shows R's path."),
  cases = list(cihybrid_case("reparam"), cihybrid_case("primitive"))
)

# ---------------------------------------------------------------------------
# 7. end_to_end.json (tolerance): full runs
# ---------------------------------------------------------------------------

e2e_run <- function(data, X, Y, H, Z, ...) {
  res <- suppressWarnings(iv_cr_test(data, X = X, Y = Y, H = H, Z = Z,
                                     verbose = FALSE, ...))
  per <- lapply(Z, function(z) {
    dg <- res$diagnostics[[z]]
    list(instrument = z, rho_hat = vec(dg$rho_hat), Sigma_rho = dg$Sigma_rho,
         mcub = vec(dg$mcub$CI_h), mcub_c = dg$mcub$c,
         lambda_hat = vec(dg$lambda_hat), sigma_b = vec(dg$sigma_b))
  })
  names(per) <- Z
  list(settings = list(...), results = cols(res$results), instruments = per)
}

sim_e2e <- e2e_run(sim, "x", "y", H_sim, c("z1", "z2", "z3"),
                   K = 20L, B_boot = 200L, B = 100L, Blarge = 1000L,
                   B_im = 199L, compute_pvalue = TRUE)
sim_e2e_prim <- e2e_run(sim, "x", "y", H_sim, "z1", route = "primitive",
                        K = 20L, B_boot = 200L, B = 100L, Blarge = 1000L,
                        B_im = 199L, compute_pvalue = FALSE)

small <- list(K = 10L, B = 50L, Blarge = 500L, tol = 0.01)
sens <- list(
  valid_range = cr_valid_range(df1, rho_bar = 0.8, cut_grid = c(0.2, 0.4, 0.6),
                               k = 1, K = small$K, Sigma_rho = S_default,
                               B = small$B, Blarge = small$Blarge,
                               tol = small$tol, verbose = FALSE),
  switch_point = cr_switch_point(df1, rho_bar = 0.8, k = 1, direction = "below",
                                 scan_grid = c(0.2, 0.4, 0.6), steps = 2L,
                                 K = small$K, Sigma_rho = S_default,
                                 B = small$B, Blarge = small$Blarge,
                                 tol = small$tol),
  pvalue = {
    A10 <- cr_Amat(seq(0, 0.8, length.out = small$K))
    pvalue_mcub_zero(cr_delta(rho1), cr_omega_delta(rho1, S_default), A10, A10,
                     cr_grid_map(A10), B = small$B, Blarge = small$Blarge,
                     tol = small$tol, alpha_grid = c(0.02, 0.05, 0.1, 0.2, 0.4),
                     refine_steps = 3L)
  }
)
sens$valid_range$table <- cols(sens$valid_range$table)

card <- NULL
if (requireNamespace("wooldridge", quietly = TRUE)) {
  H_card <- c("exper", "expersq", "black", "south", "smsa", "smsa66",
              paste0("reg66", 1:8))
  card <- e2e_run(wooldridge::card, "educ", "lwage", H_card,
                  c("fatheduc", "motheduc", "nearc2", "nearc4"),
                  k = 1, rxu_range = c(0, 0.8), compute_pvalue = FALSE)
  card$data <- list(source = "wooldridge::card",
                    wooldridge_version = as.character(packageVersion("wooldridge")),
                    controls = I(H_card))
}

write_fixture(
  "end_to_end", "tolerance",
  note = paste("Interval columns in results are rounded text, as iv_cr_test()",
               "prints them; instruments holds the unrounded MCUB intervals.",
               "Sensitivity runs use data.json residualised$z1 with",
               "Sigma_rho = bootstrap.json cov_default$S."),
  sim = sim_e2e,
  sim_primitive = sim_e2e_prim,
  sensitivity = c(list(settings = small), sens),
  card = card
)

cat("done\n")
