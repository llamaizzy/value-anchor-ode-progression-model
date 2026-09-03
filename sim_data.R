## =============================================================================
##  GENERATE OBSERVATIONS FROM GROUND TRUTH
##
##  CORRECTED VERSION -- see FIXES.md, issue #2.
##
##  The original sim_data.R generated data by calling the SAME getTraj()
##  nimbleFunction the fitting model uses, at the SAME maxSub = 0.25 and the
##  SAME nGrid = 201. Every discretisation error was therefore identical in
##  generation and in fitting, and cancelled exactly. That makes the study
##  structurally incapable of detecting integrator bias: it would report clean
##  recovery of a discretisation scheme that might be badly wrong on real data.
##
##  This version shares NOTHING with the estimator:
##
##                        generator (here)          fitter (model.R)
##    rate curve          analytic function         K=10 B-spline, quantile knots
##    representation      closed form               linear interp on 701-pt grid
##    time integration    EXACT (see below)         Crank-Nicolson, 0.25 yr steps
##
##  THE EXACT SOLVE. The governing ODE is autonomous and separable:
##
##      dmu/dt = R(mu) * exp(delta)   =>   dt = dmu / (R(mu) * exp(delta))
##
##  so the time to travel from value a to value b is a ONE-DIMENSIONAL
##  INTEGRAL, with no time-stepping at all:
##
##      t(b) - t(a) = exp(-delta) * Integral_a^b dy / R(y)
##
##  Define G(y) = Integral_{y_ref}^{y} dy'/R(y') once, on a fine value grid.
##  G is strictly increasing (R > 0), hence invertible, and the trajectory is
##  available in closed form:
##
##      mu(t) = G^{-1}( G(x0) + (t - t0) * exp(delta) )
##
##  This is both far more accurate than any stepping scheme AND completely
##  structurally different from what the estimator does -- which is exactly
##  the property the original version lacked. It also makes the true age of
##  positivity exact rather than step-interpolated, and it is vectorised over
##  all subjects at once, so it costs milliseconds for N = 1000.
## =============================================================================


## -----------------------------------------------------------------------
##  1. Build the exact solver for a given truth
## -----------------------------------------------------------------------
##
##  The solver grid deliberately EXTENDS BEYOND the fitting domain [YL, YU].
##  A fast subject (large delta) integrated forward over a long follow-up can
##  legitimately run past YU; if the generator clamped there, it would inject a
##  flat, unphysical ceiling into the "truth" and then the estimator would be
##  scored against an artefact. Extending the generator's domain keeps the
##  truth honest; whether the resulting observations still fall inside the
##  FITTING domain is a separate question, measured and reported by
##  simulate_amyloid_data() below rather than assumed.
build_truth_solver <- function(truth, y_lo = 0.15, y_hi = 3.00,
                               step_y_solve = 0.0002) {
  ygrid <- seq(y_lo, y_hi, by = step_y_solve)
  R <- truth$rate_fun(ygrid)
  stopifnot("rate must be strictly positive on the solver grid" = all(R > 0))

  ## G(y) = cumulative Integral dy/R(y), by trapezoid on a fine grid.
  inv_R <- 1 / R
  G <- c(0, cumsum((inv_R[-1] + inv_R[-length(inv_R)]) / 2 * diff(ygrid)))

  ## G is strictly increasing, so both directions are monotone interpolations.
  y_to_G <- function(y) approx(ygrid, G, xout = y, rule = 2)$y
  G_to_y <- function(g) approx(G, ygrid, xout = g, rule = 2)$y

  list(ygrid = ygrid, G = G, y_to_G = y_to_G, G_to_y = G_to_y,
       y_lo = y_lo, y_hi = y_hi, step_y_solve = step_y_solve,
       G_lo = G[1], G_hi = G[length(G)])
}


## mu at an arbitrary time, given the anchor (x0 at t0) and delta.
## Vectorised over subjects. `t` may be before or after t0.
truth_mu_at <- function(solver, x0, t0, t, delta) {
  solver$G_to_y(solver$y_to_G(x0) + (t - t0) * exp(delta))
}


## Exact age at which the trajectory crosses `thres`.
## Vectorised.
##
## max_years MATTERS AND IS NOT COSMETIC. Because R(y) > 0 everywhere, the
## trajectory crosses ANY threshold eventually, so the bare separable solve
## always returns a finite answer -- even when the crossing is eight hundred
## years away. The estimator's positivity_age_vec() searches only +/- 60 years
## and returns NA beyond that, so without the same window the "true" and
## "estimated" positivity ages would not be the same quantity and any coverage
## comparison between them would be meaningless.
##
## Concretely: under the ORIGINAL truth (rates ~0.001 SUVR/yr) an unwindowed
## solve reports that 100% of subjects have a positivity age, when in fact only
## ~23% cross within any humanly plausible span. The window is what makes that
## visible.
truth_positivity_age <- function(solver, x0, t0, delta, thres, max_years = 60) {
  g_thr <- solver$y_to_G(thres)
  g_x0  <- solver$y_to_G(x0)
  age <- t0 + (g_thr - g_x0) * exp(-delta)
  ## Not a real crossing if the threshold lies outside the solver's domain...
  age[thres < solver$y_lo | thres > solver$y_hi] <- NA_real_
  ## ...or if it lies outside the window the estimator can search.
  age[abs(age - t0) > max_years] <- NA_real_
  age
}


## -----------------------------------------------------------------------
##  2. Generate observations
## -----------------------------------------------------------------------
##
##  delta_i ~ N(0, sigma_delta), then the latent trajectory is evaluated
##  EXACTLY at every visit age, then i.i.d. Gaussian noise is added to every
##  visit -- including the first. (The true anchor x0 is never fed back as an
##  observation; the estimator only ever sees noisy values, exactly as with
##  real data.)
simulate_amyloid_data <- function(truth, design, seed = NULL, solver = NULL,
                                  fit_YL = 0.30, fit_YU = 1.70, verbose = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(solver)) solver <- build_truth_solver(truth)

  N    <- design$N
  Jmax <- design$Jmax
  sigma_delta <- design$sigma_delta_true
  sigma_eps   <- design$sigma_eps_true

  delta_true <- rnorm(N, 0, sigma_delta)

  mu_true <- matrix(NA_real_, N, Jmax)
  for (j in seq_len(Jmax)) {
    has_j <- design$J >= j
    if (!any(has_j)) next
    mu_true[has_j, j] <- truth_mu_at(
      solver,
      x0    = design$x0_true[has_j],
      t0    = design$t0_true[has_j],
      t     = design$tvisit[has_j, j],
      delta = delta_true[has_j]
    )
  }

  y <- mu_true + matrix(rnorm(N * Jmax, 0, sigma_eps), N, Jmax)
  y[is.na(mu_true)] <- NA_real_

  dat <- do.call(rbind, lapply(seq_len(N), function(i) {
    Ji <- design$J[i]
    data.frame(id = i, age = design$tvisit[i, 1:Ji], suvr = y[i, 1:Ji])
  }))

  ## ---- domain audit -------------------------------------------------------
  ## The FITTING model clamps its rate lookup to [fit_YL, fit_YU]. Any latent
  ## value outside that window is territory the estimator structurally cannot
  ## represent, so the fraction landing there is reported, not assumed away.
  obs <- dat$suvr
  lat <- as.numeric(mu_true[!is.na(mu_true)])
  frac_out <- mean(lat < fit_YL | lat > fit_YU)
  audit <- list(
    latent_range = range(lat),
    latent_q     = quantile(lat, c(0.001, 0.01, 0.99, 0.999)),
    obs_range    = range(obs),
    frac_latent_outside_fit_domain = frac_out,
    n_latent_above = sum(lat > fit_YU), n_latent_below = sum(lat < fit_YL)
  )
  if (verbose) {
    cat(sprintf("simulate_amyloid_data: N=%d, %d observations, J in [%d, %d]\n",
                N, nrow(dat), min(design$J), max(design$J)))
    cat(sprintf("  latent SUVR range   : %.3f - %.3f  (1%%-99%%: %.3f - %.3f)\n",
                audit$latent_range[1], audit$latent_range[2],
                audit$latent_q[2], audit$latent_q[3]))
    cat(sprintf("  observed SUVR range : %.3f - %.3f\n",
                audit$obs_range[1], audit$obs_range[2]))
    cat(sprintf("  latent values outside fitting domain [%.2f, %.2f]: %d above, %d below (%.3f%%)\n",
                fit_YL, fit_YU, audit$n_latent_above, audit$n_latent_below,
                100 * frac_out))
  }

  list(dat = dat, y = y, mu_true = mu_true, delta_true = delta_true,
       solver = solver, audit = audit)
}


## -----------------------------------------------------------------------
##  3. Signal-to-noise audit
## -----------------------------------------------------------------------
##
##  The diagnostic the original study most needed and did not have: how much
##  does each subject's latent trajectory actually MOVE over their own
##  follow-up, relative to the measurement noise the estimator has to see
##  through? If this is below ~1x sigma_eps for most subjects, delta_i is not
##  identifiable and no amount of sampling will recover it.
signal_to_noise_audit <- function(truth, design, sim, amy_thres = 0.75,
                                  solver = NULL, max_years = 60) {
  if (is.null(solver)) solver <- sim$solver
  N <- design$N
  first <- design$tvisit[cbind(seq_len(N), 1)]
  last  <- design$tvisit[cbind(seq_len(N), design$J)]
  span  <- last - first

  mu_first <- truth_mu_at(solver, design$x0_true, design$t0_true, first, sim$delta_true)
  mu_last  <- truth_mu_at(solver, design$x0_true, design$t0_true, last,  sim$delta_true)
  dS <- mu_last - mu_first

  se <- design$sigma_eps_true
  alpha <- truth_positivity_age(solver, design$x0_true, design$t0_true,
                                sim$delta_true, amy_thres, max_years = max_years)

  list(
    span = span, dSUVR = dS, ratio = dS / se, alpha_true = alpha,
    summary = data.frame(
      metric = c("follow-up span (yr)", "dSUVR over span", "dSUVR / sigma_eps",
                 "subjects with signal > 1x sigma_eps",
                 "subjects with signal > 2x sigma_eps",
                 "subjects with a positivity age"),
      value = c(sprintf("median %.2f  (range %.2f - %.2f)", median(span), min(span), max(span)),
                sprintf("median %.4f  (IQR %.4f - %.4f)", median(dS),
                        quantile(dS, .25), quantile(dS, .75)),
                sprintf("median %.2f x  (IQR %.2f - %.2f x)", median(dS / se),
                        quantile(dS / se, .25), quantile(dS / se, .75)),
                sprintf("%.0f%%", 100 * mean(dS > se)),
                sprintf("%.0f%%", 100 * mean(dS > 2 * se)),
                sprintf("%.0f%%", 100 * mean(!is.na(alpha)))),
      stringsAsFactors = FALSE
    )
  )
}
