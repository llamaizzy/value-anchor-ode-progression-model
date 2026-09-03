## =============================================================================
##  GROUND TRUTH: population rate curve + subject-level anchors + visit design
##
##  CORRECTED VERSION -- see FIXES.md, issue #1 (the single most important fix).
##
##  The original data_gen.R built the truth as a B-spline with
##      theta_true <- seq(-9, -4.5, length.out = K)
##  which yields R(y) = 0.0008 - 0.0033 SUVR/yr over the range where subjects
##  actually live. Over a subject's ENTIRE follow-up that is a total change of
##  0.1x - 0.5x sigma_eps: the longitudinal signal was smaller than a single
##  measurement's noise, so delta_i was unidentifiable BY CONSTRUCTION and 77%
##  of subjects never crossed the positivity threshold. No estimator could have
##  recovered that truth; the study was measuring the prior.
##
##  TWO STRUCTURAL CHANGES HERE:
##
##  1. Realistic rate magnitude. R(y) now peaks at ~0.05 SUVR/yr, in line with
##     published amyloid accumulation rates for accumulators (~0.01-0.05
##     SUVR/yr). See tests/test_snr.R for the measured signal-to-noise this
##     produces -- it is checked, not asserted.
##
##  2. The truth is an ANALYTIC FUNCTION, not a B-spline. The original drew
##     theta_true on a B-spline basis and then fitted with a B-spline basis, so
##     the estimator's function space contained the truth exactly. Using an
##     analytic logistic instead means the fitted K=10 quantile-knot basis is
##     genuinely misspecified relative to the truth, which is the honest test.
##     Nothing downstream needs theta_true -- all recovery is assessed
##     functionally, on R(y) itself.
## =============================================================================

library(ggplot2)


## -----------------------------------------------------------------------
##  1. True population rate curve R(y), as an analytic function
## -----------------------------------------------------------------------
##
##  shape = "logistic" (DEFAULT, the approved primary truth):
##      R(y) = R_lo + (R_hi - R_lo) / (1 + exp(-slope * (y - y_mid)))
##    Monotone increasing in value: accumulation accelerates as burden rises
##    and then saturates. Peaks at ~0.052 SUVR/yr.
##
##  shape = "hump" (SECONDARY, wired up but not run by default):
##      a Gaussian bump peaking at ~0.135 SUVR/yr near y = 0.63, matching the
##      known-good synthetic truth used in the parent project's own simulation
##      check (scripts/btaj_value_anchored_amyloid.R section 7). This is the
##      biologically-plateauing case -- accumulation is fastest mid-trajectory
##      and slows at high burden. It gives the estimator MORE curvature to
##      detect, so it is the easier of the two for recovering the SHAPE of
##      R(y), and is included as a robustness contrast rather than the primary.
##
##  Returns the rate on a FINE grid (default step 0.0005 -- 4x finer than the
##  fitting grid) because this grid is used by the data generator, which must
##  not share the estimator's discretisation. See sim_data.R.
simulate_true_rate <- function(YL = 0.30, YU = 1.70,
                               step_y_true = 0.0005,
                               shape = c("logistic", "hump"),
                               R_lo = 0.002, R_hi = 0.052,
                               y_mid = 1.00, slope = 4,
                               hump_peak = 0.135, hump_at = 0.63,
                               hump_sd = 0.30, hump_floor = 0.001) {
  shape <- match.arg(shape)

  rate_fun <- switch(
    shape,
    logistic = function(y) R_lo + (R_hi - R_lo) / (1 + exp(-slope * (y - y_mid))),
    hump     = function(y) hump_floor +
                 (hump_peak - hump_floor) * exp(-(y - hump_at)^2 / (2 * hump_sd^2))
  )

  ygrid <- seq(YL, YU, by = step_y_true)
  Rgrid_true <- rate_fun(ygrid)

  ## R(y) is a RATE -- it must be strictly positive everywhere, or the ODE
  ## dmu/dt = R(mu) stops being monotone and the whole value-anchored
  ## construction (a unique crossing time for any threshold) breaks down.
  stopifnot("true rate curve must be strictly positive" = all(Rgrid_true > 0))

  list(rate_fun = rate_fun, shape = shape,
       ygrid = ygrid, Rgrid_true = Rgrid_true,
       step_y_true = step_y_true, YL = YL, YU = YU)
}


## -----------------------------------------------------------------------
##  2. Subject-level anchors (x0, t0), visit counts J, and visit ages
## -----------------------------------------------------------------------
##
##  Every subject gets between 2 and 7 observations. No single-visit subjects:
##  with J = 1 the trajectory contributes no longitudinal information at all,
##  delta_i is sampled purely from its prior, and including such subjects
##  inflates apparent delta coverage while adding a wasted sampler per subject.
##  (The parent project filters to J >= 2 for exactly this reason.)
##
##  x0 is the subject's TRUE latent SUVR at their own reference time
##  t0 = mean(visit ages) -- the same reference time the estimator uses
##  (prepare_amyloid() in model.R), so x0_true is directly comparable to the
##  x_tilde posterior without any re-integration.
simulate_design <- function(N = 1000,
                            sigma_delta_true = 0.4,
                            sigma_eps_true = 0.05,
                            x0_range = c(0.6, 1.3),
                            baseline_age_range = c(40, 90),
                            visits_range = 2:7,
                            visit_prob = c(0.42, 0.26, 0.13, 0.10, 0.06, 0.03),
                            visit_gap_mean = 1.5,
                            visit_gap_sd = 0.3,
                            min_gap = 0.25) {

  stopifnot("visit_prob must match visits_range" =
              length(visit_prob) == length(visits_range))
  stopifnot("visit_prob must sum to 1" = abs(sum(visit_prob) - 1) < 1e-8)

  x0 <- runif(N, x0_range[1], x0_range[2])

  ## Sample POSITIONS, then index -- never sample(visits_range, ...) directly.
  ## R's sample() treats a length-1 first argument as 1:x, so
  ## sample(2:2, N, prob = 1) does not draw 2's: it errors, or worse, silently
  ## draws from 1:2. That matters because a J-stratified diagnostic (all
  ## subjects with the same visit count) is exactly the kind of thing this
  ## generator gets used for.
  Jvec <- visits_range[sample.int(length(visits_range), N, prob = visit_prob,
                                  replace = TRUE)]
  Jmax <- max(Jvec)

  baseline_age <- runif(N, baseline_age_range[1], baseline_age_range[2])
  tvisit <- matrix(NA_real_, N, Jmax)
  t0 <- numeric(N)

  for (i in seq_len(N)) {
    Ji <- Jvec[i]
    gaps <- pmax(rnorm(Ji - 1, visit_gap_mean, visit_gap_sd), min_gap)
    tvisit[i, 1:Ji] <- baseline_age[i] + c(0, cumsum(gaps))
    t0[i] <- mean(tvisit[i, 1:Ji])
  }

  ## Hard guarantee on the requested design, not a silent assumption.
  stopifnot("every subject must have 2 to 7 observations" =
              min(Jvec) >= 2L && max(Jvec) <= 7L)

  list(N = N, J = Jvec, Jmax = Jmax, tvisit = tvisit,
       x0_true = x0, t0_true = t0, baseline_age = baseline_age,
       sigma_delta_true = sigma_delta_true, sigma_eps_true = sigma_eps_true)
}


## -----------------------------------------------------------------------
##  3. Quick look at whatever truth is currently configured
## -----------------------------------------------------------------------
plot_true_rate <- function(truth, amy_thres = 0.75) {
  df <- data.frame(y = truth$ygrid, R = truth$Rgrid_true)
  ggplot(df, aes(y, R)) +
    geom_line(linewidth = 0.9, colour = "black") +
    geom_vline(xintercept = amy_thres, linetype = "dashed", colour = "gray40") +
    labs(x = "SUVR (y)", y = "R(y)  (SUVR / year)",
         title = sprintf("True rate curve (shape = '%s')", truth$shape),
         subtitle = sprintf("peak %.3f SUVR/yr; dashed line = positivity threshold %.2f",
                            max(truth$Rgrid_true), amy_thres)) +
    theme_minimal(base_size = 10)
}
