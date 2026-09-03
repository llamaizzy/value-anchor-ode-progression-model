## =============================================================================
##  INTEGRATOR TESTS
##
##  Three separate questions, deliberately not conflated:
##
##   A. ACCURACY  -- does Crank-Nicolson at the operational step size Delta =
##      0.25 yr actually track the true trajectory? The reference here is the
##      EXACT solution, available in closed form because the ODE is separable
##      (see sim_data.R): t(b) - t(a) = exp(-delta) * Integral_a^b dy/R(y).
##      This is what the original study could not test at all, because it
##      generated its data with the estimator's own integrator at the
##      estimator's own step size, so every discretisation error cancelled.
##
##   B. ORDER     -- does the error shrink like h^2, confirming the scheme is
##      second-order and not silently degraded by the grid interpolation?
##
##   C. SYMMETRY  -- is a forward-then-backward round trip exact? This is the
##      property the whole value-anchored construction rests on: integrating
##      out from a subject's reference point in both directions must be
##      self-consistent. Mirrors the parent project's
##      check_integrator_consistency.R, including its central design point:
##      the backward leg must start from the FORWARD LEG'S OWN endpoint, not
##      from an externally computed reference, or the test conflates symmetry
##      with truncation error.
## =============================================================================

source("tests/_helpers.R")
source("data_gen.R")
source("sim_data.R")
source("model.R")

test_init("INTEGRATOR")

## ---- setup: the fitting grid, built from the analytic truth --------------
truth  <- simulate_true_rate()
solver <- build_truth_solver(truth)

YL <- 0.30; YU <- 1.70; step_y <- 0.002
y_grid <- seq(YL, YU, by = step_y); n_ygrid <- length(y_grid)
r_grid <- truth$rate_fun(y_grid)          # the estimator's piecewise-linear R(y)

## Exact trajectory from the separable solve.
exact_mu <- function(x0, dt, delta) solver$G_to_y(solver$y_to_G(x0) + dt * exp(delta))

## Chained integration with the model's own stepper.
walk <- function(x0, dt, delta, Delta, method = 2) {
  M <- max(ceiling(abs(dt) / Delta), 1L)
  h <- dt / M
  mu <- x0
  for (m in seq_len(M)) {
    mu <- step_R(mu, h, r_grid, YL, step_y, n_ygrid, exp(delta), method)
  }
  mu
}

## ---- A. accuracy at the operational settings ----------------------------
## Spans chosen to cover realistic follow-up (up to ~9 yr) at rate multipliers
## spanning the central mass of exp(delta) for sigma_delta = 0.4.
grid <- expand.grid(x0 = c(0.55, 0.75, 0.95, 1.15, 1.35),
                    dt = c(1, 3, 6, 9),
                    delta = c(-0.8, 0, 0.8))
grid$exact  <- mapply(exact_mu, grid$x0, grid$dt, grid$delta)
grid$err_cn <- mapply(function(x0, dt, d)
  abs(walk(x0, dt, d, Delta = 0.25, method = 2) - exact_mu(x0, dt, d)),
  grid$x0, grid$dt, grid$delta)

## IN-DOMAIN vs OUT-OF-DOMAIN. The estimator's rate curve is defined only on
## [YL, YU]; lookupRate() CLAMPS outside it. A trajectory whose endpoint lies
## past YU is therefore integrating a deliberately flattened rate, and the gap
## to the exact solution there measures domain truncation, NOT integrator
## error. The two are separated rather than averaged together.
in_dom <- grid$exact <= YU & grid$exact >= YL

cat("\n  Crank-Nicolson vs EXACT solution, Delta = 0.25 yr (operational setting)\n")
cat(sprintf("    endpoints INSIDE  [%.2f, %.2f]: n = %2d, max error %.3e, median %.3e SUVR\n",
            YL, YU, sum(in_dom), max(grid$err_cn[in_dom]), median(grid$err_cn[in_dom])))
cat(sprintf("    endpoints OUTSIDE [%.2f, %.2f]: n = %2d, max error %.3e  <- domain truncation,\n",
            YL, YU, sum(!in_dom), if (any(!in_dom)) max(grid$err_cn[!in_dom]) else NA))
cat("       not integrator error: past YU the estimator integrates a clamped R(y).\n")
cat(sprintf("    worst in-domain case is %.0fx smaller than sigma_eps = 0.05\n",
            0.05 / max(grid$err_cn[in_dom])))
cat("\n    NOTE: the out-of-domain corner (x0 = 1.35, span 9 yr, delta = +0.8) does not\n")
cat("    occur under the study design -- x0 is capped at 1.3 and the realised maximum\n")
cat("    latent SUVR at N = 1000 is 1.53, with ZERO values above YU. Measured in\n")
cat("    sim_data.R's domain audit, which is printed on every data generation.\n")

check("CN error at Delta=0.25 vs exact (in-domain)", max(grid$err_cn[in_dom]), "<", 1e-4)

## Explicit Euler at the same settings, for contrast -- reported, not asserted.
err_eu <- mapply(function(x0, dt, d)
  abs(walk(x0, dt, d, Delta = 0.25, method = 0) - exact_mu(x0, dt, d)),
  grid$x0, grid$dt, grid$delta)
cat(sprintf("\n    (explicit Euler, in-domain: max error %.3e -- %.0fx worse than CN)\n",
            max(err_eu[in_dom]), max(err_eu[in_dom]) / max(grid$err_cn[in_dom])))

## ---- B. convergence order ------------------------------------------------
## The total error has TWO sources: time discretisation, which falls like h^2,
## and the piecewise-linear interpolation of R(y) on the value grid, which does
## NOT depend on h at all. Below some step size the second dominates and the
## total stops converging -- an error FLOOR, not a defect. The order is
## therefore measured in the regime where time discretisation dominates, and
## the floor is located and reported separately.
cat("\n  Convergence order (x0 = 0.75, span 8 yr, delta = 0.3)\n")
cat("    Delta      max|error|     ratio    implied order\n")
Ds <- c(1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125)
e <- sapply(Ds, function(D) abs(walk(0.75, 8, 0.3, Delta = D, method = 2) - exact_mu(0.75, 8, 0.3)))
orders <- rep(NA_real_, length(Ds))
for (i in seq_along(Ds)) {
  ratio <- if (i == 1) NA else e[i - 1] / e[i]
  orders[i] <- if (i == 1) NA else log2(ratio)
  cat(sprintf("    %7.5f   %.4e   %7s   %7s%s\n", Ds[i], e[i],
              if (is.na(ratio)) "-" else sprintf("%.2f", ratio),
              if (is.na(orders[i])) "-" else sprintf("%.2f", orders[i]),
              if (!is.na(orders[i]) && orders[i] < 1.7) "   <- interpolation floor" else ""))
}
## Only the refinements at or coarser than the operational step size are
## asserted: that is the regime the model actually runs in.
op_regime <- which(Ds <= 1.0 & Ds >= 0.125)
check("CN convergence order (Delta >= 0.125, time-error regime)",
      min(orders[op_regime], na.rm = TRUE), ">=", 1.85, fmt = "%.2f")
first_floor <- which(!is.na(orders) & orders < 1.7)[1]
cat(sprintf("\n    Interpolation error floor sets in at Delta ~ %.5f, around %.1e SUVR.\n",
            if (is.na(first_floor)) min(Ds) else Ds[first_floor], min(e)))
cat("    At the operational Delta = 0.25 the scheme is still firmly second-order,\n")
cat("    i.e. time discretisation -- not the value grid -- is the binding constraint.\n")

## ---- C. time-reversal symmetry ------------------------------------------
## Forward from (t0, x0) to get THIS METHOD'S OWN endpoint, then backward from
## that endpoint with the same method and step size; compare to x0.
round_trip <- function(x0, dt, delta, Delta, method) {
  fwd <- walk(x0, dt, delta, Delta, method)
  bwd <- walk(fwd, -dt, delta, Delta, method)
  abs(bwd - x0)
}

cat("\n  Forward-then-backward round-trip error (x0 = 0.75, span 8 yr, delta = 0.3)\n")
cat("    Delta      Euler        Heun         Crank-Nicolson\n")
rt <- sapply(c(0.25, 0.5, 1.0, 2.0), function(D)
  sapply(c(0, 1, 2), function(m) round_trip(0.75, 8, 0.3, D, m)))
for (i in seq_along(c(0.25, 0.5, 1.0, 2.0))) {
  D <- c(0.25, 0.5, 1.0, 2.0)[i]
  cat(sprintf("    %5.2f    %.3e    %.3e    %.3e\n", D, rt[1, i], rt[2, i], rt[3, i]))
}
cat("\n    Crank-Nicolson is symmetric BY CONSTRUCTION (the trapezoid equation is\n")
cat("    invariant under swapping endpoints with h -> -h), so its round-trip error\n")
cat("    sits at machine precision and does NOT grow with the step size. The\n")
cat("    explicit methods' error grows with Delta. This is why CN is the default.\n")

check("CN round-trip error at all step sizes", max(rt[3, ]), "<", 1e-10)

## Sanity: the explicit methods really are worse, i.e. the test can detect a
## failure of symmetry at all.
check("Euler round-trip is detectably worse than CN",
      max(rt[1, ]) > 1e4 * max(rt[3, ]), "TRUE")

test_summary()
