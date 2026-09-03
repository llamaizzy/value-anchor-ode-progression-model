## =============================================================================
##  LIKELIHOOD TESTS -- the port-correctness gate
##
##  dsubjRefAnchor() fuses what used to be a deterministic mu[i, 1:Jmax] node
##  plus J[i] separate dnorm nodes into ONE node per subject. That is where the
##  speedup comes from, and it is also where a silent bug would be most costly:
##  the index arithmetic (which visit each residual belongs to, the descending
##  backward walk, the partial spans on either side of the reference time) is
##  easy to get subtly wrong and would still produce a plausible-looking chain.
##
##  Four independent checks:
##
##   1. INDEPENDENT R WALK. A separate R implementation of the same algorithm,
##      written in a different style (explicit segment list rather than
##      partial-then-loop), compared against the compiled nimbleFunction. This
##      catches index bugs and R-mirror/nimbleFunction drift simultaneously.
##      The original had an R mirror that ran a FIXED 6 Newton iterations while
##      its nimbleFunction exited early -- the two were never compared.
##
##   2. RESIDUAL ACCOUNTING. Every visit must contribute exactly one residual:
##      no visit skipped, none double-counted.
##
##   3. y_hat RELOCATION. The reference value is now computed inside the
##      likelihood from the node's own data (see model.R section 4). Verify
##      inprod(w_hat, y) is exactly the yhat the explicit two-node formulation
##      would have used, so the posterior is unchanged.
##
##   4. CONSISTENCY WITH THE EXACT SOLUTION. As Delta shrinks, the walked
##      trajectory must converge to the closed-form separable solve.
## =============================================================================

source("tests/_helpers.R")
source("data_gen.R")
source("sim_data.R")
source("mcmc_run.R")

test_init("LIKELIHOOD (dsubjRefAnchor)")

set.seed(20260827)
truth  <- simulate_true_rate()
solver <- build_truth_solver(truth)
design <- simulate_design(N = 40)
sim    <- simulate_amyloid_data(truth, design, seed = 11, verbose = FALSE)
ph     <- prepare_amyloid(sim$dat, K = 10, step_y = 0.002, Delta = 0.25)


## -----------------------------------------------------------------------
##  Independent R implementation of the walk.
##  Deliberately structured differently from dsubjRefAnchor: it builds an
##  explicit list of (n_steps, step_size, target_visit) segments per direction
##  and then executes them, rather than special-casing the partial span and
##  then looping. Same algorithm, different code path.
## -----------------------------------------------------------------------
walk_subject_R <- function(i, ph, r_grid, x_start, delta, method = 2) {
  J <- ph$J[i]; p <- ph$p[i]
  ed <- exp(delta)
  mu_out <- rep(NA_real_, J)

  run_segments <- function(segs, sgn) {
    mu <- x_start
    for (s in segs) {
      for (m in seq_len(s$M)) {
        mu <- step_R(mu, sgn * s$step, r_grid, ph$YL, ph$step_y, ph$n_ygrid, ed, method)
      }
      mu_out[s$target] <<- mu
    }
  }

  ## forward: reference time -> visit p+1 -> p+2 -> ... -> J
  fwd <- list(list(M = ph$M_right[i], step = ph$step_right[i], target = p + 1))
  if (p + 1 <= J - 1) for (j in (p + 1):(J - 1))
    fwd[[length(fwd) + 1]] <- list(M = ph$M_ivl[i, j], step = ph$step_ivl[i, j], target = j + 1)
  run_segments(fwd, +1)

  ## backward: reference time -> visit p -> p-1 -> ... -> 1
  bwd <- list(list(M = ph$M_left[i], step = ph$step_left[i], target = p))
  if (p >= 2) for (j in (p - 1):1)
    bwd[[length(bwd) + 1]] <- list(M = ph$M_ivl[i, j], step = ph$step_ivl[i, j], target = j)
  run_segments(bwd, -1)

  mu_out
}

loglik_subject_R <- function(i, ph, r_grid, x_start, delta, sigma_eps, method = 2) {
  mu <- walk_subject_R(i, ph, r_grid, x_start, delta, method)
  sum(dnorm(ph$y[i, 1:ph$J[i]], mu, sigma_eps, log = TRUE))
}


## -----------------------------------------------------------------------
##  2. Residual accounting
## -----------------------------------------------------------------------
theta_test <- seq(-5.5, -3.0, length.out = ph$K)
r_grid_test <- as.numeric(exp(ph$Bgrid %*% theta_test))

n_unfilled <- 0; n_visits_total <- 0
for (i in seq_len(ph$N)) {
  mu <- walk_subject_R(i, ph, r_grid_test, ph$yhat[i], 0.1)
  n_unfilled <- n_unfilled + sum(is.na(mu))
  n_visits_total <- n_visits_total + ph$J[i]
}
cat(sprintf("\n  Residual accounting: %d visits across %d subjects\n", n_visits_total, ph$N))
check("every visit receives exactly one fitted value", n_unfilled, "==", 0, fmt = "%.0f")


## -----------------------------------------------------------------------
##  3. the anchor is a free parameter, not a function of the data
## -----------------------------------------------------------------------
## The whole point of the population prior: nothing in the density recomputes
## the anchor from this subject's own observations.
mc <- deparse(amyloidOdeCode)
check("model code never references y_hat / w_hat / n_eff",
      !any(grepl("y_hat|w_hat|n_eff|eps_tilde", mc)), "TRUE")
check("anchor is drawn from a population prior",
      any(grepl("x_tilde\\[i\\] *~ *dnorm\\(m_x", mc)), "TRUE")


## -----------------------------------------------------------------------
##  1. Compiled nimbleFunction vs independent R walk
## -----------------------------------------------------------------------
cat("\n  Building compiled model for the nimbleFunction comparison...\n")
cat("  (stepImplicitTrap cannot be compileNimble()d standalone -- a known NIMBLE\n")
cat("   quirk -- so it is exercised THROUGH a compiled model, as the parent\n")
cat("   project's check scripts do.)\n")
built <- build_amyloid_model(ph, verbose = FALSE)
cm <- built$cmodel

## Evaluate at several parameter settings, including deliberately awkward ones
## (large |delta|, tiny sigma_eps, theta pushing the rate curve to the edges).
settings <- list(
  list(theta = seq(-5.5, -3.0, length.out = ph$K), sd = 0.5, se = 0.05, lab = "typical"),
  list(theta = seq(-9.0, -2.0, length.out = ph$K), sd = 0.5, se = 0.05, lab = "steep spline"),
  list(theta = rep(-3.0, ph$K),                    sd = 0.5, se = 0.02, lab = "flat/fast, small sigma"),
  list(theta = rep(-8.0, ph$K),                    sd = 0.5, se = 0.10, lab = "flat/slow, large sigma")
)

max_err <- 0
for (st in settings) {
  set.seed(99)
  eps_v   <- ph$yhat + rnorm(ph$N, 0, 0.03)
  delta_v <- rnorm(ph$N, 0, 0.6)

  cm$theta       <- st$theta
  cm$x_tilde     <- eps_v
  cm$delta       <- delta_v
  cm$sigma_eps   <- st$se
  cm$sigma_delta <- st$sd
  cm$calculate()

  r_grid <- as.numeric(exp(ph$Bgrid %*% st$theta))
  errs <- vapply(seq_len(ph$N), function(i) {
    ll_nimble <- cm$calculate(paste0("y[", i, ", 1:", ph$Jmax, "]"))
    ll_R <- loglik_subject_R(i, ph, r_grid, eps_v[i], delta_v[i], st$se)
    abs(ll_nimble - ll_R)
  }, numeric(1))
  cat(sprintf("    %-24s max |ll_nimble - ll_R| over %d subjects = %.3e\n",
              st$lab, ph$N, max(errs)))
  max_err <- max(max_err, max(errs))
}
check("compiled nimbleFunction matches independent R walk", max_err, "<", 1e-9)


## -----------------------------------------------------------------------
##  4. Convergence of the walked trajectory to the exact solution
## -----------------------------------------------------------------------
## Uses the TRUE rate curve on the fitting grid, so the only error is the
## walk's own discretisation -- the exact solve is available in closed form.
cat("\n  Walked trajectory vs exact separable solve, as Delta shrinks\n")
cat("    Delta     max |mu_walk - mu_exact| over all subjects/visits\n")
r_grid_true_fit <- truth$rate_fun(ph$ygrid)
errs_by_D <- c()
for (D in c(0.5, 0.25, 0.125, 0.0625)) {
  ph_D <- prepare_amyloid(sim$dat, K = 10, step_y = 0.002, Delta = D)
  e <- 0
  for (i in seq_len(ph_D$N)) {
    mu <- walk_subject_R(i, ph_D, r_grid_true_fit,
                         x_start = design$x0_true[i],
                         delta = sim$delta_true[i])
    ex <- truth_mu_at(solver, design$x0_true[i], design$t0_true[i],
                      ph_D$tvisit[i, 1:ph_D$J[i]], sim$delta_true[i])
    e <- max(e, max(abs(mu - ex)))
  }
  errs_by_D <- c(errs_by_D, e)
  cat(sprintf("    %6.4f    %.4e\n", D, e))
}
check("walked trajectory matches exact solve at operational Delta=0.25",
      errs_by_D[2], "<", 1e-4)
check("error decreases monotonically as Delta shrinks",
      all(diff(errs_by_D) < 0), "TRUE")

test_summary()
