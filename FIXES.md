# What changed, and why

This study is a corrected version of an earlier one. Your original five files
(`data_gen.R`, `sim_data.R`, `model.R`, `mcmc_run.R`, `simu.R`) are in the parent
directory and are **untouched**; every file here keeps its original name, so you
can diff them one at a time:

```bash
diff ../model.R model.R
diff ../simu.R  simu.R
```

**Part I** is the four things worth understanding. **Part II** is the itemised
detail, for when you are editing the code.

---

# Part I — Key considerations

## 1. A simulated truth must contain recoverable signal

This is the one that mattered most. Everything else could have been perfect and
the study would still have failed.

The truth was built as `theta_true <- seq(-9, -4.5, length.out = K)`. Because the
B-spline basis is a partition of unity, those coefficients **are** log rates:
`exp(-9) = 0.00012`, `exp(-4.5) = 0.011` SUVR/year. Where subjects actually sat,
`R(y)` ran 0.0008–0.0033 SUVR/yr — roughly **10–40× flatter** than published
amyloid accumulation (~0.01–0.05 SUVR/yr).

Measured over the same design:

| | original truth | now |
|---|---|---|
| median change in SUVR over a subject's **entire** follow-up, in units of σ_ε | **0.09 ×** | **1.20 ×** |
| subjects with signal > 1 × σ_ε | **0 %** | 57 % |
| subjects with a positivity age within ±60 yr | **23 %** | 99.9 % |

A subject's whole trajectory moved one tenth of a single measurement's noise.
`delta_i` — the per-subject rate multiplier — was unidentifiable **by
construction**, and the positivity-age diagnostic was mostly comparing `NA`s. No
estimator, sampler, or amount of tuning could have recovered that truth.

**The general lesson:** check that the truth contains signal *before* fitting
anything. `tests/test_snr.R` now does exactly that and asserts the answer, so this
cannot silently return.

## 2. A generator must not share the estimator's numerics

The old `sim_data.R` generated data by calling **the same `getTraj()`
nimbleFunction the fitting model uses**, at the **same `maxSub = 0.25`** and the
**same `nGrid = 201`**. Every discretisation error was therefore identical in
generation and in fitting, and cancelled exactly.

That makes a study structurally incapable of detecting integrator bias — it would
report clean recovery of a scheme that might be badly wrong on real data, which is
precisely what a simulation study of an ODE model exists to check.

The generator now shares nothing with the estimator:

| | generator | fitter |
|---|---|---|
| rate curve | analytic function | K=10 B-spline, quantile knots |
| representation | closed form | linear interpolation, 701-point grid |
| time integration | **exact** | Crank–Nicolson, 0.25 yr steps |

The exact solve works because the ODE is separable: `dmu/dt = R(mu)·exp(delta)`
gives `t(b) − t(a) = exp(−delta)·∫dy/R(y)`. Precompute `G(y) = ∫dy/R(y)` once on a
fine grid and the trajectory is `mu(t) = G⁻¹(G(x0) + (t−t0)·exp(delta))` — no
time-stepping at all, and vectorised over all subjects.

The truth is also an analytic function now rather than a B-spline. Drawing the
truth on a B-spline basis and then fitting with a B-spline basis puts the truth
inside the estimator's function space; now the fitted basis is genuinely
misspecified, which is the honest test.

## 3. The anchor needs a real prior, not one centred on the data

The old specification was

```r
sigma2_ref[i] <- sigma_eps^2 * (1 / n[i] + 1)
eps[i] ~ dnorm(0, var = sigma2_ref[i])
x[i]   <- yhat[i] - eps[i]
```

`yhat[i]` is a kernel-weighted average of **subject i's own observed values**, and
those same values are then scored again in the likelihood. Each observation is
used twice — once to place the anchor, once to measure fit around it. That is not
a Bayesian prior, and the object being sampled is a pseudo-posterior with no
calibration guarantee.

There is a closed-form consequence. Because `sigma_ref` is *proportional to*
`sigma_eps`, the prior's normalising constant contributes an **extra power of
`sigma_eps`** to the marginal likelihood, shifting the divisor of the residual sum
of squares from the correct `J − 2` to `J − 1`:

```
E[sigma_eps_hat] / sigma_eps_true  =  sqrt( (J - 2) / (J - 1) )
```

At J = 2 this is zero: two observations,
two free parameters, the fit is exact and no residual remains.

**What this study does instead:**

```r
m_x   ~ dnorm(1.0, sd = 0.5)          # fixed hyperpriors, NOT data-derived
tau_x ~ T(dnorm(0, sd = 0.5), 0, )
x_tilde[i] ~ dnorm(m_x, sd = tau_x)
```

Each anchor is shrunk toward a population value learned from all subjects jointly,
rather than toward that subject's own observations. The data enter the density
exactly once, the extra power of `sigma_eps` disappears, and the estimate is
unbiased.

`yhat` is still computed in `prepare_amyloid()`, but only as an **initial value**
and for diagnostics. Initial values do not enter the posterior. (Verified:
refitting from a start containing no data at all reproduces the posterior to four
significant figures.)

**What it buys.** Over 5 replicates at N = 1000, both specifications fitted to
identical data:

| | σ_ε (true 0.050) | covered | σ_δ (true 0.400) | covered |
|---|---|---|---|---|
| old anchor | 0.0389 (0.78×) | **0/5** | 0.5363 (1.34×) | **0/5** |
| this study | 0.0494 (0.99×) | 4/5 | 0.4028 (1.01×) | **5/5** |

Per-subject shrinkage also becomes correctly calibrated. Theory says the
regression slope of the posterior mean of `delta_i` on its true value should equal
the shrinkage factor `B = σ_δ²/(σ_δ² + v)`:

| J | B (theory) | slope, old anchor | slope, this study |
|---|---|---|---|
| 2 | 0.063 | 0.159 | **0.064** |
| 4 | 0.330 | 0.515 | **0.324** |
| 7 | 0.620 | 0.811 | **0.619** |
| all | 0.207 | 0.344 | **0.207** |

The old anchor under-shrinks at every visit count — it keeps 34 % of each
subject's noisy individual signal where it should keep 21 %, so the estimated
`delta_i` come out about **1.7× more dispersed** than the data support.

**Worth being clear about what this does *not* buy.** Point-estimate accuracy
barely moves (`delta_i` correlation with truth 0.420 → 0.440; onset-age median
error 2.53 → 2.36 yr). The gain is almost entirely in **calibration**. The
clearest case is the anchor itself: identical correlation (0.989) and identical
median error (0.019) under both, but coverage **83 % → 95 %**. If you only ever
report posterior means this change buys little; if you report intervals, or treat
σ_δ as a scientific result, it is the difference between calibrated and not.

*(The head-to-head runs behind these tables are in the `sim_study_corrected/`
folder alongside this one. You do not need them to run anything here.)*

## 4. Some things cannot be fixed by modelling

`delta_i` is a multiplier on the *rate*, so it is identified only by a subject's
observed slope. For J = 2 with a ~1.5 year gap:

| | |
|---|---|
| noise on the observed slope, `√2·σ_ε/gap` | **0.0475 SUVR/yr** |
| true rate at their burden level | **0.0237 SUVR/yr** |

**The noise is twice the signal.** Recovery tracks the information content almost
exactly. Writing `SNR = rate·√Σ(t−t̄)²/σ_ε`, theory predicts the attainable
correlation is `√B` with `B = σ_δ²/(σ_δ² + 1/SNR²)` — nothing fitted:

| J | median SNR | predicted r | observed r |
|---|---|---|---|
| 2 | 0.48 | 0.243 | 0.259 |
| 4 | 1.66 | 0.573 | 0.564 |
| 7 | 3.54 | 0.791 | 0.792 |

Correlation between predicted and observed across strata: **0.985**. Nothing is
being lost to the sampler or the parameterisation — the data simply do not contain
more. This is what `figures/by_visits.png` shows.

To reach r ≥ 0.70 you would need SNR ≥ 2.7 — roughly **6 visits at 1.5-year
spacing, or 2 visits about 8 years apart**. Halving σ_ε halves the requirement.

**Practical consequence:** report the population-level quantities with confidence.
Onset age recovers at r = 0.96 even for J = 2 subjects, because it leans on the
anchor and the shared rate curve rather than on the individual slope. Report
`delta_i` for short-follow-up subjects as what it is — mostly prior.

---

# Part II — Technical detail

## Prior structure

| | Before | Now |
|---|---|---|
| anchor | `x[i] <- yhat[i] - eps[i]`, `eps[i] ~ N(0, σ_ref[i])`, `σ_ref ∝ σ_eps` | `x_tilde[i] ~ dnorm(m_x, sd = tau_x)`, hyperparameters estimated |
| `yhat`, `n_eff` | model constants entering the density | computed for **initial values and diagnostics only** |
| spline penalty | `sigma2_theta ~ dgamma(1e-3, 1e-3)` | `tau_theta ~ dgamma(1, 1)` |
| scale samplers | NIMBLE default random walk | slice (2–3× the ESS for ~10 % more runtime, same posterior) |

`dgamma(1e-3, 1e-3)` is the classic pathological variance-component prior — under
weak data it piles mass at either extreme, leaving the spline wildly wiggly or
crushed flat. It was also named as a variance while being used as a precision
(`var = lambda[k] / sigma2_theta`).

## Code structure

None of these change the model. They change how long it takes, and whether errors
stay silent.

| What | Before | Now |
|---|---|---|
| **Rate grid** | 402 scalar deterministic nodes (`logRpop[g]`, `Rpop[g]` in a loop over 201 grid points) | **one** vectorised node, `Rpop[1:nGrid] <- exp(Bgrid %*% theta)` |
| **Per-subject likelihood** | `mu[i, 1:Jmax]` deterministic node + `J[i]` separate `dnorm` nodes — ~2000 graph nodes at N = 200 | **one** custom-distribution node per subject |
| **Sub-step counts** | `ceiling(abs(dt)/maxSub)` recomputed inside the innermost loop, every likelihood evaluation | precomputed once in `prepare_amyloid()` |
| **Grid geometry** | whole `ygrid` vector passed into every rate lookup; `length()`, `ygrid[2]-ygrid[1]` recomputed each call | three scalars (`y_min`, `step_y`, `n_ygrid`) |
| **Rate lookup** | returned `c(R, dR/dy)` — allocates a vector every call | two scalar-returning functions |
| **Compilation** | model **and** MCMC recompiled inside every replicate | compiled once, reused |
| **Positivity ages** | looped over every draw × every subject in R, rebuilding the 701×K matrix product *inside* the subject loop — slower than the MCMC itself | thinned, `r_grid` hoisted to per-draw, stepping vectorised across subjects |

The innermost loop runs on the order of 10⁸–10⁹ times per fit, so the middle three
matter more than they look. Net: ~48 ms/iteration at N = 1000.

### Two guards against silent failure

- **`assert_sampler_config()`** (`mcmc_run.R`). The old code called
  `conf$removeSamplers("theta[1:8]")` and never checked it took effect. If that
  string fails to match you keep the default scalar samplers **and** add the block
  — a valid chain at double the cost that nothing would reveal. Now an error.
- **Partition-of-unity check** in `prepare_amyloid()`. That property is what makes
  `theta` interpretable directly as log-rates, and nothing verified it.

### Correctness fixes

1. **No convergence diagnostics existed.** Three chains were run and pooled with
   `as.matrix()` without ever computing R-hat or ESS — `coda` was loaded and never
   used, so "poor recovery" and "unconverged chain" were indistinguishable.
   `mcmc_diagnostics()` now reports ESS, lag-1 autocorrelation, and R-hat when more
   than one chain is run.
2. **RMSE was computed where there was no data** — over the whole `[0.30, 1.70]`
   grid, ~29 % of which had no observations, so much of the reported error was
   prior extrapolation. Now reported on both the full domain and the observed
   support, labelled. The figures shade the supported region.
3. **Two different estimators of `R(y)`**: `exp(B %*% colMeans(theta))` for the
   RMSE but `rowMeans(exp(...))` for the plot, differing by a Jensen term. The
   posterior mean of `R(y)` is used everywhere now.
4. **The true positivity age had no time window** while the estimator searched
   ±60 years, so "truth" and "estimate" were not the same quantity. Under the old
   truth this reported 100 % of subjects having a positivity age when only ~23 %
   crossed within any plausible span.
5. **Init length mismatch**: `theta` used the `K` argument while `lambda` used
   `ph$K`; these diverge silently if `unique()` collapses a knot.
6. **`p <- max(which(age <= t0))`** → `findInterval(t0, tt, all.inside = TRUE)`,
   which cannot return an out-of-range index.
7. **Only the last replicate was retained** — `simu()` returned just the final
   `df` and `last_run`, silently discarding the rest.
8. **`cmodel` was compiled and never used.**
9. **`simulate_amyloid_data(sim_build = ...)`** ignored its `truth`/`design`
   arguments for the model but used them for the returned data frame, so a stale
   `sim_build` gave silently wrong data. The NIMBLE generator is gone entirely.
10. **R mirror drift**: the R mirror ran a fixed 6 Newton iterations while the
    nimbleFunction exited early on tolerance, and the two were never compared.
    `tests/test_likelihood.R` now checks them (agreement 8e-14).

## Simulation settings

| | Before | Now |
|---|---|---|
| true rate curve | B-spline, `theta = seq(-9, -4.5)`, peak 0.012 SUVR/yr | analytic logistic, peak ~0.05 SUVR/yr |
| generator | the estimator's own `getTraj()`, same step size and grid | exact separable solve, independent grid |
| visits per subject | 2–6, unasserted | **2–7**, asserted |
| value grid | `nGrid = 201` (step 0.007) | step 0.002 (701 points) |
| replicates retained | last only | all |

### Two bugs found while building this

- **`sample()` with a length-1 range.** `sample(visits_range, N, prob = p)` breaks
  when `visits_range` has length 1, because R treats a length-1 first argument as
  `1:x`. Now `sample.int()` on positions, then index.
- **`max_intervals == 1`.** With an all-J=2 design, NIMBLE collapses the 1-wide
  matrix slice `M_ivl[i, 1:1]` to a scalar and fails `checkBasics()` against the
  distribution's `double(1)` argument type. `prepare_amyloid()` now floors
  `max_intervals` at 2; the extra column is padding no loop ever reaches.

## Known limitations

- **Single chain, one replicate by default.** Enough for direction and rough
  magnitude, not to pin coverage precisely. `CFG$nchains` and `CFG$n_rep` change it.
- **One truth shape** is run by default; the hump-shaped alternative is
  implemented but not exercised.
- **σ_δ under the corrected prior** averages 1.01× truth over five replicates,
  with individual replicates scattering both sides (0.91–1.07). More replicates
  would tighten that.
