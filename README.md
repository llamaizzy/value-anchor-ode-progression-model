# Value-anchored ODE progression model — simulation study

A simulation-recovery study for the value-anchored amyloid progression model:
simulate subjects from a known truth, fit the model, and check whether the
parameters come back.

Everything runs from this folder. Nothing outside it is needed.

---

## Run it

```bash
cd sim_study_berkeley

Rscript run_study.R tests    # verification only              (~2 min)
Rscript run_study.R          # verify, fit, and plot          (~10 min)
```

Or from inside R:

```r
setwd("<this folder>")
source("run_study.R")
```

The verification suite runs first every time. If any check fails the study
stops — recovery numbers are only meaningful if the machinery producing them is
correct.

**Requirements:** R with `nimble`, `splines2`, `coda`, `ggplot2`.
Tested on R 4.6.1 with nimble 1.4.2.

---

## What you get

```
figures/
  recovery_figures.pdf     all six figures in one document
  rate_curve.png           population rate curve R(y) vs truth
  anchor.png               anchor value x_i vs truth
  disease_age.png          age of amyloid positivity vs truth
  slope.png                rate multiplier delta_i vs truth
  by_visits.png            how recovery depends on number of visits
  traces.png               MCMC traceplots
results/
  study.rds                full posterior, truth, and design (~78 MB)
```

Every scatter is coloured by number of visits, because that is what governs how
well any individual subject can be recovered. Each carries its correlation with
truth, its 95% credible-interval coverage, and its median absolute error.

The console also prints a recovery summary: the rate curve's error over the
region where data actually exist, the variance components against their true
values, and per-subject recovery broken out by visit count.

---

## Changing the study

Everything lives in `CFG` at the top of `run_study.R`.

| Setting | Default | Notes |
|---|---|---|
| `N` | 1000 | subjects |
| `n_rep` | 1 | replicates; each extra costs ~4 min (no recompile) |
| `niter` / `nburnin` | 8000 / 3000 | MCMC iterations |
| `nchains` | 1 | set > 1 to enable R-hat — see below |
| `truth_shape` | `"logistic"` | or `"hump"` — a second truth, already implemented |
| `amy_thres` | 0.75 | positivity threshold |
| `K`, `step_y`, `Delta` | 10, 0.002, 0.25 | basis size, value grid, ODE sub-step |

Two worth knowing about:

**`nchains = 1` means no R-hat.** Convergence rests on ESS, lag-1
autocorrelation, and the traceplots — which detect poor mixing but *not* a chain
stuck in a wrong mode. The multi-chain path is written and tested; set
`nchains = 3` to use it.

**`truth_shape = "hump"`** gives a rate curve peaking mid-trajectory and slowing
at high burden — biologically the plateau case, and a useful contrast to the
default monotone rise.

---

## The files

| File | What it holds |
|---|---|
| `run_study.R` | **Start here.** Settings, verification, the fit, the figures. |
| `data_gen.R` | The truth (`simulate_true_rate`) and the design (`simulate_design`). |
| `sim_data.R` | Data generation, using an *exact* solve of the ODE. Also the signal-to-noise audit. |
| `model.R` | The integrator, plain-R mirrors of it, `prepare_amyloid()` (all deterministic precomputation), the NIMBLE distribution and model code, and the age-of-positivity functions. |
| `mcmc_run.R` | Builds and compiles once, configures samplers, checks the configuration took effect, reports convergence diagnostics. |
| `simu.R` | The driver: one replicate, the recovery summary, the printed report. |
| `plot_results.R` | The six figures. |
| `tests/` | 40 verification checks — see below. |
| `FIXES.md` | What changed from the earlier version of this study, and why. |

### The verification suite

Each test prints its measurements *and* asserts a tolerance, so you see the size
of what was checked rather than just pass/fail.

| File | What it establishes |
|---|---|
| `test_integrator.R` | The ODE integrator is accurate (4.2e-5 SUVR against an exact solve at the operational step size), second-order, and exactly symmetric under time reversal — the property the value-anchored construction depends on. |
| `test_reshape.R` | Every precomputed quantity matches a brute-force recomputation, plus the structural invariants the likelihood's indexing assumes. |
| `test_likelihood.R` | The compiled NIMBLE distribution matches an independently written R implementation to 8e-14 across four parameter regimes; every visit receives exactly one residual. |
| `test_snr.R` | The simulated truth actually contains recoverable signal — measured *before* any fitting. |

---

## Three NIMBLE quirks worked around in the code

Documented at the point of use, in case you meet them elsewhere:

1. `break` inside a `for` loop does not compile in NIMBLE 1.4.2 — the Newton
   iteration uses a flag-controlled `while` instead.
2. `compileNimble()` on the integrator *standalone* fails with a confusing type
   error. It compiles fine inside a registered distribution, so the tests
   exercise it *through* a compiled model.
3. `registerDistributions()` resolves the density function by ordinary R
   scoping, so these functions must be defined at top level.

---

## Related

`FIXES.md` cites numbers from a head-to-head comparison against the earlier
anchor specification. That comparison — both models fitted to identical data
across five replicates — lives in the `sim_study_corrected/` folder alongside
this one. You do not need it to run anything here.
