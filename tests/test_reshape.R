## =============================================================================
##  PRECOMPUTATION TESTS
##
##  prepare_amyloid() moved a large amount of arithmetic out of the MCMC's
##  innermost loop. That is only a speedup if the precomputed values are the
##  SAME values the loop would have produced. Everything it emits is checked
##  here against an independent brute-force recomputation from the raw
##  long-format data.
##
##  Also exercises the branches that actually occur in the study design:
##  J = 2 (the plain-mean branch, ~42% of subjects), J >= 3 (the kernel
##  branch), and J = 7 (the longest follow-up).
## =============================================================================

source("tests/_helpers.R")
source("data_gen.R")
source("sim_data.R")
source("model.R")

test_init("PRECOMPUTATION (prepare_amyloid)")

set.seed(20260827)
truth <- simulate_true_rate()
design <- simulate_design(N = 300)
sim <- simulate_amyloid_data(truth, design, seed = 5, verbose = FALSE)
Delta <- 0.25
ph <- prepare_amyloid(sim$dat, K = 10, step_y = 0.002, Delta = Delta)

## ---- design coverage ----------------------------------------------------
cat(sprintf("\n  Design: N = %d, J ranges %d..%d\n", ph$N, min(ph$J), max(ph$J)))
tj <- table(ph$J)
cat("  J distribution: ", paste(sprintf("J=%s:%d", names(tj), as.integer(tj)), collapse = "  "), "\n")
check("every subject has 2..7 observations",
      min(ph$J) >= 2 && max(ph$J) <= 7, "TRUE")
check("J = 2 branch is exercised", sum(ph$J == 2) > 0, "TRUE")
check("J >= 3 kernel branch is exercised", sum(ph$J >= 3) > 0, "TRUE")

## ---- brute-force recomputation, independent of prepare_amyloid ----------
sp <- split(sim$dat, factor(sim$dat$id, levels = unique(sim$dat$id)))
bf <- lapply(sp, function(d) {
  d <- d[order(d$age), ]
  tt <- d$age; yy <- d$suvr; J <- length(tt)
  t0 <- mean(tt)
  if (J <= 2) { w <- rep(1, J); h <- NA_real_ } else {
    h <- (tt[J] - tt[1]) / (J - 1)
    w <- exp(-(tt - t0)^2 / (2 * h^2))
  }
  yhat <- sum(w * yy) / sum(w)
  neff <- sum(w)^2 / sum(w^2)
  p <- findInterval(t0, tt, all.inside = TRUE)
  M <- s <- numeric(max(J - 1, 1))
  for (j in seq_len(J - 1)) {
    dt <- max(tt[j + 1] - tt[j], 1e-6)
    M[j] <- max(ceiling(dt / Delta), 1); s[j] <- dt / M[j]
  }
  dl <- max(t0 - tt[p], 1e-6); Ml <- max(ceiling(dl / Delta), 1)
  dr <- max(tt[p + 1] - t0, 1e-6); Mr <- max(ceiling(dr / Delta), 1)
  list(t0 = t0, yhat = yhat, neff = neff, p = p, J = J, w = w / sum(w),
       M = M, s = s, Ml = Ml, sl = dl / Ml, Mr = Mr, sr = dr / Mr)
})

g <- function(f) vapply(bf, function(z) z[[f]], numeric(1))
check("t0 (reference time)",        max(abs(ph$t0    - g("t0"))),   "<", 1e-12)
check("yhat (reference value)",     max(abs(ph$yhat  - g("yhat"))), "<", 1e-12)
check("n_eff",                      max(abs(ph$n     - g("neff"))), "<", 1e-12)
check("p (bracketing index)",       max(abs(ph$p     - g("p"))),    "==", 0, fmt = "%.0f")
check("M_left",                     max(abs(ph$M_left  - g("Ml"))), "==", 0, fmt = "%.0f")
check("M_right",                    max(abs(ph$M_right - g("Mr"))), "==", 0, fmt = "%.0f")
check("step_left",                  max(abs(ph$step_left  - g("sl"))), "<", 1e-12)
check("step_right",                 max(abs(ph$step_right - g("sr"))), "<", 1e-12)

## Ragged matrices: compare only the populated cells.
errM <- errS <- errW <- 0
for (i in seq_len(ph$N)) {
  Ji <- ph$J[i]
  if (Ji >= 2) {
    errM <- max(errM, max(abs(ph$M_ivl[i, 1:(Ji - 1)]    - bf[[i]]$M)))
    errS <- max(errS, max(abs(ph$step_ivl[i, 1:(Ji - 1)] - bf[[i]]$s)))
  }
  errW <- max(errW, max(abs(ph$w_hat[i, 1:Ji] - bf[[i]]$w)))
}
check("M_ivl (all populated cells)",    errM, "==", 0, fmt = "%.0f")
check("step_ivl (all populated cells)", errS, "<", 1e-12)
check("w_hat (normalised weights)",     errW, "<", 1e-12)
## w_hat and n_eff are still computed by prepare_amyloid -- they are used for
## INITIAL VALUES and diagnostics only. Under the population prior neither
## enters the model's density. See CHANGES.md, change 1.

## ---- structural invariants ----------------------------------------------
## These are what the likelihood's index arithmetic silently depends on.
check("weights sum to 1 for every subject",
      max(abs(rowSums(ph$w_hat) - 1)), "<", 1e-12)
check("padded weight cells are exactly 0",
      max(abs(vapply(seq_len(ph$N), function(i)
        if (ph$J[i] < ph$Jmax) max(abs(ph$w_hat[i, (ph$J[i] + 1):ph$Jmax])) else 0,
        numeric(1)))), "==", 0, fmt = "%.0f")
check("bracket satisfies t[p] <= t0 <= t[p+1]",
      all(vapply(seq_len(ph$N), function(i) {
        tt <- ph$tvisit[i, 1:ph$J[i]]; p <- ph$p[i]
        tt[p] <= ph$t0[i] + 1e-9 && ph$t0[i] <= tt[p + 1] + 1e-9
      }, logical(1))), "TRUE")
check("p is always in [1, J-1]",
      all(ph$p >= 1 & ph$p <= ph$J - 1), "TRUE")
check("n_eff == J exactly when J <= 2",
      max(abs(ph$n[ph$J <= 2] - ph$J[ph$J <= 2])), "<", 1e-12)
check("n_eff < J when J >= 3 (kernel downweights distant visits)",
      all(ph$n[ph$J >= 3] < ph$J[ph$J >= 3]), "TRUE")
check("sub-step sizes never exceed Delta",
      max(ph$step_ivl, ph$step_left, ph$step_right), "<=", Delta + 1e-12)
check("basis is a partition of unity",
      max(abs(rowSums(ph$Bgrid) - 1)), "<", 1e-8)
check("basis dimension equals requested K", ph$K, "==", 10, fmt = "%.0f")

test_summary()
