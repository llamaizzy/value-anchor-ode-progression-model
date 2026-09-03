## =============================================================================
##  SIGNAL-TO-NOISE AUDIT -- the test the original study most needed
##
##  A recovery study can only recover what the data contain. The original truth
##  (theta_true = seq(-9, -4.5, ...)) produced R(y) = 0.0008-0.0033 SUVR/yr
##  where subjects actually sit, so the TOTAL movement of a subject's latent
##  trajectory over their entire follow-up was 0.1x-0.5x the measurement noise.
##  delta_i was unidentifiable by construction and 77% of subjects never crossed
##  the positivity threshold, so most "true" positivity ages were NA.
##
##  No estimator can pass that test, and no amount of MCMC tuning would have
##  helped. This script measures the signal BEFORE any fitting happens, so the
##  adequacy of the truth is demonstrated rather than assumed -- and so that a
##  future change to the truth cannot silently reintroduce the same problem.
##
##  It also reports the ORIGINAL truth side by side, so the size of the fix is
##  on the record.
## =============================================================================

source("tests/_helpers.R")
source("data_gen.R")
source("sim_data.R")

test_init("SIGNAL-TO-NOISE")

set.seed(20260827)
N <- 1000
design <- simulate_design(N = N)
sigma_eps <- design$sigma_eps_true

## ---- the corrected truth ------------------------------------------------
truth <- simulate_true_rate()
sim   <- simulate_amyloid_data(truth, design, seed = 4242, verbose = FALSE)
aud   <- signal_to_noise_audit(truth, design, sim, amy_thres = 0.75)

cat(sprintf("\n  Corrected truth (shape = '%s'), N = %d, sigma_eps = %.3f\n",
            truth$shape, N, sigma_eps))
cat(sprintf("    R(y) spans %.4f - %.4f SUVR/yr over [%.2f, %.2f]\n",
            min(truth$Rgrid_true), max(truth$Rgrid_true), truth$YL, truth$YU))
for (yy in c(0.6, 0.75, 0.9, 1.05, 1.2)) {
  cat(sprintf("      R(%.2f) = %.4f SUVR/yr\n", yy, truth$rate_fun(yy)))
}
cat("\n")
print(aud$summary, row.names = FALSE)

## ---- the ORIGINAL truth, for contrast -----------------------------------
## Reconstructed exactly as data_gen.R had it before the fix: a K=8 B-spline
## with evenly-spaced interior knots and theta_true = seq(-9, -4.5, ...) plus
## the small tilt the original applied.
orig_rate_fun <- local({
  YL <- 0.30; YU <- 1.70; K <- 8; nGrid <- 201
  yg <- seq(YL, YU, length.out = nGrid)
  ik <- seq(YL, YU, length.out = K - 2)[-c(1, K - 2)]
  B <- splines::bs(yg, knots = ik, Boundary.knots = c(YL, YU), degree = 3, intercept = TRUE)
  B <- matrix(as.numeric(B), nrow = nGrid)
  th <- seq(-9, -4.5, length.out = ncol(B))
  th <- th + rev(cumsum(rev(c(0, diff(th))) * 0.15))
  Rg <- as.numeric(exp(B %*% th))
  function(y) approx(yg, Rg, xout = pmin(pmax(y, YL), YU), rule = 2)$y
})
truth_orig <- list(rate_fun = orig_rate_fun, shape = "ORIGINAL",
                   ygrid = seq(0.30, 1.70, by = 0.0005),
                   YL = 0.30, YU = 1.70)
truth_orig$Rgrid_true <- orig_rate_fun(truth_orig$ygrid)
sim_orig <- simulate_amyloid_data(truth_orig, design, seed = 4242, verbose = FALSE)
aud_orig <- signal_to_noise_audit(truth_orig, design, sim_orig, amy_thres = 0.75)

cat(sprintf("\n  ORIGINAL truth, same design, same seed\n"))
cat(sprintf("    R(y) spans %.5f - %.5f SUVR/yr\n",
            min(truth_orig$Rgrid_true), max(truth_orig$Rgrid_true)))
print(aud_orig$summary, row.names = FALSE)

cat("\n  Side by side:\n")
cmp <- data.frame(
  quantity = c("median dSUVR / sigma_eps", "% subjects signal > 1x sigma_eps",
               "% subjects with a positivity age"),
  original = c(sprintf("%.2f x", median(aud_orig$ratio)),
               sprintf("%.0f%%", 100 * mean(aud_orig$dSUVR > sigma_eps)),
               sprintf("%.0f%%", 100 * mean(!is.na(aud_orig$alpha_true)))),
  corrected = c(sprintf("%.2f x", median(aud$ratio)),
                sprintf("%.0f%%", 100 * mean(aud$dSUVR > sigma_eps)),
                sprintf("%.0f%%", 100 * mean(!is.na(aud$alpha_true)))),
  stringsAsFactors = FALSE)
print(cmp, row.names = FALSE)

## ---- assertions on the corrected truth ----------------------------------
## These are the conditions under which asking "can the parameters be
## recovered?" is a meaningful question at all.
check("median signal exceeds 1x sigma_eps", median(aud$ratio), ">", 1.0, fmt = "%.2f")
check("at least half of subjects exceed 1x sigma_eps",
      mean(aud$dSUVR > sigma_eps), ">", 0.5, fmt = "%.3f")
check("at least a quarter of subjects exceed 2x sigma_eps",
      mean(aud$dSUVR > 2 * sigma_eps), ">", 0.25, fmt = "%.3f")
## Not asserted at exactly 1: a subject who is both slow (large negative delta)
## and starts low can legitimately have no crossing within the +/- 60 year
## window the estimator searches. That is a property of the truth, not a bug --
## it just has to be RARE, or the positivity-age diagnostic loses its meaning
## the way it did under the original truth (77% missing).
n_na <- sum(is.na(aud$alpha_true))
cat(sprintf("\n    subjects with no crossing within +/- 60 yr: %d of %d (%.1f%%)\n",
            n_na, N, 100 * n_na / N))
check("at least 99% of subjects have a positivity age",
      mean(!is.na(aud$alpha_true)), ">=", 0.99, fmt = "%.3f")
check("all latent values stay inside the fitting domain",
      sim$audit$frac_latent_outside_fit_domain, "==", 0, fmt = "%.4f")
check("rate curve is strictly increasing in value (logistic truth)",
      all(diff(truth$Rgrid_true) > 0), "TRUE")

## And the diagnostic value of the contrast: confirm the original really was
## the degenerate case, so this test would have caught it.
check("ORIGINAL truth is detectably degenerate (median signal < 0.5x sigma_eps)",
      median(aud_orig$ratio), "<", 0.5, fmt = "%.2f")

test_summary()
