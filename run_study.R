## =============================================================================
##  RUN THE SIMULATION STUDY
##
##      Rscript run_study.R           verify, fit, and make figures  (~10 min)
##      Rscript run_study.R tests     verification only              (~2 min)
##
##  or from inside R:
##      setwd("<this folder>"); source("run_study.R")
##
##  Everything you would normally change lives in CFG, immediately below.
## =============================================================================

CFG <- list(
  N           = 1000,        # subjects
  n_rep       = 1,           # replicates (each extra one costs ~4 min, no recompile)
  seed        = 20260827,
  niter       = 8000,        # total MCMC iterations
  nburnin     = 3000,        #   of which discarded
  nchains     = 1,           # >1 enables R-hat; see README
  truth_shape = "hump",  # or "hump" -- see data_gen.R
  amy_thres   = 0.75,        # amyloid positivity threshold
  K           = 10,          # spline basis dimension
  step_y      = 0.002,       # value-grid resolution
  Delta       = 0.25,        # max ODE sub-step, years
  out_dir     = "results",
  fig_dir     = "figures"
)


## -----------------------------------------------------------------------
##  1. Verification -- always runs first
## -----------------------------------------------------------------------
##  The recovery numbers are only meaningful if the machinery producing them
##  is correct, so a failure here stops the study rather than being reported
##  alongside it.
args <- commandArgs(trailingOnly = TRUE)
tests_only <- length(args) > 0 && args[1] == "tests"

cat("\n===========================================================\n")
cat("  VERIFICATION SUITE\n")
cat("===========================================================\n")

all_ok <- TRUE
for (f in c("tests/test_integrator.R", "tests/test_reshape.R",
            "tests/test_likelihood.R", "tests/test_snr.R")) {
  ok <- tryCatch({
    e <- new.env(parent = globalenv()); sys.source(f, envir = e)
    df <- do.call(rbind, get(".TEST_RESULTS", envir = e)$rows)
    all(df$pass)
  }, error = function(err) {
    cat("  ERROR in ", f, ": ", conditionMessage(err), "\n", sep = ""); FALSE })
  all_ok <- all_ok && isTRUE(ok)
}

cat("\n===========================================================\n")
cat(if (all_ok) "  ALL VERIFICATION TESTS PASSED\n" else
                "  TESTS FAILED -- results would not be trustworthy\n")
cat("===========================================================\n")

if (tests_only) quit(save = "no", status = if (all_ok) 0 else 1)
if (!all_ok) stop("Verification failed; refusing to run the study.")


## -----------------------------------------------------------------------
##  2. The study
## -----------------------------------------------------------------------
source("simu.R")
source("plot_results.R")

cat(sprintf("\n===========================================================\n"))
cat(sprintf("  SIMULATION STUDY: N = %d, %d replicate(s), %d chain(s)\n",
            CFG$N, CFG$n_rep, CFG$nchains))
cat(sprintf("===========================================================\n"))

t0 <- Sys.time()
study <- simu(n_rep = CFG$n_rep, N = CFG$N, seed = CFG$seed,
              truth_shape = CFG$truth_shape, amy_thres = CFG$amy_thres,
              niter = CFG$niter, nburnin = CFG$nburnin, nchains = CFG$nchains,
              K = CFG$K, step_y = CFG$step_y, Delta = CFG$Delta,
              save_dir = NULL, verbose = TRUE)
cat(sprintf("\nTotal wall time: %.1f minutes\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))


## -----------------------------------------------------------------------
##  3. Figures
## -----------------------------------------------------------------------
make_figures(study, out_dir = CFG$fig_dir)

dir.create(CFG$out_dir, showWarnings = FALSE, recursive = TRUE)
study$built <- NULL                      # compiled handles are not serialisable
saveRDS(study, file.path(CFG$out_dir, "study.rds"))
cat(sprintf("Saved %s\n\n", file.path(CFG$out_dir, "study.rds")))
