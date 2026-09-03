## =============================================================================
##  ORCHESTRATION + RECOVERY DIAGNOSTICS
##  CORRECTED VERSION -- see FIXES.md
##
##  Fixes concentrated here:
##
##  * RMSE WHERE THERE IS DATA. The original computed rate_rmse over the whole
##    [0.30, 1.70] grid, ~29% of which had no observations at all under its
##    design -- so a large part of the reported error was pure prior
##    extrapolation, unrelated to how well the estimator did where it could
##    learn anything. Both numbers are now reported, labelled.
##
##  * ONE ESTIMATOR OF R(y), NOT TWO. The original used exp(B %*% colMeans(theta))
##    for the RMSE but rowMeans(exp(B %*% t(theta))) for the plot -- two
##    different quantities differing by a Jensen term. Everything now uses the
##    posterior mean of R(y) itself.
##
##  * ALL REPLICATES RETAINED. The original returned only the last replicate's
##    subject table and last_run, silently discarding the rest.
##
##  * POSITIVITY AGES ARE THINNED AND VECTORISED (see model.R section 6). The
##    original looped over every draw x every subject in pure R and rebuilt the
##    701 x K matrix product inside the subject loop -- slower than the MCMC.
## =============================================================================

source("data_gen.R")
source("sim_data.R")
source("mcmc_run.R")

library(ggplot2)


## =============================================================================
##  1. One replicate
## =============================================================================
##
##  `built` may be passed in to reuse an already-compiled model across
##  replicates -- valid only when the DESIGN (visit times, J) is unchanged,
##  since every constant except the observed data is a function of the times
##  alone. See mcmc_run.R.
run_one_replicate <- function(truth, design, built = NULL, seed = NULL,
                              amy_thres = 0.75,
                              YL = 0.30, YU = 1.70, K = 10, step_y = 0.002,
                              Delta = 0.25, integrator_method = 2,
                              niter = 6000, nburnin = 2000, nchains = 1,
                              thin_positivity = 200, verbose = TRUE) {

  sim <- simulate_amyloid_data(truth, design, seed = seed, verbose = verbose)
  ph  <- prepare_amyloid(sim$dat, YL = YL, YU = YU, K = K,
                         step_y = step_y, Delta = Delta)

  if (is.null(built)) {
    built <- build_amyloid_model(ph, integrator_method = integrator_method,
                                 verbose = verbose)
  } else {
    ## Reusing a compiled model: the design must match, or the baked-in
    ## constants are wrong for this data.
    stopifnot("cannot reuse a compiled model across a CHANGED design" =
                identical(dim(built$ph$y), dim(ph$y)) &&
                max(abs(built$ph$tvisit - ph$tvisit), na.rm = TRUE) < 1e-12)
    built$ph <- ph                      # y, yhat, w_hat differ; times do not
  }

  fit <- run_amyloid_mcmc(built, y_new = ph$y, niter = niter, nburnin = nburnin,
                          nchains = nchains, seed = seed, verbose = verbose)
  fit$built <- built

  rec <- recovery_summary(truth, design, sim, fit, amy_thres = amy_thres,
                          thin_positivity = thin_positivity)
  list(sim = sim, ph = ph, fit = fit, recovery = rec, built = built)
}


## =============================================================================
##  2. Recovery summary
## =============================================================================
recovery_summary <- function(truth, design, sim, fit, amy_thres = 0.75,
                             thin_positivity = 200) {
  ph  <- fit$prepped
  smp <- fit$samples

  ## ---- R(y): posterior mean of the CURVE, not exp(mean of theta) ---------
  theta_cols  <- paste0("theta[", seq_len(ph$K), "]")
  R_draws     <- exp(smp[, theta_cols, drop = FALSE] %*% t(ph$Bgrid))   # draws x nGrid
  R_mean      <- colMeans(R_draws)
  R_lo        <- apply(R_draws, 2, quantile, 0.025)
  R_hi        <- apply(R_draws, 2, quantile, 0.975)
  R_true_grid <- truth$rate_fun(ph$ygrid)

  ## Where the estimator could actually learn anything: the range of LATENT
  ## values the subjects occupied. Outside it, R(y) is prior extrapolation.
  lat <- as.numeric(sim$mu_true[!is.na(sim$mu_true)])
  sup <- quantile(lat, c(0.01, 0.99))
  in_sup <- ph$ygrid >= sup[1] & ph$ygrid <= sup[2]

  rmse_full <- sqrt(mean((R_mean - R_true_grid)^2))
  rmse_sup  <- sqrt(mean((R_mean[in_sup] - R_true_grid[in_sup])^2))
  cover_R   <- mean(R_true_grid[in_sup] >= R_lo[in_sup] & R_true_grid[in_sup] <= R_hi[in_sup])
  ## Relative error is the more interpretable scale for a rate spanning an
  ## order of magnitude across the domain.
  rel_sup   <- mean(abs(R_mean[in_sup] - R_true_grid[in_sup]) / R_true_grid[in_sup])

  ## ---- variance components ----------------------------------------------
  vc <- do.call(rbind, lapply(
    list(c("sigma_delta", design$sigma_delta_true),
         c("sigma_eps",   design$sigma_eps_true)),
    function(z) {
      d <- smp[, z[1]]; tv <- as.numeric(z[2]); ci <- quantile(d, c(0.025, 0.975))
      data.frame(parameter = z[1], truth = tv, post_mean = mean(d),
                 post_median = median(d), post_sd = sd(d),
                 ci_lo = unname(ci[1]), ci_hi = unname(ci[2]),
                 covered = tv >= ci[1] && tv <= ci[2],
                 rel_bias = (mean(d) - tv) / tv,
                 stringsAsFactors = FALSE)
    }))
  rownames(vc) <- NULL

  ## ---- subject-level delta ----------------------------------------------
  dcols <- paste0("delta[", seq_len(ph$N), "]")
  dd <- smp[, dcols, drop = FALSE]
  delta_df <- data.frame(
    id = seq_len(ph$N), J = ph$J,
    truth = sim$delta_true,
    post_mean = colMeans(dd),
    lo = apply(dd, 2, quantile, 0.025),
    hi = apply(dd, 2, quantile, 0.975),
    stringsAsFactors = FALSE)
  delta_df$covered <- delta_df$truth >= delta_df$lo & delta_df$truth <= delta_df$hi

  ## ---- anchor value x_tilde ---------------------------------------------
  xt <- reconstruct_x_tilde(smp, ph)
  x_df <- data.frame(id = seq_len(ph$N), truth = design$x0_true,
                     post_mean = colMeans(xt),
                     lo = apply(xt, 2, quantile, 0.025),
                     hi = apply(xt, 2, quantile, 0.975))
  x_df$covered <- x_df$truth >= x_df$lo & x_df$truth <= x_df$hi

  ## ---- positivity age ----------------------------------------------------
  alpha_true <- truth_positivity_age(sim$solver, design$x0_true, design$t0_true,
                                     sim$delta_true, amy_thres)
  pa <- positivity_age_posterior(smp, ph, thres = amy_thres,
                                 thin_to = thin_positivity)
  alpha_df <- data.frame(
    id = seq_len(ph$N), truth = alpha_true,
    post_median = apply(pa, 2, median, na.rm = TRUE),
    lo = apply(pa, 2, quantile, 0.025, na.rm = TRUE),
    hi = apply(pa, 2, quantile, 0.975, na.rm = TRUE),
    na_frac = apply(pa, 2, function(z) mean(is.na(z))))
  alpha_df$covered <- with(alpha_df, !is.na(truth) & !is.na(lo) &
                             truth >= lo & truth <= hi)

  list(
    rate = list(ygrid = ph$ygrid, mean = R_mean, lo = R_lo, hi = R_hi,
                truth = R_true_grid, support = sup, in_support = in_sup),
    rate_metrics = data.frame(
      metric = c("RMSE over full domain", "RMSE over observed support",
                 "mean relative error on support", "pointwise 95% coverage on support"),
      value  = c(rmse_full, rmse_sup, rel_sup, cover_R)),
    variance_components = vc,
    delta = delta_df, x_tilde = x_df, alpha = alpha_df,
    delta_coverage = mean(delta_df$covered),
    delta_cor = cor(delta_df$truth, delta_df$post_mean),
    x_coverage = mean(x_df$covered),
    alpha_coverage = mean(alpha_df$covered[!is.na(alpha_df$truth)]),
    support = sup
  )
}


## =============================================================================
##  3. Printed report
## =============================================================================
print_recovery <- function(res, diag = NULL) {
  r <- res$recovery
  cat("\n", strrep("=", 74), "\n", sep = "")
  cat("RECOVERY SUMMARY\n")
  cat(strrep("=", 74), "\n", sep = "")

  cat("\n-- population rate curve R(y) --\n")
  m <- r$rate_metrics
  for (i in seq_len(nrow(m))) cat(sprintf("  %-38s %.5f\n", m$metric[i], m$value[i]))
  cat(sprintf("  observed value support (1%%-99%% of latent): [%.3f, %.3f]\n",
              r$support[1], r$support[2]))

  cat("\n-- variance components --\n")
  print(format(r$variance_components, digits = 4), row.names = FALSE)

  cat("\n-- subject-level parameters --\n")
  cat(sprintf("  delta_i   : correlation(truth, posterior mean) = %.3f, 95%% CI coverage = %.1f%%\n",
              r$delta_cor, 100 * r$delta_coverage))
  cat(sprintf("  x_tilde_i : 95%% CI coverage = %.1f%%\n", 100 * r$x_coverage))
  cat(sprintf("  alpha_i   : 95%% CI coverage = %.1f%% (of %d subjects with a true crossing)\n",
              100 * r$alpha_coverage, sum(!is.na(r$alpha$truth))))

  ## Coverage broken out by J: subjects with more visits carry more
  ## longitudinal information, so recovery should improve with J. If it does
  ## not, that is a red flag the pooled number would hide.
  ## Correlation is the informative column here, not coverage. A subject whose
  ## delta is unidentified gets a wide interval that covers the truth almost
  ## automatically -- high coverage, zero information. Correlation between the
  ## true and estimated delta is what actually says whether the data pinned it
  ## down, and it is the column that exposes the J = 2 problem.
  cat("\n-- delta_i recovery by number of visits --\n")
  sp <- split(r$delta, r$delta$J)
  tab <- do.call(rbind, lapply(names(sp), function(k) {
    d <- sp[[k]]
    data.frame(J = k, n = nrow(d),
               coverage = sprintf("%.0f%%", 100 * mean(d$covered)),
               correlation = sprintf("%.3f",
                 if (nrow(d) > 2) cor(d$truth, d$post_mean) else NA_real_),
               mean_CI_width = sprintf("%.3f", mean(d$hi - d$lo)),
               post_sd = sprintf("%.3f", sd(d$post_mean)),
               stringsAsFactors = FALSE)
  }))
  print(tab, row.names = FALSE)
  cat(sprintf("  (true sd of delta = %.3f)\n", sd(r$delta$truth)))

  if (!is.null(diag)) {
    cat("\n-- convergence --\n")
    cat(sprintf("  %d draws, %d chain(s). R-hat %s.\n", diag$n_draws, diag$nchains,
                if (diag$rhat_available) "reported below" else
                  "NOT AVAILABLE (single chain) -- see FIXES.md on this limitation"))
    print(diag$worst, row.names = FALSE)
    if (!is.null(diag$subject_ess)) {
      s <- diag$subject_ess
      cat(sprintf("  per-subject parameters: ESS min %.0f, median %.0f, max %.0f; %d of %d below 100\n",
                  s["min"], s["median"], s["max"], s["n_below_100"], s["n_total"]))
    }
  }
  cat("\n", strrep("=", 74), "\n", sep = "")
  invisible(NULL)
}


## =============================================================================
##  4. Plots
## =============================================================================
plot_rate_curve <- function(res) {
  r <- res$recovery$rate
  df <- data.frame(y = r$ygrid, est = r$mean, lo = r$lo, hi = r$hi, truth = r$truth)
  ggplot(df, aes(y)) +
    annotate("rect", xmin = res$recovery$support[1], xmax = res$recovery$support[2],
             ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.45) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.25) +
    geom_line(aes(y = truth, colour = "Truth"), linewidth = 0.9) +
    geom_line(aes(y = est, colour = "Posterior mean"), linewidth = 0.9) +
    scale_colour_manual(values = c(Truth = "black", `Posterior mean` = "steelblue"),
                        name = NULL) +
    labs(x = "SUVR (y)", y = "R(y)  (SUVR / year)",
         title = "Population rate curve: truth vs posterior",
         subtitle = "Band: 95% pointwise credible interval. Shaded region: 1%-99% of observed latent values\n(outside it, the curve is prior extrapolation, not an estimate)") +
    theme_minimal(base_size = 10) + theme(legend.position = "bottom")
}

plot_delta_recovery <- function(res) {
  d <- res$recovery$delta
  rng <- range(c(d$truth, d$lo, d$hi))
  ggplot(d, aes(truth, post_mean, colour = covered)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey30") +
    geom_linerange(aes(ymin = lo, ymax = hi), alpha = 0.18) +
    geom_point(size = 1.1, alpha = 0.75) +
    scale_colour_manual(values = c(`TRUE` = "steelblue", `FALSE` = "firebrick"),
                        name = "95% CI covers truth") +
    coord_equal(xlim = rng, ylim = rng) +
    labs(x = expression(delta[true]), y = expression(delta[posterior~mean]),
         title = "Subject-level rate multiplier: truth vs posterior",
         subtitle = sprintf("N = %d, coverage %.1f%%, correlation %.3f",
                            nrow(d), 100 * res$recovery$delta_coverage,
                            res$recovery$delta_cor)) +
    theme_minimal(base_size = 10) + theme(legend.position = "bottom")
}

plot_delta_by_visits <- function(res) {
  d <- res$recovery$delta
  d$Jf <- factor(d$J)
  ggplot(d, aes(Jf, post_mean - truth)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30") +
    geom_boxplot(outlier.size = 0.5, fill = "steelblue", alpha = 0.35) +
    labs(x = "number of visits (J)", y = expression(delta[posterior~mean] - delta[true]),
         title = "delta recovery improves with follow-up",
         subtitle = "More visits carry more longitudinal information about the rate multiplier") +
    theme_minimal(base_size = 10)
}

plot_alpha_recovery <- function(res) {
  a <- res$recovery$alpha[!is.na(res$recovery$alpha$truth), ]
  rng <- range(c(a$truth, a$post_median), na.rm = TRUE)
  ggplot(a, aes(truth, post_median, colour = covered)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey30") +
    geom_point(size = 1.1, alpha = 0.7) +
    scale_colour_manual(values = c(`TRUE` = "steelblue", `FALSE` = "firebrick"),
                        name = "95% CI covers truth") +
    coord_equal(xlim = rng, ylim = rng) +
    labs(x = "true age of positivity", y = "posterior median",
         title = "Age of amyloid positivity: truth vs posterior",
         subtitle = sprintf("coverage %.1f%% over %d subjects with a true crossing",
                            100 * res$recovery$alpha_coverage, nrow(a))) +
    theme_minimal(base_size = 10) + theme(legend.position = "bottom")
}

## Per-subject trajectories (truth vs posterior) now live in plot_results.R as
## fig_trajectories(); see run_study.R / make_figures().

plot_traces <- function(fit, params = c("sigma_delta", "sigma_eps", "theta[1]", "theta[5]", "theta[10]", "tau_theta")) {
  params <- intersect(params, colnames(fit$samples))
  df <- do.call(rbind, lapply(params, function(p)
    data.frame(iter = seq_len(nrow(fit$samples)), value = fit$samples[, p], param = p)))
  ggplot(df, aes(iter, value)) +
    geom_line(linewidth = 0.25, colour = "steelblue") +
    facet_wrap(~param, scales = "free_y") +
    labs(x = "post-burn-in iteration", y = NULL,
         title = "Traceplots",
         subtitle = "Single chain: these show mixing, but cannot rule out a chain stuck in a wrong mode") +
    theme_minimal(base_size = 9)
}

save_all_plots <- function(res, fit,
                           file = "figures/recovery_report.pdf") {
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  pdf(file, width = 8, height = 6)
  print(plot_rate_curve(res)); print(plot_delta_recovery(res))
  print(plot_delta_by_visits(res)); print(plot_alpha_recovery(res))
  print(plot_traces(fit))
  dev.off()
  cat(sprintf("  wrote %s\n", file))
  invisible(file)
}


## =============================================================================
##  5. Multi-replicate driver
## =============================================================================
##
##  The design is held FIXED across replicates and only delta/eps/y are
##  redrawn, which is what makes the compiled model reusable: every model
##  constant except the observed data is a function of the visit times alone.
##  ALL replicates are retained, not just the last.
simu <- function(n_rep = 1, N = 1000, seed = 20260827,
                 truth_shape = "logistic", amy_thres = 0.75,
                 niter = 6000, nburnin = 2000, nchains = 1,
                 K = 10, step_y = 0.002, Delta = 0.25,
                 save_dir = NULL, verbose = TRUE) {

  set.seed(seed)
  truth  <- simulate_true_rate(shape = truth_shape)
  design <- simulate_design(N = N)

  reps <- vector("list", n_rep)
  built <- NULL
  for (it in seq_len(n_rep)) {
    if (verbose) cat(sprintf("\n===== replicate %d of %d =====\n", it, n_rep))
    r <- run_one_replicate(truth, design, built = built, seed = seed + 1000L * it,
                           amy_thres = amy_thres, K = K, step_y = step_y,
                           Delta = Delta, niter = niter, nburnin = nburnin,
                           nchains = nchains, verbose = verbose)
    built <- r$built                       # compile once, reuse thereafter
    r$diag <- mcmc_diagnostics(r$fit)
    if (verbose) print_recovery(r, r$diag)
    ## Drop the heavy compiled handles from the stored record.
    r$built <- NULL; r$fit$built <- NULL
    reps[[it]] <- r
    ## Per-replicate files exist for crash recovery on long multi-replicate
    ## runs, so they are only worth their disk cost when there is more than one
    ## replicate -- study.rds already holds everything.
    if (!is.null(save_dir) && n_rep > 1) {
      dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
      saveRDS(r, file.path(save_dir, sprintf("replicate_%03d.rds", it)))
    }
  }

  agg <- if (n_rep > 1) aggregate_replicates(reps, design) else NULL
  list(truth = truth, design = design, replicates = reps, aggregate = agg,
       built = built)
}

## Frequentist summary across replicates: coverage is a repeated-sampling
## property and needs more than one realisation.
aggregate_replicates <- function(reps, design) {
  vc <- do.call(rbind, lapply(seq_along(reps), function(i) {
    v <- reps[[i]]$recovery$variance_components; v$replicate <- i; v
  }))
  per_rep <- data.frame(
    replicate = seq_along(reps),
    rmse_support = vapply(reps, function(r) r$recovery$rate_metrics$value[2], numeric(1)),
    delta_coverage = vapply(reps, function(r) r$recovery$delta_coverage, numeric(1)),
    delta_cor = vapply(reps, function(r) r$recovery$delta_cor, numeric(1)),
    alpha_coverage = vapply(reps, function(r) r$recovery$alpha_coverage, numeric(1)))
  cover <- aggregate(covered ~ parameter, vc, mean)
  bias  <- aggregate(rel_bias ~ parameter, vc, mean)
  list(variance_components = vc, per_replicate = per_rep,
       coverage = merge(cover, bias, by = "parameter"))
}
