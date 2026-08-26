source('data_gen.R')
source('sim_data.R')  # also sources model.R
source('mcmc_run.R')

## =============================================================================
## Orchestra: generate truth -> simulate data -> fit -> compare posterior to truth.
## =============================================================================

simu <- function(outer_iter = 1,
                 # data-generation control
                 N = 200,
                 K_true = 8,
                 sigma_delta_true = 0.4,
                 sigma_eps_true = 0.05,
                 amy_thres = 0.75,
                 # fitting control (should match/bracket the truth's grid)
                 YL = 0.30, YU = 1.70, K_fit = 8, nGrid = 201, maxSub = 0.25,
                 niter = 5000, nburnin = 2000, nchains = 3, thin = 1,
                 # saving
                 save_res = FALSE, out_dir = "results", seed = NULL) {
  
  seed_list <- rep(NA_integer_, outer_iter)
  summaries <- vector("list", outer_iter)
  
  for (it in seq_len(outer_iter)) {
    
    s <- if (!is.null(seed)) seed + it - 1L else sample.int(.Machine$integer.max, 1)
    seed_list[it] <- s
    set.seed(s)
    
    ## 1. ground truth ------------------------------------------------------
    truth  <- simulate_true_rate(YL = YL, YU = YU, K = K_true, nGrid = nGrid)
    design <- simulate_design(N = N, sigma_delta_true = sigma_delta_true,
                              sigma_eps_true = sigma_eps_true)
    
    ## 2. Simulate observed data (delta -> mu -> y) --
    sim <- simulate_amyloid_data(truth, design, maxSub = maxSub)
    
    ## 3. fit -----------------------------------------------------------------
    fit <- build_and_run_amyloid(sim$dat, YL = YL, YU = YU, K = K_fit,
                                 nGrid = nGrid, maxSub = maxSub,
                                 niter = niter, nburnin = nburnin,
                                 nchains = nchains, thin = thin)
    
    ph   <- fit$prepped
    samp <- as.matrix(fit$samples)
    
    ## 4. recovery diagnostics -------------------------------------------------
    sigma_delta_est <- median(samp[, "sigma_delta"])
    sigma_eps_est   <- median(samp[, "sigma_eps"])
    
    theta_cols <- grep("^theta\\[", colnames(samp))
    Rgrid_est <- exp(as.numeric(ph$Bgrid %*% colMeans(samp[, theta_cols, drop = FALSE])))
    rate_rmse <- sqrt(mean((Rgrid_est - approx(truth$ygrid, truth$Rgrid_true,
                                               xout = ph$ygrid)$y)^2))
    
    # positivity-age recovery for a handful of subjects: true crossing age
    # from the noiseless truth vs. posterior credible interval from the fit
    check_ids <- seq_len(min(6, design$N))
    
    # True quantities for checked subjects
    true_age <- sapply(check_ids, function(i) {
      predict_positivity_age(
        design$x0_true[i],
        design$t0_true[i],
        amy_thres,
        truth$ygrid,
        truth$Rgrid_true,
        exp(sim$delta_true[i])
      )
    })
    
    # Posterior positivity-age distributions and credible intervals
    pos_check <- lapply(seq_along(check_ids), function(j) {
      i <- check_ids[j]
      
      x_col <- paste0("x[", i, "]")
      delta_col <- paste0("delta[", i, "]")
      
      ages <- sapply(seq_len(nrow(samp)), function(r) {
        positivity_age_from_draw(
          samp[r, theta_cols],
          samp[r, delta_col],
          samp[r, x_col],
          ph$t0[i],
          amy_thres,
          ph$ygrid,
          ph$Bgrid
        )
      })
      
      ci <- quantile(ages, c(0.025, 0.5, 0.975), na.rm = TRUE)
      c(true = true_age[j], ci, covered = !is.na(true_age[j]) &&
          true_age[j] >= ci[1] &&
          true_age[j] <= ci[3]
      )
    })
    
    pos_check <- do.call(rbind, pos_check)
    
    # Subject-level summary
    df <- data.frame(
      id = check_ids,
      t0 = design$t0_true[check_ids],
      x0 = design$x0_true[check_ids],
      delta = sim$delta_true[check_ids],
      mult = exp(sim$delta_true[check_ids]),
      alpha = true_age
    )
    
    # Store simulation results
    summaries[[it]] <- list(
      seed = s,
      sigma_delta_true = sigma_delta_true,
      sigma_delta_est = sigma_delta_est,
      sigma_eps_true = sigma_eps_true,
      sigma_eps_est = sigma_eps_est,
      rate_rmse = rate_rmse,
      positivity_check = pos_check,
      time = fit$time
    )
    
    if (save_res) {
      rdir <- file.path(out_dir, paste0("iter", it))
      if (!dir.exists(rdir)) dir.create(rdir, recursive = TRUE)
      saveRDS(list(truth = truth, design = design, sim = sim, fit = fit,
                   summary = summaries[[it]]), file.path(rdir, "run.rds"))
    }
    
    # keep the full objects from whichever iteration ran last, so a single
    # simu() call can be handed straight to the diagnostic plots below
    # without needing save_res = TRUE
    last_run <- list(truth = truth, design = design, sim = sim, fit = fit, samp = samp)
  }
  
  list(seed_list = seed_list, summaries = summaries, df = df, last_run = last_run)
}


## =============================================================================
## Diagnostic plots: estimated vs. truth
## =============================================================================
## All four functions below consume the pieces simu() already produces --
## truth, design, sim, fit -- plus a pooled posterior sample matrix `samp`
## (rows = draws across all chains, columns = parameters; this is exactly
## what simu() now returns as res$last_run$samp).

## -----------------------------------------------------------------------
## 1. True vs. estimated rate-vs-value curve R(y)
## -----------------------------------------------------------------------
plot_rate_curve <- function(truth, fit, samp, show_band = TRUE) {
  ph <- fit$prepped
  theta_cols  <- grep("^theta\\[", colnames(samp))
  theta_draws <- samp[, theta_cols, drop = FALSE]
  
  Rgrid_draws <- exp(ph$Bgrid %*% t(theta_draws))   # nGrid x ndraws
  Rgrid_est   <- rowMeans(Rgrid_draws)
  
  df <- data.frame(
    y = ph$ygrid,
    estimated = Rgrid_est,
    truth = approx(truth$ygrid, truth$Rgrid_true, xout = ph$ygrid)$y
  )
  
  p <- ggplot(df, aes(x = y))
  
  if (show_band) {
    band <- t(apply(Rgrid_draws, 1, quantile, probs = c(0.025, 0.975)))
    df$lo <- band[, 1]; df$hi <- band[, 2]
    p <- ggplot(df, aes(x = y)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.2)
  }
  
  p +
    geom_line(aes(y = truth, color = "Truth"), linewidth = 1) +
    geom_line(aes(y = estimated, color = "Estimated"), linewidth = 1) +
    scale_color_manual(values = c(Truth = "black", Estimated = "steelblue"), name = NULL) +
    labs(x = "SUVR (y)", y = "Rate R(y)",
         title = "Population rate curve: truth vs. posterior estimate",
         subtitle = if (show_band) "Shaded band: 95% pointwise posterior interval" else NULL) +
    theme_bw()
}

## -----------------------------------------------------------------------
## 2. True vs. estimated subject-level delta parameters
## -----------------------------------------------------------------------
plot_delta_recovery <- function(sim, samp) {
  N <- length(sim$delta_true)
  delta_cols <- paste0("delta[", seq_len(N), "]")
  stopifnot(all(delta_cols %in% colnames(samp)))
  
  draws <- samp[, delta_cols, drop = FALSE]
  est <- apply(draws, 2, median)
  lo  <- apply(draws, 2, quantile, probs = 0.025)
  hi  <- apply(draws, 2, quantile, probs = 0.975)
  
  df <- data.frame(id = seq_len(N), truth = sim$delta_true,
                   estimated = est, lo = lo, hi = hi)
  df$covered <- df$truth >= df$lo & df$truth <= df$hi
  rng <- range(c(df$truth, df$lo, df$hi))
  
  ggplot(df, aes(x = truth, y = estimated, color = covered)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    geom_linerange(aes(ymin = lo, ymax = hi), alpha = 0.3) +
    geom_point(size = 1.6) +
    scale_color_manual(values = c(`TRUE` = "steelblue", `FALSE` = "firebrick"),
                       name = "95% CI covers truth") +
    coord_equal(xlim = rng, ylim = rng) +
    labs(x = expression(delta[true]), y = expression(delta[estimated]),
         title = "Subject-level random effect: truth vs. posterior median",
         subtitle = sprintf("Coverage: %.0f%% of subjects (%d/%d)",
                            100 * mean(df$covered), sum(df$covered), N)) +
    theme_bw()
}

## -----------------------------------------------------------------------
## Helper: dense continuous trajectory via the same bidirectional trapezoid
## integrator getTraj() uses (see model.R's .trap_step_R), but evaluated on
## a fine, evenly-spaced time grid rather than only at the observed visit
## ages -- for smooth plotting.
## -----------------------------------------------------------------------
predict_traj_path <- function(x0, t0, t_min, t_max, ygrid, Rgrid, mult, maxSub = 0.25) {
  step_from <- function(t_start, t_end) {
    if (t_end == t_start) return(data.frame(t = t_start, mu = x0))
    dir <- sign(t_end - t_start)
    tCur <- t_start; muCur <- x0
    tt <- t_start; mm <- x0
    while (dir * (t_end - tCur) > 1e-9) {
      h <- dir * min(maxSub, abs(t_end - tCur))
      muCur <- .trap_step_R(muCur, h, ygrid, Rgrid, mult)
      tCur <- tCur + h
      tt <- c(tt, tCur); mm <- c(mm, muCur)
    }
    data.frame(t = tt, mu = mm)
  }
  fwd <- if (t_max > t0) step_from(t0, t_max) else data.frame(t = t0, mu = x0)
  bwd <- if (t_min < t0) step_from(t0, t_min) else data.frame(t = t0, mu = x0)
  out <- rbind(bwd[nrow(bwd):1, ], fwd[-1, , drop = FALSE])
  out[order(out$t), ]
}

## -----------------------------------------------------------------------
## 3. Subject-level trajectories: original age scale and disease-age scale
## -----------------------------------------------------------------------
## `ids`       -- which subjects to plot (a handful; each becomes one facet panel)
## `amy_thres` -- positivity threshold anchoring the disease-age scale
## `pad`       -- years to extend the continuous curve beyond the observed visits
##
## Disease-age convention: each curve is re-based to its OWN positivity-
## crossing age (truth's curve by alpha_true, the estimated curve by
## alpha_est), so the panel shows how far the recovered curve/onset is
## shifted from truth even when the two agree well on the original age
## scale. Observed points are shown against alpha_est, since that's the
## only shift available in practice (no alpha_true without truth).
## Subjects whose trajectory never crosses amy_thres (alpha = NA, for
## either truth or the estimate) are silently dropped from the disease-age
## panel only -- they still appear in the original-age panel.
plot_subject_trajectories <- function(ids, truth, design, sim, fit, samp,
                                      amy_thres = 0.75, maxSub = 0.25, pad = 2) {
  ph <- fit$prepped
  theta_cols <- grep("^theta\\[", colnames(samp))
  theta_est  <- colMeans(samp[, theta_cols, drop = FALSE])
  Rpop_est   <- as.numeric(exp(ph$Bgrid %*% theta_est))
  
  traj_list <- list(); pt_list <- list(); obs_list <- list(); lab_list <- list()
  
  for (i in ids) {
    Ji <- design$J[i]
    ages_i <- design$tvisit[i, 1:Ji]
    obs_i  <- sim$dat[sim$dat$id == i, ]
    
    ## truth: anchor, acceleration factor
    x0_true    <- design$x0_true[i]
    t0_true    <- design$t0_true[i]
    mult_true  <- exp(sim$delta_true[i])
    
    ## posterior: anchor, acceleration factor (posterior medians)
    x0_est    <- median(samp[, paste0("x[", i, "]")])
    delta_est <- median(samp[, paste0("delta[", i, "]")])
    mult_est  <- exp(delta_est)
    t0_est    <- ph$t0[i]   # anchor time is fixed by the data, not re-estimated
    
    ## continuous trajectories, original age scale
    t_min <- min(ages_i, t0_true, t0_est) - pad
    t_max <- max(ages_i, t0_true, t0_est) + pad
    
    path_true <- predict_traj_path(x0_true, t0_true, t_min, t_max,
                                   truth$ygrid, truth$Rgrid_true, mult_true, maxSub)
    path_est  <- predict_traj_path(x0_est, t0_est, t_min, t_max,
                                   ph$ygrid, Rpop_est, mult_est, maxSub)
    
    traj_list[[length(traj_list) + 1]] <- rbind(
      data.frame(id = i, curve = "Truth",     scale = "Original age", t = path_true$t, mu = path_true$mu),
      data.frame(id = i, curve = "Estimated", scale = "Original age", t = path_est$t,  mu = path_est$mu)
    )
    pt_list[[length(pt_list) + 1]] <- rbind(
      data.frame(id = i, curve = "Truth",     scale = "Original age", t = t0_true, mu = x0_true),
      data.frame(id = i, curve = "Estimated", scale = "Original age", t = t0_est,  mu = x0_est)
    )
    
    ## positivity ages, for the disease-age scale
    alpha_true <- predict_positivity_age(x0_true, t0_true, amy_thres,
                                         truth$ygrid, truth$Rgrid_true, mult_true)
    alpha_est  <- predict_positivity_age(x0_est, t0_est, amy_thres,
                                         ph$ygrid, Rpop_est, mult_est)
    
    if (!is.na(alpha_true) && !is.na(alpha_est)) {
      traj_list[[length(traj_list) + 1]] <- rbind(
        data.frame(id = i, curve = "Truth",     scale = "Disease age", t = path_true$t - alpha_true, mu = path_true$mu),
        data.frame(id = i, curve = "Estimated", scale = "Disease age", t = path_est$t - alpha_est,  mu = path_est$mu)
      )
      pt_list[[length(pt_list) + 1]] <- rbind(
        data.frame(id = i, curve = "Truth",     scale = "Disease age", t = t0_true - alpha_true, mu = x0_true),
        data.frame(id = i, curve = "Estimated", scale = "Disease age", t = t0_est - alpha_est,  mu = x0_est)
      )
      obs_list[[length(obs_list) + 1]] <- data.frame(
        id = i, age = obs_i$age, suvr = obs_i$suvr, disease_age = obs_i$age - alpha_est
      )
    } else {
      obs_list[[length(obs_list) + 1]] <- data.frame(
        id = i, age = obs_i$age, suvr = obs_i$suvr, disease_age = NA_real_
      )
    }
    
    lab_list[[length(lab_list) + 1]] <- data.frame(
      id = i, label = sprintf("acc: true=%.2f, est=%.2f", mult_true, mult_est)
    )
  }
  
  traj_df <- do.call(rbind, traj_list)
  pt_df   <- do.call(rbind, pt_list)
  obs_df  <- do.call(rbind, obs_list)
  lab_df  <- do.call(rbind, lab_list)
  
  make_panel <- function(scale_name, obs_x) {
    td <- traj_df[traj_df$scale == scale_name, ]
    pd <- pt_df[pt_df$scale == scale_name, ]
    if (nrow(td) == 0) return(NULL)
    
    # place the acceleration-factor label near the top-right corner of each panel
    lab_pos <- aggregate(cbind(t, mu) ~ id, td, max)
    lab_pos <- merge(lab_pos, lab_df, by = "id")
    
    ggplot(td, aes(x = t, y = mu, color = curve)) +
      geom_point(data = obs_df, aes(x = .data[[obs_x]], y = suvr),
                 inherit.aes = FALSE, size = 0.9, alpha = 0.6, na.rm = TRUE) +
      geom_line(linewidth = 0.7) +
      geom_point(data = pd, aes(x = t, y = mu, shape = curve), size = 2.4, na.rm = TRUE) +
      geom_text(data = lab_pos, aes(x = t, y = mu, label = label),
                inherit.aes = FALSE, hjust = 1, vjust = 1, size = 2.6) +
      scale_color_manual(values = c(Truth = "black", Estimated = "steelblue"), name = "Curve") +
      scale_shape_manual(values = c(Truth = 16, Estimated = 17), name = "Curve") +
      facet_wrap(~ id, scales = "free", labeller = label_both) +
      labs(x = if (scale_name == "Original age") "Age (years)" else "Disease age (years from positivity)",
           y = "SUVR",
           title = paste0(scale_name, " scale: truth vs. posterior estimate"),
           subtitle = "Points: observed data; triangle/circle: initial value x_i at t0") +
      theme_bw()
  }
  
  list(
    original    = make_panel("Original age", "age"),
    disease_age = make_panel("Disease age", "disease_age")
  )
}

## -----------------------------------------------------------------------
## 4. Posterior estimates of key variance components vs. truth
## -----------------------------------------------------------------------
variance_component_table <- function(design, samp) {
  pull_summary <- function(param, truth_val) {
    draws <- samp[, param]
    ci <- quantile(draws, c(0.025, 0.975))
    data.frame(
      parameter   = param,
      truth       = truth_val,
      post_mean   = mean(draws),
      post_median = median(draws),
      post_sd     = sd(draws),
      ci_lower    = unname(ci[1]),
      ci_upper    = unname(ci[2]),
      covered     = truth_val >= ci[1] && truth_val <= ci[2]
    )
  }
  
  tab <- rbind(
    pull_summary("sigma_delta", design$sigma_delta_true),
    pull_summary("sigma_eps",   design$sigma_eps_true)
  )
  rownames(tab) <- NULL
  tab
}

## Example single run:
res <- simu(outer_iter = 1, N = 200, niter = 5000, nburnin = 2000, nchains = 3)
res$summaries[[1]]
res$df

## Diagnostic plots for the run above -----------------------------------
lr <- res$last_run
 
plot_rate_curve(lr$truth, lr$fit, lr$samp)
plot_delta_recovery(lr$sim, lr$samp)

traj_plots <- plot_subject_trajectories(ids = 1:6, lr$truth, lr$design, lr$sim,
                                        lr$fit, lr$samp, amy_thres = 0.75)
traj_plots$original
traj_plots$disease_age

variance_component_table(lr$design, lr$samp)