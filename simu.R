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
    
    ## 2. forward-simulate observed data (NIMBLE simulate(): delta -> mu -> y) --
    sim <- simulate_amyloid_data(truth, design, maxSub = maxSub)
    
    ## 3. fit -----------------------------------------------------------------
    fit <- build_and_run_amyloid(sim$dat, YL = YL, YU = YU, K = K_fit,
                                 nGrid = nGrid, maxSub = maxSub,
                                 niter = niter, nburnin = nburnin,
                                 nchains = nchains, thin = thin)
    
    ph   <- fit$prepped
    samp <- as.matrix(fit$samples[[1]])
    
    ## 4. recovery diagnostics -------------------------------------------------
    sigma_delta_est <- median(samp[, "sigma_delta"])
    sigma_eps_est   <- median(samp[, "sigma_eps"])
    
    theta_cols <- grep("^theta\\[", colnames(samp))
    Rgrid_est <- exp(as.numeric(ph$Bgrid %*% colMeans(samp[, theta_cols, drop = FALSE])))
    rate_rmse <- sqrt(mean((Rgrid_est - approx(truth$ygrid, truth$Rgrid_true,
                                               xout = ph$ygrid)$y)^2))
    
    # positivity-age recovery for a handful of subjects: true crossing age
    # from the noiseless truth vs. posterior credible interval from the fit.
    # Uses the delta NIMBLE actually simulated (sim$delta_realized), not any
    # placeholder value, since simulate() draws its own fresh delta each call.
    check_ids <- seq_len(min(5, design$N))
    pos_check <- lapply(check_ids, function(i) {
      true_age <- predict_positivity_age(design$x0_true[i], design$t0_true[i],
                                         amy_thres, truth$ygrid, truth$Rgrid_true,
                                         exp(sim$delta_realized[i]))
      x_col <- paste0("x[", i, "]"); delta_col <- paste0("delta[", i, "]")
      ages <- sapply(seq_len(nrow(samp)), function(r)
        positivity_age_from_draw(samp[r, theta_cols], samp[r, delta_col],
                                 samp[r, x_col], ph$t0[i], amy_thres,
                                 ph$ygrid, ph$Bgrid))
      ci <- quantile(ages, c(0.025, 0.5, 0.975), na.rm = TRUE)
      c(true = true_age, ci, covered = (!is.na(true_age) &&
                                          true_age >= ci[1] && true_age <= ci[3]))
    })
    pos_check <- do.call(rbind, pos_check)
    
    summaries[[it]] <- list(
      seed = s,
      sigma_delta_true = sigma_delta_true, sigma_delta_est = sigma_delta_est,
      sigma_eps_true = sigma_eps_true, sigma_eps_est = sigma_eps_est,
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
  }
  
  list(seed_list = seed_list, summaries = summaries)
}

## Example single run:
res <- simu(outer_iter = 1, N = 150, niter = 3000, nburnin = 1000, nchains = 2)
res$summaries[[1]]