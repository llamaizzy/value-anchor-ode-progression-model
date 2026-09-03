## =============================================================================
##  BUILD / COMPILE / RUN the NIMBLE sampler
##  CORRECTED VERSION -- see FIXES.md
##
##  Three fixes live here:
##
##  1. BUILD AND COMPILE ONCE. The original recompiled both the model and the
##     MCMC inside every simulation replicate -- minutes of C++ compilation per
##     replicate, dominating the actual sampling. Because every precomputed
##     quantity in prepare_amyloid() except yhat is a function of the VISIT
##     TIMES alone, a fixed design means only the observed data changes between
##     replicates; and yhat now lives inside the likelihood node (see model.R
##     section 4), so a replicate is just setData() + run(). Use
##     build_amyloid_model() once, then run_amyloid_mcmc() per replicate.
##
##  2. VERIFY THE SAMPLER SWAP. The original called
##     removeSamplers("theta[1:K]") and never checked it took effect. If that
##     string fails to match, you silently keep the default scalar samplers
##     AND add the AF_slice block -- a valid chain at double the cost that
##     nobody would ever notice. assert_sampler_config() below turns that
##     silent failure into a hard error.
##
##  3. CONVERGENCE DIAGNOSTICS. The original ran 3 chains and then pooled them
##     with as.matrix() without ever computing R-hat or ESS, making "poor
##     recovery" and "unconverged chain" indistinguishable. coda was loaded and
##     never used.
##
##  NOTE ON CHAINS: this runs a SINGLE chain by default, as requested. That
##  means there is no R-hat, and the convergence evidence is ESS + lag-1
##  autocorrelation + traceplots. Those detect poor mixing but NOT a chain
##  stuck in a wrong mode -- a real limitation given this model's known ridge
##  between the level of theta and mean(delta). The multi-chain path is
##  written below and is one argument away (nchains > 1).
## =============================================================================

source("model.R")


## -----------------------------------------------------------------------
##  Constants / data / inits from a prepared design
## -----------------------------------------------------------------------
amyloid_constants <- function(ph, integrator_method = 2) {
  list(
    N = ph$N, K = ph$K, Jm = ph$Jmax, nGrid = ph$n_ygrid,
    max_intervals = ph$max_intervals,
    n_intervals = ph$n_intervals, p = ph$p,
    M_ivl = ph$M_ivl, step_ivl = ph$step_ivl,
    M_left = ph$M_left, step_left = ph$step_left,
    M_right = ph$M_right, step_right = ph$step_right,
    Bgrid = ph$Bgrid, YL = ph$YL, step_y_c = ph$step_y,
    integrator_method = integrator_method
  )
}

## Deliberately NOT at the truth -- a neutral start, so the run also tests that
## the sampler can find its way there.
amyloid_inits <- function(ph, jitter = FALSE) {
  th <- seq(-9, -2, length.out = ph$K)
  ## yhat is used HERE ONLY, as a starting value -- it never enters the
  ## density. Verified: refitting from a start containing no data at all
  ## (every anchor at 1.0) reproduces the posterior to 4 significant figures.
  out <- list(
    theta = th, lambda = rep(1, ph$K), tau_theta = 1,
    sigma_delta = 0.5, sigma_eps = 0.05,
    m_x = mean(ph$yhat), tau_x = sd(ph$yhat),
    x_tilde = ph$yhat, delta = rep(0, ph$N)
  )
  if (jitter) {
    out$theta <- th + rnorm(ph$K, 0, 0.5)
    out$sigma_delta <- runif(1, 0.2, 0.8)
    out$sigma_eps <- runif(1, 0.03, 0.08)
    out$delta <- rnorm(ph$N, 0, 0.2)
  }
  out
}


## -----------------------------------------------------------------------
##  Sampler configuration + the assertion the original lacked
## -----------------------------------------------------------------------
configure_amyloid_mcmc <- function(model, ph, verbose = TRUE) {
  conf <- configureMCMC(model, monitors = c("theta", "tau_theta", "sigma_delta",
                                            "sigma_eps", "m_x", "tau_x",
                                            "delta", "x_tilde"))

  ## theta coefficients are strongly correlated along the spline -- block them.
  conf$removeSamplers("theta")
  conf$addSampler(target = paste0("theta[1:", ph$K, "]"), type = "AF_slice")

  ## (x_tilde[i], delta[i]) trade off directly against each other: raising
  ## the anchor and lowering the rate fit the same data. Block per subject.
  for (i in seq_len(ph$N)) {
    tgt <- c(paste0("x_tilde[", i, "]"), paste0("delta[", i, "]"))
    conf$removeSamplers(tgt)
    conf$addSampler(target = tgt, type = "AF_slice")
  }

  ## Scale parameters get SLICE samplers rather than NIMBLE's default adaptive
  ## random walk. Both are global parameters constrained by all N subjects at
  ## once, so their posteriors are narrow and an RW sampler with a poorly
  ## adapted scale crawls. Measured at N = 200, 3000 post-burn-in draws:
  ##
  ##                        ESS sigma_delta   ESS sigma_eps   runtime
  ##   default RW                      55.7           279.0     25.6s
  ##   slice                          128.1           963.3     28.1s
  ##
  ## 2.3x and 3.5x more effective draws for ~10% more time. Posterior means are
  ## unchanged (sigma_eps 0.0392 vs 0.0394), i.e. this buys precision, not a
  ## different answer -- which also confirms the sigma_eps bias documented in
  ## FIXES.md is structural rather than a mixing artefact.
  conf$removeSamplers(c("sigma_delta", "sigma_eps", "m_x", "tau_x"))
  for (nm in c("sigma_delta", "sigma_eps", "m_x", "tau_x"))
    conf$addSampler(target = nm, type = "slice")

  assert_sampler_config(conf, ph, verbose = verbose)
  conf
}

## Turns the silent "removeSamplers matched nothing" failure into a hard error.
assert_sampler_config <- function(conf, ph, verbose = TRUE) {
  s <- conf$getSamplers()
  types <- vapply(s, function(z) z$name, character(1))
  targs <- lapply(s, function(z) z$target)

  is_af <- grepl("AF_slice", types)
  n_theta_block <- sum(is_af & vapply(targs, function(t)
    any(grepl("^theta\\[", t)) && length(t) == 1L && grepl(":", t[1]), logical(1)))
  ## AF_slice on theta may be reported either as one "theta[1:K]" target or as
  ## K expanded scalar targets -- accept either, but require exactly one such
  ## sampler and no non-AF_slice sampler touching theta.
  theta_samplers <- which(vapply(targs, function(t) any(grepl("^theta\\[", t)), logical(1)))
  theta_non_af <- theta_samplers[!is_af[theta_samplers]]

  subj_blocks <- sum(is_af & vapply(targs, function(t)
    length(t) == 2L && any(grepl("^x_tilde\\[", t)) && any(grepl("^delta\\[", t)),
    logical(1)))

  leftover_scalar <- which(vapply(targs, function(t)
    length(t) == 1L && grepl("^(x_tilde|delta)\\[", t[1]), logical(1)))

  if (length(theta_non_af) > 0)
    stop(sprintf("sampler config: %d non-AF_slice sampler(s) still target theta -- ",
                 length(theta_non_af)),
         "removeSamplers('theta') did not take effect. ",
         "You would be paying for both the scalar samplers and the block.")
  if (length(theta_samplers) != 1L)
    stop(sprintf("sampler config: expected exactly 1 sampler on theta, found %d.",
                 length(theta_samplers)))
  if (subj_blocks != ph$N)
    stop(sprintf("sampler config: expected %d per-subject (eps_tilde, delta) AF_slice blocks, found %d.",
                 ph$N, subj_blocks))
  if (length(leftover_scalar) > 0)
    stop(sprintf("sampler config: %d leftover scalar sampler(s) on x_tilde/delta -- ",
                 length(leftover_scalar)),
         "the per-subject removeSamplers() did not take effect.")

  ## The two scale parameters must be on slice samplers, not the defaults.
  n_scale_slice <- sum(types == "slice" & vapply(targs, function(t)
    length(t) == 1L && t[1] %in% c("sigma_delta", "sigma_eps", "m_x", "tau_x"),
    logical(1)))
  if (n_scale_slice != 4L)
    stop(sprintf("sampler config: expected slice samplers on sigma_delta, sigma_eps, m_x, tau_x; found %d.",
                 n_scale_slice))

  if (verbose) {
    cat(sprintf("  sampler config OK: 1 AF_slice block on theta[1:%d], %d per-subject AF_slice blocks, 4 scale slice samplers, %d samplers total\n",
                ph$K, subj_blocks, length(s)))
  }
  invisible(TRUE)
}


## -----------------------------------------------------------------------
##  Build + compile ONCE
## -----------------------------------------------------------------------
build_amyloid_model <- function(ph, integrator_method = 2, verbose = TRUE) {
  consts <- amyloid_constants(ph, integrator_method)
  inits  <- amyloid_inits(ph)

  if (verbose) cat("Building NIMBLE model...\n")
  t0 <- Sys.time()
  model <- nimbleModel(amyloidOdeCode, constants = consts,
                       data = list(y = ph$y), inits = inits,
                       check = FALSE, calculate = FALSE)
  t_build <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (verbose) cat("Compiling model (C++)...\n")
  t0 <- Sys.time()
  cmodel <- compileNimble(model)
  t_cmodel <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (verbose) cat("Configuring samplers...\n")
  conf <- configure_amyloid_mcmc(model, ph, verbose = verbose)
  mcmc <- buildMCMC(conf)

  if (verbose) cat("Compiling MCMC (C++)...\n")
  t0 <- Sys.time()
  cmcmc <- compileNimble(mcmc, project = model)
  t_cmcmc <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (verbose)
    cat(sprintf("  build %.1fs | compile model %.1fs | compile MCMC %.1fs (total %.1fs)\n",
                t_build, t_cmodel, t_cmcmc, t_build + t_cmodel + t_cmcmc))

  list(model = model, cmodel = cmodel, conf = conf, mcmc = mcmc, cmcmc = cmcmc,
       ph = ph, consts = consts,
       compile_time = c(build = t_build, model = t_cmodel, mcmc = t_cmcmc))
}


## -----------------------------------------------------------------------
##  Run (reusable: pass new y for a new replicate on the SAME design)
## -----------------------------------------------------------------------
##
##  Timing follows the parent project's convention: burn-in is run with
##  reset = TRUE and EXCLUDED from the reported per-iteration cost, and the
##  timed region uses cmcmc$run() directly so coda conversion stays out of it.
run_amyloid_mcmc <- function(built, y_new = NULL, niter = 6000, nburnin = 2000,
                             thin = 1, nchains = 1, reset_inits = TRUE,
                             seed = NULL, verbose = TRUE) {
  ph <- built$ph
  cm <- built$cmodel

  if (!is.null(y_new)) {
    stopifnot("y_new has the wrong shape for this compiled design" =
                all(dim(y_new) == dim(ph$y)))
    cm$y <- y_new
    built$ph$y <- y_new
  }

  chains <- vector("list", nchains)
  times  <- numeric(nchains)

  for (ch in seq_len(nchains)) {
    if (!is.null(seed)) set.seed(seed + ch - 1L)
    if (reset_inits) {
      ini <- amyloid_inits(ph, jitter = (nchains > 1))
      for (nm in names(ini)) cm[[nm]] <- ini[[nm]]
    }
    cm$calculate()

    ## burn-in: run and discard, excluded from timing
    built$cmcmc$run(nburnin, reset = TRUE, progressBar = verbose)
    t0 <- Sys.time()
    built$cmcmc$run(niter - nburnin, reset = FALSE, progressBar = verbose)
    times[ch] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    smp <- as.matrix(built$cmcmc$mvSamples)
    ## mvSamples holds burn-in + post-burn-in when reset = FALSE; keep the tail.
    if (nrow(smp) > (niter - nburnin)) smp <- smp[(nrow(smp) - (niter - nburnin) + 1):nrow(smp), , drop = FALSE]
    if (thin > 1) smp <- smp[seq(1, nrow(smp), by = thin), , drop = FALSE]
    chains[[ch]] <- smp
  }

  samples <- if (nchains == 1) chains[[1]] else do.call(rbind, chains)
  n_post <- niter - nburnin
  if (verbose)
    cat(sprintf("  sampling: %.1fs for %d post-burn-in iterations (%.1f ms/iter)\n",
                sum(times), n_post * nchains, 1000 * sum(times) / (n_post * nchains)))

  ## With one chain `chains` would be a byte-for-byte duplicate of `samples`,
  ## which doubles the size of every saved replicate (the samples matrix alone
  ## is ~80 MB at N = 1000). Only keep the per-chain split when it carries
  ## information the pooled matrix cannot -- i.e. when R-hat is computable.
  list(samples = samples, chains = if (nchains > 1) chains else NULL,
       prepped = ph,
       time = sum(times), per_iter_ms = 1000 * sum(times) / (n_post * nchains),
       niter = niter, nburnin = nburnin, nchains = nchains)
}


## One-shot convenience wrapper -- same entry point name as the original.
build_and_run_amyloid <- function(dat, YL = 0.30, YU = 1.70, K = 10,
                                  step_y = 0.002, Delta = 0.25,
                                  integrator_method = 2,
                                  niter = 6000, nburnin = 2000, nchains = 1,
                                  thin = 1, seed = NULL, verbose = TRUE) {
  ph <- prepare_amyloid(dat, YL = YL, YU = YU, K = K, step_y = step_y, Delta = Delta)
  built <- build_amyloid_model(ph, integrator_method = integrator_method, verbose = verbose)
  res <- run_amyloid_mcmc(built, niter = niter, nburnin = nburnin, thin = thin,
                          nchains = nchains, seed = seed, verbose = verbose)
  res$built <- built
  res$consts <- built$consts
  res$compile_time <- built$compile_time
  res
}


## -----------------------------------------------------------------------
##  Convergence diagnostics
## -----------------------------------------------------------------------
##
##  ESS via coda, plus lag-1 autocorrelation (which makes a stuck chain
##  obvious in a way a bare ESS number does not). R-hat is computed only when
##  more than one chain is available -- it is not defined otherwise, and the
##  report says so rather than quietly omitting it.
mcmc_diagnostics <- function(fit, params = NULL, top_n_worst = 10) {
  smp <- fit$samples
  if (is.null(params)) {
    params <- grep("^(theta|sigma_delta|sigma_eps|tau_theta)", colnames(smp), value = TRUE)
  }
  params <- intersect(params, colnames(smp))

  ess <- coda::effectiveSize(coda::as.mcmc(smp[, params, drop = FALSE]))
  ac1 <- apply(smp[, params, drop = FALSE], 2, function(x) {
    n <- length(x); if (sd(x) == 0) return(NA_real_); cor(x[-n], x[-1])
  })

  rhat <- rep(NA_real_, length(params)); names(rhat) <- params
  if (fit$nchains > 1) {
    ml <- coda::mcmc.list(lapply(fit$chains, function(c) coda::as.mcmc(c[, params, drop = FALSE])))
    rhat[] <- tryCatch(coda::gelman.diag(ml, multivariate = FALSE)$psrf[, "Point est."],
                       error = function(e) rep(NA_real_, length(params)))
  }

  tab <- data.frame(parameter = params, ESS = round(as.numeric(ess), 1),
                    ESS_pct = round(100 * as.numeric(ess) / nrow(smp), 1),
                    lag1_autocorr = round(as.numeric(ac1), 3),
                    Rhat = round(as.numeric(rhat), 4),
                    stringsAsFactors = FALSE)
  tab <- tab[order(tab$ESS), ]
  rownames(tab) <- NULL

  ## Per-subject parameters are too numerous to print; summarise them.
  subj <- grep("^(delta|eps_tilde)\\[", colnames(smp), value = TRUE)
  subj_ess <- if (length(subj)) as.numeric(coda::effectiveSize(coda::as.mcmc(smp[, subj, drop = FALSE]))) else numeric(0)

  list(table = tab,
       n_draws = nrow(smp), nchains = fit$nchains,
       rhat_available = fit$nchains > 1,
       subject_ess = if (length(subj_ess))
         c(min = min(subj_ess), median = median(subj_ess), max = max(subj_ess),
           n_below_100 = sum(subj_ess < 100), n_total = length(subj_ess)) else NULL,
       worst = head(tab, top_n_worst))
}
