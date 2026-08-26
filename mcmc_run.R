## =============================================================================
## Builds constants/data/inits from a long-format data.frame, configures and
## compiles the NIMBLE sampler for amyloidOdeCode, and runs it.
## =============================================================================

source('model.R')

build_and_run_amyloid <- function(dat, YL = 0.30, YU = 1.70, K = 8, nGrid = 201,
                                  maxSub = 0.25, niter = 5000, nburnin = 2000,
                                  nchains = 3, thin = 1) {
  
  ph <- prepare_amyloid(dat, YL = YL, YU = YU, K = K, nGrid = nGrid)
  
  consts <- list(N = ph$N, K = ph$K, Jmax = ph$Jmax, nGrid = nGrid,
                 J = ph$J, p = ph$p, t0 = ph$t0, yhat = ph$yhat, n = ph$n,
                 tvisit = ph$tvisit, ygrid = ph$ygrid, Bgrid = ph$Bgrid,
                 maxSub = maxSub)
  dat_nimble <- list(y = ph$y)
  
  inits <- function() list(
    theta = seq(-9, -2, length.out = K),
    lambda = rgamma(ph$K, 1, 1),
    sigma2_theta = 1,
    sigma_delta = 0.5,
    sigma_eps = 0.05,
    delta = rep(0, ph$N),
    eps = rep(0, ph$N)
  )
  
  model <- nimbleModel(amyloidOdeCode, constants = consts, data = dat_nimble,
                       inits = inits(), calculate = FALSE)
  cmodel <- compileNimble(model)
  
  conf <- configureMCMC(model, monitors = c("theta", "sigma2_theta", "sigma_delta",
                                            "sigma_eps", "delta", "x"))
  # spline coefficients and (eps, delta) are strongly correlated -- block them
  conf$removeSamplers(paste0("theta[1:", ph$K, "]"))
  conf$addSampler(target = paste0("theta[1:", ph$K, "]"), type = "AF_slice")
  for (i in seq_len(ph$N)) {
    eps_i   <- paste0("eps[", i, "]")
    delta_i <- paste0("delta[", i, "]")
    conf$removeSamplers(c(eps_i, delta_i))
    conf$addSampler(target = c(eps_i, delta_i), type = "AF_slice")
  }
  
  mcmc <- buildMCMC(conf)
  cmcmc <- compileNimble(mcmc, project = model)
  
  t_start <- Sys.time()
  samples <- runMCMC(cmcmc, niter = niter, nburnin = nburnin,
                     nchains = nchains, thin = thin, inits = inits,
                     samplesAsCodaMCMC = TRUE)
  elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  
  list(samples = samples, prepped = ph, consts = consts, data = dat_nimble,
       time = elapsed)
}