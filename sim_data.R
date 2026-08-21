## =============================================================================
## Generate simulated data from ground truth
## =============================================================================

source('model.R')     # getTraj nimbleFunction, reused unchanged from the fitting model

## -----------------------------------------------------------------------
## The generative (simulation) model.
##   - theta, sigma_delta, sigma_eps, Bgrid, ygrid, x0, t0, tvisit, p, J
##     are passed in as CONSTANTS
##   - delta[i] is the only stochastic node NIMBLE actually has to simulate.
##   - mu and y are downstream of delta (mu deterministically, y stochastically).
## -----------------------------------------------------------------------

amyloidSimCode <- nimbleCode({
  for (g in 1:nGrid) {
    logRpop[g] <- inprod(Bgrid[g, 1:K], theta[1:K])
    Rpop[g] <- exp(logRpop[g])
  }
  for (i in 1:N) {
    delta[i] ~ dnorm(0, sd = sigma_delta)
    mult[i] <- exp(delta[i])
    mu[i, 1:Jmax] <- getTraj(
      x0 = x0[i], t0 = t0[i],
      tvisit = tvisit[i, 1:Jmax], p = p[i], Jn = J[i],
      ygrid = ygrid[1:nGrid], Rgrid = Rpop[1:nGrid], mult = mult[i],
      maxSub = maxSub
    )
    for (j in 1:J[i]) {
      y[i, j] ~ dnorm(mu[i, j], sd = sigma_eps)
    }
  }
})

## -----------------------------------------------------------------------
## Step 1: build + compile the model.
## No `data` or `inits` are supplied for delta/y -- they're left NA until
## simulate() fills them in. Everything else is a constant.
## -----------------------------------------------------------------------

build_amyloid_sim_model <- function(truth, design, maxSub = 0.25) {
  
  consts <- list(
    N = design$N, K = truth$K, Jmax = design$Jmax, nGrid = length(truth$ygrid),
    J = design$J, p = design$p,
    x0 = design$x0_true, t0 = design$t0_true, tvisit = design$tvisit,
    ygrid = truth$ygrid, Bgrid = truth$Bgrid, theta = truth$theta_true,
    sigma_delta = design$sigma_delta_true, sigma_eps = design$sigma_eps_true,
    maxSub = maxSub
  )
  
  simModel <- nimbleModel(amyloidSimCode, constants = consts, check = FALSE, calculate = FALSE)
  cSimModel <- compileNimble(simModel)
  
  cSimModel$calculate(c('logRpop', 'Rpop'))
  
  # Step 3: identify every node downstream of the one true stochastic root in dependency order
  # this is what makes simulate() fill in delta, then mu, then y correctly instead of returning NAs
  nodesToSim <- simModel$getDependencies("delta", self = TRUE, downstream = TRUE)
  
  list(model = simModel, cmodel = cSimModel, nodesToSim = nodesToSim,
       design = design, truth = truth)
}

## -----------------------------------------------------------------------
## Step 4: simulate. Cheap to call repeatedly on an already-compiled model
## (e.g. once per replicate in a simulation study) -- only delta/mu/y get
## redrawn each time, the constants (design, truth) stay fixed.
## -----------------------------------------------------------------------

simulate_amyloid_data <- function(truth, design, maxSub = 0.25, seed = NULL,
                                  sim_build = NULL) {
  if (is.null(sim_build)) sim_build <- build_amyloid_sim_model(truth, design, maxSub)
  if (!is.null(seed)) set.seed(seed)
  
  cm <- sim_build$cmodel
  cm$simulate(sim_build$nodesToSim)
  
  y <- cm$y
  mu_true <- cm$mu
  delta_realized <- cm$delta   # the delta actually drawn this call
  N <- design$N
  dat <- do.call(rbind, lapply(seq_len(N), function(i) {
    Ji <- design$J[i]
    data.frame(id = i, age = design$tvisit[i, 1:Ji], suvr = y[i, 1:Ji])
  }))
  
  list(dat = dat, y = y, mu_true = mu_true, delta_realized = delta_realized,
       sim_build = sim_build)
}