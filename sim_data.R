## =============================================================================
## Forward-simulates OBSERVED data (id, age, suvr) from the ground truth
## produced by data_gen.R, using the same getTraj() ODE integrator that the
## fitting model uses
## =============================================================================

source('model.R')     # for the getTraj / integrateInterval nimbleFunctions

## -----------------------------------------------------------------------
## Forward-simulate one subject's trajectory at its true visit ages, then
## add iid measurement noise sigma_eps_true.
## -----------------------------------------------------------------------

simulate_amyloid_data <- function(truth, design, maxSub = 0.25, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  N <- design$N
  Jmax <- design$Jmax
  y <- matrix(NA_real_, N, Jmax)
  mu_true <- matrix(NA_real_, N, Jmax)   # noiseless trajectories, kept for validation
  
  for (i in seq_len(N)) {
    Ji <- design$J[i]
    tv <- design$tvisit[i, 1:Ji]
    mult_i <- exp(design$delta_true[i])
    
    # anchor sits at the visit closest to (but not after) t0_true[i]
    p_i <- max(which(tv <= design$t0_true[i]))
    if (!is.finite(p_i)) p_i <- 0
    
    mu_i <- getTraj(
      x0 = design$x0_true[i], t0 = design$t0_true[i],
      tvisit = tv, p = p_i, Jn = Ji,
      ygrid = truth$ygrid, Rgrid = truth$Rgrid_true, mult = mult_i,
      maxSub = maxSub
    )
    mu_true[i, 1:Ji] <- mu_i
    y[i, 1:Ji] <- mu_i + rnorm(Ji, 0, design$sigma_eps_true)
  }
  
  # long format matching for dataset
  dat <- do.call(rbind, lapply(seq_len(N), function(i) {
    Ji <- design$J[i]
    data.frame(id = i, age = design$tvisit[i, 1:Ji], suvr = y[i, 1:Ji])
  }))
  
  list(dat = dat, y = y, mu_true = mu_true)
}