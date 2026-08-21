## =============================================================================
## Value-anchored ODE model for individual amyloid (SUVR) progression
## =============================================================================

library(nimble)
library(splines)
library(coda)

## -----------------------------------------------------------------------
## 1. NIMBLE helper functions: grid lookup, one implicit-trapezoid step,
##    sub-stepped integration over an interval, full subject trajectory.
## -----------------------------------------------------------------------

# Given y, return rate R(y) and its slope via linear interpolation on a fixed grid
# Rgrid defined on an evenly-spaced grid ygrid.
rateAndDeriv <- nimbleFunction(
  run = function(y = double(0), ygrid = double(1), Rgrid = double(1),
                 mult = double(0)) {
    returnType(double(1))
    nG <- length(ygrid)
    yLo <- ygrid[1]
    yHi <- ygrid[nG]
    yq <- y
    # Clamp curr y into [yLo, yHi]
    if (yq < yLo) yq <- yLo
    if (yq > yHi) yq <- yHi
    
    gridStep <- ygrid[2] - ygrid[1]
    idx <- floor((yq - yLo) / gridStep) + 1
    if (idx < 1) idx <- 1
    if (idx > (nG - 1)) idx <- nG - 1
    
    frac <- (yq - ygrid[idx]) / gridStep
    R0 <- Rgrid[idx] * mult             # Rgrid = population rate R_h(y) pre-evaluated at each grid point
    R1 <- Rgrid[idx + 1] * mult
    out <- numeric(2)
    out[1] <- R0 + frac * (R1 - R0)          # R(y)
    out[2] <- (R1 - R0) / gridStep           # dR/dy on this sub-interval
    return(out)
  }
)

# One trapezoidal integration step of size h, solved by a
# early stopping of fixed Newton iterations starting from the Euler guess.
trapStep <- nimbleFunction(
  run = function(muPrev = double(0), h = double(0),
                 ygrid = double(1), Rgrid = double(1), mult = double(0)) {
    returnType(double(0))
    Rd0 <- rateAndDeriv(muPrev, ygrid, Rgrid, mult)  # evaluate rate at starting point
    Rprev <- Rd0[1]
    muNew <- muPrev + h * Rprev              # take one Euler step as initial guess
    iter <- 1
    step <- 1.0 # dummy to enter loop
    while (iter <= 6 & abs(step) > 1e-8) { # early stopping (6 iterations max)
      Rd1 <- rateAndDeriv(muNew, ygrid, Rgrid, mult)
      Fval <- muNew - muPrev - 0.5 * h * (Rprev + Rd1[1])
      Fprime <- 1 - 0.5 * h * Rd1[2]
      if (abs(Fprime) < 1e-8) Fprime <- 1e-8  # guard a (near) flat rate curve (division by near 0)
      step <- Fval / Fprime
      muNew <- muNew - step
      iter <- iter + 1
    }
    return(muNew) # return converged muNew
  }
)

# Integrate from (mu0, t0) to t1: find biomarker value at time 1 from t0
# default maxSub: 0.25 years
integrateInterval <- nimbleFunction(
  run = function(mu0 = double(0), t0 = double(0), t1 = double(0),
                 ygrid = double(1), Rgrid = double(1), mult = double(0),
                 maxSub = double(0)) {
    returnType(double(0))
    dt <- t1 - t0  # interval between t1 and t0
    nSteps <- ceiling(abs(dt) / maxSub) # number of sub-steps needed so doesn't exceed maxSub
    if (nSteps < 1) nSteps <- 1 # guard prevent zero step loop when t0 == t1
    h <- dt / nSteps         # sub-step size (negative if integrate backwards)
    muCur <- mu0
    # Chain Crank-Nicolson steps together
    for (s in 1:nSteps) {
      muCur <- trapStep(muCur, h, ygrid, Rgrid, mult)
    }
    return(muCur)
  }
)

# Full subject trajectory: bi-directional integration from anchor (t0, x0)
getTraj <- nimbleFunction(
  run = function(x0 = double(0), t0 = double(0),
                 tvisit = double(1), p = double(0), Jn = double(0),
                 ygrid = double(1), Rgrid = double(1), mult = double(0),
                 maxSub = double(0)) {
    returnType(double(1))
    Jmax <- length(tvisit) # tvisit padded to Jmax (longest follow-up in dataset)
    muOut <- numeric(Jmax) # initialize with zeroes to fill in
    pI <- round(p) # convert double value to usable integers
    JnI <- round(Jn)
    
    # forward: t0 -> visit p+1 -> visit p+2 -> ... -> visit Jn
    muCur <- x0
    tCur <- t0
    if (pI < JnI) { # skip if anchor is at/after last visit
      for (j in (pI + 1):JnI) {
        muCur <- integrateInterval(muCur, tCur, tvisit[j], ygrid, Rgrid, mult, maxSub)
        muOut[j] <- muCur
        tCur <- tvisit[j]
      }
    }
    
    # backward: t0 -> visit p -> visit p-1 -> ... -> visit 1
    muCur <- x0
    tCur <- t0
    if (pI >= 1) {
      for (jj in 1:pI) {
        j <- pI - jj + 1 # reverse sequence
        muCur <- integrateInterval(muCur, tCur, tvisit[j], ygrid, Rgrid, mult, maxSub)
        muOut[j] <- muCur
        tCur <- tvisit[j]
      }
    }
    return(muOut) # predicted trajectory value at subject's all actual visit ages 
  }
)

## -----------------------------------------------------------------------
## 2. Preprocessing: per-subject reference point (t-tilde, y-hat, n),
##    bracket index p, and the fixed B-spline design matrix on the rate grid.
##    All computed once, outside the model
##
##    dat must be a long-format data.frame with columns: id, age, suvr
## -----------------------------------------------------------------------

prepare_amyloid <- function(dat, YL = 0.30, YU = 1.70, K = 8, nGrid = 201) {
  dat <- dat[order(dat$id, dat$age), ]
  ids <- unique(dat$id)
  N <- length(ids)
  splitDat <- split(dat, factor(dat$id, levels = ids)) # split one chunk per subject
  Jvec <- vapply(splitDat, nrow, integer(1)) # count each subject's number of visits
  Jmax <- max(Jvec) # longest follow-up across everyone
  
  tvisit <- matrix(NA_real_, N, Jmax) # visit ages 
  yobs   <- matrix(NA_real_, N, Jmax)
  t0     <- numeric(N)   # t-tilde_i (anchor time)
  yhat   <- numeric(N)   # y-hat_i (anchor value)
  neff   <- numeric(N)   # n_i
  pidx   <- integer(N)   # bracket index: pivot/split point for integration
  
  
  for (i in seq_len(N)) {
    di <- splitDat[[i]]
    Ji <- nrow(di)
    tvisit[i, 1:Ji] <- di$age
    yobs[i, 1:Ji]   <- di$suvr
    
    t0[i] <- mean(di$age)
    
    if (Ji <= 2) {
      yhat[i] <- mean(di$suvr)
      neff[i] <- Ji
    } else {
      h_i <- mean(diff(di$age))                     # mean inter-visit interval
      w <- exp(-(di$age - t0[i])^2 / (2 * h_i^2))   # Gausian kernel weights centered at t0[i]
      yhat[i] <- sum(w * di$suvr) / sum(w)          # weighted average yhat_i
      neff[i] <- sum(w)^2 / sum(w^2)                # effective sample size 
    }
    
    # last visit with age <= t0
    pidx[i] <- max(which(di$age <= t0[i]))
  }
  
  # B-spline basis on a fixed evaluation grid (over all subjects)
  interiorProbs <- seq(0.05, 0.95, length.out = max(K - 4, 1)) # evenly spaced from 5th-95th percentiles of observed values
  interiorKnots <- unique(as.numeric(quantile(dat$suvr, probs = interiorProbs, na.rm = TRUE)))
  ygrid <- seq(YL, YU, length.out = nGrid)
  Bgrid <- splines::bs(ygrid, knots = interiorKnots, Boundary.knots = c(YL, YU),
                       degree = 3, intercept = TRUE)
  Bgrid <- matrix(as.numeric(Bgrid), nrow = nGrid)
  Kactual <- ncol(Bgrid)
  
  list(
    N = N, Jmax = Jmax, J = Jvec, ids = ids,
    tvisit = tvisit, y = yobs,
    t0 = t0, yhat = yhat, n = neff, p = pidx,
    ygrid = ygrid, Bgrid = Bgrid, K = Kactual,
    interiorKnots = interiorKnots, YL = YL, YU = YU
  )
}

## -----------------------------------------------------------------------
## 3. NIMBLE model code for amyloid
## -----------------------------------------------------------------------

amyloidOdeCode <- nimbleCode({
  
  ## population log-rate spline
  theta[1] ~ dnorm(-9, sd = 0.1)
  for (k in 2:K) {
    theta[k] ~ dnorm(theta[k - 1], var = lambda[k] / sigma2_theta)
    lambda[k] ~ dexp(1)
  }
  sigma2_theta ~ dgamma(1e-3, 1e-3)
  
  sigma_delta ~ T(dnorm(0, sd = 1), 0, )
  sigma_eps   ~ T(dnorm(0, sd = 1), 0, )
  
  ## population rate curve evaluated once per iteration on the fixed grid
  for (g in 1:nGrid) {
    logRpop[g] <- inprod(Bgrid[g, 1:K], theta[1:K])
    Rpop[g] <- exp(logRpop[g])
  }
  
  for (i in 1:N) {
    
    ## subject random intercept
    delta[i] ~ dnorm(0, sd = sigma_delta)
    mult[i] <- exp(delta[i]) # multiplicative scalar on rate curve
    
    ## reference-point measurement model
    sigma2_ref[i] <- sigma_eps^2 * (1 / n[i] + 1) # σ²_infl = σ²_ε in paper
    eps[i] ~ dnorm(0, var = sigma2_ref[i])
    x[i] <- yhat[i] - eps[i]
    
    ## latent trajectory via symmetric outward integration from anchor
    mu[i, 1:Jmax] <- getTraj(
      x0 = x[i], t0 = t0[i],
      tvisit = tvisit[i, 1:Jmax], p = p[i], Jn = J[i],
      ygrid = ygrid[1:nGrid], Rgrid = Rpop[1:nGrid], mult = mult[i],
      maxSub = maxSub
    )
    
    ## observation model
    for (j in 1:J[i]) {
      y[i, j] ~ dnorm(mu[i, j], sd = sigma_eps)
    }
  }
})

## -----------------------------------------------------------------------
## 4. Derive Age of positivity - re-runs the same bi-directional integration
##    logic until the trajectory crosses amy_thres.
## -----------------------------------------------------------------------

.rate_lookup_R <- function(y, ygrid, Rgrid, mult) {
  nG <- length(ygrid)
  yq <- pmin(pmax(y, ygrid[1]), ygrid[nG])
  gridStep <- ygrid[2] - ygrid[1]
  idx <- pmin(pmax(floor((yq - ygrid[1]) / gridStep) + 1, 1), nG - 1)
  frac <- (yq - ygrid[idx]) / gridStep
  R0 <- Rgrid[idx] * mult
  R1 <- Rgrid[idx + 1] * mult
  list(R = R0 + frac * (R1 - R0), Rp = (R1 - R0) / gridStep)
}

.trap_step_R <- function(muPrev, h, ygrid, Rgrid, mult) {
  Rd0 <- .rate_lookup_R(muPrev, ygrid, Rgrid, mult)
  muNew <- muPrev + h * Rd0$R
  for (iter in 1:6) {
    Rd1 <- .rate_lookup_R(muNew, ygrid, Rgrid, mult)
    Fval <- muNew - muPrev - 0.5 * h * (Rd0$R + Rd1$R)
    Fprime <- 1 - 0.5 * h * Rd1$Rp
    if (abs(Fprime) < 1e-8) Fprime <- 1e-8
    muNew <- muNew - Fval / Fprime
  }
  muNew
}

# Integrate from (mu0, t0) until the trajectory first crosses amy_thres
# returns NA if it never crosses within threshold
predict_positivity_age <- function(x0, t0, amy_thres, ygrid, Rgrid, mult,
                                   maxSub = 0.25, maxYears = 60) {
  dir <- if (x0 < amy_thres) 1 else -1
  muPrev <- x0
  tPrev <- t0
  h <- dir * maxSub
  for (step in seq_len(ceiling(maxYears / maxSub))) {
    muNew <- .trap_step_R(muPrev, h, ygrid, Rgrid, mult)
    tNew <- tPrev + h
    crossed <- (dir == 1 && muPrev < amy_thres && muNew >= amy_thres) ||
      (dir == -1 && muPrev > amy_thres && muNew <= amy_thres)
    if (crossed) { # linear interpolation within last step to estimate age of positivity
      frac <- (amy_thres - muPrev) / (muNew - muPrev)
      return(tPrev + frac * (tNew - tPrev))
    }
    muPrev <- muNew
    tPrev <- tNew
  }
  NA_real_
}

# Convenience wrapper: rebuild the population rate grid (not saved from MCMC) and predict that
# subject's positivity age from that draw.
positivity_age_from_draw <- function(theta_draw, delta_draw, x0, t0, amy_thres,
                                     ygrid, Bgrid, maxSub = 0.25, maxYears = 60) {
  Rgrid <- exp(as.numeric(Bgrid %*% theta_draw))
  predict_positivity_age(x0, t0, amy_thres, ygrid, Rgrid, exp(delta_draw),
                         maxSub = maxSub, maxYears = maxYears)
}