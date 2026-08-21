## =============================================================================
## Generates the GROUND TRUTH: population rate curve, subject-level random
## effects, and visit design. 
## =============================================================================

library(nimble)
library(splines2)
library(ggplot2)
library(tidyr)
library(dplyr)
library(mgcv)
library(MASS)

## -----------------------------------------------------------------------
## 1. True population log-rate spline + evaluation grid (independent of any observed data)
## -----------------------------------------------------------------------

simulate_true_rate <- function(YL = 0.30, YU = 1.70, K = 8, nGrid = 201,
                               shape = c("sigmoid", "custom"), theta_true = NULL) {
  shape <- match.arg(shape)
  ygrid <- seq(YL, YU, length.out = nGrid)
  interiorKnots <- seq(YL, YU, length.out = K - 2)[-c(1, K - 2)]
  Bgrid <- splines::bs(ygrid, knots = interiorKnots, Boundary.knots = c(YL, YU),
                       degree = 3, intercept = TRUE)
  Bgrid <- matrix(as.numeric(Bgrid), nrow = nGrid)
  Kactual <- ncol(Bgrid)
  
  if (is.null(theta_true)) {
    if (shape == "sigmoid") {
      # slow -> fast -> slow log-rate profile, loosely matching the RW1 prior
      theta_true <- seq(-9, -4.5, length.out = Kactual)
      theta_true <- theta_true + rev(cumsum(rev(c(0, diff(theta_true))) * 0.15))
    } else {
      stop("Supply theta_true explicitly when shape = 'custom'.")
    }
  }
  stopifnot(length(theta_true) == Kactual)
  
  Rgrid_true <- as.numeric(exp(Bgrid %*% theta_true))
  
  list(theta_true = theta_true, ygrid = ygrid, Bgrid = Bgrid, K = Kactual,
       Rgrid_true = Rgrid_true, YL = YL, YU = YU)
}

## -----------------------------------------------------------------------
## 2. Subject-level random effects and visit design
## -----------------------------------------------------------------------

simulate_design <- function(N = 200,
                            sigma_delta_true = 0.4,
                            sigma_eps_true = 0.05,
                            x0_range = c(0.9, 1.3),   # anchor SUVR value range
                            t0_range = c(50, 80),      # anchor age range (yrs)
                            visits_range = 2:6,
                            visit_prob = c(0.44, .27, .13, .10, .06),
                            visit_gap_mean = 1.5,       # years between visits
                            visit_gap_sd = 0.3) {
  
  delta <- rnorm(N, 0, sigma_delta_true)          # subject random intercept (log-scale mult.)
  x0 <- runif(N, x0_range[1], x0_range[2])        # anchor value
  t0 <- runif(N, t0_range[1], t0_range[2])        # anchor age
  
  Jvec <- sample(visits_range, N, prob = visit_prob, replace = TRUE)
  Jmax <- max(Jvec)
  
  tvisit <- matrix(NA_real_, N, Jmax)
  for (i in seq_len(N)) {
    Ji <- Jvec[i]
    # place anchor near the middle of the visit sequence, then space visits
    # forward/backward from it by ~visit_gap_mean years (matches the bidirectional integration)
    gaps <- rnorm(Ji - 1, visit_gap_mean, visit_gap_sd)
    gaps <- pmax(gaps, 0.25)
    offsets <- cumsum(c(0, gaps)) - cumsum(c(0, gaps))[ceiling(Ji / 2)]
    tvisit[i, 1:Ji] <- sort(t0[i] + offsets)
  }
  
  list(N = N, J = Jvec, Jmax = Jmax, tvisit = tvisit,
       delta_true = delta, x0_true = x0, t0_true = t0,
       sigma_delta_true = sigma_delta_true, sigma_eps_true = sigma_eps_true)
}