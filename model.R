## =============================================================================
##  Value-anchored ODE model for individual amyloid (SUVR) progression
##  CORRECTED / OPTIMISED VERSION -- see FIXES.md
##
##  Integrator semantics are ported from the parent project's validated
##  implementation (Claude DPM, scripts/btaj_value_anchored_amyloid.R sections 3
##  and its R mirrors) rather than reinvented. Three NIMBLE-specific facts that
##  project already established the hard way, and that this file respects:
##
##   * NIMBLE 1.4.2 does NOT compile `break` inside a `for` loop. The Newton
##     iteration below therefore uses a flag-controlled `while`.
##   * compileNimble() on stepImplicitTrap() STANDALONE fails with "Problem
##     with type of arg1 in sizeBinaryCwise". It compiles fine inside a
##     registered distribution. Hence tests/test_rmirror_vs_nimble.R exercises
##     it THROUGH a compiled model, never standalone.
##   * registerDistributions() resolves the density by ordinary R scoping, so
##     these functions must be defined at top level (they are).
## =============================================================================

library(nimble)
library(splines2)
library(coda)


## =============================================================================
##  1. RATE LOOKUP + INTEGRATOR (nimbleFunctions)
## =============================================================================
##
##  FIX (see FIXES.md issue #5): the original passed the whole `ygrid` vector
##  (201 doubles) into every rate evaluation and then recomputed
##  `length(ygrid)`, `ygrid[1]` and `ygrid[2]-ygrid[1]` inside the innermost
##  loop -- a loop executed on the order of 10^8-10^9 times per fit. The grid
##  is uniform by construction, so it is fully described by three SCALARS
##  (y_min, step_y, n_ygrid). Only the rate values `r_grid` need to be a vector.
##
##  The original also returned c(R, dR/dy) as a length-2 vector, allocating on
##  every call; splitting into two scalar-returning functions removes that.

lookupRate <- nimbleFunction(
  run = function(y = double(0), r_grid = double(1), y_min = double(0),
                 step_y = double(0), n_ygrid = double(0)) {
    returnType(double(0))
    declare(lo, integer(0))
    declare(hi, integer(0))
    pos <- (y - y_min) / step_y + 1          # fractional 1-based grid index
    if (pos < 1) pos <- 1                     # clamp: R is defined only on the grid
    if (pos > n_ygrid) pos <- n_ygrid
    lo <- floor(pos)
    hi <- lo + 1
    if (hi > n_ygrid) hi <- n_ygrid
    w <- pos - lo
    return((1 - w) * r_grid[lo] + w * r_grid[hi])
  }
)

## Derivative of the SAME piecewise-linear interpolant lookupRate evaluates --
## not an analytic B-spline derivative. Since R(y) is only ever accessed
## through the grid, the function the ODE actually integrates IS the linear
## interpolant, and Newton needs the derivative of the function it is really
## using: the local secant slope. Exact for the interpolant, not an
## approximation to it.
lookupRateDeriv <- nimbleFunction(
  run = function(y = double(0), r_grid = double(1), y_min = double(0),
                 step_y = double(0), n_ygrid = double(0)) {
    returnType(double(0))
    declare(lo, integer(0))
    declare(hi, integer(0))
    pos <- (y - y_min) / step_y + 1
    if (pos < 1) pos <- 1
    if (pos > n_ygrid) pos <- n_ygrid
    lo <- floor(pos)
    hi <- lo + 1
    if (hi > n_ygrid) hi <- n_ygrid
    if (hi == lo) return(0)
    return((r_grid[hi] - r_grid[lo]) / step_y)
  }
)

## Implicit trapezoid (Crank-Nicolson):
##     mu_next - mu_current - (h/2)*(R(mu_current) + R(mu_next)) = 0
## solved for mu_next by Newton, starting from the explicit-Euler guess.
##
## This equation is symmetric under swapping (mu_current, mu_next) with
## h -> -h, which is why the forward walk from a subject's reference point and
## the backward walk agree to machine precision -- the property the whole
## value-anchored construction depends on, and the reason this is the default
## rather than the cheaper explicit methods. Verified in
## tests/test_integrator_symmetry.R.
##
## FALLBACK: if Newton stalls (|F'| ~ 0 on a locally flat rate curve, or the
## iterate escapes a generous admissible range), fall back to fixed-point
## iteration, which needs no derivative and so cannot fail the same way. This
## guarantees the function never returns NaN into the MCMC.
stepImplicitTrap <- nimbleFunction(
  run = function(mu_current = double(0), h = double(0), r_grid = double(1),
                 y_min = double(0), step_y = double(0), n_ygrid = double(0),
                 exp_delta = double(0),
                 max_newton_iter = double(0, default = 10),
                 newton_tol = double(0, default = 1e-10)) {
    returnType(double(0))
    declare(iter, integer(0))
    declare(converged, integer(0))
    declare(k, integer(0))

    R_current <- lookupRate(mu_current, r_grid, y_min, step_y, n_ygrid) * exp_delta
    mu_euler  <- mu_current + h * R_current
    mu_next   <- mu_euler

    y_lo_bound <- y_min - 1
    y_hi_bound <- y_min + (n_ygrid - 1) * step_y + 1

    iter <- 0
    converged <- 0
    while (iter < max_newton_iter & converged == 0) {   # `break` does not compile in NIMBLE 1.4.2
      R_next <- lookupRate(mu_next, r_grid, y_min, step_y, n_ygrid) * exp_delta
      F_val  <- mu_next - mu_current - (h / 2) * (R_current + R_next)
      if (abs(F_val) < newton_tol) {
        converged <- 1
      } else {
        Rp_next <- lookupRateDeriv(mu_next, r_grid, y_min, step_y, n_ygrid) * exp_delta
        Fp_val  <- 1 - (h / 2) * Rp_next
        if (abs(Fp_val) < 1e-10 | mu_next < y_lo_bound | mu_next > y_hi_bound) {
          iter <- max_newton_iter          # bail out to the fixed-point fallback
        } else {
          mu_next <- mu_next - F_val / Fp_val
          iter <- iter + 1
        }
      }
    }

    if (converged == 0) {
      mu_next <- mu_euler
      for (k in 1:max_newton_iter) {
        R_next_fp <- lookupRate(mu_next, r_grid, y_min, step_y, n_ygrid) * exp_delta
        mu_next <- mu_current + (h / 2) * (R_current + R_next_fp)
      }
    }
    return(mu_next)
  }
)

## One integration step of dmu/dt = R(mu)*exp(delta).
##   method < 0.5   : explicit Euler        -- O(h^2) local, NOT time-symmetric
##   0.5 <= m < 1.5 : explicit Heun         -- O(h^3) local, NOT time-symmetric
##   method >= 1.5  : Crank-Nicolson        -- time-symmetric. THE DEFAULT.
## Every caller threads `method` from the single `integrator_method` constant,
## so switching schemes for a speed comparison is a one-line change.
stepValue <- nimbleFunction(
  run = function(mu = double(0), h = double(0), r_grid = double(1),
                 y_min = double(0), step_y = double(0), n_ygrid = double(0),
                 exp_delta = double(0), method = double(0),
                 max_newton_iter = double(0, default = 10),
                 newton_tol = double(0, default = 1e-10)) {
    returnType(double(0))
    if (method < 0.5) {
      r0 <- lookupRate(mu, r_grid, y_min, step_y, n_ygrid) * exp_delta
      return(mu + h * r0)
    }
    if (method < 1.5) {
      r0 <- lookupRate(mu, r_grid, y_min, step_y, n_ygrid) * exp_delta
      mu_pred <- mu + h * r0
      r1 <- lookupRate(mu_pred, r_grid, y_min, step_y, n_ygrid) * exp_delta
      return(mu + 0.5 * h * (r0 + r1))
    }
    return(stepImplicitTrap(mu, h, r_grid, y_min, step_y, n_ygrid, exp_delta,
                            max_newton_iter, newton_tol))
  }
)


## =============================================================================
##  2. PLAIN-R MIRRORS
## =============================================================================
##
##  Deliberately SEPARATE code from the nimbleFunctions above, not a wrapper --
##  so that a bug in one shows up as a disagreement rather than being silently
##  shared. tests/test_rmirror_vs_nimble.R checks they agree.
##
##  FIX: the original's R mirror ran a FIXED 6 Newton iterations while its
##  nimbleFunction exited early at |step| < 1e-8, so the two were not
##  guaranteed to agree and nothing checked them. These mirror the NIMBLE
##  convergence logic exactly.

lookup_rate_R <- function(y, r_grid, y_min, step_y, n_ygrid) {
  pos <- min(max((y - y_min) / step_y + 1, 1), n_ygrid)
  lo <- floor(pos); hi <- min(lo + 1, n_ygrid); w <- pos - lo
  (1 - w) * r_grid[lo] + w * r_grid[hi]
}

lookup_rate_deriv_R <- function(y, r_grid, y_min, step_y, n_ygrid) {
  pos <- min(max((y - y_min) / step_y + 1, 1), n_ygrid)
  lo <- floor(pos); hi <- min(lo + 1, n_ygrid)
  if (hi == lo) return(0)
  (r_grid[hi] - r_grid[lo]) / step_y
}

step_implicit_trap_R <- function(mu_current, h, r_grid, y_min, step_y, n_ygrid,
                                 exp_delta, max_newton_iter = 10, newton_tol = 1e-10) {
  R_current <- lookup_rate_R(mu_current, r_grid, y_min, step_y, n_ygrid) * exp_delta
  mu_euler <- mu_current + h * R_current
  mu_next  <- mu_euler
  y_lo_bound <- y_min - 1
  y_hi_bound <- y_min + (n_ygrid - 1) * step_y + 1
  converged <- FALSE
  for (iter in seq_len(max_newton_iter)) {
    R_next <- lookup_rate_R(mu_next, r_grid, y_min, step_y, n_ygrid) * exp_delta
    F_val  <- mu_next - mu_current - (h / 2) * (R_current + R_next)
    if (abs(F_val) < newton_tol) { converged <- TRUE; break }
    Rp_next <- lookup_rate_deriv_R(mu_next, r_grid, y_min, step_y, n_ygrid) * exp_delta
    Fp_val  <- 1 - (h / 2) * Rp_next
    if (abs(Fp_val) < 1e-10 || mu_next < y_lo_bound || mu_next > y_hi_bound) break
    mu_next <- mu_next - F_val / Fp_val
  }
  if (!converged) {
    mu_next <- mu_euler
    for (k in seq_len(max_newton_iter)) {
      R_fp <- lookup_rate_R(mu_next, r_grid, y_min, step_y, n_ygrid) * exp_delta
      mu_next <- mu_current + (h / 2) * (R_current + R_fp)
    }
  }
  mu_next
}

step_R <- function(mu, h, r_grid, y_min, step_y, n_ygrid, exp_delta, method,
                   max_newton_iter = 10, newton_tol = 1e-10) {
  if (method < 0.5) {
    return(mu + h * lookup_rate_R(mu, r_grid, y_min, step_y, n_ygrid) * exp_delta)
  }
  if (method < 1.5) {
    r0 <- lookup_rate_R(mu, r_grid, y_min, step_y, n_ygrid) * exp_delta
    r1 <- lookup_rate_R(mu + h * r0, r_grid, y_min, step_y, n_ygrid) * exp_delta
    return(mu + 0.5 * h * (r0 + r1))
  }
  step_implicit_trap_R(mu, h, r_grid, y_min, step_y, n_ygrid, exp_delta,
                       max_newton_iter, newton_tol)
}


## =============================================================================
##  3. PREPROCESSING -- everything deterministic, computed ONCE
## =============================================================================
##
##  FIX (FIXES.md issue #5): the original recomputed the sub-step count
##  ceiling(|dt|/maxSub) and the step size dt/nSteps INSIDE the model, on every
##  likelihood evaluation, for every interval, of every subject, on every MCMC
##  iteration. None of it depends on any parameter -- it is a function of the
##  visit times alone. It is all precomputed here instead.
##
##  Per subject this builds:
##    t0      -- reference time t_tilde = mean(visit ages)
##    yhat    -- kernel-weighted reference value (plain mean when J <= 2)
##    w_hat   -- the NORMALISED kernel weights, so yhat = sum(w_hat * y).
##               Passed into the likelihood so the model needs no separate
##               yhat constant -- see section 5.
##    n_eff   -- weighted effective sample size behind yhat
##    p       -- bracketing visit index: t[p] <= t_tilde <= t[p+1].
##               FIX: computed with findInterval(all.inside = TRUE), which
##               cannot return an out-of-range index, rather than the
##               original's max(which(age <= t0)).
##    M_ivl / step_ivl        -- sub-stepping for each whole visit-to-visit gap
##    M_left / step_left      -- partial span from t_tilde back to visit p
##    M_right / step_right    -- partial span from t_tilde out to visit p+1
prepare_amyloid <- function(dat, YL = 0.30, YU = 1.70, K = 10, step_y = 0.002,
                            Delta = 0.25, degree = 3) {

  dat <- dat[order(dat$id, dat$age), ]
  ids <- unique(dat$id)
  N <- length(ids)
  splitDat <- split(dat, factor(dat$id, levels = ids))
  Jvec <- vapply(splitDat, nrow, integer(1))
  Jmax <- max(Jvec)
  ## Floor of 2, not 1. If every subject has exactly 2 visits then Jmax - 1 = 1,
  ## and NIMBLE collapses the 1-wide slice M_ivl[i, 1:1] to a scalar, which
  ## fails checkBasics() against the distribution's double(1) argument type.
  ## The extra column is zero padding that no loop ever reaches (the interval
  ## loops are bounded by n_intervals[i] = J_i - 1), so this is free.
  ## Only ever binds for a degenerate all-J=2 design, e.g. the J-stratified
  ## diagnostic in FIXES.md; the study design has max_intervals = 6.
  max_intervals <- max(Jmax - 1L, 2L)

  tvisit <- matrix(NA_real_, N, Jmax)
  yobs   <- matrix(0,        N, Jmax)      # 0-padded: padded weights are 0 too
  t0 <- yhat <- neff <- bandwidth <- numeric(N)
  w_hat <- matrix(0, N, Jmax)
  pidx <- n_intervals <- integer(N)
  M_ivl <- step_ivl <- matrix(0, N, max_intervals)
  M_left <- step_left <- M_right <- step_right <- numeric(N)

  for (i in seq_len(N)) {
    di <- splitDat[[i]]
    Ji <- nrow(di)
    tt <- di$age; yy <- di$suvr
    tvisit[i, 1:Ji] <- tt
    yobs[i, 1:Ji]   <- yy
    n_intervals[i]  <- Ji - 1L

    t0[i] <- mean(tt)

    ## Reference value and its effective sample size.
    ## For J <= 2 the kernel is skipped: with one point it is undefined, and
    ## with two points t_tilde is the exact midpoint so the weights are
    ## symmetric for ANY bandwidth -- the weighted result is identically the
    ## plain mean. n_eff is then exactly J.
    if (Ji <= 2L) {
      w <- rep(1, Ji)
      bandwidth[i] <- NA_real_
    } else {
      h_i <- (tt[Ji] - tt[1]) / (Ji - 1)          # subject's own mean inter-visit interval
      bandwidth[i] <- h_i
      w <- exp(-(tt - t0[i])^2 / (2 * h_i^2))
    }
    yhat[i] <- sum(w * yy) / sum(w)
    neff[i] <- sum(w)^2 / sum(w^2)
    w_hat[i, 1:Ji] <- w / sum(w)

    ## Bracketing pair. all.inside = TRUE clamps to [1, Ji-1], so an exact
    ## endpoint coincidence yields a valid bracket with a zero-length partial
    ## span rather than an out-of-range index.
    p <- findInterval(t0[i], tt, all.inside = TRUE)
    pidx[i] <- p

    ## Whole visit-to-visit intervals.
    for (j in seq_len(Ji - 1L)) {
      dt <- max(tt[j + 1] - tt[j], 1e-6)          # guard duplicate timestamps
      Mj <- max(ceiling(dt / Delta), 1L)
      M_ivl[i, j] <- Mj
      step_ivl[i, j] <- dt / Mj
    }

    ## Partial spans out of the reference time, same M/step convention.
    dtl <- max(t0[i] - tt[p], 1e-6)
    Ml <- max(ceiling(dtl / Delta), 1L); M_left[i]  <- Ml; step_left[i]  <- dtl / Ml
    dtr <- max(tt[p + 1] - t0[i], 1e-6)
    Mr <- max(ceiling(dtr / Delta), 1L); M_right[i] <- Mr; step_right[i] <- dtr / Mr
  }

  ## ---- B-spline basis on the value grid ----------------------------------
  ## Interior knot count is not hand-picked: a basis with `degree` and
  ## intercept = TRUE has dimension n_interior + degree + 1, so
  ## n_interior = K - degree - 1 is what holds the basis dimension at K.
  ## Knots at data quantiles concentrate resolution where observations are.
  n_interior <- max(K - degree - 1L, 1L)
  interiorKnots <- as.numeric(quantile(dat$suvr,
                                       probs = seq(0.05, 0.95, length.out = n_interior),
                                       na.rm = TRUE))
  interiorKnots <- unique(pmin(pmax(interiorKnots, YL + 1e-6), YU - 1e-6))
  ygrid <- seq(YL, YU, by = step_y)
  n_ygrid <- length(ygrid)
  Bgrid <- splines2::bSpline(ygrid, knots = interiorKnots, degree = degree,
                             intercept = TRUE, Boundary.knots = c(YL, YU))
  Bgrid <- matrix(as.numeric(Bgrid), nrow = n_ygrid)
  Kactual <- ncol(Bgrid)

  ## Partition of unity -- the property that makes theta interpretable
  ## directly as log-rates. The original never checked this.
  stopifnot("B-spline basis is not a partition of unity" =
              all(abs(rowSums(Bgrid) - 1) < 1e-8))

  list(N = N, Jmax = Jmax, J = Jvec, ids = ids,
       tvisit = tvisit, y = yobs,
       t0 = t0, yhat = yhat, w_hat = w_hat, n = neff, bandwidth = bandwidth,
       p = pidx, n_intervals = n_intervals, max_intervals = max_intervals,
       M_ivl = M_ivl, step_ivl = step_ivl,
       M_left = M_left, step_left = step_left,
       M_right = M_right, step_right = step_right,
       ygrid = ygrid, n_ygrid = n_ygrid, step_y = step_y,
       Bgrid = Bgrid, K = Kactual, interiorKnots = interiorKnots,
       YL = YL, YU = YU, Delta = Delta)
}


## =============================================================================
## =============================================================================
##  4. CUSTOM DISTRIBUTION -- one node per subject
## =============================================================================
##
##  Fuses the trajectory integration and the likelihood accumulation into ONE
##  node per subject. The original built an N x Jmax deterministic mu node plus
##  J[i] separate dnorm nodes -- roughly 2000 graph nodes at N = 200, all
##  re-traversed on every one of the hundreds of log-density evaluations an
##  AF_slice sweep over theta performs per iteration.
##
##  The anchor value x_start arrives as a PLAIN PARAMETER. It is sampled from a
##  population prior in the model code below -- it is NOT rebuilt from this
##  subject's own observations. See CHANGES.md, change 1.
##
##  The walk: one partial sub-interval from the reference time out to the
##  nearest visit on each side, then the ordinary visit-to-visit intervals
##  further out. The backward direction is an ASCENDING k-loop with a computed
##  descending index, because NIMBLE's `:` operator does not reliably count down.
dsubjRefAnchor <- nimbleFunction(
  run = function(x = double(1), n_intervals = double(0), p_hi = double(0),
                 M_ivl = double(1), step_ivl = double(1),
                 M_left = double(0), step_left = double(0),
                 M_right = double(0), step_right = double(0),
                 x_start = double(0), delta = double(0), sigma_eps = double(0),
                 r_grid = double(1), y_min = double(0), step_y = double(0),
                 n_ygrid = double(0), method = double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    declare(nJ, integer(0)); declare(pH, integer(0)); declare(Mj, integer(0))
    declare(Ml, integer(0)); declare(Mr, integer(0)); declare(j, integer(0))

    nJ <- round(n_intervals)
    pH <- round(p_hi)
    ll <- 0
    exp_delta <- exp(delta)

    if (nJ >= 1) {
      ## ---- forward: partial to visit p+1, then whole intervals outward ----
      Mr <- round(M_right)
      mu_f <- x_start
      for (m in 1:Mr) {
        mu_f <- stepValue(mu_f, step_right, r_grid, y_min, step_y, n_ygrid, exp_delta, method)
      }
      ll <- ll + dnorm(x[pH + 1], mean = mu_f, sd = sigma_eps, log = TRUE)
      if (pH < nJ) {
        for (j in (pH + 1):nJ) {
          Mj <- round(M_ivl[j])
          for (m in 1:Mj) {
            mu_f <- stepValue(mu_f, step_ivl[j], r_grid, y_min, step_y, n_ygrid, exp_delta, method)
          }
          ll <- ll + dnorm(x[j + 1], mean = mu_f, sd = sigma_eps, log = TRUE)
        }
      }

      ## ---- backward: partial to visit p, then whole intervals inward ------
      Ml <- round(M_left)
      mu_b <- x_start
      for (m in 1:Ml) {
        mu_b <- stepValue(mu_b, -step_left, r_grid, y_min, step_y, n_ygrid, exp_delta, method)
      }
      ll <- ll + dnorm(x[pH], mean = mu_b, sd = sigma_eps, log = TRUE)
      if (pH >= 2) {
        for (k in 1:(pH - 1)) {
          j <- pH - k
          Mj <- round(M_ivl[j])
          for (m in 1:Mj) {
            mu_b <- stepValue(mu_b, -step_ivl[j], r_grid, y_min, step_y, n_ygrid, exp_delta, method)
          }
          ll <- ll + dnorm(x[j], mean = mu_b, sd = sigma_eps, log = TRUE)
        }
      }
    }

    if (log) return(ll) else return(exp(ll))
  }
)

## Random-draw companion NIMBLE requires for registerDistributions().
## Never invoked: y is always supplied as observed data.
rsubjRefAnchor <- nimbleFunction(
  run = function(n = integer(0), n_intervals = double(0), p_hi = double(0),
                 M_ivl = double(1), step_ivl = double(1),
                 M_left = double(0), step_left = double(0),
                 M_right = double(0), step_right = double(0),
                 x_start = double(0), delta = double(0), sigma_eps = double(0),
                 r_grid = double(1), y_min = double(0), step_y = double(0),
                 n_ygrid = double(0), method = double(0)) {
    returnType(double(1))
    out <- numeric(length(M_ivl) + 1)
    return(out)
  }
)

registerDistributions(list(
  dsubjRefAnchor = list(
    BUGSdist = paste0("dsubjRefAnchor(n_intervals, p_hi, M_ivl, step_ivl, ",
                      "M_left, step_left, M_right, step_right, x_start, ",
                      "delta, sigma_eps, r_grid, y_min, step_y, n_ygrid, method)"),
    types = c("value = double(1)", "M_ivl = double(1)", "step_ivl = double(1)",
              "r_grid = double(1)"),
    discrete = FALSE, pqAvail = FALSE
  )
), verbose = FALSE)


## =============================================================================
##  5. MODEL
## =============================================================================
amyloidOdeCode <- nimbleCode({

  ## ---- population log-rate spline, random-walk prior with local scales ----
  ## dgamma(1, 1), NOT dgamma(1e-3, 1e-3): the latter is the classic
  ## pathological variance-component prior and piles mass at either extreme
  ## under weak data. (The original also named it as a variance while using it
  ## as a precision.)
  tau_theta ~ dgamma(1, 1)
  theta[1] ~ dnorm(-9, sd = 2)
  for (k in 2:K) {
    lambda[k] ~ dexp(1)
    theta[k] ~ dnorm(theta[k - 1], tau = tau_theta / lambda[k])
  }

  ## ONE vectorised node, not 402 scalar ones.
  Rpop[1:nGrid] <- exp(Bgrid[1:nGrid, 1:K] %*% theta[1:K])

  sigma_delta ~ T(dnorm(0, sd = 1), 0, )
  sigma_eps   ~ T(dnorm(0, sd = 0.1), 0, )

  ## ---- POPULATION PRIOR ON THE ANCHOR (CHANGES.md, change 1) --------------
  ## m_x and tau_x are estimated, so each subject's anchor is shrunk toward a
  ## population value learned from ALL subjects jointly -- not toward that
  ## subject's own observations. The data therefore enter the density exactly
  ## once, which is what makes this a posterior rather than a pseudo-posterior.
  ##
  ## The hyperpriors are fixed constants, deliberately NOT derived from the
  ## data. With N in the hundreds or thousands both hyperparameters are
  ## strongly identified, so weakly-informative is enough.
  m_x   ~ dnorm(1.0, sd = 0.5)
  tau_x ~ T(dnorm(0, sd = 0.5), 0, )

  for (i in 1:N) {
    x_tilde[i] ~ dnorm(m_x, sd = tau_x)
    delta[i]   ~ dnorm(0, sd = sigma_delta)

    y[i, 1:Jm] ~ dsubjRefAnchor(
      n_intervals = n_intervals[i], p_hi = p[i],
      M_ivl = M_ivl[i, 1:max_intervals], step_ivl = step_ivl[i, 1:max_intervals],
      M_left = M_left[i], step_left = step_left[i],
      M_right = M_right[i], step_right = step_right[i],
      x_start = x_tilde[i], delta = delta[i], sigma_eps = sigma_eps,
      r_grid = Rpop[1:nGrid], y_min = YL, step_y = step_y_c,
      n_ygrid = nGrid, method = integrator_method
    )
  }
})


## x_tilde is now a SAMPLED node, so it is monitored directly. This helper is
## kept so downstream code that used to reconstruct it keeps working.
reconstruct_x_tilde <- function(samples, ph) {
  cols <- paste0("x_tilde[", seq_len(ph$N), "]")
  stopifnot("x_tilde was not monitored" = all(cols %in% colnames(samples)))
  samples[, cols, drop = FALSE]
}

## =============================================================================
##  6. AGE OF POSITIVITY
## =============================================================================
##
##  FIX (FIXES.md issue #6): the original called this per posterior draw PER
##  SUBJECT, and rebuilt the 701 x K rate-grid matrix product inside the
##  subject loop even though it only varies by draw. With 15,000 draws x 6
##  subjects that is ~90,000 redundant matrix products and ~10^8 R-level
##  stepping calls -- longer than the MCMC itself.
##
##  Here the r_grid is built ONCE per draw and shared across all subjects, and
##  the stepping is vectorised ACROSS SUBJECTS (all subjects advance together),
##  so the cost is (n_draws x n_steps) vector operations rather than
##  (n_draws x n_subjects x n_steps) scalar ones.

## Vectorised over subjects: one r_grid, many (x0, t0, delta).
positivity_age_vec <- function(x0, t0, delta, thres, r_grid, y_min, step_y,
                               n_ygrid, Delta = 0.25, method = 2, max_years = 60) {
  n <- length(x0)
  out <- rep(NA_real_, n)
  exp_delta <- exp(delta)
  dir <- ifelse(x0 < thres, 1, -1)
  mu <- x0
  t  <- t0
  active <- rep(TRUE, n)
  nstep <- ceiling(max_years / Delta)

  ## Vectorised rate lookup over the currently-active subjects.
  lk <- function(y) {
    pos <- pmin(pmax((y - y_min) / step_y + 1, 1), n_ygrid)
    lo <- floor(pos); hi <- pmin(lo + 1, n_ygrid); w <- pos - lo
    (1 - w) * r_grid[lo] + w * r_grid[hi]
  }

  for (s in seq_len(nstep)) {
    if (!any(active)) break
    idx <- which(active)
    h <- dir[idx] * Delta
    ## Heun step, vectorised -- the crossing search only needs the trajectory
    ## to locate a sign change, and is then refined by linear interpolation.
    r0 <- lk(mu[idx]) * exp_delta[idx]
    r1 <- lk(mu[idx] + h * r0) * exp_delta[idx]
    mu_new <- mu[idx] + 0.5 * h * (r0 + r1)
    t_new  <- t[idx] + h

    crossed <- (dir[idx] ==  1 & mu[idx] <  thres & mu_new >= thres) |
               (dir[idx] == -1 & mu[idx] >  thres & mu_new <= thres)
    if (any(crossed)) {
      ci <- idx[crossed]
      frac <- (thres - mu[ci]) / (mu_new[crossed] - mu[ci])
      out[ci] <- t[ci] + frac * (t_new[crossed] - t[ci])
      active[ci] <- FALSE
    }
    mu[idx] <- mu_new
    t[idx]  <- t_new
  }
  out
}

## Posterior distribution of the positivity age, thinned. Returns a
## draws x subjects matrix.
positivity_age_posterior <- function(samples, ph, thres, subjects = NULL,
                                     thin_to = 200, method = 2) {
  if (is.null(subjects)) subjects <- seq_len(ph$N)
  nd <- nrow(samples)
  idx <- if (nd > thin_to) round(seq(1, nd, length.out = thin_to)) else seq_len(nd)

  theta_cols <- paste0("theta[", seq_len(ph$K), "]")
  delta_cols <- paste0("delta[", subjects, "]")
  x_tilde <- reconstruct_x_tilde(samples, ph)[, subjects, drop = FALSE]

  out <- matrix(NA_real_, length(idx), length(subjects))
  for (r in seq_along(idx)) {
    d <- idx[r]
    r_grid <- as.numeric(exp(ph$Bgrid %*% samples[d, theta_cols]))   # ONCE per draw
    out[r, ] <- positivity_age_vec(
      x0 = x_tilde[d, ], t0 = ph$t0[subjects],
      delta = samples[d, delta_cols], thres = thres,
      r_grid = r_grid, y_min = ph$YL, step_y = ph$step_y,
      n_ygrid = ph$n_ygrid, Delta = ph$Delta, method = method
    )
  }
  colnames(out) <- paste0("subject_", subjects)
  out
}
