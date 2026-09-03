## =============================================================================
##  RECOVERY FIGURES
##
##  Seven figures comparing what the model recovered against the known truth.
##  Called automatically by run_study.R; you can also call make_figures() on a
##  saved study object -- but the trajectory figure additionally needs
##  truth_mu_at() (sim_data.R) and step_R() (model.R) in scope; make_figures()
##  quietly drops it when they are absent.
##
##  Points are coloured by number of visits (J) throughout, because that is the
##  quantity that governs how much any individual subject can be recovered --
##  see FIXES.md section 4.
## =============================================================================

library(ggplot2)

## One-hue sequential ramp for J (an ordered quantity), ink for truth.
SEQ_J <- c("#b7d3f6", "#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281")
EST   <- "#2a78d6"
INK   <- "#0b0b0b"; INK2 <- "#52514e"; GRIDC <- "#e4e4e0"

theme_study <- function(base = 10) {
  theme_minimal(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = GRIDC, linewidth = 0.3),
          axis.text = element_text(colour = INK2),
          axis.title = element_text(colour = INK2),
          plot.title = element_text(colour = INK, face = "bold"),
          plot.subtitle = element_text(colour = INK2),
          legend.position = "bottom",
          legend.text = element_text(colour = INK2))
}


## -----------------------------------------------------------------------
##  Pull everything the figures need out of one or more replicates
## -----------------------------------------------------------------------
gather_results <- function(study) {
  reps <- study$replicates
  rate <- do.call(rbind, lapply(seq_along(reps), function(i) {
    r <- reps[[i]]$recovery$rate
    data.frame(rep = i, y = r$ygrid, est = r$mean, lo = r$lo, hi = r$hi,
               truth = r$truth)
  }))
  subj <- do.call(rbind, lapply(seq_along(reps), function(i) {
    rc <- reps[[i]]$recovery; ph <- reps[[i]]$ph
    data.frame(rep = i, id = seq_len(ph$N), J = ph$J,
               x_truth = rc$x_tilde$truth, x_est = rc$x_tilde$post_mean,
               x_lo = rc$x_tilde$lo, x_hi = rc$x_tilde$hi,
               d_truth = rc$delta$truth, d_est = rc$delta$post_mean,
               d_lo = rc$delta$lo, d_hi = rc$delta$hi,
               a_truth = rc$alpha$truth, a_est = rc$alpha$post_median,
               a_lo = rc$alpha$lo, a_hi = rc$alpha$hi)
  }))
  subj$Jf <- factor(subj$J)
  list(rate = rate, subj = subj, n_rep = length(reps))
}

lab_n <- function(G) if (G$n_rep == 1)
  sprintf("%d subjects", length(unique(G$subj$id))) else
  sprintf("%d replicates, %d subjects each", G$n_rep, length(unique(G$subj$id)))


## -----------------------------------------------------------------------
##  1. Population rate curve
## -----------------------------------------------------------------------
fig_rate_curve <- function(G) {
  sup <- quantile(G$subj$x_truth, c(0.01, 0.99))
  p <- ggplot(G$rate, aes(y)) +
    annotate("rect", xmin = sup[1], xmax = sup[2], ymin = -Inf, ymax = Inf,
             fill = GRIDC, alpha = 0.55)
  if (G$n_rep == 1) {
    p <- p + geom_ribbon(aes(ymin = lo, ymax = hi), fill = EST, alpha = 0.25)
  }
  p + geom_line(aes(y = est, group = rep, colour = "Posterior mean"),
                linewidth = 0.6, alpha = if (G$n_rep > 1) 0.8 else 1) +
    geom_line(aes(y = truth, colour = "Truth"), linewidth = 0.9) +
    scale_colour_manual(values = c(Truth = INK, `Posterior mean` = EST), name = NULL) +
    labs(x = "SUVR (y)", y = "R(y)   (SUVR / year)",
         title = "Population rate curve: recovered vs truth",
         subtitle = paste0(
           if (G$n_rep == 1) "Band: 95% pointwise credible interval. " else "One line per replicate. ",
           "Shaded region: 1%-99% of true anchor values;\noutside it the curve is prior extrapolation, not an estimate. ",
           lab_n(G))) +
    theme_study()
}


## -----------------------------------------------------------------------
##  2-4. Per-subject recovery against truth
## -----------------------------------------------------------------------
fig_recovery <- function(G, xv, yv, lo, hi, title, axlab, unit = "") {
  d <- G$subj[is.finite(G$subj[[xv]]) & is.finite(G$subj[[yv]]), ]
  rng <- range(c(d[[xv]], d[[yv]]))
  lab <- sprintf("r = %.3f    95%% CI coverage %.0f%%    median |error| %.3f%s",
                 cor(d[[xv]], d[[yv]]),
                 100 * mean(d[[xv]] >= d[[lo]] & d[[xv]] <= d[[hi]]),
                 median(abs(d[[yv]] - d[[xv]])), unit)
  ggplot(d, aes(.data[[xv]], .data[[yv]])) +
    geom_abline(slope = 1, intercept = 0, colour = INK, linewidth = 0.4,
                linetype = "dashed") +
    geom_point(aes(colour = Jf), size = 0.6, alpha = 0.35) +
    annotate("text", x = rng[1], y = rng[2], label = lab, hjust = 0, vjust = 1,
             size = 3, colour = INK2) +
    scale_colour_manual(values = SEQ_J, name = "visits (J)") +
    guides(colour = guide_legend(override.aes = list(size = 2.6, alpha = 1), nrow = 1)) +
    coord_equal(xlim = rng, ylim = rng) +
    labs(x = paste("true", axlab), y = paste("posterior", axlab), title = title,
         subtitle = paste0("Dashed line = perfect recovery. ", lab_n(G))) +
    theme_study()
}


## -----------------------------------------------------------------------
##  5. How recovery depends on follow-up
## -----------------------------------------------------------------------
fig_by_visits <- function(G) {
  d <- G$subj
  tab <- do.call(rbind, lapply(split(d, d$J), function(z) data.frame(
    J = z$J[1], n = nrow(z),
    delta = if (nrow(z) > 2) cor(z$d_truth, z$d_est) else NA_real_,
    onset = if (nrow(z) > 2) cor(z$a_truth, z$a_est, use = "complete.obs") else NA_real_)))
  long <- rbind(data.frame(J = tab$J, n = tab$n, r = tab$delta,
                           what = "delta_i  (rate multiplier)"),
                data.frame(J = tab$J, n = tab$n, r = tab$onset,
                           what = "disease age (onset)"))
  ggplot(long, aes(factor(J), r, group = what, colour = what)) +
    geom_line(linewidth = 0.7) + geom_point(size = 2.2) +
    geom_text(aes(label = sprintf("n=%d", n)), y = 0.02, colour = INK2,
              size = 2.5, show.legend = FALSE) +
    scale_colour_manual(values = c("#2a78d6", "#eb6834"), name = NULL) +
    ylim(0, 1) +
    labs(x = "number of visits (J)", y = "correlation with truth",
         title = "Recovery improves with follow-up -- but not equally",
         subtitle = paste0("delta_i needs an individual slope, which short follow-up cannot supply.\n",
                           "Disease age leans on the anchor and the shared rate curve, so it survives. ",
                           lab_n(G))) +
    theme_study()
}


## -----------------------------------------------------------------------
##  6. Traceplots (single replicate)
## -----------------------------------------------------------------------
fig_traces <- function(fit, params = c("sigma_delta", "sigma_eps", "m_x", "tau_x",
                                       "theta[1]", "theta[5]", "theta[10]", "tau_theta")) {
  params <- intersect(params, colnames(fit$samples))
  df <- do.call(rbind, lapply(params, function(p)
    data.frame(iter = seq_len(nrow(fit$samples)), value = fit$samples[, p], param = p)))
  ggplot(df, aes(iter, value)) +
    geom_line(linewidth = 0.22, colour = EST) +
    facet_wrap(~param, scales = "free_y") +
    labs(x = "post-burn-in iteration", y = NULL, title = "Traceplots",
         subtitle = "Single chain: these show mixing, but cannot rule out a chain stuck in a wrong mode") +
    theme_study(9)
}


## -----------------------------------------------------------------------
##  7. Individual trajectories: truth vs posterior mean (single replicate)
## -----------------------------------------------------------------------
##  The truth curve comes from the EXACT separable solve (truth_mu_at), not
##  from re-running the estimator's own integrator, so truth and estimate are
##  genuinely independent here. The estimate integrates the posterior-mean
##  rate curve with the model's own stepper (step_R) outward from the
##  reference time t0, using the posterior-mean anchor and delta.
fig_trajectories <- function(rp, design, ids = NULL, pad = 2) {
  ph <- rp$ph; rc <- rp$recovery; sim <- rp$sim
  R_post <- rc$rate$mean

  ## Default: one subject at each distinct follow-up length, up to six, so the
  ## panel shows both the well-observed and the barely-observed cases.
  if (is.null(ids))
    ids <- head(vapply(sort(unique(ph$J)),
                       function(j) which(ph$J == j)[1], integer(1)), 6)

  traj <- do.call(rbind, lapply(ids, function(i) {
    ages <- ph$tvisit[i, 1:ph$J[i]]
    tt   <- seq(min(ages) - pad, max(ages) + pad, length.out = 120)

    mu_true <- truth_mu_at(sim$solver, design$x0_true[i], design$t0_true[i],
                           tt, sim$delta_true[i])

    x_est <- rc$x_tilde$post_mean[i]; d_est <- rc$delta$post_mean[i]
    mu_est <- vapply(tt, function(tg) {
      dt <- tg - ph$t0[i]; M <- max(ceiling(abs(dt) / ph$Delta), 1L); h <- dt / M
      mu <- x_est
      for (m in seq_len(M))
        mu <- step_R(mu, h, R_post, ph$YL, ph$step_y, ph$n_ygrid, exp(d_est), 2)
      mu
    }, numeric(1))

    rbind(data.frame(id = i, t = tt, mu = mu_true, curve = "Truth"),
          data.frame(id = i, t = tt, mu = mu_est,  curve = "Posterior mean"))
  }))

  obs <- do.call(rbind, lapply(ids, function(i)
    data.frame(id = i, t = ph$tvisit[i, 1:ph$J[i]], y = ph$y[i, 1:ph$J[i]])))

  ## Anchor x_i at the reference time t0
  anc <- do.call(rbind, lapply(ids, function(i) rbind(
    data.frame(id = i, t = design$t0_true[i], mu = design$x0_true[i], curve = "Truth"),
    data.frame(id = i, t = ph$t0[i], mu = rc$x_tilde$post_mean[i], curve = "Posterior mean"))))

  ## Per-panel note: acceleration factor exp(delta), true vs recovered.
  acc_lab <- data.frame(id = ids, txt = sprintf(
    "exp(delta):  true %.2f   est %.2f",
    exp(sim$delta_true[ids]), exp(rc$delta$post_mean[ids])))

  ## Facet label carries the follow-up length, the quantity that governs how
  ## much any one subject can be recovered.
  id_lab <- sprintf("subject %d  (J = %d)", ids, ph$J[ids])
  as_id  <- function(d) { d$id <- factor(d$id, levels = ids, labels = id_lab); d }
  traj <- as_id(traj); obs <- as_id(obs); anc <- as_id(anc); acc_lab <- as_id(acc_lab)

  ggplot(traj, aes(t, mu, colour = curve)) +
    geom_line(linewidth = 0.7) +
    geom_point(data = obs, aes(t, y), inherit.aes = FALSE, size = 1.3,
               colour = INK2, alpha = 0.85) +
    geom_point(data = anc, shape = 23, size = 2.3, stroke = 0.7, fill = "white") +
    ## bottom-right: trajectories rise left-to-right, so that corner stays clear
    geom_text(data = acc_lab, aes(x = Inf, y = -Inf, label = txt), inherit.aes = FALSE,
              hjust = 1.04, vjust = -0.8, size = 2.7, colour = INK2) +
    scale_colour_manual(values = c(Truth = INK, `Posterior mean` = EST), name = NULL) +
    facet_wrap(~id, scales = "free") +
    labs(x = "age (years)", y = "SUVR",
         title = "Individual trajectories: truth vs posterior mean",
         subtitle = paste0("One subject per distinct follow-up length. Diamonds mark the anchor ",
                           "x_i (where integration begins);\npoints are the observed (noisy) values the estimator saw.")) +
    theme_study(9)
}


## -----------------------------------------------------------------------
##  Write everything
## -----------------------------------------------------------------------
make_figures <- function(study, out_dir = "figures") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  G <- gather_results(study)
  r1 <- study$replicates[[1]]

  figs <- list(
    rate_curve  = fig_rate_curve(G),
    anchor      = fig_recovery(G, "x_truth", "x_est", "x_lo", "x_hi",
                               "Anchor value x_i: posterior mean vs truth", "anchor SUVR"),
    disease_age = fig_recovery(G, "a_truth", "a_est", "a_lo", "a_hi",
                               "Disease age (age of amyloid positivity) vs truth",
                               "onset age (yr)", " yr"),
    slope       = fig_recovery(G, "d_truth", "d_est", "d_lo", "d_hi",
                               "Slope multiplier delta_i: posterior mean vs truth",
                               "delta (log rate multiplier)"),
    by_visits   = fig_by_visits(G),
    traces      = fig_traces(r1$fit),
    trajectories = fig_trajectories(r1, study$design)
  )

  pdf(file.path(out_dir, "recovery_figures.pdf"), width = 9, height = 5.4)
  for (f in figs) print(f)
  invisible(dev.off())
  for (nm in names(figs))
    ggsave(file.path(out_dir, paste0(nm, ".png")), figs[[nm]],
           width = 9, height = 5.4, dpi = 130)

  cat(sprintf("\nWrote %s/recovery_figures.pdf and %d PNGs:\n  %s\n",
              out_dir, length(figs), paste(names(figs), collapse = ", ")))
  invisible(figs)
}
