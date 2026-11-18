# =============================================================================
# ANALYSIS FUNCTIONS — describe agreement and estimate conditional repeatability
# =============================================================================
# These functions assume repeated biological measurements, not independent cells
# from the same culture. They do not establish persistence mechanisms or GWAS
# power. Read the README for censoring, confounding and independence limitations.
# The exported ICCs are within-plate quantities, not genetic heritability.
# =============================================================================
# Check data identity and design before transformation. kind determines the
# expected count column. On success, return the unchanged data invisibly; the
# function is a guard, not a cleaning operation.
validate_observations <- function(d, kind) {
  value <- if (kind == "recovery") "recovery_cfu_ml" else "initial_cfu_ml"
# Identifier fields must not be blank. Missing numeric counts are allowed
  # and are audited later, but missing keys make biological matching ambiguous.
  ids <- c("isolate_id", "date", "replicate", if (kind == "recovery") "plate")
  if (any(vapply(d[ids], function(z) anyNA(z) || any(trimws(as.character(z)) == ""), logical(1)))) {
    stop("Missing identifier in ", kind)
  }
  if (any(!grepl("^R[0-9]+$", d$replicate))) stop("Replicates must be R1, R2, ...")
# A sample/date/repeat key identifies one observation. Do not let joins or
  # models silently count technical replicates as extra biological evidence.
  if (anyDuplicated(d[c("isolate_id", "date", "replicate")])) {
    stop("Duplicate isolate/date/replicate in ", kind, "; resolve technical replicates first.")
  }
  if (any(d[[value]] < 0, na.rm = TRUE)) stop("Negative CFU values in ", kind)
  if (kind == "recovery") {
# This template assumes each isolate belongs to one plate. A crossed design
    # is legitimate scientifically but needs a different model and identifier plan.
    assignment <- unique(d[c("isolate_id", "plate")])
    if (anyDuplicated(assignment$isolate_id)) stop("Isolate IDs must be globally unique and assigned to one plate.")
    if (anyDuplicated(d[c("isolate_id", "replicate")])) {
      stop("Multiple dates for one isolate/replicate; use unambiguous biological repeat labels.")
    }
  }
  invisible(d)
}

# Preserve all rows and add an exclusion reason plus a log10 column.
# zero_policy="stop" demands clarification before modelling; "exclude" is an
# explicit preliminary analysis of positive observations, not a censoring model.
prepare_recovery <- function(d, zero_policy = "stop") {
  validate_observations(d, "recovery")
  if (!zero_policy %in% c("stop", "exclude")) stop("Unknown zero_policy.")
  if (any(d$recovery_cfu_ml == 0, na.rm = TRUE) && zero_policy == "stop") {
    stop("Zero counts present. Resolve detection limits or explicitly choose preliminary exclusion.")
  }
# NA counts and measured zeros have different meanings and separate labels.
  # For both, leave log10_recovery as NA rather than generating -Inf or adding 1.
  d$exclusion_reason <- ifelse(is.na(d$recovery_cfu_ml), "missing_count",
                               ifelse(d$recovery_cfu_ml == 0, "zero_count", NA_character_))
  d$log10_recovery <- NA_real_
# Here ok marks rows WITH an exclusion reason. The complement !ok is the
  # set eligible for log transformation. Original CFU values remain unchanged.
  ok <- !is.na(d$exclusion_reason)
  if (any(ok)) warning(sum(ok), " missing/zero observations excluded from log-scale summaries and models.")
  d$log10_recovery[!ok] <- log10(d$recovery_cfu_ml[!ok])
# Order labels numerically: R2 must precede R10. Factor levels preserve this
  # intended order in plots and pairwise comparisons.
  reps <- unique(d$replicate)
  reps <- reps[order(as.integer(sub("^R", "", reps)))]
  d$replicate <- factor(d$replicate, levels = reps)
  d$plate <- factor(d$plate, levels = unique(d$plate))
  d
}

# Convert the long table temporarily to one row per isolate, one column per
# repeat. Then generate every two-repeat combination rather than hard-coding
# R1/R2/R3. Return only pairs with two finite log-scale measurements.
make_pairs <- function(d) {
  wide <- tidyr::pivot_wider(
    dplyr::select(d, plate, isolate_id, replicate, log10_recovery),
    names_from = replicate, values_from = log10_recovery)
  reps <- levels(d$replicate)
  if (length(reps) < 2L) stop("At least two repeat labels are required.")
# combn enumerates repeat pairs once. x is the earlier-numbered repeat and y
  # the later-numbered repeat, so y-x has a consistent signed interpretation.
  pairs <- combn(reps, 2L, simplify = FALSE)
  dplyr::bind_rows(lapply(pairs, function(p) {
    data.frame(plate = wide$plate, isolate_id = wide$isolate_id,
               comparison = paste(p, collapse = " vs "),
               x = wide[[p[1]]], y = wide[[p[2]]])
  })) |>
    dplyr::filter(is.finite(x), is.finite(y))
}

# Extract variance components from an lme4 random-intercept model.
# isolate_key must name the grouping factor in that fitted model.
# Formula: ICC = between-isolate variance / total conditional variance.
# A high ICC means consistent differences are large relative to residual
# fluctuations; it says nothing directly about which differences are genetic.
get_repeatability <- function(model) {
  vi <- as.numeric(lme4::VarCorr(model)$isolate_key[1, 1])
  ve <- stats::sigma(model)^2
  data.frame(isolate_variance = vi, residual_variance = ve,
             repeatability = vi / (vi + ve))
}

# Fit the same eligible observations with two fixed-effect specifications.
# Returns fitted models plus their data for diagnostics; no adjusted GWAS
# phenotype is constructed automatically.
fit_models <- function(d) {
  m <- d[is.finite(d$log10_recovery), , drop = FALSE]
  m$isolate_key <- factor(m$isolate_id)
  m$plate <- droplevels(m$plate)
# A run is provisionally one plate on one date. Change this definition if the
  # real experiment includes sessions within a day or a shared paired-plate run.
  m$run_id <- interaction(m$plate, m$date, drop = TRUE)
  if (nlevels(m$plate) < 2L) stop("This workflow requires at least two plates.")
  if (nlevels(m$isolate_key) >= nrow(m)) stop("No repeated measurements remain.")
  per_plate <- tapply(as.character(m$date), m$plate, function(z) length(unique(z)))
  if (any(per_plate < 2L)) stop("Each plate must have at least two measured dates.")
# Primary model: subtract plate means before estimating isolate variance.
  # Random intercepts represent consistent isolate deviations; residuals combine
  # unexplained biological fluctuations and measurement variation. REML estimates
  # variance components while allowing for fitted fixed effects.
  a <- lme4::lmer(log10_recovery ~ plate + (1 | isolate_key), data = m, REML = TRUE)
# Diagnostic model: give each plate/date a separate mean. The "0 +" syntax
  # fits all run means directly instead of an intercept plus reference contrasts.
  # This ICC is conditional on observed run means and may not describe agreement
  # on a future, unseen day. A small ICC change does not prove no batch effects.
  b <- lme4::lmer(log10_recovery ~ 0 + run_id + (1 | isolate_key), data = m, REML = TRUE)
  list(before = a, adjusted = b, data = m)
}

# Estimate a model-based percentile interval for the PRIMARY model ICC.
# Inputs: fitted model, number of simulations and reproducible random seed.
# Output: interval/status table, all draws, and bootstrap refitting messages.
# This is not resampling raw rows; rows from one isolate are not independent.
bootstrap_icc <- function(model, nsim, seed) {
  if (length(nsim) != 1L || !is.finite(nsim) || nsim < 2L || nsim != as.integer(nsim)) {
    stop("bootstrap_n must be an integer >= 2.")
  }
# type="parametric" simulates outcomes from the fitted Gaussian model.
  # use.u=FALSE simulates new isolate effects as well as new residual errors.
  # Every simulated dataset is refitted and its ICC is extracted.
  b <- lme4::bootMer(model, FUN = function(fit) get_repeatability(fit)$repeatability,
                     nsim = nsim, seed = seed, type = "parametric", use.u = FALSE)
  vals <- as.numeric(b$t[, 1])
# Retain the full draw vector, including failed/non-finite entries, for audit.
  # Only finite entries enter quantiles. A computed interval should not be
  # trusted automatically when many refits fail or convergence is problematic.
  good <- is.finite(vals)
  ci <- if (sum(good) >= 2L) stats::quantile(vals[good], c(.025, .975), names = FALSE) else c(NA, NA)
  list(summary = data.frame(
    lower_95 = ci[1], upper_95 = ci[2], requested = nsim,
    finite_estimates = sum(good), failed_refits = if (is.null(attr(b, "bootFail"))) 0L else attr(b, "bootFail")),
    draws = data.frame(draw = seq_along(vals), repeatability = vals),
    messages = attr(b, "boot.all.msgs"))
}

# Record numerical fit diagnostics independently of the scientific estimate.
# A finite ICC does not guarantee successful optimisation or adequate assumptions.
model_status <- function(model, label) {
  msgs <- model@optinfo$conv$lme4$messages
  data.frame(model = label, singular = lme4::isSingular(model),
             optimizer_code = paste(model@optinfo$conv$opt, collapse = "; "),
             convergence_messages = if (is.null(msgs)) "none" else paste(msgs, collapse = "; "))
}

# Save three types of diagnostic: trajectories, paired agreement and model
# residual checks. These are descriptive diagnostics, not automatic tests.
# out is the configured LOCAL output folder; plots can reveal private IDs/data.
save_figures <- function(d, pairs, models, out) {
# Grey trajectories show individual isolates; red medians show common shifts.
  # Opposing isolate changes can cancel at the median, so inspect both layers.
  p1 <- ggplot2::ggplot(d, ggplot2::aes(replicate, log10_recovery)) +
    ggplot2::geom_line(ggplot2::aes(group = isolate_id), alpha = .15, na.rm = TRUE) +
    ggplot2::geom_point(alpha = .2, size = .7, na.rm = TRUE) +
    ggplot2::stat_summary(ggplot2::aes(group = 1), fun = median, geom = "line", colour = "red", na.rm = TRUE) +
    ggplot2::facet_wrap(~plate, ncol = 3) + ggplot2::theme_bw() +
    ggplot2::labs(x = "Biological repeat", y = "Log10 recovery CFU/ml",
                  caption = "Lines connect isolates; red lines show plate medians. Nonpositive/missing counts omitted.")
  nr <- ceiling(nlevels(d$plate) / 3)
  ggplot2::ggsave(file.path(out, "recovery_across_repeats.png"), p1, width = 11, height = max(4, nr * 2.8), dpi = 180)
# Use separate files for comparisons so panels remain legible as plates grow.
  # The equality line y=x represents perfect paired agreement; correlation
  # alone would not detect a shared shift away from that line.
  for (cmp in unique(pairs$comparison)) {
    sub <- pairs[pairs$comparison == cmp, ]
    p2 <- ggplot2::ggplot(sub, ggplot2::aes(x, y)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
      ggplot2::geom_point(alpha = .4, size = .8) + ggplot2::facet_wrap(~plate, ncol = 3) +
      ggplot2::coord_equal() + ggplot2::theme_bw() +
      ggplot2::labs(title = cmp, x = "Earlier repeat log10 CFU/ml", y = "Later repeat log10 CFU/ml")
    ggplot2::ggsave(file.path(out, paste0(gsub(" ", "_", cmp), ".png")), p2,
                    width = 10, height = max(4, nr * 2.8), dpi = 180)
  }
# A residual-versus-fitted plot checks structure and changing spread.
  # A normal quantile plot checks strong distributional departures in residuals.
  # Neither validates the biological phenotype or detects every design problem.
  # Close the PDF device explicitly so the output is complete on disk.
  grDevices::pdf(file.path(out, "model_diagnostics.pdf"), width = 9, height = 4.5)
  for (name in c("before", "adjusted")) {
    f <- models[[name]]
    graphics::par(mfrow = c(1, 2))
    graphics::plot(stats::fitted(f), stats::residuals(f), pch = 16, cex = .4,
                   xlab = "Fitted log10 recovery", ylab = "Residual", main = name)
    graphics::abline(h = 0, lty = 2)
    stats::qqnorm(stats::residuals(f), pch = 16, cex = .4)
    stats::qqline(stats::residuals(f))
  }
  grDevices::dev.off()
}
