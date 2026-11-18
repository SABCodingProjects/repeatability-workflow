#!/usr/bin/env Rscript
# =============================================================================
# MAIN WORKFLOW — read this script first
# =============================================================================
# Purpose: organise the complete analysis; the detailed operations live in R/.
# Run with: Rscript run_analysis.R, from the folder containing this script.
# Read config.R before changing inputs or exclusion rules.
#
# Scientific sequence:
#   import -> validate -> audit inoculum matching -> describe repeat agreement
#   -> estimate variance components -> bootstrap uncertainty -> save diagnostics.
# The script estimates repeatability. It does not perform GWAS or establish
# that experimental variation has no effect on genetic association results.
#
# Important R notation used below:
#   <- assigns an object; |> passes an object into the next function.
#   package::function identifies the package explicitly.
#   source() runs another R file; it makes that file's functions/settings available.
#   NULL means an optional input is absent, not that its values are zero.
# =============================================================================
# Run from the project root. Private inputs and all generated outputs stay local.
# 1. Check dependencies without installing anything during a scientific run.
# This keeps software installation separate from analysis execution.
packages <- c("readxl", "dplyr", "tidyr", "purrr", "ggplot2", "lme4")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install dependencies first: Rscript scripts/install_dependencies.R")
if (utils::packageVersion("dplyr") < "1.1.0") stop("Update dplyr to version 1.1.0 or later.")
suppressPackageStartupMessages(library(dplyr))
# Load reusable functions first, then the settings that determine this run.
source("R/import.R")
source("R/analysis.R")
source("config.R")
# 2. Prepare a local output folder. Refuse to overwrite prior results by default.
# CSV outputs may contain private sample IDs and must not be committed to Git.
out <- config$output_dir
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE)) &&
    !isTRUE(config$overwrite_outputs)) stop("Output directory is not empty. Choose a new directory or enable overwrite_outputs.")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
# Small helper: export a table without artificial row-number columns.
# Empty fields represent missing values in the exported CSV, not measured zeros.
write_table <- function(d, name) utils::write.csv(d, file.path(out, name), row.names = FALSE, na = "")

# 3. Import one of three supported layouts. The example mode is wholly synthetic.
# All import routes return one row per isolate/date/biological-repeat observation.
inoculum <- NULL
if (config$input_mode == "example") {
  recovery <- read_long("data/example/recovery.csv", "recovery")
  inoculum <- read_long("data/example/inoculum.csv", "inoculum")
} else if (config$input_mode == "long_csv") {
  recovery <- read_long(config$recovery_file, "recovery")
  if (!is.null(config$inoculum_file)) inoculum <- read_long(config$inoculum_file, "inoculum")
} else if (config$input_mode == "excel_blocks") {
  recovery <- read_blocks(config$workbook, config$recovery_sheet, "recovery")
  if (!is.null(config$inoculum_sheet)) inoculum <- read_blocks(config$workbook, config$inoculum_sheet, "inoculum")
} else stop("Unknown input_mode.")

# 4. Reject ambiguous identifiers or duplicate observations before modelling.
# A duplicate might be a technical replicate that requires separate aggregation.
# It must not silently become an extra independent biological observation.
validate_observations(recovery, "recovery")
write_table(recovery, "recovery_long.csv")
if (!is.null(inoculum)) {
  validate_observations(inoculum, "inoculum")
  write_table(inoculum, "inoculum_long.csv")
# Match only measurements sharing all three keys. Block number is deliberately
  # absent: the sixth inoculum block need not describe the sixth recovery block.
  keys <- c("isolate_id", "date", "replicate")
# anti_join returns inoculum records with no exact recovery match.
  # These are an audit list for checking the source, not an instruction to delete.
  write_table(anti_join(inoculum, recovery, by = keys), "unmatched_inoculum.csv")
  # Retain only the matched initial counts; never substitute a run mean.
# left_join keeps every recovery observation, including unsampled isolates.
  # The many-to-one check prevents duplicate inoculum keys from multiplying rows.
  joined <- left_join(recovery, select(inoculum, all_of(keys), initial_cfu_ml),
                      by = keys, relationship = "many-to-one")
  write_table(joined, "recovery_with_matched_inoculum.csv")
  # Descriptive summaries only: no survival ratios are used in the models.
  write_table(inoculum |> group_by(date, replicate) |>
                summarise(n_observed = sum(!is.na(initial_cfu_ml)),
                          median_initial_cfu_ml = median(initial_cfu_ml, na.rm = TRUE),
                          .groups = "drop"), "inoculum_summary.csv")
}
# 5. Create log10 values under the explicitly chosen zero-handling policy.
# The full imported counts stay intact; exclusions are exported for review.
# No automatic outlier removal or positive pseudocount is applied.
rec <- prepare_recovery(recovery, config$zero_policy)
write_table(filter(rec, !is.na(exclusion_reason)), "excluded_observations.csv")
# Count usable observations for each isolate. Missing/zero exclusions can leave
# an unbalanced design, even if the original spreadsheet had three repeat rows.
write_table(rec |> group_by(plate, isolate_id) |>
              summarise(n_rows = n(), n_finite = sum(is.finite(log10_recovery)),
                        .groups = "drop"), "isolate_completeness.csv")
# Summarise runs descriptively. A stable median does not guarantee that the
# individual isolates agree: upward and downward changes can cancel.
write_table(rec |> group_by(plate, date, replicate) |>
              summarise(n_rows = n(), n_finite = sum(is.finite(log10_recovery)),
                        median_log10 = median(log10_recovery, na.rm = TRUE),
                        .groups = "drop"), "run_summary.csv")
# 6. Compare every available pair of repeats within the same isolate.
# These paired comparisons share isolates and are not independent sample groups.
pairs <- make_pairs(rec)
# On the log10 scale, y - x = log10(later / earlier).
# Median signed difference describes direction; absolute difference describes size.
# A twofold difference is log10(2); a tenfold difference is 1.
# The tiny tolerance avoids classifying exact boundary values incorrectly due to
# floating-point rounding. Neither margin is a universal quality threshold.
write_table(pairs |> group_by(plate, comparison) |>
              summarise(n = n(), median_shift_log10 = median(y - x),
                        median_absolute_difference = median(abs(y - x)),
                        percent_within_2fold = 100 * mean(abs(y - x) <= log10(2) + 1e-12),
                        percent_within_10fold = 100 * mean(abs(y - x) <= 1 + 1e-12),
                        .groups = "drop"), "agreement_summary.csv")
# 7. Fit two diagnostic models with the same finite positive observations.
# The first accounts for plate means; the second for plate-by-date means.
# Compare ICCs descriptively; do not use REML likelihoods to test these different
# fixed-effect specifications with a naive likelihood-ratio test.
models <- fit_models(rec)
results <- bind_rows(before_run_adjustment = get_repeatability(models$before),
                     after_run_adjustment = get_repeatability(models$adjusted), .id = "model")
write_table(results, "repeatability.csv")
# Convergence and singularity need review even when a model returns numbers.
# A singular fit commonly indicates a variance estimate close to a boundary.
status <- bind_rows(model_status(models$before, "before_run_adjustment"),
                    model_status(models$adjusted, "after_run_adjustment"))
write_table(status, "model_status.csv")
if (any(status$singular) || any(status$convergence_messages != "none") ||
    any(!status$optimizer_code %in% c("0", ""))) warning("Inspect model_status.csv before interpreting estimates.")
# 8. Quantify uncertainty in the primary, plate-adjusted ICC.
# Simulated datasets come from the fitted model, so the interval is conditional
# on that model being an adequate description of the experiment.
boot <- bootstrap_icc(models$before, config$bootstrap_n, config$seed)
write_table(boot$summary, "bootstrap_ci.csv")
write_table(boot$draws, "bootstrap_draws.csv")
utils::capture.output(boot$messages, file = file.path(out, "bootstrap_messages.txt"))
if (boot$summary$failed_refits > 0 || boot$summary$finite_estimates < config$bootstrap_n) {
  warning("Bootstrap refits failed or produced non-finite estimates; inspect the bootstrap diagnostics.")
}
# 9. Approximate the reliability of hypothetical averages of k independent
# repeats. Use the primary ICC, not the run-adjusted ICC, as the planning input.
# Formula: R_k = k * ICC / (1 + (k - 1) * ICC).
# Equal residual variance and independent errors are essential assumptions;
# systematic bias does not disappear by averaging. These are log-scale averages,
# corresponding to geometric means of original counts, not arithmetic means.
k <- config$reliability_repeats
if (any(!is.finite(k) | k < 1 | k != as.integer(k))) stop("reliability_repeats must be positive integers.")
icc <- results$repeatability[1]
write_table(data.frame(biological_repeats = k,
                        estimated_reliability = k * icc / (1 + (k - 1) * icc)),
            "average_reliability.csv")
# 10. Export plots, fitted-model summaries and the computational context.
# Inspect residual patterns and quantile plots before reporting estimates.
# Settings and sessionInfo document what ran; package versions are recorded
# but are not frozen by a dependency lockfile in this template.
save_figures(rec, pairs, models, out)
utils::capture.output(summary(models$before), summary(models$adjusted),
                       file = file.path(out, "model_summaries.txt"))
utils::capture.output(dput(config), file = file.path(out, "analysis_config.txt"))
utils::capture.output(sessionInfo(), file = file.path(out, "session_info.txt"))
message("Analysis complete. Inspect exclusions, model diagnostics and bootstrap status before reporting results.")
