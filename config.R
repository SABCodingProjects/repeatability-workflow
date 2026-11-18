# =============================================================================
# USER SETTINGS — this is the main file to edit for a new analysis
# =============================================================================
# Keep the generic, synthetic defaults in the public repository. Use local copies
# of real data in data/private/ and review tracked config edits before committing.
# No sample data are modified by these settings; results go to output_dir.
#
# input_mode:
#   example      = bundled synthetic CSVs; ignores your private input paths.
#   long_csv     = recovery_file and optional inoculum_file.
#   excel_blocks = workbook plus recovery_sheet and optional inoculum_sheet.
#
# zero_policy:
#   stop    = refuse a log-scale fit while measured zeros are unresolved.
#   exclude = retain raw zeros but omit them from log-scale summaries/models;
#             exports an audit list. This is a preliminary positive-only analysis.
# Neither option is a formal treatment of left-censored observations.
#
# bootstrap_n controls Monte Carlo precision and runtime, not study sample size.
# seed makes the simulated bootstrap repeatable in the same software environment.
# reliability_repeats lists hypothetical numbers of independent repeats to average.
# overwrite_outputs is FALSE to protect results from an earlier analysis.
# =============================================================================
# Run from the project root: Rscript run_analysis.R
config <- list(
  input_mode = "example",                 # "example", "long_csv", "excel_blocks"
  recovery_file = "data/private/recovery.csv",
  inoculum_file = NULL,                    # optional long CSV; never guessed
  workbook = "data/private/assay.xlsx",
  recovery_sheet = 2,
  inoculum_sheet = 3,                      # NULL skips inoculum import
  output_dir = "outputs",
  zero_policy = "exclude",                # "stop" or "exclude"; never adds a pseudocount
  bootstrap_n = 500L,                     # 1000+ may be useful for final reporting
  seed = 104729L,
  reliability_repeats = c(1L, 2L, 3L, 4L),
  overwrite_outputs = FALSE
)
