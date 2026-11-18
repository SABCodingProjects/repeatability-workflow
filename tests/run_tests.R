# =============================================================================
# SYNTHETIC CHECKS — verify data identity and pipeline behaviour
# =============================================================================
# Run: Rscript tests/run_tests.R from the project root.
# Tests exercise failure modes that could change a scientific conclusion:
# wrong block headers, row/column misalignment, dropped unnamed measurements,
# duplicate observations, invalid numeric cells, zeros, and false inoculum joins.
# An end-to-end check then runs the models and exports in a temporary folder.
# No actual biological data or original study results are used here.
# stopifnot raises an error on a failed condition, producing a failing CI job.
# =============================================================================
# Meaningful integration checks. No real research data are required.
library(dplyr)
source("R/import.R")
source("R/analysis.R")
# This helper succeeds only when invalid input raises an error. It ensures
# that validation guards reject problematic data rather than silently cleaning it.
expect_error <- function(expr) {
  result <- tryCatch({ force(expr); FALSE }, error = function(e) TRUE)
  stopifnot(result)
}
# Build two spreadsheet-like blocks with deliberately DIFFERENT isolate IDs.
# If the second block inherits the first header, these tests must fail.
# Different counts also check the row-wise ordering of isolate/repeat assignments.
raw <- rbind(
 c("Plate 1", "Date", "Repeat", "SYN_A", "SYN_B", "Average"),
 c(NA, "2020-01-01", "R1", "10", "20", "15"),
 c(NA, "2020-01-02", "R2", "30", "40", "35"),
 rep(NA, 6),
 c("Plate 2", "Date", "Repeat", "SYN_C", "SYN_D", "Average"),
 c(NA, "2020-01-03", "R1", "50", "60", "55"),
 c(NA, "2020-01-04", "R2", "70", "80", "75")
)
d <- parse_blocks(raw, "recovery")
# Assert both identities and values: a correct row count alone would not catch
# a valid-looking table whose measurements were attached to the wrong isolate.
stopifnot(nrow(d) == 8L, identical(d$isolate_id, c("SYN_A", "SYN_B", "SYN_A", "SYN_B", "SYN_C", "SYN_D", "SYN_C", "SYN_D")),
          identical(d$recovery_cfu_ml, c(10,20,30,40,50,60,70,80)),
          all(d$plate[d$isolate_id == "SYN_C"] == "Plate 2"))
stopifnot(parse_dates("43831") == as.Date("2020-01-01"))
# Remove a header while leaving its counts: the importer must not drop them.
# Then test a censored value string, which needs a scientific decision first.
bad <- raw; bad[1, 5] <- NA
expect_error(parse_blocks(bad, "recovery"))
bad <- raw; bad[2, 4] <- "<5"
expect_error(parse_blocks(bad, "recovery"))
expect_error(validate_observations(bind_rows(d, d[1, ]), "recovery"))
# Test both explicit zero policies. Exclusion must create NA for the transformed
# value, never -Inf, while retaining the original measurement for audit.
z <- d; z$recovery_cfu_ml[1] <- 0
expect_error(prepare_recovery(z, "stop"))
z <- suppressWarnings(prepare_recovery(z, "exclude"))
stopifnot(sum(is.na(z$log10_recovery)) == 1L, !any(is.infinite(z$log10_recovery)))
# Inoculum should match only on exact identifiers, date and replicate.
i <- d[1:2, c("date", "replicate", "isolate_id")]
i$initial_cfu_ml <- c(1000, 2000)
i$date[1] <- i$date[1] + 1
stopifnot(nrow(anti_join(i, d, by = c("isolate_id", "date", "replicate"))) == 1L)
# End-to-end model/bootstrap/output integration check in a temporary directory.
# Use the synthetic example mode and a temporary output folder. A tiny
# bootstrap count checks execution only; its interval must never be reported
# as a scientifically adequate confidence interval.
source("config.R")
config$bootstrap_n <- 5L
config$output_dir <- tempfile("synthetic-analysis-")
config$overwrite_outputs <- FALSE
# run_analysis.R sources config.R, so evaluate a temporary entry point with that
# source statement replaced by an explicit configuration assignment.
# Load the entry-point code but prevent it from reloading default settings
# over this test configuration. The remaining workflow executes unchanged.
entry <- readLines("run_analysis.R")
config_line <- which(entry == 'source("config.R")')
stopifnot(length(config_line) == 1L)
entry[config_line] <- "# Configuration supplied by integration test"
suppressWarnings(eval(parse(text = entry), envir = environment()))
stopifnot(file.exists(file.path(config$output_dir, "repeatability.csv")),
          file.exists(file.path(config$output_dir, "model_diagnostics.pdf")))
# Check basic mathematical bounds and finite estimates, not a hard-coded ICC.
# A fixed expected estimate would be brittle across legitimate software versions.
# These checks do not establish that modelling assumptions hold for real data.
r <- read.csv(file.path(config$output_dir, "repeatability.csv"))
stopifnot(nrow(r) == 2L, all(is.finite(r$repeatability)),
          all(r$repeatability >= 0 & r$repeatability <= 1))
# Delete only the temporary test output folder created above. Actual user
# inputs and the normal outputs/ directory are not touched.
unlink(config$output_dir, recursive = TRUE)
message("All synthetic checks passed.")
