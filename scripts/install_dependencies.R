# =============================================================================
# ONE-TIME DEPENDENCY SETUP
# =============================================================================
# Run from the project root: Rscript scripts/install_dependencies.R
# This installs missing packages from CRAN into your configured R library.
# It does not read research data or run models. You may need build tools if a
# platform-specific binary is unavailable. Already-installed packages are not
# upgraded automatically; see the README for minimum supported versions.
# =============================================================================
packages <- c("readxl", "dplyr", "tidyr", "purrr", "ggplot2", "lme4")
# requireNamespace checks availability without attaching every package.
# vapply returns one TRUE/FALSE per package using a specified logical result type.
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
# Install only when something is missing. Analysis scripts never do this
# automatically, making software changes a separate, deliberate setup step.
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
message("Dependencies available. Run Rscript run_analysis.R from the project root.")
