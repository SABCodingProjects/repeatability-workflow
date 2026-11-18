# Recovery phenotype repeatability workflow

A focused R workflow for evaluating whether repeated, positive CFU/ml recovery measurements distinguish isolates consistently before downstream genetic association analysis.

This repository contains reusable code and **wholly synthetic data**. It contains no real isolate identifiers, collection names, original measurements, original experiment dates, unpublished estimates, original figures, personal paths or correspondence. Synthetic examples demonstrate the workflow; their results are not biological findings.

## What researchers can investigate

| Research question | Analysis provided | What the result supports |
|---|---|---|
| Do isolate measurements agree across independent biological repeats? | Matched repeat plots, median paired differences and fold-agreement summaries | Describing reproducibility and identifying substantial disagreement |
| Does an experimental run shift many isolates together? | Plate/date summaries, trajectories and a model allowing each plate/date run its own mean | Investigating shared shifts; not establishing their cause |
| Are consistent isolate differences large relative to within-isolate fluctuations? | Plate-adjusted random-intercept model and variance components | Estimating conditional within-plate repeatability |
| Does allowing run shifts materially change repeatability? | Second random-intercept model with fixed plate/date run effects | Assessing sensitivity to common run shifts |
| How uncertain is the primary repeatability estimate? | Model-based parametric bootstrap percentile interval | Sampling uncertainty under the assumed model |
| What might another independent repeat add? | Reliability of averages under an equal-error model | Planning precision gains, not calculating GWAS power |
| Do sampled starting counts genuinely match recovery observations? | Exact-key matching and unmatched-record export | Auditing sample provenance before any normalisation |

The workflow is primarily an **estimation and diagnostic workflow**, not a collection of significance tests. It does not automatically compute an omnibus batch-effect p-value, classify runs as failed, remove statistical outliers, establish persistence biologically, estimate genetic heritability, perform GWAS, or identify causal variants.

## Appropriate study design

The supplied models assume:

- Each isolate is assigned to one plate and has a globally unique identifier.
- Recovery is measured repeatedly on the same isolate, ideally using independently grown cultures and independent assay runs.
- There are at least two plates, at least two measured dates per plate and repeated measurements remaining after exclusions.
- Plate/date defines the experimental unit for the second model's fixed run means. If the real run spans paired plates, or plates are processed in multiple sessions on one date, adapt the run definition to the actual experiment.
- Conditional on the fitted effects, residual errors are independent with a common variance, and random isolate intercepts follow an approximately normal distribution on the log scale.

Technical plating repeats of one culture are not equivalent to independent biological repeats. Resolve them using a scientifically justified counting/aggregation approach before this workflow; preserve dilution and plated-volume information separately. Do not average biological and technical levels indiscriminately.

This is a template for nested repeated-measure studies. Crossed designs, isolates reassigned across plates, repeated measures over biological time, strong heteroscedasticity, substantial censoring, and non-independent cultures may require different models.

## Project contents

| Path | Purpose |
|---|---|
| `run_analysis.R` | Entry point for import, validation, summaries, models and exported outputs |
| `config.R` | Public configuration template; defaults to synthetic inputs |
| `R/import.R` | Strict long-CSV and stacked-Excel import functions |
| `R/analysis.R` | Validation, paired comparisons, model fitting, bootstrap and plots |
| `scripts/install_dependencies.R` | Install missing R dependencies from CRAN |
| `tests/run_tests.R` | Synthetic import, matching, zero-handling and end-to-end checks |
| `data/example/` | Synthetic CSV inputs, safe to distribute |
| `docs/GITHUB_VSCODE_GUIDE.md` | Step-by-step publication and update instructions |
| `.github/workflows/r-checks.yml` | GitHub Actions checks using synthetic inputs only |
| `.gitignore` | Keeps private inputs and normal output paths out of version control |

## Requirements and quick start

Use R 4.1 or later. A recent stable R release is recommended. Packages: `readxl`, `dplyr` **1.1.0 or later**, `tidyr`, `purrr`, `ggplot2`, and `lme4`. The dplyr minimum is required for join relationship validation. `lme4` installation may need platform-specific compilation tools if a binary is unavailable.

Open the extracted project folder in VS Code. In its integrated terminal, run:

```sh
Rscript scripts/install_dependencies.R
Rscript tests/run_tests.R
Rscript run_analysis.R
```

The analysis reads synthetic examples by default and writes to `outputs/`. One synthetic observation is intentionally zero, demonstrating the exclusion audit. The integration test uses five bootstrap draws **only to exercise the code**; that is not an inferentially adequate bootstrap run.

The default analysis uses 500 draws. It may take minutes, depending on dataset size and hardware. Increase the draw count for final reporting if computationally practical and check interval stability across seeds or larger runs. Package versions are recorded in `session_info.txt`; this release does not claim a frozen dependency environment. For a publication, consider adding an `renv` lockfile after a successful local run.

Run commands from the project root. The scripts deliberately do not call `setwd()` or install packages during analysis. Existing nonempty output directories are protected from accidental overwriting; select a fresh `output_dir` or explicitly enable `overwrite_outputs`.

## Your own data

Place real inputs in `data/private/`, which is ignored by Git. Do not move them into `data/example/`. Edit `config.R` locally, but restore its generic defaults before committing; `.gitignore` does not hide edits to tracked configuration files.

### Option A — long CSV files

Set:

```r
input_mode = "long_csv"
recovery_file = "data/private/recovery.csv"
inoculum_file = "data/private/inoculum.csv" # or NULL
```

Recovery columns:

| Column | Type | Meaning |
|---|---|---|
| `isolate_id` | Text | Globally unique isolate identifier |
| `plate` | Text | Plate assignment; one per isolate |
| `date` | ISO date | Actual assay date, `YYYY-MM-DD` |
| `replicate` | Text | Biological repeat label: `R1`, `R2`, ... |
| `recovery_cfu_ml` | Numeric | Recovery concentration in CFU/ml, not raw colony numbers |

Optional inoculum columns are `isolate_id`, `date`, `replicate`, `initial_cfu_ml`. Blanks or `NA` denote missing counts; they are not zeros. Identifiers and dates cannot be missing. Numeric values must already use the same CFU/ml units; dilution conversion is outside this workflow.

Duplicate isolate/date/repeat records stop the workflow. Multiple dates for one recovery isolate/repeat also stop it: give truly different biological repeats unambiguous labels. Do not rename technical repeats just to bypass this validation.

### Option B — stacked Excel blocks

Set:

```r
input_mode = "excel_blocks"
workbook = "data/private/assay.xlsx"
recovery_sheet = 2
inoculum_sheet = 3 # or NULL
```

The workbook is read with `col_names = FALSE` and text cell types so the first block is handled exactly like later blocks.

Recovery sheet layout:

| Plate 1 | Date | Repeat | SYN_A | SYN_B | Average |
|---|---|---|---|---|---|
| | 2020-01-01 | R1 | 10 | 20 | 15 |
| | 2020-01-02 | R2 | 15 | 25 | 20 |
| | | | | | |
| Plate 2 | Date | Repeat | SYN_C | SYN_D | Average |
| | 2020-01-03 | R1 | 30 | 40 | 35 |
| | 2020-01-04 | R2 | 35 | 45 | 40 |

The plate label belongs in the first cell of each header row. Each block can have its own isolate IDs and number of isolates. Blank separator rows are allowed. Inoculum blocks omit the plate column and begin with `Date`, `Repeat`, then their own sampled-isolate IDs. Blocks are **not** paired by ordinal position.

Header matching for Date/Repeat is case-insensitive; replicate labels are uppercase `R` followed by digits. `Average` and `Mean` columns are ignored. Blank headers containing measurements cause an error rather than silent data loss. Notes or other nonblank rows inside blocks are rejected: remove them from the input region or adapt the importer explicitly. The first header must not be preceded by nonblank title rows. Dates must be ISO strings or Excel serial dates using the 1900 date system; 1904-system workbooks must be converted first. No ambiguous day/month strings are guessed. Formulas must have usable cached values or be converted to values before import. Hidden rows/columns are not automatically excluded.

### Inoculum matching

The key is **isolate ID + date + replicate**. Matched initial counts are exported without altering recovery values. Unmatched inoculum records are exported for investigation. Unsampled isolates keep missing initial counts. Inoculum sampling is often nonrandom; the sampled subset should not be assumed representative of every isolate.

A matching key verifies labels, not biology: independently confirm that initial and surviving counts refer to the same culture/treatment experiment. This workflow deliberately does not calculate survival fractions or replace missing starting counts with a plate mean. If normalisation is justified, define and validate that phenotype before adapting the model.

## Zero counts, missing data and assay limits

`zero_policy = "stop"` prevents modelling when zeros occur. `zero_policy = "exclude"` permits a **preliminary positive-observation analysis** and records all zero/missing exclusions in `excluded_observations.csv`. It does not discard entire isolates or add 1 to counts. The distributed configuration uses exclusion to demonstrate this behaviour on synthetic data.

A no-colony observation usually provides information about a detection bound, not a precisely measured zero concentration. Positive values may also be imposed floors or ceilings. The code cannot infer these from repeated numbers. Check assay documentation, plated volume and dilution factors. Strings such as `<5` are intentionally rejected; converting them to 5 without recording censoring changes their meaning.

If censoring is substantial, use a suitable censored or count-based model and appropriate assay metadata rather than relying on this Gaussian positive-only model. Even a small excluded fraction warrants a sensitivity check when exclusions are concentrated by isolate, genotype or run. Missingness may also be informative.

## Statistical methods and interpretation

### Paired agreement

For a pair of repeats, define `d_i = log10(later_i) - log10(earlier_i)` for the same isolate. This equals `log10(later_i / earlier_i)`.

- Median paired difference: typical directional shift; exponentiating it gives the corresponding fold ratio.
- Median absolute paired difference: typical disagreement regardless of direction.
- Within-twofold percentage: fraction with `abs(d_i) <= log10(2)`.
- Within-tenfold percentage: fraction with `abs(d_i) <= 1`.

The margins are descriptive and configurable in the code, not universal quality thresholds. A zero median paired shift does not establish agreement. Shared detection floors can make agreement look better than it is. Pairwise summaries are not independent observations because one isolate contributes to several comparisons; do not use their row count as an independent sample size.

### Primary model

```r
log10_recovery ~ plate + (1 | isolate_key)
```

On the log scale, `y_ij = mu + plate_effect + u_i + error_ij`. The isolate random intercept has variance `sigma_isolate^2`; the residual has variance `sigma_residual^2`. Both components are estimated with restricted maximum likelihood (REML).

Conditional within-plate repeatability is:

`ICC = sigma_isolate^2 / (sigma_isolate^2 + sigma_residual^2)`.

It describes consistent isolate differences relative to unexplained within-isolate fluctuations **after plate means are accounted for**. It is not genetic heritability, per-observation accuracy, a fraction of correctly identified isolates, or the probability of GWAS success. Residual variation can contain genuine biological stochasticity as well as measurement error.

Plate is fixed to avoid counting broad plate differences as reproducible isolate variation. This also means between-plate genetic or biological differences are not evaluated by this ICC. If isolate genetic backgrounds cluster by plate, plate and biology may be confounded.

### Model allowing run shifts

```r
log10_recovery ~ 0 + run_id + (1 | isolate_key)
```

`run_id` is plate × date. The model estimates a separate mean for every observed plate/date combination. `0 +` is a parameterisation choice, not a zero biological baseline. Compare variance estimates and ICC values descriptively. A small change means common mean shifts explain relatively little under these models; it does not prove that batch effects are absent. REML likelihoods from models with different fixed effects should **not** be compared with a naive likelihood-ratio test.

Run-adjusted ICC is conditional on the fitted run means and may look optimistic for performance on an unseen day. The workflow does not export adjusted isolate scores as ready-made GWAS phenotypes.

### Bootstrap uncertainty

The primary model is bootstrapped with `bootMer(type = "parametric", use.u = FALSE)`, generating new random isolate effects and residuals under the fitted model. Each simulated dataset is refitted; percentile limits are taken from finite ICC estimates. Requested draws, finite estimates, failed refits and messages are retained.

Inspect those diagnostics. A warning-free finite estimate can still come from a singular simulated fit; review the message export. If many draws fail or models are unstable, investigate before treating a percentile interval as reliable. A narrow interval is conditional on the model assumptions; it does not protect against systematic bias, wrong labels or unsuitable censoring treatment.

### Reliability of averages

For `k` independent repeats with equal residual variance:

`R_k = k * ICC / (1 + (k - 1) * ICC)`.

This is a planning calculation for averages of **log-transformed recovery**, equivalent to geometric rather than arithmetic means on the original scale. Correlated errors, shared biases, unequal precision and informative exclusions invalidate a simple interpretation. With unbalanced data, it describes a hypothetical complete `k`-repeat average, not the exact reliability of every observed isolate mean.

The primary model's residual also contains unmodelled run shifts. Consequently, treat this calculation as approximate unless its independence assumptions are defensible. This workflow does not automatically decide that an extra repeat is unnecessary; the scientific value of the expected gain must be weighed against cost and study goals.

## Output reference

All outputs are local and ignored by Git by default.

| File | Use |
|---|---|
| `recovery_long.csv` | Audit imported recovery measurements |
| `inoculum_long.csv` | Audit optional starting-count data |
| `unmatched_inoculum.csv` | Investigate labels/dates; never force a join |
| `recovery_with_matched_inoculum.csv` | Recovery observations with only valid exact-key matches |
| `inoculum_summary.csv` | Descriptive starting-count summaries by date/repeat; block composition may differ |
| `excluded_observations.csv` | Missing/zero observations excluded from log-scale analysis |
| `isolate_completeness.csv` | Remaining measurements per isolate |
| `run_summary.csv` | Plate/date/repeat recovery summaries |
| `agreement_summary.csv` | Paired shifts and fold agreement |
| `repeatability.csv` | Variance components and both ICC estimates |
| `model_status.csv` | Singular-fit and convergence information |
| `bootstrap_ci.csv`, `bootstrap_draws.csv`, `bootstrap_messages.txt` | Primary ICC uncertainty and bootstrap audit |
| `average_reliability.csv` | Planning values for repeat averages |
| `recovery_across_repeats.png` | Individual trajectories and plate medians |
| `R*_vs_R*.png` | Separate, readable agreement plots for each repeat pair |
| `model_diagnostics.pdf` | Residual-versus-fitted and normal quantile plots for both models |
| `model_summaries.txt` | Model coefficient and fit summaries |
| `analysis_config.txt`, `session_info.txt` | Settings and software environment for reproduction |

Plot residuals against fitted values for changing spread or structure; inspect quantile plots for strong departures. The diagnostics are not an automatic pass/fail test. Missing/zero exclusions and strong ceiling effects must be considered even when residual plots look acceptable.

## Turning results into a research decision

A defensible recommendation combines agreement, variance components, confidence intervals, biological controls, experimental provenance and model diagnostics. For example:

> “There is a reproducible within-plate isolate signal. The estimated gain from another independent repeat appears modest under the fitted model. We suggest exploratory genetic analysis after resolving assay-limit handling and considering plate allocation.”

Do not replace that with “variability will not affect GWAS.” A reproducible phenotype can lack genetic signal; a noisy one can still reveal strong effects. Formal GWAS power depends on effect sizes, allele frequencies, sample size, relatedness and genetic architecture. Sensitivity analyses using different valid phenotype summaries or independent repeats can complement, but not replace, validation. Do not select whichever repeat yields the most attractive association results.

## Privacy and publication

Only the synthetic data and generic source code are intended for publication. Outputs can contain original IDs, dates, counts, model coefficients and local paths; treat them as research data. Keep real data in `data/private/` and output under `outputs/`. Custom output folders need their own ignore rule. `.gitignore` does not remove previously tracked files, encrypt data or check spreadsheet metadata. Review the staged diff and tracked file list before publishing. Do not publish original workbooks, identifiable screenshots, correspondence, private configuration, credentials or mapping keys.

The code is provided under the MIT license. Confirm institutional/employer ownership and third-party obligations before applying that license to your published repository. No original dataset is licensed or distributed here.

## Validation status

Synthetic fixtures and a GitHub Actions workflow are included. The supplied source was statically reviewed and the example tables and distribution contents were checked during preparation. **R was unavailable in the preparation environment, so the R tests and model workflow were not executed there.** Run the tests and example locally and confirm the GitHub Actions result before relying on or describing this release as tested. No original scientific estimates are reproduced by the synthetic example.

## References

- [lme4 package and mixed-model implementation](https://lme4.github.io/lme4/)
- [bootMer reference and simulation options](https://lme4.github.io/lme4/reference/bootMer.html)
- [VS Code source control quickstart](https://code.visualstudio.com/docs/sourcecontrol/quickstart)
- [Git documentation](https://git-scm.com/doc)

## Suggested reporting information

Record assay conditions, independent versus technical repeat definitions, starting-count sampling, units, detection bounds, missingness/exclusions, plate allocation, model formulas, variance components, ICC and interval method, bootstrap diagnostics, residual checks, software versions and any sensitivity analyses. Cite this repository using its eventual public URL and a tagged release or commit; no DOI or repository address has been invented in this template.
