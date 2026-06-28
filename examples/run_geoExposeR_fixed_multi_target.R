#' =============================================================================
#' Example: Fixed-Effects Models with Multiple Targets (incl. Categorical)
#' =============================================================================
#'
#' This script demonstrates a FIXED-effects analysis (model_type = "fixed", no
#' random effects) run over MULTIPLE health outcomes at once, where some of the
#' outcomes are CATEGORICAL.
#'
#' Because a single perform_sensitivity_analysis() call applies one response
#' distribution (model_family) to all of its targets, outcomes are grouped by
#' type and analysed in three fixed-effects calls:
#'
#'   1. Continuous outcomes         -> linear fixed effects
#'        (model_family = NULL, stats::lm)
#'   2. Binary categorical outcomes -> logistic fixed effects
#'        (model_family = "binomial", stats::glm)
#'   3. Multi-level categorical     -> multinomial fixed effects
#'        (model_family = "multinomial", nnet::multinom)
#'
#' Categorical outcomes are handled directly:
#'   - A two-level character column such as "high"/"normal" (as it is read from
#'     a text file) is coerced to a factor automatically before fitting, and a
#'     0/1 indicator works as-is.
#'   - A column with three or more categories (e.g. "normal"/"elevated"/"stage2")
#'     is modeled with multinomial logistic regression. Multinomial results gain
#'     a `y.level` column naming the outcome category (relative to the reference
#'     level) that each estimate refers to.
#'
#' Data come from examples/input_data/ (the Dominion arsenic example), which
#' ships a continuous outcome (SBP). We derive the additional continuous and
#' categorical outcomes used here for demonstration purposes.
#'
#' Author: Dr. Sayantan Majumdar (sayantan.majumdar@dri.edu)
#' =============================================================================

# ==============================================================================
# Step 1: Load Required Packages
# ==============================================================================
if (file.exists("DESCRIPTION")) {
  devtools::load_all()
} else {
  library(geoExposeR)
}

# ==============================================================================
# Step 2: Define File Paths
# ==============================================================================
pkg_dir <- getwd()

prob_data_path <- file.path(pkg_dir, "examples", "input_data", "prob_dom.csv")
conc_data_path <- file.path(pkg_dir, "examples", "input_data", "conc_dom.csv")
demographics_data_path <- file.path(
  pkg_dir, "examples", "input_data", "Demographics_dom.txt"
)

output_dir <- file.path(pkg_dir, "examples", "output_fixed_multi_target")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

stopifnot(
  "Probability data file not found" = file.exists(prob_data_path),
  "Concentration data file not found" = file.exists(conc_data_path),
  "Demographics data file not found" = file.exists(demographics_data_path)
)

# ==============================================================================
# Step 3: Build Multiple Continuous and Categorical Outcomes
# ==============================================================================
# The shipped demographics file contains a single continuous outcome (SBP).
# We augment it with:
#   Continuous outcomes (linear fixed effects):
#     - SBP            : systolic blood pressure (mmHg), as shipped
#     - DBP            : diastolic blood pressure (mmHg), synthesised from SBP
#     - MAP_est        : mean arterial pressure (mmHg), MAP = DBP + 1/3(SBP-DBP)
#   Categorical outcomes (logistic fixed effects):
#     - bp_category    : a TWO-LEVEL CHARACTER label, "high" vs "normal"
#                        (high if SBP >= 140 OR DBP >= 90). Stored as text on
#                        disk to show that character categorical outcomes work.
#     - hypertension   : a 0/1 indicator using the ACC/AHA 2017 thresholds
#                        (1 if SBP >= 130 OR DBP >= 80), the numeric form.
demographics <- data.table::fread(demographics_data_path)

set.seed(12345)

# Derived diastolic blood pressure (mmHg). The shipped file ships only SBP, so
# we synthesise a plausible DBP (~half of SBP plus a small offset, giving a
# realistic pulse pressure of ~40-55 mmHg) purely so the example can use the
# textbook MAP formula and DBP-aware categories below.
demographics$DBP <- round(
  0.5 * demographics$SBP + stats::rnorm(nrow(demographics), 15, 4)
)

# Estimated mean arterial pressure (mmHg) using the standard formula:
#   MAP = DBP + (1/3) * (SBP - DBP) = (2 * DBP + SBP) / 3
demographics$MAP_est <- round(
  demographics$DBP + (1 / 3) * (demographics$SBP - demographics$DBP) +
    stats::rnorm(nrow(demographics), 0, 4)
)

# Two-level CHARACTER label ("high"/"normal"); "high" follows the stage-2
# thresholds (SBP >= 140 OR DBP >= 90). Stored as text on disk to show that
# character categorical outcomes are handled.
demographics$bp_category <- ifelse(
  demographics$SBP >= 140 | demographics$DBP >= 90, "high", "normal"
)

# 0/1 indicator using the ACC/AHA 2017 hypertension thresholds
#   (SBP >= 130 OR DBP >= 80).
demographics$hypertension <- as.integer(
  demographics$SBP >= 130 | demographics$DBP >= 80
)

# Multi-level (3-category) DBP-aware label stored as text:
#   stage2   : SBP >= 140 OR DBP >= 90
#   elevated : SBP >= 120 OR DBP >= 80 (but not stage2)
#   normal   : otherwise
demographics$bp_stage <- ifelse(
  demographics$SBP >= 140 | demographics$DBP >= 90, "stage2",
  ifelse(
    demographics$SBP >= 120 | demographics$DBP >= 80, "elevated", "normal"
  )
)

augmented_demographics_path <- file.path(
  output_dir, "Demographics_dom_multi_target.txt"
)
data.table::fwrite(demographics, augmented_demographics_path, sep = "\t")

cat("Outcome columns written:",
    paste(
      c("SBP", "DBP", "MAP_est", "bp_category", "hypertension", "bp_stage"),
      collapse = ", "
    ),
    "\n")
cat("  DBP/MAP_est are continuous; bp_category is character",
    "('high'/'normal');\n")
cat("  hypertension is 0/1; bp_stage is 3-level",
    "('normal'/'elevated'/'stage2').\n\n")

# ==============================================================================
# Step 4: Shared Analysis Parameters
# ==============================================================================
exposure_level_col <- "ExposureLevel"
cat_labels <- c("<5", "5-10", "10+")
drop_cat_label_ref <- c("<5")
ndraws <- 5
seed <- 12345

# A FIXED-effects formula: no random-effects terms are allowed. age_decade
# enters as an ordinary fixed covariate (via as.factor()).
fixed_formula <- "~ as.factor(ExposureLevel) + as.factor(age_decade)"

# Helper wrapping the constant arguments so the two calls below differ only in
# their target set and model_family.
run_fixed <- function(targets, model_family, results_subdir) {
  perform_sensitivity_analysis(
    ndraws = ndraws,
    model_prob_csv = prob_data_path,
    conc_data_csv = conc_data_path,
    birth_data_txt = augmented_demographics_path,
    regression_formula = fixed_formula,
    output_dir = file.path(output_dir, results_subdir),
    targets = targets,
    model_type = "fixed",
    model_family = model_family,
    cat_label = cat_labels,
    drop_cat_label_ref = drop_cat_label_ref,
    id_col = "id",
    prob_cols = c("RFCprobC1", "RFCprobC2", "RFCprobC3"),
    pop_well_col = "pub_water_pop",
    record_id_col = "id",
    exposure_level_col = exposure_level_col,
    conc_raw_col = "weighted_as20062008",
    default_sdlog = 1.0,
    pwell_col = "Dom_Pop_Percent",
    conc_cutoffs = c(5, 10),
    seed = seed,
    apply_imputation_fallback = TRUE
  )
}

# ==============================================================================
# Step 5: Continuous Outcomes - Linear Fixed Effects (stats::lm)
# ==============================================================================
cat("=====================================================\n")
cat("Linear fixed-effects models for CONTINUOUS outcomes\n")
cat("  targets      = SBP, DBP, MAP_est\n")
cat("  model_type   = fixed,  model_family = NULL (Gaussian/lm)\n")
cat("-----------------------------------------------------\n")

continuous_results <- tryCatch(
  run_fixed(
    targets = c("SBP", "DBP", "MAP_est"),
    model_family = NULL,
    results_subdir = "continuous"
  ),
  error = function(e) {
    cat("Continuous analysis failed:", conditionMessage(e), "\n")
    NULL
  }
)

# ==============================================================================
# Step 6: Categorical Outcomes - Logistic Fixed Effects (stats::glm)
# ==============================================================================
cat("\n=====================================================\n")
cat("Logistic fixed-effects models for CATEGORICAL outcomes\n")
cat("  targets      = bp_category (character), hypertension (0/1)\n")
cat("  model_type   = fixed,  model_family = binomial (glm)\n")
cat("-----------------------------------------------------\n")

categorical_results <- tryCatch(
  run_fixed(
    targets = c("bp_category", "hypertension"),
    model_family = "binomial",
    results_subdir = "categorical"
  ),
  error = function(e) {
    cat("Categorical analysis failed:", conditionMessage(e), "\n")
    NULL
  }
)

# ==============================================================================
# Step 6b: Multi-Level Categorical Outcome - Multinomial Fixed Effects
# ==============================================================================
cat("\n=====================================================\n")
cat("Multinomial fixed-effects model for a MULTI-LEVEL outcome\n")
cat("  target       = bp_stage ('normal'/'elevated'/'stage2')\n")
cat("  model_type   = fixed,  model_family = multinomial (nnet::multinom)\n")
cat("-----------------------------------------------------\n")

multinomial_results <- tryCatch(
  run_fixed(
    targets = "bp_stage",
    model_family = "multinomial",
    results_subdir = "multinomial"
  ),
  error = function(e) {
    cat("Multinomial analysis failed:", conditionMessage(e), "\n")
    NULL
  }
)

# ==============================================================================
# Step 7: Combine and Display Results
# ==============================================================================
all_results <- c(continuous_results, categorical_results, multinomial_results)

cat("\n=====================================================\n")
cat("Pooled fixed-effects results across all targets\n")
cat("=====================================================\n")
for (target in names(all_results)) {
  scale_note <- if (target %in% c("SBP", "DBP", "MAP_est")) {
    "difference in the outcome vs reference exposure (<5)"
  } else if (target == "bp_stage") {
    paste(
      "log-odds ratio for each outcome category (y.level) vs the reference",
      "category, comparing each exposure level to the reference exposure (<5)"
    )
  } else {
    "log-odds ratio vs reference exposure (<5); exp(q.mi) = odds ratio"
  }
  cat("\nTarget:", target, "\n")
  cat("  Interpretation of q.mi:", scale_note, "\n")
  print(all_results[[target]])
}

cat("\nOutput written under:", output_dir, "\n")
cat("\n=== Fixed-effects multi-target example completed ===\n")
