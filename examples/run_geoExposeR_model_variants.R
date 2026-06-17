#' =============================================================================
#' Example: Fixed-, Random-, and Mixed-Effects Models (Linear & Non-Linear)
#' =============================================================================
#'
#' This script demonstrates the flexible regression model options in geoExposeR.
#' The same exposure-outcome question is fitted with every combination of:
#'
#'   effects structure (model_type):
#'     - "fixed"  : no random effects                -> stats::lm  / stats::glm
#'     - "random" : random-intercept only            -> lme4::lmer / lme4::glmer
#'     - "mixed"  : fixed + random effects (default) -> lme4::lmer / lme4::glmer
#'
#'   response type (model_family):
#'     - NULL          : linear / Gaussian model (continuous outcome)
#'     - "binomial"    : non-linear logistic model (binary outcome)
#'     - "poisson"     : non-linear count model (count outcome)
#'     - "multinomial" : multi-level categorical outcome (3+ categories)
#'
#' Engine selection (model_type x model_family):
#'
#'              linear (NULL)  non-linear (family)  multinomial
#'   "fixed"    stats::lm      stats::glm           nnet::multinom
#'   "random"   lme4::lmer     lme4::glmer          mclogit::mblogit
#'   "mixed"    lme4::lmer     lme4::glmer          mclogit::mblogit
#'
#' Data come from examples/input_data/ (the Dominion arsenic example). The
#' continuous outcome is systolic blood pressure (SBP); for the other variants
#' we derive a binary outcome (high_SBP = SBP >= 140), a count outcome
#' (sbp_events), and a 3-level categorical outcome (bp_stage) so the same
#' workflow can showcase generalized and multinomial models.
#'
#' NOTE ON THE RANDOM-EFFECTS GROUPING FACTOR
#' ------------------------------------------------------------------------------
#' The random-/mixed-effects variants below group by `age_decade` purely to
#' exercise the random-effects machinery in a compact way. Every `(1 | ...)`
#' term here is a random INTERCEPT (a group-specific baseline) with NO random
#' slope, so the arsenic effect is assumed constant across groups. Grouping on
#' `age_decade` treats age as a *contextual* variable (one baseline per age
#' band); with only 5 levels the between-group variance is estimated from very
#' few groups and is shown for illustration, not as a recommended analysis. A
#' real geographic study would group by county, `(1 | id)`, which has many
#' levels with multiple participants each (see run_geoExposeR_example.R). To
#' simply adjust for age, use the fixed effect `+ as.factor(age_decade)`.
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

output_dir <- file.path(pkg_dir, "examples", "output_model_variants")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

stopifnot(
  "Probability data file not found" = file.exists(prob_data_path),
  "Concentration data file not found" = file.exists(conc_data_path),
  "Demographics data file not found" = file.exists(demographics_data_path)
)

# ==============================================================================
# Step 3: Derive Binary and Count Outcomes for the Non-Linear Variants
# ==============================================================================
# The shipped demographics file only contains a continuous outcome (SBP).
# Generalized (non-linear) models need a binary or count response, so we
# augment the file with:
#   - high_SBP    : 1 if SBP >= 140 (stage-2 hypertension), else 0  (binomial)
#   - sbp_events  : a small non-negative count derived from SBP       (poisson)
#   - bp_stage    : 3-level category "normal"/"elevated"/"stage2" (multinomial)
demographics <- data.table::fread(demographics_data_path)
demographics$high_SBP <- as.integer(demographics$SBP >= 140)
demographics$sbp_events <- pmax(0L, round((demographics$SBP - 110) / 15))
demographics$bp_stage <- as.character(cut(
  demographics$SBP,
  breaks = c(-Inf, 120, 140, Inf),
  labels = c("normal", "elevated", "stage2")
))

augmented_demographics_path <- file.path(
  output_dir, "Demographics_dom_augmented.txt"
)
data.table::fwrite(
  demographics, augmented_demographics_path, sep = "\t"
)

# ==============================================================================
# Step 4: Shared Analysis Parameters
# ==============================================================================
exposure_level_col <- "ExposureLevel"
cat_labels <- c("<5", "5-10", "10+")
drop_cat_label_ref <- c("<5")
ndraws <- 5
seed <- 12345

# A helper that wraps perform_sensitivity_analysis with the parameters that are
# constant across every variant, so each call below differs only in the model
# specification (formula, target, model_type, model_family).
run_variant <- function(label, regression_formula, target,
                        model_type, model_family = NULL) {
  cat("\n=====================================================\n")
  family_label <- if (is.null(model_family)) {
    "NULL (linear)"
  } else if (inherits(model_family, "family")) {
    model_family$family
  } else {
    as.character(model_family)
  }
  cat("Variant:", label, "\n")
  cat("  model_type   =", model_type, "\n")
  cat("  model_family =", family_label, "\n")
  cat("  formula      =", regression_formula, "\n")
  cat("-----------------------------------------------------\n")

  variant_dir <- file.path(output_dir, gsub("[^A-Za-z0-9]+", "_", label))
  if (!dir.exists(variant_dir)) {
    dir.create(variant_dir, recursive = TRUE)
  }

  result <- tryCatch(
    perform_sensitivity_analysis(
      ndraws = ndraws,
      model_prob_csv = prob_data_path,
      conc_data_csv = conc_data_path,
      birth_data_txt = augmented_demographics_path,
      regression_formula = regression_formula,
      output_dir = variant_dir,
      targets = target,
      model_type = model_type,
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
    ),
    error = function(e) {
      cat("  Variant failed:", conditionMessage(e), "\n")
      NULL
    }
  )

  if (!is.null(result) && target %in% names(result)) {
    print(result[[target]])
  }
  invisible(result)
}

# ==============================================================================
# Step 5: Linear Variants (continuous outcome: SBP)
# ==============================================================================

# 1. Linear FIXED-effects model: ordinary least squares, no random effects.
#    age_decade enters as a fixed covariate (not a random effect).
run_variant(
  label = "Linear fixed effects",
  regression_formula =
    "~ as.factor(ExposureLevel) + as.factor(age_decade)",
  target = "SBP",
  model_type = "fixed",
  model_family = NULL
)

# 2. Linear RANDOM-effects model: random intercept for age_decade only.
run_variant(
  label = "Linear random effects",
  regression_formula = "~ as.factor(ExposureLevel) + (1 | age_decade)",
  target = "SBP",
  model_type = "random",
  model_family = NULL
)

# 3. Linear MIXED-effects model (the package default): fixed covariate plus a
#    random intercept.
run_variant(
  label = "Linear mixed effects",
  regression_formula =
    "~ as.factor(ExposureLevel) + as.factor(age_decade) + (1 | age_decade)",
  target = "SBP",
  model_type = "mixed",
  model_family = NULL
)

# ==============================================================================
# Step 6: Non-Linear Variants (generalized models)
# ==============================================================================

# 4. Non-linear FIXED-effects model: logistic regression (binary outcome) via
#    stats::glm with a binomial family.
run_variant(
  label = "Logistic fixed effects",
  regression_formula =
    "~ as.factor(ExposureLevel) + as.factor(age_decade)",
  target = "high_SBP",
  model_type = "fixed",
  model_family = "binomial"
)

# 5. Non-linear RANDOM-effects model: logistic mixed model with a random
#    intercept via lme4::glmer.
run_variant(
  label = "Logistic random effects",
  regression_formula = "~ as.factor(ExposureLevel) + (1 | age_decade)",
  target = "high_SBP",
  model_type = "random",
  model_family = "binomial"
)

# 6. Non-linear MIXED-effects model: Poisson count model with a random
#    intercept via lme4::glmer. A family object is passed instead of a name to
#    show that both forms are accepted.
run_variant(
  label = "Poisson mixed effects",
  regression_formula =
    "~ as.factor(ExposureLevel) + as.factor(age_decade) + (1 | age_decade)",
  target = "sbp_events",
  model_type = "mixed",
  model_family = poisson()
)

# ==============================================================================
# Step 7: Multinomial Variants (multi-level categorical outcome: bp_stage)
# ==============================================================================

# 7. Multinomial FIXED-effects model via nnet::multinom. The pooled results
#    include a y.level column for each non-reference outcome category.
run_variant(
  label = "Multinomial fixed effects",
  regression_formula =
    "~ as.factor(ExposureLevel) + as.factor(age_decade)",
  target = "bp_stage",
  model_type = "fixed",
  model_family = "multinomial"
)

# 8. Multinomial RANDOM-effects model: random intercept for age_decade via
#    mclogit::mblogit (random-effects multinomial logistic regression).
run_variant(
  label = "Multinomial random effects",
  regression_formula = "~ as.factor(ExposureLevel) + (1 | age_decade)",
  target = "bp_stage",
  model_type = "random",
  model_family = "multinomial"
)

# 9. Multinomial MIXED-effects model: fixed covariate plus a random intercept
#    via mclogit::mblogit.
run_variant(
  label = "Multinomial mixed effects",
  regression_formula =
    "~ as.factor(ExposureLevel) + as.factor(age_decade) + (1 | age_decade)",
  target = "bp_stage",
  model_type = "mixed",
  model_family = "multinomial"
)

# ==============================================================================
# Step 8: Notes on Interpretation
# ==============================================================================
cat("\n=====================================================\n")
cat("All model variants completed.\n")
cat("Output written under:", output_dir, "\n")
cat("\nInterpretation reminders:\n")
cat("- Linear models: q.mi is a difference in the outcome (e.g. mmHg of SBP)\n")
cat("  relative to the reference exposure category (<5).\n")
cat("- Logistic models: q.mi is a log-odds ratio; exp(q.mi) gives the odds\n")
cat("  ratio relative to the reference category.\n")
cat("- Poisson models: q.mi is a log-rate ratio; exp(q.mi) gives the rate\n")
cat("  ratio relative to the reference category.\n")
cat("- Multinomial models: one estimate per outcome category (y.level); q.mi\n")
cat("  is a log-odds ratio for that category vs the reference category.\n")
cat("- Random effects here are random INTERCEPTS grouped by age_decade (a\n")
cat("  contextual variable, no random slope); a geographic analysis would\n")
cat("  group by county, (1 | id). See run_geoExposeR_example.R.\n")
cat("=====================================================\n")
