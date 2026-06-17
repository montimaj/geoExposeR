#' =============================================================================
#' Example: Running geoExposeR with Dominion Data
#' =============================================================================
#' 
#' This script demonstrates how to use the geoExposeR package to perform
#' a sensitivity analysis for arsenic exposure effects on systolic blood
#' pressure (SBP) using data from the examples/input_data/ folder.
#' 
#' The data files include:
#'   - prob_dom.csv: Modeled probability data with multinomial probabilities of arsenic levels
#'       Columns: id (FIPS code), pub_water_pop (public water population),
#'       RFCprobC1/C2/C3 (P(As<5), P(5<=As<10), P(As>=10) in µg/L),
#'       RFC_maxprob, RFC_maxprob_num, Dom_Pop_Percent (% on private wells)
#'   - conc_dom.csv: Measured concentration data from public water systems
#'       Columns: id (FIPS code), weighted_as20062008 (pop-weighted mean As in µg/L),
#'       county_pop_served20062008, well, Dom_Pop_Percent (% on private wells)
#'   - Demographics_dom.txt: Individual-level health outcomes
#'       (multiple participants per county)
#'       Columns: id (FIPS code), SBP (systolic blood pressure in mmHg),
#'       age_decade (age group: 30-39, 40-49, etc.)
#'
#' See examples/input_data/README.md for full column documentation.
#' 
#' Author: Dr. Sayantan Majumdar (sayantan.majumdar@dri.edu)
#' =============================================================================

# ==============================================================================
# Step 1: Load Required Packages
# ==============================================================================
# Use devtools::load_all() when running from source, or library() if installed
if (file.exists("DESCRIPTION")) {
  devtools::load_all()
} else {
  library(geoExposeR)
}

# ==============================================================================
# Step 2: Define File Paths
# ==============================================================================

# Get the package root directory (assumes running from geoExposeR folder)
# If running from a different directory, update this path accordingly
pkg_dir <- getwd()

# Define paths to the data files in the examples/input_data/ folder
prob_data_path <- file.path(pkg_dir, "examples", "input_data", "prob_dom.csv")
conc_data_path <- file.path(pkg_dir, "examples", "input_data", "conc_dom.csv")
demographics_data_path <- file.path(pkg_dir, "examples", "input_data", "Demographics_dom.txt")

# Create output directory for results
output_dir <- file.path(pkg_dir, "examples", "output_results")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Verify files exist
stopifnot(
  "Probability data file not found" = file.exists(prob_data_path),
  "Concentration data file not found" = file.exists(conc_data_path),
  "Demographics data file not found" = file.exists(demographics_data_path)
)

cat("Data files found:\n")
cat("  Probability data:", prob_data_path, "\n")
cat("  Concentration data:", conc_data_path, "\n")
cat("  Demographics data:", demographics_data_path, "\n\n")

# ==============================================================================
# Step 3: Concentration Data - Automatic Lognormal Conversion
# ==============================================================================
#
# The geoExposeR package now supports automatic conversion of raw
# concentration values to lognormal parameters. Simply specify:
#   - conc_raw_col: The column containing raw concentrations
#   - default_sdlog: Default sdlog value (optional, defaults to 1.0)
#
# The package will automatically calculate:
#   - meanlog = log(concentration + 0.1)
#   - sdlog = default_sdlog (or from conc_sd_col if provided)
#
# This eliminates the need for manual preprocessing!
# ==============================================================================

# Read concentration data to show available columns
conc_raw <- data.table::fread(conc_data_path)
cat("Concentration data columns:", paste(names(conc_raw), collapse = ", "), "\n")
cat("Using 'weighted_as20062008' for automatic lognormal conversion.\n\n")

# ==============================================================================
# Step 4: Define Analysis Parameters
# ==============================================================================

# Number of imputed datasets to generate
# Use a small number (5-10) for testing, larger (20-100) for publication
ndraws <- 5

# Target health outcome(s)
# SBP = Systolic Blood Pressure
targets <- c("SBP")

# Define the regression formula
# geoExposeR uses mixed-effects models (lme4::lmer), so the formula needs at
# least one random-effects term.
#
# PRIMARY EXAMPLE: a geographic random INTERCEPT by county, (1 | id).
# ---------------------------------------------------------------------------
# This is the structure most exposure studies want. Arsenic exposure is
# assigned at the county level, while SBP is measured on individuals nested
# within counties (the Demographics file has multiple participants per county).
# The (1 | id) term lets each county have its own baseline SBP (a random
# intercept), separating between-county variation from individual (residual)
# variation. The arsenic exposure coefficients are then interpreted as the
# average SBP difference between exposure categories, net of county baselines.
# There is NO random slope here, so the exposure effect is assumed constant
# across counties.
exposure_level_col <- "ExposureLevel"
regression_formula <- paste0(
  "~ as.factor(", exposure_level_col, ") + ",
  "(1 | id)"
)

# ALTERNATIVE (for demonstration only): a random intercept by age group,
# (1 | age_decade).
# ---------------------------------------------------------------------------
# This does NOT model individual age. It treats age_decade as a *contextual
# grouping variable*: each age band (30-39, 40-49, ...) gets its own random
# intercept (baseline SBP shift), with NO random slope, so the arsenic effect
# is still assumed constant across age groups. Because age_decade has only a
# handful of levels (5 here), the between-group variance is estimated from very
# few groups and should be interpreted with caution -- prefer a fixed effect,
# e.g. "+ as.factor(age_decade)", if you simply want to adjust for age. To use
# this contextual demonstration instead, swap the formula above for:
#   regression_formula <- paste0(
#     "~ as.factor(", exposure_level_col, ") + (1 | age_decade)"
#   )

# Exposure category labels (matching 3-category probability model)
cat_labels <- c("<5", "5-10", "10+")

# Reference category for regression (typically lowest exposure)
drop_cat_label_ref <- c("<5")

# Random seed for reproducibility
seed <- 12345

cat("Analysis parameters:\n")
cat("  Number of draws:", ndraws, "\n")
cat("  Target outcome(s):", paste(targets, collapse = ", "), "\n")
cat("  Regression formula:", regression_formula, "\n")
cat("  Exposure categories:", paste(cat_labels, collapse = ", "), "\n")
cat("  Reference category:", drop_cat_label_ref, "\n")
cat("  Random seed:", seed, "\n\n")

# ==============================================================================
# Step 5: Run the Sensitivity Analysis
# ==============================================================================

cat("Starting sensitivity analysis...\n")
cat("This may take a few minutes depending on the number of draws.\n\n")

# Run the main analysis function using AUTOMATIC LOGNORMAL CONVERSION
# This approach uses the original concentration file directly - no preprocessing needed!
results <- tryCatch({
  perform_sensitivity_analysis(
    ndraws = ndraws,
    model_prob_csv = prob_data_path,             # Probability data file
    conc_data_csv = conc_data_path,              # Concentration data file (no preprocessing!)
    birth_data_txt = demographics_data_path,     # Demographics file
    regression_formula = regression_formula,
    output_dir = output_dir,
    targets = targets,
    cat_label = cat_labels,
    drop_cat_label_ref = drop_cat_label_ref,
    # Custom column names matching the original data files
    id_col = "id",                            # ID column in prob/conc data
    prob_cols = c("RFCprobC1", "RFCprobC2", "RFCprobC3"),  # Probability columns
    pop_well_col = "pub_water_pop",              # Population column in prob data
    record_id_col = "id",                     # ID column in demographics
    exposure_level_col = exposure_level_col,
    # Use raw concentration column for automatic lognormal conversion
    conc_raw_col = "weighted_as20062008",     # Raw concentration column
    default_sdlog = 1.0,                     # Default sdlog value
    pwell_col = "Dom_Pop_Percent",           # Private well percentage
    conc_cutoffs = c(5, 10),
    # Imputation settings
    seed = seed,
    mice_m = 1,
    mice_maxit = 5,
    mice_method = "pmm",
    mice_covs = NULL,  # No additional covariates for MICE
    impute_vars = NULL,  # No additional variables to impute
    apply_imputation_fallback = TRUE
  )
}, error = function(e) {
  cat("\nError encountered during analysis:\n")
  cat(conditionMessage(e), "\n")
  cat("\nThis may be due to data format issues or missing values.\n")
  cat("Please check the data files for issues.\n")
  return(NULL)
})

# ==============================================================================
# Step 6: Examine and Interpret Results
# ==============================================================================
# The results table contains the following columns (pooled via Rubin's Rules):
#   term      - Regression term (e.g., "5-10", "10+")
#   q.mi      - Pooled effect estimate (mean across imputations)
#   se.mi     - Pooled standard error (accounts for within- and between-imputation variance)
#   statistic - Test statistic (q.mi / se.mi)
#   conf.low  - Lower bound of 95% confidence interval
#   conf.high - Upper bound of 95% confidence interval
#   p.value   - Two-sided p-value for the test statistic
# ==============================================================================

if (!is.null(results)) {
  cat("\n========================================\n")
  cat("Analysis completed successfully!\n")
  cat("========================================\n\n")
  
  # Display results for each target
  for (target in targets) {
    cat("Results for", target, ":\n")
    cat("----------------------------------------\n")
    
    if (target %in% names(results)) {
      result_df <- results[[target]]
      
      # Print summary table
      print(result_df)
      
      cat("\nInterpretation:\n")
      cat("- The 'q.mi' column shows the pooled effect estimate (Rubin's Rules)\n")
      cat("- The 'se.mi' shows the pooled standard error\n")
      cat("- The 'p.value' indicates statistical significance\n")
      cat("- Positive q.mi values for 5-10 or 10+ indicate higher SBP\n")
      cat("  compared to the reference category (<5)\n\n")
    } else {
      cat("No results found for this target.\n\n")
    }
  }
  
  # List output files
  cat("Output files saved to:", output_dir, "\n")
  cat("Files:\n")
  output_files <- list.files(output_dir, full.names = FALSE)
  for (f in output_files) {
    cat("  -", f, "\n")
  }
  
} else {
  cat("\nAnalysis did not complete. Please review the error messages above.\n")
}

cat("\n=== Example script completed ===\n")
