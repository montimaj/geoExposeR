#' =============================================================================
#' Example: Running geoExposeR with Dominion Data
#' =============================================================================
#' 
#' This script demonstrates how to use the geoExposeR package to perform
#' a sensitivity analysis for arsenic exposure effects on systolic blood
#' pressure (SBP) using data from the examples/input_data/ folder.
#' 
#' The data files include:
#'   - USGS_dom.csv: USGS model-based multinomial probabilities of arsenic levels
#'   - EPA_dom.csv: EPA arsenic data with lognormal parameters
#'   - Demographics_dom.txt: Demographics data with health outcomes
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
usgs_data_path <- file.path(pkg_dir, "examples", "input_data", "USGS_dom.csv")
epa_data_path <- file.path(pkg_dir, "examples", "input_data", "EPA_dom.csv")
demographics_data_path <- file.path(pkg_dir, "examples", "input_data", "Demographics_dom.txt")

# Create output directory for results
output_dir <- file.path(pkg_dir, "examples", "output_results")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Verify files exist
stopifnot(
  "USGS data file not found" = file.exists(usgs_data_path),
  "EPA data file not found" = file.exists(epa_data_path),
  "Demographics data file not found" = file.exists(demographics_data_path)
)

cat("Data files found:\n")
cat("  USGS data:", usgs_data_path, "\n")
cat("  EPA data:", epa_data_path, "\n")
cat("  Demographics data:", demographics_data_path, "\n\n")

# ==============================================================================
# Step 3: EPA Data - Automatic Lognormal Conversion
# ==============================================================================
# 
# The geoExposeR package now supports automatic conversion of raw arsenic
# concentration values to lognormal parameters. Simply specify:
#   - epa_as_conc_col: The column containing raw arsenic concentrations
#   - epa_default_sdlog: Default sdlog value (optional, defaults to 1.0)
#
# The package will automatically calculate:
#   - meanlog = log(concentration + 0.1)
#   - sdlog = epa_default_sdlog (or from epa_as_sd_col if provided)
#
# This eliminates the need for manual preprocessing!
# ==============================================================================

# Read EPA data to show available columns
epa_raw <- data.table::fread(epa_data_path)
cat("EPA data columns:", paste(names(epa_raw), collapse = ", "), "\n")
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
# This models SBP as a function of arsenic level and age decade
# geoExposeR uses mixed-effects models (lme4::lmer), so we need a random effect
# Here we use age_decade as a random effect for demonstration purposes
# In real studies, you would typically use geographic grouping (e.g., county, state)
as_level_col <- "AsLevel"
regression_formula <- paste0(
  "~ as.factor(", as_level_col, ") + ",
  "(1 | age_decade)"
)

# Arsenic category labels (matching USGS 3-category model)
as_cat_labels <- c("As<5", "As5-10", "As10+")

# Reference category for regression (typically lowest exposure)
drop_as_cat_label_reg <- c("As<5")

# Random seed for reproducibility
seed <- 12345

cat("Analysis parameters:\n")
cat("  Number of draws:", ndraws, "\n")
cat("  Target outcome(s):", paste(targets, collapse = ", "), "\n")
cat("  Regression formula:", regression_formula, "\n")
cat("  Arsenic categories:", paste(as_cat_labels, collapse = ", "), "\n")
cat("  Reference category:", drop_as_cat_label_reg, "\n")
cat("  Random seed:", seed, "\n\n")

# ==============================================================================
# Step 5: Run the Sensitivity Analysis
# ==============================================================================

cat("Starting sensitivity analysis...\n")
cat("This may take a few minutes depending on the number of draws.\n\n")

# Run the main analysis function using AUTOMATIC LOGNORMAL CONVERSION
# This approach uses the original EPA file directly - no preprocessing needed!
results <- tryCatch({
  perform_sensitivity_analysis(
    ndraws = ndraws,
    as_usgs_prob_csv = usgs_data_path,           # Original USGS file
    as_epa_prob_csv = epa_data_path,             # Original EPA file (no preprocessing!)
    birth_data_txt = demographics_data_path,     # Original demographics file
    regression_formula = regression_formula,
    output_dir = output_dir,
    targets = targets,
    as_cat_label = as_cat_labels,
    drop_as_cat_label_reg = drop_as_cat_label_reg,
    # Custom column names matching the original data files
    id_col = "id",                            # ID column in USGS/EPA data
    as_conc_cols = c("RFCprobC1", "RFCprobC2", "RFCprobC3"),  # USGS probability cols
    pop_well_col = "pub_water_pop",              # Population column in USGS
    record_id_col = "id",                     # ID column in demographics
    as_level_col = as_level_col,
    # Use raw concentration column for automatic lognormal conversion
    epa_as_conc_col = "weighted_as20062008",     # Raw concentration column
    epa_default_sdlog = 1.0,                     # Default sdlog value
    epa_pwell_col = "Dom_Pop_Percent",           # Private well percentage
    epa_cutoffs = c(5, 10),
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
      cat("- The 'estimate' column shows the pooled effect size\n")
      cat("- The 'std.error' shows the pooled standard error\n")
      cat("- The 'p.value' indicates statistical significance\n")
      cat("- Positive estimates for As5-10 or As10+ indicate higher SBP\n")
      cat("  compared to the reference category (As<5)\n\n")
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
