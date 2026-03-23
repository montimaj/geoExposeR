#' =============================================================================
#' Example: Running geoExposeR for Uranium Exposure Analysis
#' =============================================================================
#' 
#' This script demonstrates how to adapt geoExposeR for analyzing uranium
#' exposure effects on kidney function (estimated Glomerular Filtration Rate).
#' 
#' BACKGROUND:
#' Uranium is a naturally occurring radioactive element found in groundwater,
#' particularly in areas with certain geological formations. Chronic uranium
#' exposure is associated with nephrotoxicity (kidney damage), measured by
#' decreased eGFR (estimated Glomerular Filtration Rate).
#' 
#' ADAPTATION NOTES:
#' - Cutoffs are customized for uranium: <15, 15-30, >=30 µg/L
#'
#' DATA FILES:
#'   - prob_uranium.csv: Multinomial probabilities of uranium concentration categories
#'   - conc_uranium.csv: Uranium concentration data with lognormal parameters
#'   - Health_uranium.txt: Kidney function data (eGFR measurements)
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

# Data files for uranium analysis
prob_data_path <- file.path(pkg_dir, "examples", "input_data", "prob_uranium.csv")
conc_data_path <- file.path(pkg_dir, "examples", "input_data", "conc_uranium.csv")
health_data_path <- file.path(pkg_dir, "examples", "input_data", "Health_uranium.txt")

# Output directory
output_dir <- file.path(pkg_dir, "examples", "output_uranium")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Verify files exist
stopifnot(
  "Probability data file not found" = file.exists(prob_data_path),
  "Concentration data file not found" = file.exists(conc_data_path),
  "Health data file not found" = file.exists(health_data_path)
)

cat("=== Uranium Exposure Analysis ===\n\n")
cat("Data files:\n")
cat("  Probability data:", prob_data_path, "\n")
cat("  Concentration data:", conc_data_path, "\n")
cat("  Health outcomes:", health_data_path, "\n\n")

# ==============================================================================
# Step 3: Define Uranium-Specific Parameters
# ==============================================================================

# Uranium concentration cutoffs (µg/L)
# Category 1: <15 µg/L (below half of MCL)
# Category 2: 15-30 µg/L (approaching MCL)
# Category 3: >=30 µg/L (exceeds MCL)
uranium_cutoffs <- c(15, 30)

# Number of imputed datasets
ndraws <- 5

# Target health outcome
# eGFR = estimated Glomerular Filtration Rate (mL/min/1.73m²)
# Lower eGFR indicates worse kidney function
targets <- c("eGFR")

# Column containing uranium concentration in concentration data
conc_col <- "weighted_uranium"

# Regression formula
# Model eGFR as function of uranium level
# Using sex as random effect (categorical variable)
exposure_level_col <- "ExposureLevel"  # This column name is used internally by geoExposeR
regression_formula <- paste0(
  "~ as.factor(", exposure_level_col, ") + ",
  "(1 | sex)"
)

# Category labels for uranium levels (matching 3 categories from cutoffs)
cat_labels <- c("U<15", "U_15-30", "U_30+")

# Reference category for regression (lowest exposure)
drop_cat_label_ref <- c("U<15")

cat("Analysis parameters:\n")
cat("  Uranium cutoffs:", paste(uranium_cutoffs, collapse = ", "), "µg/L\n")
cat("  Categories: <15 (low), 15-30 (moderate), >=30 (exceeds MCL)\n")
cat("  Target outcome: eGFR (kidney function)\n")
cat("  Number of imputations:", ndraws, "\n\n")

# ==============================================================================
# Step 4: Run the Analysis
# ==============================================================================

cat("Running sensitivity analysis...\n\n")

results <- perform_sensitivity_analysis(
  # Number of imputed datasets
  ndraws = ndraws,
  
  # Data file paths
  model_prob_csv = prob_data_path,
  conc_data_csv = conc_data_path,
  birth_data_txt = health_data_path,
  
  # Concentration data configuration - auto-convert concentration to lognormal
  conc_raw_col = conc_col,
  default_sdlog = 0.9,  # Uranium typically shows moderate variance
  
  # Custom cutoffs for uranium
  conc_cutoffs = uranium_cutoffs,
  
  # Category labels for results display
  cat_label = cat_labels,
  drop_cat_label_ref = drop_cat_label_ref,
  
  # Column name mappings for uranium data
  id_col = "id",                                    # ID column in prob/conc data
  prob_cols = c("RFCprobC1", "RFCprobC2", "RFCprobC3"),  # Probability columns
  pop_well_col = "pub_water_pop",                   # Population column
  record_id_col = "id",                             # ID column in health data
  pwell_col = "Dom_Pop_Percent",                # Private well percentage
  
  # Regression configuration
  regression_formula = regression_formula,
  exposure_level_col = exposure_level_col,
  targets = targets,
  
  # Output configuration
  output_dir = output_dir,
  
  # Identifier formatting
  format_ids = TRUE,
  
  # Imputation settings
  seed = 12345,
  mice_m = 1,
  mice_maxit = 5,
  mice_method = "pmm",
  apply_imputation_fallback = TRUE
)

# ==============================================================================
# Step 5: Display Results
# ==============================================================================

cat("\n=== Analysis Complete ===\n\n")
cat("Pooled regression results for eGFR (kidney function):\n\n")

# Access results by target name (not 'pooled_results')
if (!is.null(results) && length(results) > 0 && targets[1] %in% names(results)) {
  print(results[[targets[1]]])
  
  cat("\n\nInterpretation:\n")
  cat("- 'U_15-30' = Effect of uranium 15-30 µg/L vs <15 µg/L\n")
  cat("- 'U_30+' = Effect of uranium >=30 µg/L vs <15 µg/L\n")
  cat("- NEGATIVE coefficients indicate decreased kidney function (worse outcome)\n")
  cat("- eGFR <60 indicates chronic kidney disease (CKD)\n")
} else {
  cat("No pooled results generated. Check the output directory for details.\n")
}

cat("\nResults saved to:", output_dir, "\n")

# ==============================================================================
# Step 6: Summary Statistics
# ==============================================================================

cat("\n=== Summary Statistics ===\n")
cat("Exposure categories based on uranium cutoffs:\n")
cat("  Category 1 (Low): <", uranium_cutoffs[1], " µg/L\n")
cat("  Category 2 (Moderate): ", uranium_cutoffs[1], "-", uranium_cutoffs[2], " µg/L\n")
cat("  Category 3 (High): >=", uranium_cutoffs[2], " µg/L (exceeds EPA MCL)\n\n")
cat("EPA Maximum Contaminant Level (MCL) for Uranium: 30 µg/L\n")
cat("eGFR Categories:\n")
cat("  >=90: Normal kidney function\n")
cat("  60-89: Mildly decreased\n")
cat("  30-59: Moderately decreased (CKD Stage 3)\n")
cat("  15-29: Severely decreased (CKD Stage 4)\n")
cat("  <15: Kidney failure (CKD Stage 5)\n")
