#' =============================================================================
#' Example: Running geoExposeR for Nitrate Exposure Analysis
#' =============================================================================
#' 
#' This script demonstrates how to adapt geoExposeR for analyzing nitrate
#' exposure effects on infant methemoglobinemia (blue baby syndrome risk).
#' 
#' BACKGROUND:
#' Nitrate contamination in groundwater is a significant public health concern,
#' particularly in agricultural areas. High nitrate levels (>10 mg/L as NO3-N)
#' can cause methemoglobinemia in infants by interfering with oxygen transport.
#' 
#' ADAPTATION NOTES:
#' - The geoExposeR framework uses "EPA" and "USGS" terminology internally,
#'   but these can represent ANY two probability data sources
#' - "EPA" data = concentration measurements from monitored water systems
#' - "USGS" data = model-predicted probability distributions (any source)
#' - Cutoffs are customized for nitrate: <5, 5-10, >=10 mg/L NO3-N
#' 
#' DATA FILES:
#'   - USGS_nitrate.csv: Multinomial probabilities of nitrate concentration categories
#'   - EPA_nitrate.csv: Nitrate concentration data with lognormal parameters
#'   - Health_nitrate.txt: Infant health data (methemoglobin levels)
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

# Data files for nitrate analysis
usgs_data_path <- file.path(pkg_dir, "examples", "input_data", "USGS_nitrate.csv")
epa_data_path <- file.path(pkg_dir, "examples", "input_data", "EPA_nitrate.csv")
health_data_path <- file.path(pkg_dir, "examples", "input_data", "Health_nitrate.txt")

# Output directory
output_dir <- file.path(pkg_dir, "examples", "output_nitrate")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Verify files exist
stopifnot(
  "USGS data file not found" = file.exists(usgs_data_path),
  "EPA data file not found" = file.exists(epa_data_path),
  "Health data file not found" = file.exists(health_data_path)
)

cat("=== Nitrate Exposure Analysis ===\n\n")
cat("Data files:\n")
cat("  Probability data (USGS):", usgs_data_path, "\n")
cat("  Concentration data (EPA):", epa_data_path, "\n")
cat("  Health outcomes:", health_data_path, "\n\n")

# ==============================================================================
# Step 3: Define Nitrate-Specific Parameters
# ==============================================================================

# Nitrate concentration cutoffs (mg/L as NO3-N)
# Category 1: <5 mg/L (safe)
# Category 2: 5-10 mg/L (elevated)
# Category 3: >=10 mg/L (exceeds MCL)
nitrate_cutoffs <- c(5, 10)

# Number of imputed datasets
ndraws <- 5

# Target health outcome
# methemoglobin_pct = percentage of hemoglobin converted to methemoglobin
targets <- c("methemoglobin_pct")

# Column containing nitrate concentration in EPA data
conc_col <- "weighted_nitrate"

# Regression formula
# Model methemoglobin as function of nitrate level
# Using infant_age_group as random effect (categorical variable)
as_level_col <- "AsLevel"  # This column name is used internally by geoExposeR
regression_formula <- paste0(
  "~ as.factor(", as_level_col, ") + ",
  "(1 | infant_age_group)"
)

# Category labels for nitrate levels (matching 3 categories from cutoffs)
as_cat_labels <- c("NO3<5", "NO3_5-10", "NO3_10+")

# Reference category for regression (lowest exposure)
drop_as_cat_label_reg <- c("NO3<5")

cat("Analysis parameters:\n")
cat("  Nitrate cutoffs:", paste(nitrate_cutoffs, collapse = ", "), "mg/L NO3-N\n")
cat("  Categories: <5 (safe), 5-10 (elevated), >=10 (exceeds MCL)\n")
cat("  Target outcome: methemoglobin_pct\n")
cat("  Number of imputations:", ndraws, "\n\n")

# ==============================================================================
# Step 4: Run the Analysis
# ==============================================================================

cat("Running sensitivity analysis...\n\n")

results <- perform_sensitivity_analysis(
  # Number of imputed datasets
  ndraws = ndraws,
  
  # Data file paths (using EPA/USGS terminology for any data sources)
  as_usgs_prob_csv = usgs_data_path,
  as_epa_prob_csv = epa_data_path,
  birth_data_txt = health_data_path,
  
  # EPA data configuration - auto-convert concentration to lognormal
  epa_as_conc_col = conc_col,
  epa_default_sdlog = 0.8,  # Typically lower variance for nitrate
  
  # Custom cutoffs for nitrate
  epa_cutoffs = nitrate_cutoffs,
  
  # Category labels for results display
  as_cat_label = as_cat_labels,
  drop_as_cat_label_reg = drop_as_cat_label_reg,
  
  # Column name mappings (adapting package parameters for nitrate data)
  id_col = "id",                                    # ID column in USGS/EPA data
  as_conc_cols = c("RFCprobC1", "RFCprobC2", "RFCprobC3"),  # Probability columns
  pop_well_col = "pub_water_pop",                   # Population column
  record_id_col = "id",                             # ID column in health data
  epa_pwell_col = "Dom_Pop_Percent",                # Private well percentage
  
  # Regression configuration
  regression_formula = regression_formula,
  as_level_col = as_level_col,
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
cat("Pooled regression results for methemoglobin percentage:\n\n")

# Access results by target name (not 'pooled_results')
if (!is.null(results) && length(results) > 0 && targets[1] %in% names(results)) {
  print(results[[targets[1]]])
  
  cat("\n\nInterpretation:\n")
  cat("- 'NO3_5-10' = Effect of nitrate 5-10 mg/L vs <5 mg/L\n")
  cat("- 'NO3_10+' = Effect of nitrate >=10 mg/L vs <5 mg/L\n")
  cat("- POSITIVE coefficients indicate higher methemoglobin levels (worse outcome)\n")
  cat("- Methemoglobin >2% is considered elevated in infants\n")
} else {
  cat("No pooled results generated. Check the output directory for details.\n")
}

cat("\nResults saved to:", output_dir, "\n")

# ==============================================================================
# Step 6: Summary Statistics
# ==============================================================================

cat("\n=== Summary Statistics ===\n")
cat("Exposure categories based on nitrate cutoffs:\n")
cat("  Category 1 (Low): <", nitrate_cutoffs[1], " mg/L NO3-N\n")
cat("  Category 2 (Elevated): ", nitrate_cutoffs[1], "-", nitrate_cutoffs[2], " mg/L NO3-N\n")
cat("  Category 3 (High): >=", nitrate_cutoffs[2], " mg/L NO3-N (exceeds EPA MCL)\n\n")
cat("EPA Maximum Contaminant Level (MCL) for Nitrate: 10 mg/L as NO3-N\n")
