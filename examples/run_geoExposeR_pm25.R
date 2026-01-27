#' =============================================================================
#' Example: Running geoExposeR for PM2.5 Air Quality Analysis
#' =============================================================================
#' 
#' This script demonstrates how to adapt geoExposeR for analyzing PM2.5
#' (fine particulate matter) exposure effects on respiratory function.
#' 
#' BACKGROUND:
#' PM2.5 refers to particulate matter with diameter less than 2.5 micrometers.
#' These fine particles can penetrate deep into lungs and bloodstream, causing
#' respiratory and cardiovascular health effects. The EPA annual standard is
#' 9 µg/m³ (updated in 2024), while WHO recommends <5 µg/m³.
#'
#' ADAPTATION NOTES:
#' - The geoExposeR framework uses "EPA" and "USGS" terminology internally,
#'   but these represent ANY two probability data sources:
#'   - "USGS" data = Model-predicted probabilities (satellite/hybrid models)
#'   - "EPA" data = Ground monitoring station measurements
#' - This example shows air quality data (not groundwater) demonstrating
#'   the package's broader applicability
#' - Cutoffs are customized for PM2.5: <9, 9-35, >=35 µg/m³
#' 
#' DATA FILES:
#'   - Model_pm25.csv: Satellite/model-based PM2.5 probability estimates
#'   - EPA_pm25.csv: Ground monitoring station measurements
#'   - Health_pm25.txt: Respiratory function data (FEV1 % predicted)
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

# Data files for PM2.5 analysis
# Note: We use "USGS" parameter for model predictions and "EPA" for ground monitors
model_data_path <- file.path(pkg_dir, "examples", "input_data", "Model_pm25.csv")
monitor_data_path <- file.path(pkg_dir, "examples", "input_data", "EPA_pm25.csv")
health_data_path <- file.path(pkg_dir, "examples", "input_data", "Health_pm25.txt")

# Output directory
output_dir <- file.path(pkg_dir, "examples", "output_pm25")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Verify files exist
stopifnot(
  "Model data file not found" = file.exists(model_data_path),
  "Monitor data file not found" = file.exists(monitor_data_path),
  "Health data file not found" = file.exists(health_data_path)
)

cat("=== PM2.5 Air Quality Exposure Analysis ===\n\n")
cat("Data files:\n")
cat("  Satellite/model predictions:", model_data_path, "\n")
cat("  Ground monitoring data:", monitor_data_path, "\n")
cat("  Health outcomes:", health_data_path, "\n\n")

# ==============================================================================
# Step 3: Define PM2.5-Specific Parameters
# ==============================================================================

# PM2.5 concentration cutoffs (µg/m³)
# Category 1: <9 µg/m³ (meets EPA annual standard, updated 2024)
# Category 2: 9-35 µg/m³ (exceeds annual standard, below 24-hr standard)
# Category 3: >=35 µg/m³ (exceeds EPA 24-hour standard)
pm25_cutoffs <- c(9, 35)

# Number of imputed datasets
ndraws <- 5

# Target health outcome
# FEV1_pct_predicted = Forced Expiratory Volume in 1 second (% of predicted)
# Lower values indicate worse respiratory function
targets <- c("FEV1_pct_predicted")

# Column containing PM2.5 concentration in monitoring data
conc_col <- "weighted_pm25"

# Regression formula
# Model FEV1 as function of PM2.5 level
# Using smoking_status as random effect (categorical variable)
as_level_col <- "AsLevel"  # This column name is used internally by geoExposeR
regression_formula <- paste0(
  "~ as.factor(", as_level_col, ") + ",
  "(1 | smoking_status)"
)

# Category labels for PM2.5 levels (matching 3 categories from cutoffs)
as_cat_labels <- c("PM25<9", "PM25_9-35", "PM25_35+")

# Reference category for regression (lowest exposure)
drop_as_cat_label_reg <- c("PM25<9")

cat("Analysis parameters:\n")
cat("  PM2.5 cutoffs:", paste(pm25_cutoffs, collapse = ", "), "µg/m³\n")
cat("  Categories: <9 (meets annual standard), 9-35 (moderate), >=35 (exceeds 24-hr)\n")
cat("  Target outcome: FEV1_pct_predicted (lung function)\n")
cat("  Number of imputations:", ndraws, "\n\n")

# ==============================================================================
# Step 4: Run the Analysis
# ==============================================================================

cat("Running sensitivity analysis...\n\n")

results <- perform_sensitivity_analysis(
  # Number of imputed datasets
  ndraws = ndraws,
  
  # Data file paths
  # as_usgs_prob_csv = model predictions (satellite/hybrid)
  # as_epa_prob_csv = ground monitor measurements
  as_usgs_prob_csv = model_data_path,
  as_epa_prob_csv = monitor_data_path,
  birth_data_txt = health_data_path,
  
  # EPA (monitor) data configuration - auto-convert concentration to lognormal
  epa_as_conc_col = conc_col,
  epa_default_sdlog = 0.5,  # PM2.5 typically shows lower variance than groundwater contaminants
  
  # Custom cutoffs for PM2.5
  epa_cutoffs = pm25_cutoffs,
  
  # Category labels for results display
  as_cat_label = as_cat_labels,
  drop_as_cat_label_reg = drop_as_cat_label_reg,
  
  # Column name mappings (adapting package parameters for PM2.5 data)
  id_col = "id",                                    # ID column in model/monitor data
  as_conc_cols = c("RFCprobC1", "RFCprobC2", "RFCprobC3"),  # Probability columns
  pop_well_col = "pub_water_pop",                   # Population column
  record_id_col = "id",                             # ID column in health data
  epa_pwell_col = "Dom_Pop_Percent",                # Weighting percentage
  
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
cat("Pooled regression results for FEV1 % predicted (lung function):\n\n")

# Access results by target name (not 'pooled_results')
if (!is.null(results) && length(results) > 0 && targets[1] %in% names(results)) {
  print(results[[targets[1]]])
  
  cat("\n\nInterpretation:\n")
  cat("- 'PM25_9-35' = Effect of PM2.5 9-35 µg/m³ vs <9 µg/m³\n")
  cat("- 'PM25_35+' = Effect of PM2.5 >=35 µg/m³ vs <9 µg/m³\n")
  cat("- NEGATIVE coefficients indicate decreased lung function (worse outcome)\n")
  cat("- FEV1 <80% predicted often indicates obstructive lung disease\n")
} else {
  cat("No pooled results generated. Check the output directory for details.\n")
}

cat("\nResults saved to:", output_dir, "\n")

# ==============================================================================
# Step 6: Summary Statistics
# ==============================================================================

cat("\n=== Summary Statistics ===\n")
cat("Exposure categories based on PM2.5 cutoffs:\n")
cat("  Category 1 (Good): <", pm25_cutoffs[1], " µg/m³ (meets EPA annual standard)\n")
cat("  Category 2 (Moderate): ", pm25_cutoffs[1], "-", pm25_cutoffs[2], " µg/m³\n")
cat("  Category 3 (High): >=", pm25_cutoffs[2], " µg/m³ (unhealthy)\n\n")
cat("Reference Standards:\n")
cat("  EPA Annual Standard: 9 µg/m³ (updated 2024)\n")
cat("  WHO Guideline: 5 µg/m³ (annual mean)\n")
cat("  EPA 24-hour Standard: 35 µg/m³\n\n")
cat("FEV1 Interpretation (% predicted):\n")
cat("  >=80%: Normal lung function\n")
cat("  60-79%: Mild obstruction\n")
cat("  40-59%: Moderate obstruction\n")
cat("  <40%: Severe obstruction\n")

# ==============================================================================
# Step 7: Adaptation Notes
# ==============================================================================

cat("\n=== Adaptation Notes for Air Quality Studies ===\n")
cat("This example demonstrates geoExposeR's broader applicability beyond\n")
cat("groundwater contaminants. Key adaptations include:\n\n")
cat("1. DATA SOURCES:\n")
cat("   - 'USGS' parameter → Satellite/hybrid model predictions\n")
cat("   - 'EPA' parameter → Ground monitoring station measurements\n\n")
cat("2. POPULATION WEIGHTING:\n")
cat("   - Dom_Pop_Percent can represent rural/urban population fraction\n")
cat("   - Or monitoring station density weighting\n\n")
cat("3. CUTOFFS:\n")
cat("   - Customize epa_cutoffs for your pollutant's health thresholds\n")
cat("   - Consider regulatory standards (EPA, WHO) as category boundaries\n\n")
cat("4. HEALTH OUTCOMES:\n")
cat("   - Respiratory: FEV1, asthma incidence, COPD exacerbations\n")
cat("   - Cardiovascular: blood pressure, heart rate variability\n")
cat("   - Mortality: all-cause, cause-specific\n")
