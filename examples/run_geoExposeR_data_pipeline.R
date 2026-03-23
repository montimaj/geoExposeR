#' =============================================================================
#' Example: Running geoExposeR with Data Ingestion Pipeline
#' =============================================================================
#' 
#' This script demonstrates a complete data ingestion pipeline using all
#' the data preparation functions in geoExposeR:
#' 
#'   - prepare_prob_data(): Extract probabilities from GeoTIFF rasters
#'   - prepare_conc_data(): Convert concentration data to lognormal parameters
#'   - format_health_data(): Standardize health outcome data
#'   - add_well_percentages(): Merge private well population data
#'   - subset_by_geography(): Filter data by geographic region
#'   - validate_prepared_inputs(): Comprehensive input validation
#' 
#' This example creates mock data in various formats to demonstrate each
#' function's capabilities.
#' 
#' Author: Dr. Sayantan Majumdar (sayantan.majumdar@dri.edu)
#' =============================================================================

# ==============================================================================
# Step 1: Load Required Packages
# ==============================================================================

# Load geoExposeR - use devtools::load_all() if running from source
if (requireNamespace("devtools", quietly = TRUE) && 
    file.exists(file.path(getwd(), "DESCRIPTION"))) {
  devtools::load_all()
} else {
  library(geoExposeR)
}

# Load other packages
library(terra)
library(sf)
library(data.table)

cat("=== geoExposeR Data Ingestion Pipeline Example ===\n\n")

# ==============================================================================
# Step 2: Set Up Environment
# ==============================================================================

pkg_dir <- getwd()
output_dir <- file.path(pkg_dir, "examples", "output_pipeline")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Load original data for reference (from examples/input_data/)
prob_orig <- fread(file.path(pkg_dir, "examples", "input_data", "prob_dom.csv"))
conc_orig <- fread(file.path(pkg_dir, "examples", "input_data", "conc_dom.csv"))
demo_orig <- fread(file.path(pkg_dir, "examples", "input_data", "Demographics_dom.txt"))

cat("Loaded original data:\n")
cat("  Probability data:", nrow(prob_orig), "records\n")
cat("  Concentration data:", nrow(conc_orig), "records\n")
cat("  Demographics:", nrow(demo_orig), "records\n\n")

# ==============================================================================
# Step 3: Demonstrate prepare_conc_data() with Different Input Types
# ==============================================================================

cat("=== Testing prepare_conc_data() ===\n\n")

# 3a. Raw concentration data with population weighting
cat("3a. Processing raw concentration data with population weighting...\n")

conc_raw <- data.frame(
  county_fips = conc_orig$id,
  concentration = ifelse(is.na(conc_orig$weighted_as20062008),
                        runif(nrow(conc_orig), 0.5, 15),  # Fill NAs with random
                        conc_orig$weighted_as20062008),
  pop_served = ifelse(is.na(conc_orig$county_pop_served20062008),
                      sample(1000:50000, nrow(conc_orig), replace = TRUE),
                      conc_orig$county_pop_served20062008),
  pwell_pct = conc_orig$Dom_Pop_Percent
)

conc_processed_raw <- prepare_conc_data(
  conc_data = conc_raw,
  id_col = "county_fips",
  concentration_col = "concentration",
  population_col = "pop_served",
  private_well_pct_col = "pwell_pct",
  input_type = "raw"
)

cat("  Result columns:", paste(names(conc_processed_raw), collapse = ", "), "\n")
cat("  Number of geographic units:", nrow(conc_processed_raw), "\n\n")

# 3b. Pre-aggregated data (already has mean/sd)
cat("3b. Processing pre-aggregated concentration data...\n")

conc_aggregated <- data.frame(
  geoid = conc_orig$id,
  mean_conc = runif(nrow(conc_orig), 1, 10),
  sd_conc = runif(nrow(conc_orig), 0.5, 3),
  pwell = conc_orig$Dom_Pop_Percent
)

conc_processed_agg <- prepare_conc_data(
  conc_data = conc_aggregated,
  id_col = "geoid",
  mean_col = "mean_conc",
  sd_col = "sd_conc",
  private_well_pct_col = "pwell",
  input_type = "summary"
)

cat("  Result has meanlog and sdlog columns:",
    all(c("meanlog", "sdlog") %in% names(conc_processed_agg)), "\n\n")

# ==============================================================================
# Step 4: Demonstrate format_health_data() with Various Formats
# ==============================================================================

cat("=== Testing format_health_data() ===\n\n")

# 4a. Numeric FIPS codes (need padding)
cat("4a. Formatting health data with numeric FIPS codes...\n")

health_numeric <- data.frame(
  county_id = demo_orig$id,
  blood_pressure = demo_orig$SBP,
  age_group = demo_orig$age_decade
)

health_formatted <- format_health_data(
  health_data = health_numeric,
  id_col = "county_id",
  outcome_cols = c("blood_pressure"),
  covariate_cols = c("age_group"),
  geoid_digits = 5,
  validate = TRUE
)

cat("  Original county_id type:", class(health_numeric$county_id), "\n")
cat("  Formatted FIPS example:", paste(head(health_formatted$FIPS, 3), collapse = ", "), "\n\n")

# 4b. Health data with multiple outcomes
cat("4b. Formatting health data with multiple outcomes...\n")

health_multi <- data.frame(
  fips = sprintf("%05d", demo_orig$id),
  systolic_bp = demo_orig$SBP,
  diastolic_bp = demo_orig$SBP - 40 + rnorm(nrow(demo_orig), 0, 5),
  heart_rate = 70 + rnorm(nrow(demo_orig), 0, 10),
  age = demo_orig$age_decade
)

health_multi_formatted <- format_health_data(
  health_data = health_multi,
  id_col = "fips",
  outcome_cols = c("systolic_bp", "diastolic_bp", "heart_rate"),
  covariate_cols = c("age"),
  validate = FALSE  # Skip validation for demo
)

cat("  Formatted columns:", paste(names(health_multi_formatted), collapse = ", "), "\n\n")

# ==============================================================================
# Step 5: Demonstrate add_well_percentages()
# ==============================================================================

cat("=== Testing add_well_percentages() ===\n\n")

# Create base dataset
base_data <- data.frame(
  GEOID = sprintf("%05d", prob_orig$id),
  prob_low = prob_orig$RFCprobC1,
  prob_med = prob_orig$RFCprobC2,
  prob_high = prob_orig$RFCprobC3
)

# Create well percentage dataset (simulating Census data)
well_data <- data.frame(
  fips_code = sprintf("%05d", prob_orig$id),
  total_pop = sample(10000:500000, nrow(prob_orig), replace = TRUE),
  well_pop = sample(500:50000, nrow(prob_orig), replace = TRUE)
)
well_data$well_pct <- well_data$well_pop / well_data$total_pop * 100

# Merge using add_well_percentages()
merged_with_wells <- add_well_percentages(
  target_data = base_data,
  well_data = well_data,
  target_id_col = "GEOID",
  well_id_col = "fips_code",
  well_pct_col = "well_pct"
)

cat("  Original columns:", ncol(base_data), "\n")
cat("  After adding wells:", ncol(merged_with_wells), "\n")
cat("  New column 'PWELL_private_pct' added:", "PWELL_private_pct" %in% names(merged_with_wells), "\n\n")

# ==============================================================================
# Step 6: Demonstrate subset_by_geography()
# ==============================================================================

cat("=== Testing subset_by_geography() ===\n\n")

# Create dataset with state information (using FIPS prefix to denote state)
full_data <- data.frame(
  FIPS = c(
    sprintf("06%03d", 1:20),   # California (state FIPS 06)
    sprintf("32%03d", 1:20),   # Nevada (state FIPS 32)
    sprintf("04%03d", 1:20),   # Arizona (state FIPS 04)
    sprintf("41%03d", 1:20),   # Oregon (state FIPS 41)
    sprintf("53%03d", 1:20)    # Washington (state FIPS 53)
  ),
  value = rnorm(100)
)

# 6a. Subset by state (CA and NV using state FIPS codes)
cat("6a. Subsetting by state FIPS codes (CA=06, NV=32)...\n")
subset_states <- subset_by_geography(
  data = full_data,
  id_col = "FIPS",
  state_fips = c("06", "32")
)
cat("  Original rows:", nrow(full_data), "\n")
cat("  After subsetting:", nrow(subset_states), "\n\n")

# 6b. Subset by specific FIPS list
cat("6b. Subsetting by specific FIPS list...\n")
subset_fips <- subset_by_geography(
  data = full_data,
  id_col = "FIPS",
  geoid_list = c("06001", "06002", "06003", "32001", "32002")
)
cat("  Matching FIPS codes:", nrow(subset_fips), "\n\n")

# ==============================================================================
# Step 7: Demonstrate validate_prepared_inputs()
# ==============================================================================

cat("=== Testing validate_prepared_inputs() ===\n\n")

# Create properly formatted files for validation
temp_dir <- tempdir()

# Probability file
prob_valid <- data.frame(
  GEOID10 = sprintf("%05d", prob_orig$id),
  prob_C1 = prob_orig$RFCprobC1,
  prob_C2 = prob_orig$RFCprobC2,
  prob_C3 = prob_orig$RFCprobC3,
  Wells_2010 = ifelse(is.na(prob_orig$pub_water_pop), 1000, prob_orig$pub_water_pop),
  PWELL = prob_orig$Dom_Pop_Percent / 100
)
prob_path <- file.path(temp_dir, "prob_valid.csv")
fwrite(prob_valid, prob_path)

# Concentration file
conc_valid <- data.frame(
  GEOID10 = sprintf("%05d", conc_orig$id),
  conc_meanlog = log(pmax(0.1, ifelse(is.na(conc_orig$weighted_as20062008),
                                         1, conc_orig$weighted_as20062008))),
  conc_sdlog = 1.0,
  PWELL_private_pct = conc_orig$Dom_Pop_Percent / 100
)
conc_path <- file.path(temp_dir, "conc_valid.csv")
fwrite(conc_valid, conc_path)

# Health file
health_valid <- data.frame(
  FIPS = sprintf("%05d", demo_orig$id),
  SBP = demo_orig$SBP,
  age_decade = demo_orig$age_decade
)
health_path <- file.path(temp_dir, "health_valid.csv")
fwrite(health_valid, health_path)

# Run validation
cat("Validating prepared input files...\n")
validation_result <- validate_prepared_inputs(
  prob_file = prob_path,
  conc_file = conc_path,
  health_file = health_path
)

cat("\nValidation Results:\n")
cat("  Valid:", validation_result$valid, "\n")
cat("  Issues:", length(validation_result$issues), "\n")
cat("  Warnings:", length(validation_result$warnings), "\n")

if (length(validation_result$warnings) > 0) {
  cat("\n  Warnings:\n")
  for (w in validation_result$warnings) {
    cat("    -", w, "\n")
  }
}

cat("\n  Summary:\n")
cat("    Probability records:", validation_result$summary$prob$n_rows, "\n")
cat("    Concentration records:", validation_result$summary$conc$n_rows, "\n")
cat("    Health records:", validation_result$summary$health$n_rows, "\n\n")

# ==============================================================================
# Step 8: Run Analysis with Validated Data
# ==============================================================================

cat("=== Running Analysis with Validated Data ===\n\n")

# Merge probability and concentration data for analysis
prob_for_analysis <- fread(prob_path)
setnames(prob_for_analysis, "GEOID10", "id")
prob_for_analysis[, id := as.integer(id)]
setnames(prob_for_analysis, "prob_C1", "RFCprobC1")
setnames(prob_for_analysis, "prob_C2", "RFCprobC2")
setnames(prob_for_analysis, "prob_C3", "RFCprobC3")

conc_for_analysis <- fread(conc_path)
setnames(conc_for_analysis, "GEOID10", "id")
conc_for_analysis[, id := as.integer(id)]

merged <- merge(prob_for_analysis, conc_for_analysis[, .(id, conc_meanlog, conc_sdlog)],
                by = "id", all.x = TRUE)
merged[, Dom_Pop_Percent := PWELL * 100]

merged_path <- file.path(output_dir, "merged_data.csv")
fwrite(merged, merged_path)

# Prepare health data
health_for_analysis <- fread(health_path)
setnames(health_for_analysis, "FIPS", "id")
health_for_analysis[, id := as.integer(id)]
health_analysis_path <- file.path(output_dir, "health_data.txt")
fwrite(health_for_analysis, health_analysis_path, sep = "\t")

# Run analysis
cat("Running sensitivity analysis...\n")

results <- tryCatch({
  perform_sensitivity_analysis(
    ndraws = 5,
    model_prob_csv = merged_path,
    conc_data_csv = merged_path,
    birth_data_txt = health_analysis_path,
    regression_formula = "~ as.factor(ExposureLevel) + (1 | age_decade)",
    output_dir = output_dir,
    targets = c("SBP"),
    cat_label = c("<5", "5-10", "10+"),
    drop_cat_label_ref = c("<5"),
    id_col = "id",
    prob_cols = c("RFCprobC1", "RFCprobC2", "RFCprobC3"),
    pop_well_col = "Wells_2010",
    record_id_col = "id",
    exposure_level_col = "ExposureLevel",
    conc_mean_col = "conc_meanlog",
    conc_sd_col = "conc_sdlog",
    pwell_col = "Dom_Pop_Percent",
    conc_cutoffs = c(5, 10),
    seed = 42,
    mice_m = 1,
    mice_maxit = 5,
    apply_imputation_fallback = TRUE
  )
}, error = function(e) {
  cat("Error:", conditionMessage(e), "\n")
  return(NULL)
})

# ==============================================================================
# Step 9: Display Results Summary
# ==============================================================================

if (!is.null(results)) {
  cat("\n========================================\n")
  cat("Pipeline Analysis Complete!\n")
  cat("========================================\n\n")
  
  cat("SBP Analysis Results:\n")
  print(results$SBP)
  
  cat("\n\nOutput files saved to:", output_dir, "\n")
  for (f in list.files(output_dir)) {
    cat("  -", f, "\n")
  }
} else {
  cat("\nAnalysis encountered an error.\n")
}

# Cleanup
unlink(c(prob_path, conc_path, health_path))

cat("\n=== Data Ingestion Pipeline Example Complete ===\n")
