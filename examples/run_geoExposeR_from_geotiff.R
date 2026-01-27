#' =============================================================================
#' Example: Running geoExposeR from GeoTIFF Rasters
#' =============================================================================
#' 
#' This script demonstrates the complete workflow for running geoExposeR
#' starting from raw GeoTIFF rasters containing USGS arsenic probability data.
#' 
#' Data Citation:
#'   Lombard, M.A., 2021, Data used to model and map arsenic concentration
#'   exceedances in private wells throughout the conterminous United States
#'   for human health studies: U.S. Geological Survey data release,
#'   https://doi.org/10.5066/P90RBJXS.
#'   Ravalli, F., Yu, Y., Bostick, B. C., Chillrud, S. N., Schilling, K., Basu, A., 
#'   Navas-Acien, A., Nigra, A. E., 2022, Sociodemographic inequalities in uranium and other metals 
#'   in community water systems across the USA, 2006-11: a cross-sectional study: a cross-sectional study. 
#'   The Lancet Planetary Health, 6(4), e320–e330.
#'   https://doi.org/10.1016/S2542-5196(22)00043-2
#' 
#' The workflow includes:
#'   1. Creating mock GeoTIFF rasters (simulating USGS probability surfaces)
#'   2. Extracting probabilities by geographic unit using prepare_usgs_data()
#'   3. Preparing EPA data using prepare_epa_data()
#'   4. Formatting health data using format_health_data()
#'   5. Running the sensitivity analysis
#' 
#' This example uses the same underlying data as run_geoExposeR_example.R
#' but demonstrates the GeoTIFF → CSV conversion workflow.
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

# Load spatial packages
library(terra)
library(sf)

cat("=== geoExposeR GeoTIFF Workflow Example ===\n\n")

# ==============================================================================
# Step 2: Set Up Directories and Paths
# ==============================================================================

pkg_dir <- getwd()
output_dir <- file.path(pkg_dir, "examples", "output_geotiff")
raster_dir <- file.path(output_dir, "rasters")

# Create output directories
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(raster_dir)) dir.create(raster_dir, recursive = TRUE)

# Original data paths (from examples/input_data/)
usgs_csv_path <- file.path(pkg_dir, "examples", "input_data", "USGS_dom.csv")
epa_csv_path <- file.path(pkg_dir, "examples", "input_data", "EPA_dom.csv")
demographics_path <- file.path(pkg_dir, "examples", "input_data", "Demographics_dom.txt")

cat("Output directory:", output_dir, "\n")
cat("Raster directory:", raster_dir, "\n\n")

# ==============================================================================
# Step 3: Create Mock GeoTIFF Rasters from CSV Data
# ==============================================================================
#
# In a real workflow, you would have actual USGS GeoTIFF probability surfaces.
# Here we create mock rasters that contain the same probability values as the
# CSV data, demonstrating the workflow.
#
# USGS typically provides rasters with exceedance probabilities:
#   - prob_gt5.tif: P(As > 5 µg/L)
#   - prob_gt10.tif: P(As > 10 µg/L)
# ==============================================================================

cat("Step 3: Creating mock GeoTIFF rasters...\n")

# Read the original USGS data to get probability values
usgs_data <- data.table::fread(usgs_csv_path)
cat("  Loaded USGS data with", nrow(usgs_data), "geographic units\n")

# Create a simple raster grid covering a hypothetical region
# Each pixel represents one geographic unit from the CSV data
n_units <- nrow(usgs_data)
ncols <- ceiling(sqrt(n_units))
nrows <- ceiling(n_units / ncols)

# Create rasters for P(As > 5) and P(As > 10)
# From the 3-category model:
#   P(As > 5) = RFCprobC2 + RFCprobC3
#   P(As > 10) = RFCprobC3

# Pad data to fill raster grid
prob_gt5_values <- c(
  usgs_data$RFCprobC2 + usgs_data$RFCprobC3,
  rep(NA, nrows * ncols - n_units)
)
prob_gt10_values <- c(
  usgs_data$RFCprobC3,
  rep(NA, nrows * ncols - n_units)
)

# Create rasters with arbitrary extent
r_gt5 <- terra::rast(
  nrows = nrows, ncols = ncols,
  xmin = -120, xmax = -110,
  ymin = 35, ymax = 45
)
terra::values(r_gt5) <- prob_gt5_values
terra::crs(r_gt5) <- "EPSG:4326"

r_gt10 <- terra::rast(
  nrows = nrows, ncols = ncols,
  xmin = -120, xmax = -110,
  ymin = 35, ymax = 45
)
terra::values(r_gt10) <- prob_gt10_values
terra::crs(r_gt10) <- "EPSG:4326"

# Save rasters as GeoTIFF
gt5_path <- file.path(raster_dir, "prob_gt5.tif")
gt10_path <- file.path(raster_dir, "prob_gt10.tif")
terra::writeRaster(r_gt5, gt5_path, overwrite = TRUE)
terra::writeRaster(r_gt10, gt10_path, overwrite = TRUE)

cat("  Created GeoTIFF rasters:\n")
cat("    - prob_gt5.tif: P(As > 5 µg/L)\n")
cat("    - prob_gt10.tif: P(As > 10 µg/L)\n\n")

# ==============================================================================
# Step 4: Create Geographic Boundaries (Mock County Polygons)
# ==============================================================================
#
# In a real workflow, you would load actual county/tract boundaries from a
# shapefile or Census TIGER files. Here we create mock polygons matching
# the raster pixels.
# ==============================================================================

cat("Step 4: Creating geographic boundaries...\n")

# Create a polygon for each geographic unit (matching raster pixels)
# Calculate pixel dimensions
pixel_width <- (terra::xmax(r_gt5) - terra::xmin(r_gt5)) / ncols
pixel_height <- (terra::ymax(r_gt5) - terra::ymin(r_gt5)) / nrows

polygons <- list()
valid_idx <- 0

for (i in 1:nrows) {
  for (j in 1:ncols) {
    idx <- (i - 1) * ncols + j
    if (idx <= n_units) {
      valid_idx <- valid_idx + 1
      
      # Calculate polygon corners
      xmin <- terra::xmin(r_gt5) + (j - 1) * pixel_width
      xmax <- xmin + pixel_width
      ymax <- terra::ymax(r_gt5) - (i - 1) * pixel_height
      ymin <- ymax - pixel_height
      
      # Create polygon
      coords <- matrix(c(
        xmin, ymin,
        xmax, ymin,
        xmax, ymax,
        xmin, ymax,
        xmin, ymin
      ), ncol = 2, byrow = TRUE)
      
      polygons[[valid_idx]] <- sf::st_polygon(list(coords))
    }
  }
}

# Create sf object with GEOID
boundaries <- sf::st_as_sf(
  data.frame(GEOID = as.character(usgs_data$id)),
  geometry = sf::st_sfc(polygons, crs = 4326)
)

# Add population data
boundaries$private_well_pop <- usgs_data$pub_water_pop

cat("  Created", nrow(boundaries), "polygon boundaries\n\n")

# ==============================================================================
# Step 5: Extract USGS Probabilities Using prepare_usgs_data()
# ==============================================================================

cat("Step 5: Extracting USGS probabilities from GeoTIFFs...\n")

# Define raster paths
prob_rasters <- list(
  gt5 = gt5_path,
  gt10 = gt10_path
)

# Prepare population data
pop_data <- data.frame(
  GEOID = as.character(usgs_data$id),
  private_well_pop = usgs_data$pub_water_pop
)

# Extract and convert to multinomial probabilities
usgs_processed_path <- file.path(output_dir, "usgs_from_geotiff.csv")

usgs_processed <- prepare_usgs_data(
  prob_rasters = prob_rasters,
  boundaries = boundaries,
  id_col = "GEOID",
  pop_data = pop_data,
  pop_id_col = "GEOID",
  pop_well_col = "private_well_pop",
  extraction_method = "centroid",  # Using centroid for mock pixel data
  cutoffs = c(5, 10),
  output_file = usgs_processed_path
)

cat("  USGS data extracted and saved to:", usgs_processed_path, "\n")
cat("  Columns:", paste(names(usgs_processed), collapse = ", "), "\n\n")

# ==============================================================================
# Step 6: Prepare EPA Data Using prepare_epa_data()
# ==============================================================================

cat("Step 6: Preparing EPA data...\n")

# Read raw EPA data
epa_raw <- data.table::fread(epa_csv_path)

# Use prepare_epa_data() to convert to lognormal parameters
# The function outputs columns: GEOID10, EPA_AS_meanlog, EPA_AS_sdlog, PWELL_private_pct
epa_processed <- prepare_epa_data(
  epa_data = epa_raw,
  id_col = "id",
  concentration_col = "weighted_as20062008",
  population_col = "county_pop_served20062008",
  private_well_pct_col = "Dom_Pop_Percent",
  input_type = "raw"
)

# Add percentage conversion for compatibility
epa_processed[, Dom_Pop_Percent := PWELL_private_pct]  # Keep original percentage

epa_processed_path <- file.path(output_dir, "epa_processed.csv")
data.table::fwrite(epa_processed, epa_processed_path)

cat("  EPA data processed and saved to:", epa_processed_path, "\n\n")

# ==============================================================================
# Step 7: Format Health Data Using format_health_data()
# ==============================================================================

cat("Step 7: Formatting health data...\n")

# Read raw health data
health_raw <- data.table::fread(demographics_path)

# Use format_health_data() to standardize
health_processed <- format_health_data(
  health_data = health_raw,
  id_col = "id",
  outcome_cols = c("SBP"),
  covariate_cols = c("age_decade"),
  geoid_digits = 5,
  validate = TRUE
)

health_processed_path <- file.path(output_dir, "health_processed.csv")
data.table::fwrite(health_processed, health_processed_path)

cat("  Health data formatted and saved to:", health_processed_path, "\n\n")

# ==============================================================================
# Step 8: Validate All Prepared Inputs
# ==============================================================================

cat("Step 8: Validating prepared inputs...\n")

validation <- validate_prepared_inputs(
  usgs_file = usgs_processed_path,
  epa_file = epa_processed_path,
  health_file = health_processed_path
)

if (validation$valid) {
  cat("  ✓ All inputs validated successfully\n\n")
} else {
  cat("  ✗ Validation issues found:\n")
  for (issue in validation$issues) {
    cat("    -", issue, "\n")
  }
  cat("\n")
}

# ==============================================================================
# Step 9: Run the Sensitivity Analysis
# ==============================================================================

cat("Step 9: Running sensitivity analysis...\n")

# Merge USGS and EPA data for the analysis
# The perform_sensitivity_analysis expects specific column formats
usgs_for_analysis <- data.table::fread(usgs_processed_path)
epa_for_analysis <- data.table::fread(epa_processed_path)

# Add EPA columns to USGS data
# Note: epa_for_analysis may have different column names based on the prepare_epa_data output
epa_cols <- c("GEOID10", "EPA_AS_meanlog", "EPA_AS_sdlog")
if ("PWELL_private_pct" %in% names(epa_for_analysis)) {
  epa_cols <- c(epa_cols, "PWELL_private_pct")
}
if ("Dom_Pop_Percent" %in% names(epa_for_analysis)) {
  epa_cols <- c(epa_cols, "Dom_Pop_Percent")
}

merged_data <- merge(
  usgs_for_analysis,
  epa_for_analysis[, ..epa_cols],
  by = "GEOID10",
  all.x = TRUE
)

# Rename for compatibility with original example
data.table::setnames(merged_data, "GEOID10", "id")
data.table::setnames(merged_data, "RFC3_C1v2", "RFCprobC1")
data.table::setnames(merged_data, "RFC3_C2v2", "RFCprobC2")
data.table::setnames(merged_data, "RFC3_C3v2", "RFCprobC3")
# Dom_Pop_Percent should already exist from prepare_epa_data
if (!"Dom_Pop_Percent" %in% names(merged_data) && "PWELL_private_pct" %in% names(merged_data)) {
  merged_data[, Dom_Pop_Percent := PWELL_private_pct]  # Already in percentage (0-100)
}

# Save merged USGS file
usgs_merged_path <- file.path(output_dir, "usgs_merged.csv")
data.table::fwrite(merged_data, usgs_merged_path)

# Health data needs id column
health_for_analysis <- data.table::fread(health_processed_path)
data.table::setnames(health_for_analysis, "FIPS", "id")
health_for_analysis[, id := as.integer(id)]

health_merged_path <- file.path(output_dir, "health_merged.txt")
data.table::fwrite(health_for_analysis, health_merged_path, sep = "\t")

# Run analysis
results <- tryCatch({
  perform_sensitivity_analysis(
    ndraws = 5,
    as_usgs_prob_csv = usgs_merged_path,
    as_epa_prob_csv = usgs_merged_path,  # EPA columns merged into USGS
    birth_data_txt = health_merged_path,
    regression_formula = "~ as.factor(AsLevel) + (1 | age_decade)",
    output_dir = output_dir,
    targets = c("SBP"),
    as_cat_label = c("As<5", "As5-10", "As10+"),
    drop_as_cat_label_reg = c("As<5"),
    id_col = "id",
    as_conc_cols = c("RFCprobC1", "RFCprobC2", "RFCprobC3"),
    pop_well_col = "Wells_2010",
    record_id_col = "id",
    as_level_col = "AsLevel",
    epa_as_mean_col = "EPA_AS_meanlog",
    epa_as_sd_col = "EPA_AS_sdlog",
    epa_pwell_col = "Dom_Pop_Percent",
    epa_cutoffs = c(5, 10),
    seed = 12345,
    mice_m = 1,
    mice_maxit = 5,
    mice_method = "pmm",
    apply_imputation_fallback = TRUE
  )
}, error = function(e) {
  cat("\nError during analysis:\n")
  cat(conditionMessage(e), "\n")
  return(NULL)
})

# ==============================================================================
# Step 10: Display Results
# ==============================================================================

if (!is.null(results)) {
  cat("\n========================================\n")
  cat("Analysis Complete (GeoTIFF Workflow)\n")
  cat("========================================\n\n")
  
  cat("Results for SBP:\n")
  print(results$SBP)
  
  cat("\nOutput files saved to:", output_dir, "\n")
  cat("Files:\n")
  for (f in list.files(output_dir)) {
    cat("  -", f, "\n")
  }
} else {
  cat("\nAnalysis did not complete successfully.\n")
}

cat("\n=== GeoTIFF Workflow Example Complete ===\n")
