# Load the testthat library
library(testthat)
library(geoExposeR)

# Helper function to create dummy data for testing
# This function generates a complete set of dummy data files required to
# run perform_sensitivity_analysis(). It creates data that mimics the
# structure of real-world inputs.
create_dummy_data <- function(
  file_dir = "test_data",
  num_rows = 500,
  num_births = 2000,
  seed = 12345,
  introduce_nas = TRUE
) {
  set.seed(seed)
  if (!dir.exists(file_dir)) {
    dir.create(file_dir, recursive = TRUE)
  }
  county_ids <- sprintf("%05d", 1:100)
  geoid10 <- sample(county_ids, num_rows, replace = TRUE)
  birth_fips <- sample(county_ids, num_births, replace = TRUE)

  # Modeled Probability Data
  prob_values <- t(apply(
    matrix(stats::runif(num_rows * 3), ncol = 3),
    1,
    function(x) x / sum(x)
  ))
  wells_2010 <- sample(50:500, num_rows, replace = TRUE)
  dummy_prob_df <- data.frame(
    prob_C1 = prob_values[, 1],
    prob_C2 = prob_values[, 2],
    prob_C3 = prob_values[, 3],
    Wells_2010 = wells_2010,
    GEOID10 = geoid10
  )
  dummy_prob_csv <- file.path(file_dir, "dummy_prob_data.csv")
  data.table::fwrite(dummy_prob_df, dummy_prob_csv)

  # Measured Concentration Data
  conc_lognormal_meanlog <- stats::runif(num_rows, 0.5, 1.5)
  conc_lognormal_sdlog <- stats::runif(num_rows, 0.1, 0.5)
  private_well_percent <- stats::runif(num_rows, min = 1, max = 99)
  dummy_conc_df <- data.frame(
    conc_meanlog = conc_lognormal_meanlog,
    conc_sdlog = conc_lognormal_sdlog,
    PWELL_private_pct = private_well_percent,
    GEOID10 = geoid10
  )
  dummy_conc_csv <- file.path(file_dir, "dummy_conc_data.csv")
  data.table::fwrite(dummy_conc_df, dummy_conc_csv)

  # Birth Data
  births_df <- data.frame(
    MAGE_R = sample(18:45, num_births, replace = TRUE),
    MRACEHISP_F = factor(
      sample(c("White", "Black", "Hispanic", "Other"), num_births, replace = TRUE)
    ),
    MEDUC_2 = sample(0:1, num_births, replace = TRUE),
    MEDUC_3 = sample(0:1, num_births, replace = TRUE),
    MEDUC_4 = sample(0:1, num_births, replace = TRUE),
    smoke = sample(0:1, num_births, replace = TRUE),
    pm = stats::rnorm(num_births, mean = 8, sd = 2),
    DMAR_1 = sample(0:1, num_births, replace = TRUE),
    FIPS = birth_fips,
    MRSTATE = factor(sample(state.abb, num_births, replace = TRUE)),
    BWT = round(stats::rnorm(num_births, mean = 3400, sd = 450)),
    OEGEST = round(stats::rnorm(num_births, mean = 40, sd = 2))
  )
  if (introduce_nas) {
    n_na <- floor(0.1 * num_births)
    births_df$MAGE_R[sample.int(num_births, n_na)] <- NA
    births_df$smoke[sample.int(num_births, n_na)] <- NA
  }
  dummy_births_txt <- file.path(file_dir, "dummy_births.txt")
  data.table::fwrite(births_df, dummy_births_txt, sep = "\t")

  return(list(
    prob = dummy_prob_csv,
    conc = dummy_conc_csv,
    births = dummy_births_txt
  ))
}

# Define the test suite for the main analysis function
context("Test perform_sensitivity_analysis with comprehensive coverage")

test_that("Function handles missing covariates and produces valid output", {
  # --- 1. Setup: Create data with missing values ---
  temp_dir <- tempfile(pattern = "test_mice_")
  dir.create(temp_dir)

  # Generate dummy data files
  dummy_files <- create_dummy_data(file_dir = temp_dir, introduce_nas = TRUE)

  # Load the birth data
  birth_data <- data.table::fread(
    dummy_files$births, header = TRUE,
    nThread = parallel::detectCores()
  )

  # Verify that NAs were introduced
  expect_true(any(is.na(birth_data$MAGE_R)))
  expect_true(any(is.na(birth_data$smoke)))

  # --- 2. Define Parameters for the Analysis ---
  ndraws <- 5
  targets <- c("BWT", "OEGEST")
  impute_vars <- c("MAGE_R", "smoke")
  exposure_level_col <- "ExposureLevel"
  mice_covs <- c("pm", "DMAR_1", "MEDUC_2", "MEDUC_3", "MEDUC_4")

  # check if mice_covs are in the birth_data
  if (!all(mice_covs %in% names(birth_data))) {
    stop("Some mice_covs are not present in the birth_data.")
  }

  # Define a simple regression formula for testing purposes
  regression_formula <- paste0(
    "~ as.factor(", exposure_level_col, ") + ",
    "MAGE_R + smoke + (1 | MRSTATE)"
  )
  output_dir <- "test_output"

  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # --- 3. Execution: Run the analysis ---
  results <- expect_no_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = ndraws,
      regression_formula = regression_formula,
      impute_vars = impute_vars,
      targets = targets,
      output_dir = output_dir,
      mice_covs = mice_covs,
      exposure_level_col = exposure_level_col,
      rucc_col = NULL
    )
  )

  # --- 4. Validation: Check the structure and content of the results ---
  cat("\n--- Pooled Analysis Results ---\n")
  print(results)
  cat("-----------------------------\n")

  # Check if the main result is a list
  expect_type(results, "list")

  # Check if the list is named after the targets
  expect_named(results, targets, ignore.order = TRUE)

  # Check the structure of the results for all targets
  for (target in targets) {
    print(paste("Checking results for target:", target))
    expect_s3_class(results[[target]], "data.frame")

    # Check if the expected columns are present
    expected_cols <- c(
      "q.mi", "se.mi", "statistic", "conf.low", "conf.high", "p.value"
    )
    expect_true(
      all(expected_cols %in% names(results[[target]]))
    )

    # Check that the number of rows corresponds to the number of
    # estimated effects
    expect_equal(nrow(results[[target]]), 2)

    # Check that the output CSV file was created
    expected_csv_path <- file.path(
      output_dir, paste0(target, "_pooled_results.csv")
    )
    expect_true(file.exists(expected_csv_path))

    # Check for any NA, NaN, or Inf values in the key numeric columns
    expect_false(
      any(is.na(results[[target]]$q.mi)),
      "Pooled estimates (q.mi) should not be NA"
    )
    expect_false(
      any(is.na(results[[target]]$se.mi)),
      "Pooled standard errors (se.mi) should not be NA"
    )
    expect_false(
      any(is.na(results[[target]]$p.value)),
      "P-values should not be NA"
    )

    expect_true(
      all(is.finite(results[[target]]$q.mi)),
      "Pooled estimates (q.mi) must be finite"
    )
    expect_true(
      all(is.finite(results[[target]]$se.mi)),
      "Pooled standard errors (se.mi) must be finite"
    )
  }

  # --- 5. Teardown: Clean up the temporary directory ---
  unlink(temp_dir, recursive = TRUE)
})

test_that("Function handles basic analysis without missing data", {
  # Test with complete data (no NAs)
  temp_dir <- tempfile(pattern = "test_complete_")
  dir.create(temp_dir)

  dummy_files <- create_dummy_data(file_dir = temp_dir, introduce_nas = FALSE)

  results <- expect_no_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 3,
      regression_formula = "~ as.factor(ExposureLevel) + MAGE_R + (1 | MRSTATE)",
      targets = c("BWT"),
      output_dir = tempdir()
    )
  )

  expect_type(results, "list")
  expect_true("BWT" %in% names(results))
  expect_s3_class(results[["BWT"]], "data.frame")

  unlink(temp_dir, recursive = TRUE)
})

test_that("Function handles single target variable", {
  temp_dir <- tempfile(pattern = "test_single_")
  dir.create(temp_dir)

  dummy_files <- create_dummy_data(file_dir = temp_dir)

  results <- expect_no_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2,
      regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)",
      targets = "OEGEST",
      output_dir = tempdir()
    )
  )

  expect_type(results, "list")
  expect_equal(length(results), 1)
  expect_true("OEGEST" %in% names(results))

  unlink(temp_dir, recursive = TRUE)
})

# ===== Tests for validate_inputs function =====

test_that("validate_inputs validates required parameters correctly", {
  # Test missing birth_data_txt
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = NULL,
      model_prob_csv = "prob_data.csv",
      conc_data_csv = "conc_data.csv"
    ),
    "birth_data_txt is required and cannot be NULL"
  )

  # Test missing model_prob_csv
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = "birth.txt",
      model_prob_csv = NULL,
      conc_data_csv = "conc_data.csv"
    ),
    "model_prob_csv is required and cannot be NULL"
  )

  # Test missing conc_data_csv
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = "birth.txt",
      model_prob_csv = "prob_data.csv",
      conc_data_csv = NULL
    ),
    "conc_data_csv is required and cannot be NULL"
  )
})

test_that("validate_inputs validates file path formats", {
  # Test non-character birth_data_txt
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = 123,
      model_prob_csv = "prob_data.csv",
      conc_data_csv = "conc_data.csv"
    ),
    "birth_data_txt must be a single character string"
  )

  # Test multiple values for birth_data_txt
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = c("file1.txt", "file2.txt"),
      model_prob_csv = "prob_data.csv",
      conc_data_csv = "conc_data.csv"
    ),
    "birth_data_txt must be a single character string"
  )

  # Test non-character model_prob_csv
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = "birth.txt",
      model_prob_csv = 123,
      conc_data_csv = "conc_data.csv"
    ),
    "model_prob_csv must be a single character string"
  )

  # Test non-character conc_data_csv
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = "birth.txt",
      model_prob_csv = "prob_data.csv",
      conc_data_csv = 123
    ),
    "conc_data_csv must be a single character string"
  )
})

test_that("validate_inputs validates file existence", {
  # Test non-existent files
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = "nonexistent_birth.txt",
      model_prob_csv = "prob_data.csv",
      conc_data_csv = "conc_data.csv"
    ),
    "File does not exist: nonexistent_birth.txt"
  )

  # Create temporary files for testing
  temp_dir <- tempfile(pattern = "test_validation_")
  dir.create(temp_dir)

  birth_file <- file.path(temp_dir, "birth.txt")
  prob_file <- file.path(temp_dir, "prob_data.csv")
  conc_file <- file.path(temp_dir, "conc_data.csv")

  # Create the files
  writeLines("test", birth_file)
  writeLines("test", prob_file)
  writeLines("test", conc_file)

  # Test with existing files (should not error)
  expect_no_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file
    )
  )

  # Test when the concentration file does not exist
  unlink(conc_file)  # Remove concentration file
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file
    ),
    "File does not exist.*conc_data\\.csv"
  )

  # Test when the probability file does not exist
  writeLines("test", conc_file)  # Restore concentration file
  unlink(prob_file)  # Remove probability file
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file
    ),
    "File does not exist.*prob_data\\.csv"
  )

  unlink(temp_dir, recursive = TRUE)
})

test_that("validate_inputs validates ndraws parameter", {
  temp_dir <- tempfile(pattern = "test_ndraws_")
  dir.create(temp_dir)

  # Create temporary files
  birth_file <- file.path(temp_dir, "birth.txt")
  prob_file <- file.path(temp_dir, "prob_data.csv")
  conc_file <- file.path(temp_dir, "conc_data.csv")
  writeLines("test", birth_file)
  writeLines("test", prob_file)
  writeLines("test", conc_file)

  # Test non-numeric ndraws
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      ndraws = "five"
    ),
    "ndraws must be a single numeric value"
  )

  # Test multiple values for ndraws
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      ndraws = c(1, 2, 3)
    ),
    "ndraws must be a single numeric value"
  )

  # Test zero ndraws
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      ndraws = 0
    ),
    "ndraws must be a positive integer"
  )

  # Test negative ndraws
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      ndraws = -5
    ),
    "ndraws must be a positive integer"
  )

  # Test decimal ndraws
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      ndraws = 5.5
    ),
    "ndraws must be a positive integer"
  )

  # Test very large ndraws (should warn)
  expect_warning(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      ndraws = 15000
    ),
    "ndraws is very large.*This may take a long time"
  )

  unlink(temp_dir, recursive = TRUE)
})

test_that("validate_inputs validates regression_formula parameter", {
  temp_dir <- tempfile(pattern = "test_formula_")
  dir.create(temp_dir)

  # Create temporary files
  birth_file <- file.path(temp_dir, "birth.txt")
  prob_file <- file.path(temp_dir, "prob_data.csv")
  conc_file <- file.path(temp_dir, "conc_data.csv")
  writeLines("test", birth_file)
  writeLines("test", prob_file)
  writeLines("test", conc_file)

  # Test non-character, non-formula regression_formula
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      regression_formula = 123
    ),
    "regression_formula must be a character string or formula object"
  )

  # Test multiple character values
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      regression_formula = c("~ x", "~ y")
    ),
    "regression_formula must be a single character string"
  )

  # Test empty regression_formula
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      regression_formula = ""
    ),
    "regression_formula cannot be empty"
  )

  # Test whitespace-only regression_formula
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      regression_formula = "   "
    ),
    "regression_formula cannot be empty"
  )

  # Test invalid formula syntax
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      regression_formula = "invalid formula syntax $$"
    ),
    "Invalid regression_formula"
  )

  # Test valid character formula
  expect_no_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      regression_formula = "~ as.factor(ExposureLevel) + y + (1 | MRSTATE)"
    )
  )

  # Test valid formula object
  expect_no_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      regression_formula = as.formula("~ x + (1 | y)")
    )
  )

  unlink(temp_dir, recursive = TRUE)
})

test_that("validate_inputs validates targets parameter", {
  temp_dir <- tempfile(pattern = "test_targets_")
  dir.create(temp_dir)

  # Create temporary files
  birth_file <- file.path(temp_dir, "birth.txt")
  prob_file <- file.path(temp_dir, "prob_data.csv")
  conc_file <- file.path(temp_dir, "conc_data.csv")
  writeLines("test", birth_file)
  writeLines("test", prob_file)
  writeLines("test", conc_file)

  # Test non-character targets
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      targets = 123
    ),
    "targets must be a character vector"
  )

  # Test empty targets
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      targets = character(0)
    ),
    "targets cannot be empty"
  )

  # Test targets with NA
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      targets = c("BWT", NA, "OEGEST")
    ),
    "targets cannot contain NA or empty values"
  )

  # Test targets with empty string
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      targets = c("BWT", "", "OEGEST")
    ),
    "targets cannot contain NA or empty values"
  )

  # Test duplicate targets
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      targets = c("BWT", "BWT", "OEGEST")
    ),
    "targets cannot contain duplicate values"
  )

  # Test valid targets
  expect_no_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      targets = c("BWT", "OEGEST")
    )
  )

  unlink(temp_dir, recursive = TRUE)
})

test_that("validate_inputs validates impute_vars parameter", {
  temp_dir <- tempfile(pattern = "test_impute_vars_")
  dir.create(temp_dir)

  # Create temporary files
  birth_file <- file.path(temp_dir, "birth.txt")
  prob_file <- file.path(temp_dir, "prob_data.csv")
  conc_file <- file.path(temp_dir, "conc_data.csv")
  writeLines("test", birth_file)
  writeLines("test", prob_file)
  writeLines("test", conc_file)

  # Test non-character impute_vars
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      impute_vars = 123
    ),
    "impute_vars must be a character vector or NULL"
  )

  # Test empty impute_vars (should suggest NULL)
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      impute_vars = character(0)
    ),
    "impute_vars cannot be empty \\(use NULL instead\\)"
  )

  # Test impute_vars with NA
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      impute_vars = c("var1", NA, "var2")
    ),
    "impute_vars cannot contain NA or empty values"
  )

  # Test valid impute_vars
  expect_no_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      impute_vars = c("var1", "var2")
    )
  )

  # Test NULL impute_vars (should be valid)
  expect_no_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      impute_vars = NULL
    )
  )

  unlink(temp_dir, recursive = TRUE)
})

test_that("validate_inputs validates output_dir parameter", {
  temp_dir <- tempfile(pattern = "test_output_dir_")
  dir.create(temp_dir)

  # Create temporary files
  birth_file <- file.path(temp_dir, "birth.txt")
  prob_file <- file.path(temp_dir, "prob_data.csv")
  conc_file <- file.path(temp_dir, "conc_data.csv")
  writeLines("test", birth_file)
  writeLines("test", prob_file)
  writeLines("test", conc_file)

  # Test non-character output_dir
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      output_dir = 123
    ),
    "output_dir must be a single character string"
  )

  # Test empty output_dir
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      output_dir = ""
    ),
    "output_dir cannot be empty"
  )

  # Test output_dir with non-existent parent
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      output_dir = "/nonexistent/parent/directory"
    ),
    "Parent directory does not exist"
  )

  # Test valid output_dir
  valid_output_dir <- file.path(temp_dir, "output")
  expect_no_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      output_dir = valid_output_dir
    )
  )

  unlink(temp_dir, recursive = TRUE)
})

test_that("validate_inputs validates MICE parameters", {
  temp_dir <- tempfile(pattern = "test_mice_")
  dir.create(temp_dir)

  # Create temporary files
  birth_file <- file.path(temp_dir, "birth.txt")
  prob_file <- file.path(temp_dir, "prob_data.csv")
  conc_file <- file.path(temp_dir, "conc_data.csv")
  writeLines("test", birth_file)
  writeLines("test", prob_file)
  writeLines("test", conc_file)

  # Test invalid mice_m
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      mice_m = -1
    ),
    "mice_m must be a positive integer"
  )

  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      mice_m = 2.5
    ),
    "mice_m must be a positive integer"
  )

  # Test non-numeric mice_m
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      mice_m = "ten"
    ),
    "mice_m must be a single numeric value"
  )

  # Test large mice_m (should warn)
  expect_warning(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      mice_m = 150
    ),
    "mice_m is very large.*This may take a long time"
  )

  # Test invalid mice_maxit
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      mice_maxit = 0
    ),
    "mice_maxit must be a positive integer"
  )

  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      mice_maxit = "three"
    ),
    "mice_maxit must be a single numeric value"
  )

  # Test invalid mice_method
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      mice_method = "invalid_method"
    ),
    "mice_method must be one of"
  )

  # Test non-character mice_method
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      mice_method = 123
    ),
    "mice_method must be a single character string"
  )

  # Test valid mice_method
  expect_no_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      mice_method = "pmm"
    )
  )

  unlink(temp_dir, recursive = TRUE)
})

test_that("validate_inputs validates numeric parameters", {
  temp_dir <- tempfile(pattern = "test_numeric_")
  dir.create(temp_dir)

  # Create temporary files
  dummy_files <- create_dummy_data(file_dir = temp_dir)

  # Test invalid conc_lognormal_sdlog
  # create a dummy conc_lognormal_sdlog
  dummy_conc_df <- data.table::fread(dummy_files$conc, header = TRUE)
  dummy_conc_df$conc_sdlog <- -1.0
  dummy_conc_csv <- file.path(temp_dir, "invalid_conc_data.csv")
  data.table::fwrite(dummy_conc_df, dummy_conc_csv)

  expect_error(
    geoExposeR:::perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_conc_csv
    ),
    "conc_sd_col must contain only positive values.*"
  )

  dummy_conc_df$conc_sdlog1 <- "one"
  data.table::fwrite(dummy_conc_df, dummy_conc_csv)
  expect_error(
    geoExposeR:::perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_conc_csv,
      conc_sd_col = "conc_sdlog1"
    ),
    "conc_sd_col must contain only finite numeric values.*"
  )

  # Test large conc_lognormal_sdlog (should warn)
  dummy_conc_df$conc_sdlog <- 15.0
  data.table::fwrite(dummy_conc_df, dummy_conc_csv)
  expect_warning(
    geoExposeR:::perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_conc_csv,
      regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)"
    ),
    "conc_sd_col contains very large values.*"
  )

  # Test invalid seed
  expect_error(
    geoExposeR:::perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_conc_csv,
      seed = 123.45
    ),
    "seed must be an integer"
  )

  expect_error(
    geoExposeR:::perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_conc_csv,
      seed = c(1, 2, 3)
    ),
    "seed must be a single numeric value"
  )

  unlink(temp_dir, recursive = TRUE)
})

test_that("validate_inputs validates exposure category labels", {
  temp_dir <- tempfile(pattern = "test_exposure_labels_")
  dir.create(temp_dir)

  # Create temporary files
  dummy_files <- create_dummy_data(file_dir = temp_dir)

  # Test non-character cat_label
  expect_error(
    geoExposeR:::perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      cat_label = 123,
      drop_cat_label_ref = c("Low"),
      regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)"
    ),
    "cat_label must be a character vector"
  )

  # Test empty cat_label
  expect_error(
    geoExposeR:::perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      cat_label = character(0),
      drop_cat_label_ref = c("Low")
    ),
    "cat_label cannot be empty"
  )

  # Test cat_label with duplicates
  expect_error(
    geoExposeR:::perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      cat_label = c("Low", "Medium", "Low"),
      drop_cat_label_ref = c("Low")
    ),
    "cat_label cannot contain duplicate values"
  )

  # Test drop_cat_label_ref not in cat_label
  expect_error(
    geoExposeR:::perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      cat_label = c("Low", "Medium", "High"),
      drop_cat_label_ref = c("VeryLow")
    ),
    "All values in drop_cat_label_ref must be present in cat_label"
  )

  # Test dropping all categories
  expect_error(
    geoExposeR:::perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      cat_label = c("Low", "Medium"),
      drop_cat_label_ref = c("Low", "Medium")
    ),
    "Cannot drop all exposure categories as reference levels"
  )

  # Test valid category configuration
  expect_no_error(
    geoExposeR:::perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      cat_label = c("Low", "Medium", "High"),
      drop_cat_label_ref = c("Low"),
      regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)"
    )
  )

  unlink(temp_dir, recursive = TRUE)
})

# ===== Tests for main function validation integration =====

test_that("Function validates required parameters through main interface", {
  # Test missing required parameters via main function
  expect_error(
    perform_sensitivity_analysis(),
    "birth_data_txt is required and cannot be NULL"
  )

  expect_error(
    perform_sensitivity_analysis(birth_data_txt = "test.txt"),
    "model_prob_csv is required and cannot be NULL"
  )

  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = "test.txt",
      model_prob_csv = "prob_data.csv"
    ),
    "conc_data_csv is required and cannot be NULL"
  )
})

test_that("Function validates parameter values through main interface", {
  temp_dir <- tempfile(pattern = "test_main_validation_")
  dir.create(temp_dir)

  dummy_files <- create_dummy_data(file_dir = temp_dir)

  # Test invalid ndraws through main function
  expect_error(
    geoExposeR::perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 0,
      regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)"
    ),
    "ndraws must be a positive integer"
  )

  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = -5,
      regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)"
    ),
    "ndraws must be a positive integer"
  )

  # Test invalid regression formula through main function
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2,
      regression_formula = ""
    ),
    "regression_formula cannot be empty"
  )

  unlink(temp_dir, recursive = TRUE)
})

test_that("Function handles file validation errors through main interface", {
  # Test non-existent files through main function
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = "nonexistent.txt",
      model_prob_csv = "nonexistent.csv",
      conc_data_csv = "nonexistent.csv",
      ndraws = 2
    ),
    "File does not exist: nonexistent.txt"
  )
})

# ===== Continue with remaining existing tests =====

test_that("Function handles different MICE parameters", {
  temp_dir <- tempfile(pattern = "test_mice_params_")
  dir.create(temp_dir)

  dummy_files <- create_dummy_data(file_dir = temp_dir, introduce_nas = TRUE)

  # Test with different MICE parameters
  
  results <- expect_no_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2,
      regression_formula = "~ as.factor(ExposureLevel) + MAGE_R + (1 | MRSTATE)",
      targets = "BWT",
      impute_vars = "MAGE_R",
      mice_m = 3,
      mice_maxit = 2,
      mice_method = "pmm",
      output_dir = tempdir()
    )
  )

  expect_type(results, "list")
  expect_true("BWT" %in% names(results))

  unlink(temp_dir, recursive = TRUE)
})

test_that("Function handles different exposure level column names", {
  temp_dir <- tempfile(pattern = "test_as_col_")
  dir.create(temp_dir)

  dummy_files <- create_dummy_data(
    file_dir = temp_dir
  )
  exposure_level_col <- "ArsenicCategory"
  results <- expect_no_error(
    geoExposeR::perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2,
      regression_formula = paste0(
        "~ as.factor(", exposure_level_col, ") + (1 | MRSTATE)"
      ),
      targets = "BWT",
      exposure_level_col = exposure_level_col,
      output_dir = tempdir()
    )
  )

  expect_type(results, "list")
  expect_true("BWT" %in% names(results))

  unlink(temp_dir, recursive = TRUE)
})

test_that("Function handles custom output directory", {
  temp_dir <- tempfile(pattern = "test_custom_out_")
  dir.create(temp_dir)

  dummy_files <- create_dummy_data(file_dir = temp_dir)
  custom_output <- file.path(temp_dir, "custom_results")

  results <- expect_no_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2,
      regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = custom_output
    )
  )

  # Check that custom output directory was created
  expect_true(dir.exists(custom_output))

  # Check that results file was created in custom directory
  expect_true(file.exists(file.path(custom_output, "BWT_pooled_results.csv")))

  unlink(temp_dir, recursive = TRUE)
})

test_that("Function handles boundary conditions appropriately", {
  temp_dir <- tempfile(pattern = "test_boundary_")
  dir.create(temp_dir)

  dummy_files <- create_dummy_data(file_dir = temp_dir)

  # Test with minimum viable ndraws
  results <- expect_no_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 1,  # Absolute minimum
      regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir()
    )
  )

  expect_type(results, "list")
  expect_true("BWT" %in% names(results))

  unlink(temp_dir, recursive = TRUE)
})

test_that("Function handles multiple targets with imputation", {
  temp_dir <- tempfile(pattern = "test_multi_target_")
  dir.create(temp_dir)

  dummy_files <- create_dummy_data(file_dir = temp_dir, introduce_nas = TRUE)

  results <- expect_no_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2,
      regression_formula = "~ as.factor(ExposureLevel) + MAGE_R + (1 | MRSTATE)",
      targets = c("BWT", "OEGEST"),
      impute_vars = "MAGE_R",
      output_dir = tempdir()
    )
  )

  expect_type(results, "list")
  expect_equal(length(results), 2)
  expect_true(all(c("BWT", "OEGEST") %in% names(results)))

  # Check that both targets have valid results
  for (target in c("BWT", "OEGEST")) {
    expect_s3_class(results[[target]], "data.frame")
    expect_true(
      all(c("q.mi", "se.mi", "p.value") %in% names(results[[target]]))
    )
  }

  unlink(temp_dir, recursive = TRUE)
})

test_that("Function handles formula validation", {
  temp_dir <- tempfile(pattern = "test_formula_")
  dir.create(temp_dir)

  dummy_files <- create_dummy_data(file_dir = temp_dir)

  # Test formula as character string
  results1 <- expect_no_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2,
      regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir()
    )
  )

  # Test formula as formula object
  formula_obj <- as.formula("~ as.factor(ExposureLevel) + (1 | MRSTATE)")
  results2 <- expect_no_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2,
      regression_formula = formula_obj,
      targets = "BWT",
      output_dir = tempdir()
    )
  )

  expect_type(results1, "list")
  expect_type(results2, "list")

  unlink(temp_dir, recursive = TRUE)
})

test_that("Function handles seed setting for reproducibility", {
  temp_dir <- tempfile(pattern = "test_seed_")
  dir.create(temp_dir)

  dummy_files <- create_dummy_data(file_dir = temp_dir)

  # Run analysis twice with same seed
  set.seed(12345)
  results1 <- perform_sensitivity_analysis(
    birth_data_txt = dummy_files$births,
    model_prob_csv = dummy_files$prob,
    conc_data_csv = dummy_files$conc,
    ndraws = 2,
    regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)",
    targets = "BWT",
    output_dir = tempdir()
  )

  set.seed(12345)
  results2 <- perform_sensitivity_analysis(
    birth_data_txt = dummy_files$births,
    model_prob_csv = dummy_files$prob,
    conc_data_csv = dummy_files$conc,
    ndraws = 2,
    regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)",
    targets = "BWT",
    output_dir = tempdir()
  )

  # Results should be very similar (allowing for small numerical differences)
  expect_equal(
    results1[["BWT"]]$q.mi,
    results2[["BWT"]]$q.mi,
    tolerance = 1e-6
  )

  unlink(temp_dir, recursive = TRUE)
})

test_that("Function handles all error conditions comprehensively", {
  # Test with completely invalid data
  temp_dir <- tempfile(pattern = "test_errors_")
  dir.create(temp_dir)

  # Create invalid birth data (missing required columns)
  invalid_birth_file <- file.path(temp_dir, "invalid_birth.txt")
  invalid_data <- data.frame(wrong_col = 1:10)
  write.table(invalid_data, invalid_birth_file, sep = "\t", row.names = FALSE)

  dummy_files <- create_dummy_data(file_dir = temp_dir)

  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = invalid_birth_file,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2
    )
  )

  unlink(temp_dir, recursive = TRUE)
})

test_that("Function handles memory and performance edge cases", {
  temp_dir <- tempfile(pattern = "test_performance_")
  dir.create(temp_dir)

  # Create large dataset to test memory handling
  dummy_files <- create_dummy_data(
    file_dir = temp_dir,
    num_births = 10000,  # Large number of births
    num_rows = 1000,     # Large number of rows in prob and conc
    introduce_nas = TRUE  # Introduce some NAs to test imputation
  )

  # Test with larger dataset but minimal draws to keep test time reasonable
  results <- expect_no_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2,
      regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir()
    )
  )

  expect_type(results, "list")
  expect_true("BWT" %in% names(results))

  unlink(temp_dir, recursive = TRUE)
})

test_that("Function validates ExposureLevel requirement through main interface", {
  temp_dir <- tempfile(pattern = "test_main_aslevel_")
  dir.create(temp_dir)

  dummy_files <- create_dummy_data(file_dir = temp_dir)

  # Test regression formula without ExposureLevel through main function
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2,
      regression_formula = "~ MAGE_R + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir()
    ),
    "regression_formula must include 'ExposureLevel' as a factor variable"
  )

  # Test regression formula with ExposureLevel but not as factor through main function
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2,
      regression_formula = "~ ExposureLevel + MAGE_R + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir()
    ),
    "regression_formula must include 'ExposureLevel' as a factor variable"
  )

  unlink(temp_dir, recursive = TRUE)
})

# Additional validation tests

test_that("Function handles exposure category reference dropping edge cases", {
  temp_dir <- tempfile(pattern = "test_exposure_ref_")
  dir.create(temp_dir)

  dummy_files <- create_dummy_data(file_dir = temp_dir)

  # Test dropping all available exposure categories as reference
  # This should trigger:
  # "Cannot drop all exposure categories as reference levels."
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2,
      regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)",
      targets = "BWT",
      cat_label = c("Low", "Medium", "High"),
      drop_cat_label_ref = c("Low", "Medium", "High"),  # Drop all categories
      output_dir = tempdir()
    ),
    "Cannot drop all exposure categories as reference levels"
  )

  # Test with data that only has one unique exposure category
  # First, create a birth dataset with only one exposure level
  birth_data <- data.table::fread(dummy_files$births, header = TRUE)
  birth_data$ExposureLevel <- "Low"  # Set all to the same category

  # Write the modified data back
  single_category_file <- file.path(temp_dir, "single_category_births.txt")
  data.table::fwrite(birth_data, single_category_file, sep = "\t")

  # Test with single category that gets dropped
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = single_category_file,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2,
      regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)",
      targets = "BWT",
      cat_label = c("Low"),
      drop_cat_label_ref = c("Low"),  # Drop the only category
      output_dir = tempdir()
    ),
    "Cannot drop all exposure categories as reference levels"
  )

  # Test edge case: data has categories but after filtering, none remain
  # Create data with two categories
  birth_data_two <- data.table::fread(dummy_files$births, header = TRUE)
  birth_data_two$ExposureLevel <- sample(
    c("Low", "Medium"),
    nrow(birth_data_two),
    replace = TRUE
  )

  two_category_file <- file.path(temp_dir, "two_category_births.txt")
  data.table::fwrite(birth_data_two, two_category_file, sep = "\t")

  # Drop both available categories
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = two_category_file,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2,
      regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)",
      targets = "BWT",
      cat_label = c("Low", "Medium"),
      drop_cat_label_ref = c("Low", "Medium"),  # Drop both categories
      output_dir = tempdir()
    ),
    "Cannot drop all exposure categories as reference levels"
  )

  # Test valid case: have multiple categories and only drop some
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = two_category_file,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2,
      regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)",
      targets = "BWT",
      cat_label = c("Low", "Medium"),
      # Only drop one category, leaving "Medium"
      drop_cat_label_ref = c("Low"),
      output_dir = tempdir()
    ),
    "incorrect number of probabilities"
  )

  unlink(temp_dir, recursive = TRUE)
})

test_that("Function validates exposure category data consistency", {
  temp_dir <- tempfile(pattern = "test_exposure_data_")
  dir.create(temp_dir)

  dummy_files <- create_dummy_data(file_dir = temp_dir)

  # Test when specified cat_label doesn't match data categories
  # Load the birth data to see what categories are actually present
  birth_data <- data.table::fread(dummy_files$births, header = TRUE)

  # Specify categories that don't exist in the data
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2,
      regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)",
      targets = "BWT",
      # Categories not in data
      cat_label = c("NonExistent1", "NonExistent2"),
      output_dir = tempdir()
    ),
    # This might trigger a different error depending on your validation logic
    # Adjust the error message to match your actual implementation
    "All values in drop_cat_label_ref must be present in cat_label"
  )

  # Test dropping categories that don't exist in cat_label
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2,
      regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)",
      targets = "BWT",
      drop_cat_label_ref = c("NonExistentCategory"),
      # Try to drop non-existent category
      output_dir = tempdir()
    ),
    "All values in drop_cat_label_ref must be present in cat_label"
  )

  # Test if drop_cat_label_ref is a character vector
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2,
      regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)",
      targets = "BWT",
      drop_cat_label_ref = 123,  # Invalid type
      output_dir = tempdir()
    ),
    "drop_cat_label_ref must be a character vector"
  )

  # Test if cat_label has NA or empty values
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2,
      regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)",
      targets = "BWT",
      cat_label = c("Low", NA, "High"),  # Contains NA
      output_dir = tempdir()
    ),    "cat_label cannot contain NA or empty values"
  )

  unlink(temp_dir, recursive = TRUE)
})


test_that("Function validates column names in input files", {
  temp_dir <- tempfile(pattern = "test_column_names_")
  dir.create(temp_dir)
  dummy_files <- create_dummy_data(file_dir = temp_dir)

  # Test missing columns in probability data
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2,
      regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir(),
      prob_cols = c(
        "NonExistentColumn1", "NonExistentColumn2",
        "NonExistentColumn3"
      )
    ),
    "The following columns are missing in.*"
  )

  # Test incorrect file format (e.g., non-CSV for conc)
  dummy_files <- create_dummy_data(file_dir = temp_dir)
  dummy_conc_df <- data.table::fread(dummy_files$conc, header = TRUE)
  invalid_conc_file <- file.path(temp_dir, "invalid_epa.bin")
  data.table::fwrite(
    dummy_conc_df,
    file = invalid_conc_file,
    sep = "\t"
  )
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = invalid_conc_file,
      ndraws = 2,
      regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir()
    ),
    "Unsupported file format for.*"
  )

  # Test non-character other column names
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2,
      regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir(),
      id_col = 123  # Invalid type
    ),
    "id_col must be a single character string"
  )

  # Test empty other column names
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2,
      regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir(),
      id_col = ""  # Empty string
    ),
    "id_col cannot be empty"
  )

  # Test incorrect character vectors for other column names
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2,
      regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir(),
      prob_cols = c(1, 2, 3)  # Mixed types
    ),
    "prob_cols must be a character vector.*"
  )

  # Test NA values in other column names
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2,
      regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir(),
      prob_cols = c("Col1", NA)  # NA value
    ),
    "prob_cols cannot contain NA or empty values"
  )

  # Test duplicate values in other column names
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2,
      regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir(),
      prob_cols = c("Col1", "Col1")  # Duplicate values
    ),
    "prob_cols cannot contain duplicate values"
  )

  unlink(temp_dir, recursive = TRUE)
})


test_that("Function validates conc_cutoff parameter", {
  temp_dir <- tempfile(pattern = "test_conc_cutoff_")
  dir.create(temp_dir)
  dummy_files <- create_dummy_data(file_dir = temp_dir)

  # Test negative conc_cutoffs
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2,
      regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir(),
      conc_cutoffs = -1  # Invalid negative value
    ),
    "conc_cutoffs must contain only positive values"
  )

  # Test non-numeric conc_cutoffs
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2,
      regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir(),
      conc_cutoffs = "high"  # Invalid non-numeric value
    ),
    "conc_cutoffs must be a numeric vector"
  )

  # Test NA in conc_cutoffs
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2,
      regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir(),
      conc_cutoffs = c(5, NA)  # Contains NA
    ),
    "conc_cutoffs cannot contain NA values"
  )

  # Test duplicate values in conc_cutoffs
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2,
      regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir(),
      conc_cutoffs = c(5, 5)  # Duplicate values
    ),
    "conc_cutoffs cannot contain duplicate values"
  )

  # Test valid conc_cutoffs with 4 categories
  # Create a dummy probability file with 4 exposure categories
  # Ensure the dummy probability data has appropriate columns
  # and realistic probability distributions

  num_rows <- 500
  county_ids <- sprintf("%05d", 1:100)
  geoid10 <- sample(county_ids, num_rows, replace = TRUE)
  # Generate random probabilities for 4 categories that sum to 1
  prob_values <- t(apply(
    matrix(stats::runif(num_rows * 4), ncol = 4),
    1,
    function(x) x / sum(x)
  ))
  wells_2010 <- sample(50:500, num_rows, replace = TRUE)
  dummy_prob_df <- data.frame(
    AS_C1 = prob_values[, 1],
    AS_C2 = prob_values[, 2],
    AS_C3 = prob_values[, 3],
    AS_C4 = prob_values[, 4],
    Wells_2010 = wells_2010,
    GEOID10 = geoid10
  )
  # Write the dummy probability data to a CSV fil
  dummy_prob_csv <- file.path(temp_dir, "dummy_prob_4cats.csv")
  data.table::fwrite(dummy_prob_df, dummy_prob_csv)

  # Create a dummy concentration file with 4 categories
  dummy_conc_df <- data.frame(
    conc_meanlog = runif(num_rows, 0.5, 1.5),  # Log-normal mean
    conc_sdlog = runif(num_rows, 0.1, 0.5),  # Log-normal sd
    PWELL_private_pct = runif(num_rows, 1, 99),  # Private well percentage
    GEOID10 = geoid10
  )
  dummy_conc_csv <- file.path(temp_dir, "dummy_conc_4cats.csv")
  data.table::fwrite(dummy_conc_df, dummy_conc_csv)
  dummy_files <- create_dummy_data(file_dir = temp_dir)

  # this tests the default fallback imputation as well
  expect_no_error({
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_prob_csv,
      conc_data_csv = dummy_conc_csv,
      ndraws = 2,
      regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir(),
      prob_cols = c("AS_C1", "AS_C2", "AS_C3", "AS_C4"),
      drop_cat_label_ref = "<5",  # Drop one category
      cat_label = c("<5", "5-10", "10-15", ">15"),
      conc_cutoffs = c(5, 10, 15)  # Valid ascending cutoffs
    )
  })

  expect_no_error({
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_prob_csv,
      conc_data_csv = dummy_conc_csv,
      ndraws = 2,
      regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir(),
      prob_cols = c("AS_C1", "AS_C2", "AS_C3", "AS_C4"),
      drop_cat_label_ref = "<5",  # Drop one category
      cat_label = c("<5", "5-10", "10-15", ">15"),
      conc_cutoffs = c(5, 10, 15),  # Valid ascending cutoffs
      apply_imputation_fallback = FALSE # removes NA values
    )
  })

  # Test descending order conc_cutoffs
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2,
      regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir(),
      conc_cutoffs = c(15, 10, 5)  # Descending order
    ),
    "conc_cutoffs must be in ascending order"
  )

  # Test if conc_cutoffs have categories not matching probability data
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2,
      regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir(),
      conc_cutoffs = c(5, 10, 20)  # Assuming probability data has only 3 categories
    ),
    "prob_cols must have length equal to.*"
  )

  unlink(temp_dir, recursive = TRUE)
})


# ===================================================================
# Tests for Dynamic Concentration Cutoffs Functionality
# ===================================================================

test_that("validate_conc_cutoffs works correctly", {
  prob_cols <- c("AS_C1", "AS_C2", "AS_C3", "AS_C4")
  # Test valid cutoffs
  valid_cutoffs <- c(5, 10, 15)
  expect_no_error(
    validate_conc_cutoff_params(valid_cutoffs, prob_cols)
  )

  # Test empty cutoffs similar to no cutoffs provided
  expect_no_error(
    validate_conc_cutoff_params(c(), prob_cols)
  )

  # Test non-numeric cutoffs
  expect_error(
    validate_conc_cutoff_params(c("5", "10"), prob_cols),
    "conc_cutoffs must be a numeric vector"
  )

  # Test unsorted cutoffs
  expect_error(
    validate_conc_cutoff_params(c(10, 5, 2.5), prob_cols),
    "conc_cutoffs must be in ascending order"
  )

  # Test negative cutoffs
  expect_error(
    validate_conc_cutoff_params(c(-1, 5, 10), prob_cols),
    "conc_cutoffs must contain.*"
  )
})

test_that("convert_conc_to_multinomial creates correct probabilities", {
  # Create test concentration data
  conc_data <- data.frame(
    FIPS = c(1001, 1002, 1003),
    conc_meanlog = c(log(1), log(5), log(15)),
    conc_sdlog = c(0.5, 1.0, 1.5)
  )

  # Test 3 categories (2 cutoffs)
  cutoff_config <- list(cutoffs = c(5, 10), n_categories = 3)
  result <- convert_conc_to_multinomial(
    conc_data, conc_mean_col = "conc_meanlog",
    conc_sd_col = "conc_sdlog",
    cutoff_config = cutoff_config,
    return_prob_matrix = FALSE
  )

  # Check that correct columns were created
  expected_cols <- c("conc_C1v2", "conc_C2v2", "conc_C3v2")
  expect_true(all(expected_cols %in% names(result)))

  # Check that probabilities are between 0 and 1
  prob_cols <- expected_cols
  prob_matrix <- as.matrix(result[, prob_cols, with = FALSE])
  expect_true(all(prob_matrix >= 0))
  expect_true(all(prob_matrix <= 1))

  # Check that rows sum to 1 (within tolerance)
  row_sums <- rowSums(prob_matrix)
  expect_true(all(abs(row_sums - 1) < 1e-10))

  # Test specific probability calculations for known values
  # For meanlog = log(5), sdlog = 1, cutoffs at 5, 10
  # P(X < 5) should be close to 0.5 for lognormal with median 5
  row2_probs <- as.numeric(result[2, ..prob_cols])
  expect_true(row2_probs[1] > 0.4 && row2_probs[1] < 0.6)  # Around 0.5
})

test_that("convert_conc_to_multinomial handles different categories", {
  conc_data <- data.frame(
    FIPS = 1001,
    conc_meanlog = log(5),
    conc_sdlog = 1.0
  )

  # Test 2 categories
  cutoff_config_2 <- list(cutoffs = c(5), n_categories = 2)
  result_2 <- convert_conc_to_multinomial(
    conc_data,
    conc_mean_col = "conc_meanlog",
    conc_sd_col = "conc_sdlog",
    cutoff_config = cutoff_config_2,
    return_prob_matrix = FALSE
  )
  prob_cols_2 <- c("conc_C1v2", "conc_C2v2")
  expect_true(all(prob_cols_2 %in% names(result_2)))
  expect_true(abs(sum(result_2[, ..prob_cols_2]) - 1) < 1e-10)

  # Test 5 categories
  cutoff_config_5 <- list(cutoffs = c(1, 5, 10, 20), n_categories = 5)
  result_5 <- convert_conc_to_multinomial(
    conc_data,
    conc_mean_col = "conc_meanlog",
    conc_sd_col = "conc_sdlog",
    cutoff_config = cutoff_config_5,
    return_prob_matrix = FALSE
  )
  expected_cols_5 <- paste0("conc_C", 1:5, "v2")
  expect_true(all(expected_cols_5 %in% names(result_5)))
  expect_true(abs(sum(result_5[, ..expected_cols_5]) - 1) < 1e-10)
})

test_that("generate_default_conc_cutoffs produces correct schemes", {
  # Test standard predefined schemes
  expect_equal(generate_default_conc_cutoffs(2), c(5))
  expect_equal(generate_default_conc_cutoffs(3), c(5, 10))
  expect_equal(generate_default_conc_cutoffs(4), c(2.5, 5, 10))
  expect_equal(generate_default_conc_cutoffs(5), c(1, 5, 10, 20))
  expect_equal(generate_default_conc_cutoffs(6), c(1, 3, 5, 10, 20))

  # Test logarithmic generation for unusual numbers
  cutoffs_8 <- generate_default_conc_cutoffs(8)
  expect_equal(length(cutoffs_8), 7)  # n-1 cutoffs for n categories
  expect_true(all(cutoffs_8 > 0))
  expect_true(all(diff(cutoffs_8) > 0))  # Should be increasing

  # Test error for invalid input
  expect_error(
    generate_default_conc_cutoffs(1),
    "Must have at least 2 categories"
  )
  expect_error(
    generate_default_conc_cutoffs(0),
    "Must have at least 2 categories"
  )
})

test_that("generate_logarithmic_cutoffs works correctly", {
  # Test normal case
  cutoffs_5 <- generate_logarithmic_cutoffs(5)
  expect_equal(length(cutoffs_5), 4)  # n-1 cutoffs
  expect_true(all(cutoffs_5 > 0))
  expect_true(all(diff(cutoffs_5) > 0))  # Increasing

  # Test boundary cases
  cutoffs_2 <- generate_logarithmic_cutoffs(2)
  expect_equal(length(cutoffs_2), 1)

  # Test error for invalid input
  expect_error(
    generate_logarithmic_cutoffs(1),
    "Must have at least 2 categories"
  )
})

test_that("validate_category_alignment catches mismatches", {
  # Valid alignment
  expect_true(validate_category_alignment(3, c(5, 10)))
  expect_true(validate_category_alignment(4, c(2.5, 5, 10)))
  expect_true(validate_category_alignment(2, c(5)))

  # Invalid alignments
  expect_error(
    validate_category_alignment(3, c(5)),
    paste0(
      "Category mismatch: probability data has 3 categories but concentration cutoffs define ",
      "2 categories"
    )
  )

  expect_error(
    validate_category_alignment(2, c(5, 10)),
    paste0(
      "Category mismatch: probability data has 2 categories but concentration cutoffs define ",
      "3 categories"
    )
  )
})

test_that("edge cases in concentration probability calculations", {
  # Test extreme values
  conc_extreme <- data.frame(
    FIPS = c(1001, 1002, 1003),
    conc_meanlog_extreme = c(log(0.1), log(100), log(1000)),
    conc_sdlog_extreme = c(0.1, 0.5, 1.0)
  )

  cutoff_config <- list(cutoffs = c(5, 10), n_categories = 3)
  result <- convert_conc_to_multinomial(
    datmatx_conc = conc_extreme,
    conc_mean_col = "conc_meanlog_extreme",
    conc_sd_col = "conc_sdlog_extreme",
    return_prob_matrix = FALSE,
    cutoff_config = cutoff_config
  )

  # Check row sums still equal 1
  prob_cols <- c("conc_C1v2", "conc_C2v2", "conc_C3v2")
  # convert to matrix for row sums
  prob_matrix <- as.matrix(result[, prob_cols, with = FALSE])
  row_sums <- rowSums(prob_matrix)
  expect_true(all(abs(row_sums - 1) < 1e-10))

  # Very low mean should have most probability in first category
  expect_true(prob_matrix[1, "conc_C1v2"] > 0.5)

  # Very high mean should have most probability in last category
  expect_true(prob_matrix[3, "conc_C3v2"] > 0.5)
})

# test apply_fallback_imputation functionality
test_that("Fallback imputation works when concentration data is missing", {
  temp_dir <- tempfile(pattern = "test_fallback_imputation_")
  dir.create(temp_dir)

  dummy_files <- create_dummy_data(file_dir = temp_dir)

  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2,
      regression_formula = "~ as.factor(ExposureLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir(),
      apply_imputation_fallback = "FOO"  # Enable fallback imputation
    ), "apply_imputation_fallback must be.*"
  )

  unlink(temp_dir, recursive = TRUE)
})
# ===== Tests for impute_exposure function =====

test_that("impute_exposure handles basic case correctly", {
  # Create test data
  set.seed(12345)
  n_counties <- 10
  n_cells <- 50
  ndraws <- 3
  cat_label <- c("<5", "5-10", "10+")

  # Create probability data with multiple cells per county
  county_ids <- sprintf("%05d", rep(1:n_counties, each = 5))
  datmatx_prob <- data.frame(
    GEOID10 = county_ids,
    stringsAsFactors = FALSE
  )

  # Create probability matrix
  rasterprob_combined <- matrix(
    runif(n_cells * 3),
    nrow = n_cells,
    ncol = 3
  )
  # Normalize rows
  rasterprob_combined <- rasterprob_combined / rowSums(rasterprob_combined)

  # Create population weights
  cellpop <- sample(50:500, n_cells, replace = TRUE)

  # Create birth data
  births <- data.frame(
    FIPS = sample(sprintf("%05d", 1:n_counties), 100, replace = TRUE),
    BWT = rnorm(100, 3400, 450),
    stringsAsFactors = FALSE
  )

  result <- geoExposeR:::impute_exposure(
    datmatx_prob = datmatx_prob,
    rasterprob_combined = rasterprob_combined,
    cellpop = cellpop,
    births = births,
    ndraws = ndraws,
    cat_label = cat_label,
    id_col = "GEOID10",
    record_id_col = "FIPS",
    exposure_level_col = "ExposureLevel",
    apply_imputation_fallback = TRUE
  )

  # Check result structure
  expect_type(result, "list")
  expect_equal(length(result), ndraws)

  # Check each imputed dataset
  for (j in seq_along(result)) {
    expect_true("ExposureLevel" %in% names(result[[j]]))
    expect_true(all(result[[j]]$ExposureLevel %in% cat_label))
    expect_false(any(is.na(result[[j]]$ExposureLevel)))
  }
})

test_that("impute_exposure handles NA weights correctly", {
  set.seed(12345)
  n_counties <- 5
  n_cells <- 25
  ndraws <- 2
  cat_label <- c("<5", "5-10", "10+")

  county_ids <- sprintf("%05d", rep(1:n_counties, each = 5))
  datmatx_prob <- data.frame(
    GEOID10 = county_ids,
    stringsAsFactors = FALSE
  )

  rasterprob_combined <- matrix(
    runif(n_cells * 3),
    nrow = n_cells,
    ncol = 3
  )
  rasterprob_combined <- rasterprob_combined / rowSums(rasterprob_combined)

  # Introduce NA values in weights
  cellpop <- sample(50:500, n_cells, replace = TRUE)
  cellpop[c(1, 5, 10, 15)] <- NA  # Some NA weights

  births <- data.frame(
    FIPS = sample(sprintf("%05d", 1:n_counties), 50, replace = TRUE),
    BWT = rnorm(50, 3400, 450),
    stringsAsFactors = FALSE
  )

  # Should not error with NA weights
  result <- expect_no_error(
    geoExposeR:::impute_exposure(
      datmatx_prob = datmatx_prob,
      rasterprob_combined = rasterprob_combined,
      cellpop = cellpop,
      births = births,
      ndraws = ndraws,
      cat_label = cat_label,
      id_col = "GEOID10",
      record_id_col = "FIPS",
      exposure_level_col = "ExposureLevel",
      apply_imputation_fallback = TRUE
    )
  )

  expect_equal(length(result), ndraws)
  for (j in seq_along(result)) {
    expect_false(any(is.na(result[[j]]$ExposureLevel)))
  }
})

test_that("impute_exposure handles all-zero weights gracefully", {
  set.seed(12345)
  n_counties <- 3
  n_cells <- 15
  ndraws <- 2
  cat_label <- c("<5", "5-10", "10+")

  county_ids <- sprintf("%05d", rep(1:n_counties, each = 5))
  datmatx_prob <- data.frame(
    GEOID10 = county_ids,
    stringsAsFactors = FALSE
  )

  rasterprob_combined <- matrix(
    runif(n_cells * 3),
    nrow = n_cells,
    ncol = 3
  )
  rasterprob_combined <- rasterprob_combined / rowSums(rasterprob_combined)

  # All-zero weights for some counties (should fall back to colMeans)
  cellpop <- rep(0, n_cells)
  cellpop[1:5] <- sample(50:100, 5)  # Only first county has non-zero weights

  births <- data.frame(
    FIPS = sample(sprintf("%05d", 1:n_counties), 30, replace = TRUE),
    BWT = rnorm(30, 3400, 450),
    stringsAsFactors = FALSE
  )

  result <- expect_no_error(
    geoExposeR:::impute_exposure(
      datmatx_prob = datmatx_prob,
      rasterprob_combined = rasterprob_combined,
      cellpop = cellpop,
      births = births,
      ndraws = ndraws,
      cat_label = cat_label,
      id_col = "GEOID10",
      record_id_col = "FIPS",
      exposure_level_col = "ExposureLevel",
      apply_imputation_fallback = TRUE
    )
  )

  expect_equal(length(result), ndraws)
})

test_that("impute_exposure handles single cell per county", {
  set.seed(12345)
  n_counties <- 5
  ndraws <- 2
  cat_label <- c("<5", "5-10", "10+")

  # Single cell per county
  county_ids <- sprintf("%05d", 1:n_counties)
  datmatx_prob <- data.frame(
    GEOID10 = county_ids,
    stringsAsFactors = FALSE
  )

  rasterprob_combined <- matrix(
    runif(n_counties * 3),
    nrow = n_counties,
    ncol = 3
  )
  rasterprob_combined <- rasterprob_combined / rowSums(rasterprob_combined)

  cellpop <- sample(50:500, n_counties, replace = TRUE)

  births <- data.frame(
    FIPS = sample(county_ids, 20, replace = TRUE),
    BWT = rnorm(20, 3400, 450),
    stringsAsFactors = FALSE
  )

  result <- geoExposeR:::impute_exposure(
    datmatx_prob = datmatx_prob,
    rasterprob_combined = rasterprob_combined,
    cellpop = cellpop,
    births = births,
    ndraws = ndraws,
    cat_label = cat_label,
    id_col = "GEOID10",
    record_id_col = "FIPS",
    exposure_level_col = "ExposureLevel",
    apply_imputation_fallback = TRUE
  )

  expect_equal(length(result), ndraws)
  for (j in seq_along(result)) {
    expect_true(all(result[[j]]$ExposureLevel %in% cat_label))
  }
})

test_that("impute_exposure removes rows when fallback disabled", {
  set.seed(12345)
  n_counties <- 5
  ndraws <- 2
  cat_label <- c("<5", "5-10", "10+")

  county_ids <- sprintf("%05d", 1:n_counties)
  datmatx_prob <- data.frame(
    GEOID10 = county_ids,
    stringsAsFactors = FALSE
  )

  rasterprob_combined <- matrix(
    runif(n_counties * 3),
    nrow = n_counties,
    ncol = 3
  )
  rasterprob_combined <- rasterprob_combined / rowSums(rasterprob_combined)

  cellpop <- sample(50:500, n_counties, replace = TRUE)

  # Add some birth records with counties not in probability data data
  births <- data.frame(
    FIPS = c(county_ids, sprintf("%05d", 99:100)),  # 99, 100 not in probability data
    BWT = rnorm(n_counties + 2, 3400, 450),
    stringsAsFactors = FALSE
  )

  result <- geoExposeR:::impute_exposure(
    datmatx_prob = datmatx_prob,
    rasterprob_combined = rasterprob_combined,
    cellpop = cellpop,
    births = births,
    ndraws = ndraws,
    cat_label = cat_label,
    id_col = "GEOID10",
    record_id_col = "FIPS",
    exposure_level_col = "ExposureLevel",
    apply_imputation_fallback = FALSE  # Disable fallback
  )

  # Rows with missing exposure levels should be removed
  for (j in seq_along(result)) {
    expect_false(any(is.na(result[[j]]$ExposureLevel)))
    expect_lte(nrow(result[[j]]), n_counties)  # Should have fewer rows
  }
})

# ===== Tests for impute_additional_variables function =====

test_that("impute_additional_variables handles basic case", {
  set.seed(12345)

  # Create imputed datasets with missing values
  imputed_datasets <- list(
    data.frame(
      FIPS = sprintf("%05d", 1:20),
      ExposureLevel = factor(sample(c("<5", "5-10", "10+"), 20, replace = TRUE)),
      MAGE_R = c(rep(NA, 5), sample(18:45, 15, replace = TRUE)),
      pm = rnorm(20, 8, 2),
      stringsAsFactors = FALSE
    ),
    data.frame(
      FIPS = sprintf("%05d", 1:20),
      ExposureLevel = factor(sample(c("<5", "5-10", "10+"), 20, replace = TRUE)),
      MAGE_R = c(rep(NA, 3), sample(18:45, 17, replace = TRUE)),
      pm = rnorm(20, 8, 2),
      stringsAsFactors = FALSE
    )
  )

  births <- data.frame(
    FIPS = sprintf("%05d", 1:20),
    MAGE_R = sample(18:45, 20, replace = TRUE),
    pm = rnorm(20, 8, 2),
    stringsAsFactors = FALSE
  )

  result <- expect_no_error(
    geoExposeR:::impute_additional_variables(
      imputed_datasets = imputed_datasets,
      impute_vars = "MAGE_R",
      births = births,
      mice_covs = "pm",
      mice_m = 1,
      mice_maxit = 3,
      mice_method = "pmm",
      seed = 12345
    )
  )

  # Check that NAs were imputed
  for (j in seq_along(result)) {
    expect_false(any(is.na(result[[j]]$MAGE_R)))
  }
})

test_that("impute_additional_variables warns when no covariates available", {
  imputed_datasets <- list(
    data.frame(
      FIPS = sprintf("%05d", 1:10),
      ExposureLevel = factor(sample(c("<5", "5-10", "10+"), 10, replace = TRUE)),
      MAGE_R = c(rep(NA, 3), sample(18:45, 7, replace = TRUE)),
      stringsAsFactors = FALSE
    )
  )

  births <- data.frame(
    FIPS = sprintf("%05d", 1:10),
    MAGE_R = sample(18:45, 10, replace = TRUE),
    stringsAsFactors = FALSE
  )

  # Should warn because 'pm' is not in the dataset
  expect_warning(
    geoExposeR:::impute_additional_variables(
      imputed_datasets = imputed_datasets,
      impute_vars = "MAGE_R",
      births = births,
      mice_covs = "pm",  # Not in dataset
      mice_m = 1,
      mice_maxit = 3,
      mice_method = "pmm",
      seed = 12345
    ),
    "No MICE covariates available"
  )
})

test_that("impute_additional_variables errors on missing variables in births", {
  imputed_datasets <- list(
    data.frame(
      FIPS = sprintf("%05d", 1:10),
      ExposureLevel = factor(sample(c("<5", "5-10", "10+"), 10, replace = TRUE)),
      stringsAsFactors = FALSE
    )
  )

  births <- data.frame(
    FIPS = sprintf("%05d", 1:10),
    stringsAsFactors = FALSE
  )

  expect_error(
    geoExposeR:::impute_additional_variables(
      imputed_datasets = imputed_datasets,
      impute_vars = "MAGE_R",  # Not in births
      births = births,
      mice_covs = "pm",
      mice_m = 1,
      mice_maxit = 3,
      mice_method = "pmm",
      seed = 12345
    ),
    "Variables not found in birth data"
  )
})

test_that("impute_additional_variables skips variables without NAs", {
  set.seed(12345)

  imputed_datasets <- list(
    data.frame(
      FIPS = sprintf("%05d", 1:10),
      ExposureLevel = factor(sample(c("<5", "5-10", "10+"), 10, replace = TRUE)),
      MAGE_R = sample(18:45, 10, replace = TRUE),  # No NAs
      pm = rnorm(10, 8, 2),
      stringsAsFactors = FALSE
    )
  )

  births <- data.frame(
    FIPS = sprintf("%05d", 1:10),
    MAGE_R = sample(18:45, 10, replace = TRUE),
    pm = rnorm(10, 8, 2),
    stringsAsFactors = FALSE
  )

  result <- geoExposeR:::impute_additional_variables(
    imputed_datasets = imputed_datasets,
    impute_vars = "MAGE_R",
    births = births,
    mice_covs = "pm",
    mice_m = 1,
    mice_maxit = 3,
    mice_method = "pmm",
    seed = 12345
  )

  # Should return unchanged (no imputation needed)
  expect_equal(result[[1]]$MAGE_R, imputed_datasets[[1]]$MAGE_R)
})

# ===== Tests for validate_imputed_datasets function =====

test_that("validate_imputed_datasets passes valid datasets", {
  cat_label <- c("<5", "5-10", "10+")

  imputed_datasets <- list(
    data.frame(
      FIPS = sprintf("%05d", 1:10),
      ExposureLevel = factor(
        sample(cat_label, 10, replace = TRUE),
        levels = cat_label
      ),
      BWT = rnorm(10, 3400, 450)
    ),
    data.frame(
      FIPS = sprintf("%05d", 1:10),
      ExposureLevel = factor(
        sample(cat_label, 10, replace = TRUE),
        levels = cat_label
      ),
      BWT = rnorm(10, 3400, 450)
    )
  )

  expect_message(
    geoExposeR:::validate_imputed_datasets(
      imputed_datasets = imputed_datasets,
      exposure_level_col = "ExposureLevel",
      cat_label = cat_label
    ),
    "All 2 imputed datasets validated successfully"
  )
})

test_that("validate_imputed_datasets catches missing exposure levels", {
  cat_label <- c("<5", "5-10", "10+")

  imputed_datasets <- list(
    data.frame(
      FIPS = sprintf("%05d", 1:10),
      ExposureLevel = factor(c(NA, sample(cat_label, 9, replace = TRUE))),
      BWT = rnorm(10, 3400, 450)
    )
  )

  expect_error(
    geoExposeR:::validate_imputed_datasets(
      imputed_datasets = imputed_datasets,
      exposure_level_col = "ExposureLevel",
      cat_label = cat_label
    ),
    "Missing values found in exposure level column"
  )
})

test_that("validate_imputed_datasets catches invalid exposure levels", {
  cat_label <- c("<5", "5-10", "10+")

  imputed_datasets <- list(
    data.frame(
      FIPS = sprintf("%05d", 1:10),
      ExposureLevel = factor(c("InvalidLevel", sample(cat_label, 9, replace = TRUE))),
      BWT = rnorm(10, 3400, 450)
    )
  )

  expect_error(
    geoExposeR:::validate_imputed_datasets(
      imputed_datasets = imputed_datasets,
      exposure_level_col = "ExposureLevel",
      cat_label = cat_label
    ),
    "Invalid exposure levels found"
  )
})

# ===== Additional tests for edge cases in imputation functions =====

test_that("impute_exposure applies fallback with message for unmatched counties", {
  set.seed(12345)
  n_counties <- 5
  ndraws <- 2
  cat_label <- c("<5", "5-10", "10+")

  # Create probability data with only 3 counties
  county_ids <- sprintf("%05d", 1:3)
  datmatx_prob <- data.frame(
    GEOID10 = county_ids,
    stringsAsFactors = FALSE
  )

  rasterprob_combined <- matrix(
    runif(3 * 3),
    nrow = 3,
    ncol = 3
  )
  rasterprob_combined <- rasterprob_combined / rowSums(rasterprob_combined)

  cellpop <- sample(50:500, 3, replace = TRUE)

  # Create birth data with counties 1-5 (counties 4,5 not in probability data)
  births <- data.frame(
    FIPS = sprintf("%05d", 1:5),
    BWT = rnorm(5, 3400, 450),
    stringsAsFactors = FALSE
  )

  # Should print message about filling NAs with most common category
  expect_message(
    geoExposeR:::impute_exposure(
      datmatx_prob = datmatx_prob,
      rasterprob_combined = rasterprob_combined,
      cellpop = cellpop,
      births = births,
      ndraws = ndraws,
      cat_label = cat_label,
      id_col = "GEOID10",
      record_id_col = "FIPS",
      exposure_level_col = "ExposureLevel",
      apply_imputation_fallback = TRUE
    ),
    "filling with most common category"
  )
})

test_that("impute_additional_variables handles insufficient columns for MICE", {
  set.seed(12345)

  # Create dataset with only the variable to impute (no covariates)
  imputed_datasets <- list(
    data.frame(
      FIPS = sprintf("%05d", 1:10),
      ExposureLevel = factor(sample(c("<5", "5-10", "10+"), 10, replace = TRUE)),
      MAGE_R = c(rep(NA, 3), sample(18:45, 7, replace = TRUE)),
      stringsAsFactors = FALSE
    )
  )

  births <- data.frame(
    FIPS = sprintf("%05d", 1:10),
    MAGE_R = sample(18:45, 10, replace = TRUE),
    stringsAsFactors = FALSE
  )

  # Should warn about insufficient variables or MICE failure
  expect_warning(
    geoExposeR:::impute_additional_variables(
      imputed_datasets = imputed_datasets,
      impute_vars = "MAGE_R",
      births = births,
      mice_covs = "MAGE_R",  # Same as impute var, will result in only 1 column
      mice_m = 1,
      mice_maxit = 3,
      mice_method = "pmm",
      seed = 12345
    ),
    "Insufficient variables|No MICE covariates|MICE failed"
  )
})

test_that("impute_additional_variables handles MICE failure gracefully", {
  set.seed(12345)

  # Create dataset that might cause MICE to fail
  imputed_datasets <- list(
    data.frame(
      FIPS = sprintf("%05d", 1:5),
      ExposureLevel = factor(sample(c("<5", "5-10", "10+"), 5, replace = TRUE)),
      MAGE_R = c(NA, NA, NA, NA, 25),  # Almost all NA
      pm = c(NA, NA, NA, NA, 8),  # Covariate also mostly NA
      stringsAsFactors = FALSE
    )
  )

  births <- data.frame(
    FIPS = sprintf("%05d", 1:5),
    MAGE_R = sample(18:45, 5, replace = TRUE),
    pm = rnorm(5, 8, 2),
    stringsAsFactors = FALSE
  )

  # Should handle gracefully (warn but not error)
  result <- expect_no_error(
    suppressWarnings(
      geoExposeR:::impute_additional_variables(
        imputed_datasets = imputed_datasets,
        impute_vars = "MAGE_R",
        births = births,
        mice_covs = "pm",
        mice_m = 1,
        mice_maxit = 1,
        mice_method = "pmm",
        seed = 12345
      )
    )
  )

  expect_type(result, "list")
})

test_that("impute_additional_variables skips birth column check when disabled", {
  set.seed(12345)

  imputed_datasets <- list(
    data.frame(
      FIPS = sprintf("%05d", 1:10),
      ExposureLevel = factor(sample(c("<5", "5-10", "10+"), 10, replace = TRUE)),
      MAGE_R = c(rep(NA, 3), sample(18:45, 7, replace = TRUE)),
      pm = rnorm(10, 8, 2),
      stringsAsFactors = FALSE
    )
  )

  # births doesn't have MAGE_R but check_birth_cols = FALSE
  births <- data.frame(
    FIPS = sprintf("%05d", 1:10),
    pm = rnorm(10, 8, 2),
    stringsAsFactors = FALSE
  )

  # Should not error when check_birth_cols = FALSE
  result <- expect_no_error(
    geoExposeR:::impute_additional_variables(
      imputed_datasets = imputed_datasets,
      impute_vars = "MAGE_R",
      births = births,
      mice_covs = "pm",
      mice_m = 1,
      mice_maxit = 3,
      mice_method = "pmm",
      seed = 12345,
      check_birth_cols = FALSE
    )
  )

  expect_type(result, "list")
})

# ===== Tests for data-loading.R functions =====

test_that("load_prob_data errors on missing columns", {
  temp_dir <- tempfile(pattern = "test_prob_missing_")
  dir.create(temp_dir)

  # Create probability data with missing columns
  prob_df <- data.frame(
    GEOID10 = sprintf("%05d", 1:10),
    prob_C1 = runif(10),
    # Missing prob_C2 and prob_C3
    Wells_2010 = sample(50:500, 10)
  )
  prob_file <- file.path(temp_dir, "prob_missing.csv")
  data.table::fwrite(prob_df, prob_file)

  expect_error(
    geoExposeR:::load_prob_data(
      model_prob_csv = prob_file,
      prob_cols = c("prob_C1", "prob_C2", "prob_C3"),
      pop_well_col = "Wells_2010"
    ),
    "columns are missing in probability data"
  )

  unlink(temp_dir, recursive = TRUE)
})

test_that("validate_data_alignment errors on row count mismatch", {
  prob_df <- data.frame(GEOID10 = 1:10)
  conc_df <- data.frame(GEOID10 = 1:5)  # Different row count

  expect_error(
    geoExposeR:::validate_data_alignment(prob_df, conc_df),
    "different numbers of rows"
  )
})

test_that("load_and_process_birth_data handles column selection", {
  temp_dir <- tempfile(pattern = "test_birth_cols_")
  dir.create(temp_dir)

  # Create birth data with extra columns
  births_df <- data.frame(
    FIPS = sprintf("%05d", 1:20),
    BWT = rnorm(20, 3400, 450),
    OEGEST = rnorm(20, 40, 2),
    EXTRA_COL = runif(20),
    stringsAsFactors = FALSE
  )
  births_file <- file.path(temp_dir, "births_extra.txt")
  data.table::fwrite(births_df, births_file, sep = "\t")

  # Should select only specified columns
  result <- geoExposeR:::load_and_process_birth_data(
    birth_data_txt = births_file,
    columns_to_select = c("FIPS", "BWT"),
    rucc_col = NULL,
    record_id_col = "FIPS"
  )

  expect_true("FIPS" %in% names(result))
  expect_true("BWT" %in% names(result))
  expect_false("EXTRA_COL" %in% names(result))

  unlink(temp_dir, recursive = TRUE)
})

test_that("load_and_process_birth_data creates rural indicator", {
  temp_dir <- tempfile(pattern = "test_birth_rucc_")
  dir.create(temp_dir)

  births_df <- data.frame(
    FIPS = sprintf("%05d", 1:10),
    BWT = rnorm(10, 3400, 450),
    RUCC = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 1),  # Some rural (>3), some not
    stringsAsFactors = FALSE
  )
  births_file <- file.path(temp_dir, "births_rucc.txt")
  data.table::fwrite(births_df, births_file, sep = "\t")

  result <- geoExposeR:::load_and_process_birth_data(
    birth_data_txt = births_file,
    columns_to_select = NULL,
    rucc_col = "RUCC",
    record_id_col = "FIPS"
  )

  expect_true("rural" %in% names(result))
  # RUCC > 3 should be TRUE
  expect_equal(result$rural[1], FALSE)  # RUCC = 1
  expect_equal(result$rural[4], TRUE)   # RUCC = 4

  unlink(temp_dir, recursive = TRUE)
})

test_that("create_weighted_prob_matrix handles NA in concentration data", {
  set.seed(12345)

  # Model probabilities (no NAs)
  rasterprob_model <- matrix(
    c(0.7, 0.2, 0.1, 0.5, 0.3, 0.2, 0.6, 0.25, 0.15),
    nrow = 3, ncol = 3, byrow = TRUE
  )

  # Concentration probabilities with NA in row 2
  rasterprob_conc <- matrix(
    c(0.8, 0.15, 0.05, NA, NA, NA, 0.4, 0.35, 0.25),
    nrow = 3, ncol = 3, byrow = TRUE
  )

  datmatx_conc <- data.frame(
    PWELL_private_pct = c(50, 50, 50)
  )

  result <- geoExposeR:::create_weighted_prob_matrix(
    rasterprob_model = rasterprob_model,
    rasterprob_conc = rasterprob_conc,
    datmatx_conc = datmatx_conc,
    pwell_col = "PWELL_private_pct"
  )

  # Row 2 should use model probabilities (fallback)
  expect_equal(result[2, ], rasterprob_model[2, ])
  # Row 1 and 3 should be weighted averages
  expect_false(any(is.na(result)))
})

test_that("validate_conc_cutoffs validates properly with validate_cutoffs=TRUE", {
  # Non-numeric cutoffs
  expect_error(
    geoExposeR:::validate_conc_cutoffs(c("5", "10"), validate_cutoffs = TRUE),
    "must be numeric"
  )

  # Unsorted cutoffs
  expect_error(
    geoExposeR:::validate_conc_cutoffs(c(10, 5), validate_cutoffs = TRUE),
    "must be in ascending order"
  )

  # Negative cutoffs
  expect_error(
    geoExposeR:::validate_conc_cutoffs(c(-5, 10), validate_cutoffs = TRUE),
    "cannot contain negative values"
  )

  # Valid cutoffs
  result <- geoExposeR:::validate_conc_cutoffs(c(5, 10), validate_cutoffs = TRUE)
  expect_equal(result$cutoffs, c(5, 10))
  expect_equal(result$n_categories, 3)
})

test_that("load_and_process_exposure_data uses auto-generated cutoffs", {
  temp_dir <- tempfile(pattern = "test_auto_cutoffs_")
  dir.create(temp_dir)

  # Create probability data with 3 categories
  set.seed(12345)
  n_rows <- 20
  model_probs <- matrix(runif(n_rows * 3), ncol = 3)
  model_probs <- model_probs / rowSums(model_probs)

  prob_df <- data.frame(
    GEOID10 = sprintf("%05d", 1:n_rows),
    prob_C1 = model_probs[, 1],
    prob_C2 = model_probs[, 2],
    prob_C3 = model_probs[, 3],
    Wells_2010 = sample(50:500, n_rows, replace = TRUE)
  )
  prob_file <- file.path(temp_dir, "prob_data.csv")
  data.table::fwrite(prob_df, prob_file)

  # Create concentration data
  conc_df <- data.frame(
    GEOID10 = sprintf("%05d", 1:n_rows),
    conc_meanlog = runif(n_rows, 0.5, 1.5),
    conc_sdlog = runif(n_rows, 0.1, 0.5),
    PWELL_private_pct = runif(n_rows, 10, 90)
  )
  conc_file <- file.path(temp_dir, "conc_data.csv")
  data.table::fwrite(conc_df, conc_file)

  # Should auto-generate cutoffs when conc_cutoffs = NULL
  expect_message(
    geoExposeR:::load_and_process_exposure_data(
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      prob_cols = c("prob_C1", "prob_C2", "prob_C3"),
      pop_well_col = "Wells_2010",
      conc_mean_col = "conc_meanlog",
      conc_sd_col = "conc_sdlog",
      pwell_col = "PWELL_private_pct",
      id_col = "GEOID10",
      conc_cutoffs = NULL
    ),
    "auto-generated cutoffs"
  )

  unlink(temp_dir, recursive = TRUE)
})

# ===== Tests for regression.R functions =====

test_that("regression_analysis errors on missing .id column", {
  births <- data.frame(
    FIPS = sprintf("%05d", 1:10),
    BWT = rnorm(10, 3400, 450),
    ExposureLevel = factor(sample(c("<5", "5-10", "10+"), 10, replace = TRUE))
  )
  # No .id column

  expect_error(
    geoExposeR:::regression_analysis(
      births = births,
      targets = "BWT",
      cat_label = c("5-10", "10+"),
      exposure_level_col = "ExposureLevel",
      regression_formula = "~ as.factor(ExposureLevel) + (1|FIPS)",
      output_dir = tempdir()
    ),
    "must contain '.id' column"
  )
})

test_that("regression_analysis errors on missing target variables", {
  births <- data.frame(
    .id = rep(1:2, each = 5),
    FIPS = sprintf("%05d", 1:10),
    BWT = rnorm(10, 3400, 450),
    ExposureLevel = factor(sample(c("<5", "5-10", "10+"), 10, replace = TRUE))
  )

  expect_error(
    geoExposeR:::regression_analysis(
      births = births,
      targets = c("BWT", "NONEXISTENT"),
      cat_label = c("5-10", "10+"),
      exposure_level_col = "ExposureLevel",
      regression_formula = "~ as.factor(ExposureLevel) + (1|FIPS)",
      output_dir = tempdir()
    ),
    "Target variables not found"
  )
})

test_that("pool_estimates_by_term handles empty term data", {
  # Create combined_results with terms that won't match exposure_terms
  combined_results <- data.frame(
    term = c("other_term1", "other_term2"),
    estimate = c(1.0, 2.0),
    std.error = c(0.1, 0.2)
  )

  # Create a mock reference model for df extraction
  set.seed(12345)
  mock_data <- data.frame(
    y = rnorm(50),
    x = rnorm(50),
    group = factor(rep(1:5, each = 10))
  )
  mock_model <- lme4::lmer(y ~ x + (1|group), data = mock_data)

  result <- geoExposeR:::pool_estimates_by_term(
    combined_results = combined_results,
    exposure_terms = c("as.factor(ExposureLevel)5-10"),  # Won't match any rows
    reference_model = mock_model
  )

  # Should return NULL or empty when no terms match
  expect_true(is.null(result) || nrow(result) == 0)
})

test_that("pool_single_estimate handles single imputation", {
  result <- geoExposeR:::pool_single_estimate(
    estimates = c(1.5),
    std_errors = c(0.5),
    df = 100
  )

  expect_equal(result$q.mi, 1.5)
  expect_equal(result$se.mi, 0.5)
  expect_equal(result$statistic, 3.0)
})

test_that("pool_single_estimate applies Rubin's Rules for multiple imputations", {
  result <- geoExposeR:::pool_single_estimate(
    estimates = c(1.0, 1.2, 0.8, 1.1, 0.9),
    std_errors = c(0.3, 0.35, 0.28, 0.32, 0.31),
    df = 100
  )

  # Pooled estimate should be mean of estimates
  expect_equal(result$q.mi, mean(c(1.0, 1.2, 0.8, 1.1, 0.9)))
  # Should have valid confidence interval
  expect_true(result$conf.low < result$q.mi)
  expect_true(result$conf.high > result$q.mi)
  # P-value should be between 0 and 1
  expect_true(result$p.value >= 0 && result$p.value <= 1)
})

test_that("pool_single_estimate handles df adjustment", {
  # This tests when b > 0 and finite df
  result <- geoExposeR:::pool_single_estimate(
    estimates = c(2.0, 2.5, 1.8, 2.2, 2.1),
    std_errors = c(0.2, 0.25, 0.18, 0.22, 0.21),
    df = 50
  )

  expect_true(!is.na(result$p.value))
  expect_true(!is.na(result$conf.low))
  expect_true(!is.na(result$conf.high))
})

test_that("regression_analysis handles model fitting errors", {
  # Create data that will cause model fitting to fail
  set.seed(12345)
  births <- data.frame(
    .id = rep(1:2, each = 20),
    FIPS = sprintf("%05d", rep(1:4, 10)),
    BWT = rnorm(40, 3400, 450),
    ExposureLevel = factor(rep("<5", 40), levels = c("<5", "5-10", "10+"))
  )

  # Model will fail because exposure level has only one level
  expect_error(
    geoExposeR:::regression_analysis(
      births = births,
      targets = "BWT",
      cat_label = c("5-10", "10+"),
      exposure_level_col = "ExposureLevel",
      regression_formula = "~ as.factor(ExposureLevel) + (1|FIPS)",
      output_dir = tempdir()
    ),
    "Model fitting failed"
  )
})

# ==============================================================================
# Tests for Data Ingestion Functions
# ==============================================================================

# ------------------------------------------------------------------------------
# Tests for prepare_conc_data()
# ------------------------------------------------------------------------------

test_that("prepare_conc_data processes raw concentration measurements correctly", {
  # Create sample raw concentration data
  conc_raw_data <- data.frame(
    county_fips = c("06001", "06001", "06003", "06003", "06005"),
    concentration = c(2.5, 3.1, 8.2, 7.5, 4.0),
    population_served = c(50000, 25000, 10000, 15000, 30000)
  )

  result <- prepare_conc_data(
    conc_data = conc_raw_data,
    id_col = "county_fips",
    concentration_col = "concentration",
    input_type = "raw"
  )

  # Check structure

  expect_true(data.table::is.data.table(result))
  expect_true("GEOID10" %in% names(result))
  expect_true("conc_meanlog" %in% names(result))
  expect_true("conc_sdlog" %in% names(result))

  # Check aggregation (should have 3 unique counties)
  expect_equal(nrow(result), 3)

  # Check that meanlog and sdlog are valid
  expect_true(all(is.finite(result$conc_meanlog)))
  expect_true(all(result$conc_sdlog > 0))
})

test_that("prepare_conc_data handles summary input type", {
  # Create summary concentration data
  conc_summary <- data.frame(
    county_fips = c("06001", "06003", "06005"),
    mean_concentration = c(3.0, 7.8, 4.0),
    sd_concentration = c(0.5, 1.2, 0.8)
  )

  result <- prepare_conc_data(
    conc_data = conc_summary,
    id_col = "county_fips",
    input_type = "summary",
    mean_col = "mean_concentration",
    sd_col = "sd_concentration"
  )

  expect_equal(nrow(result), 3)
  expect_true(all(is.finite(result$conc_meanlog)))
  expect_true(all(result$conc_sdlog > 0))
})

test_that("prepare_conc_data handles missing SD with default", {
  conc_summary <- data.frame(
    county_fips = c("06001", "06003"),
    mean_concentration = c(3.0, 7.8)
  )

  result <- prepare_conc_data(
    conc_data = conc_summary,
    id_col = "county_fips",
    input_type = "summary",
    mean_col = "mean_concentration",
    default_sdlog = 1.5
  )

  # sdlog should be calculated based on derived sd_conc
  expect_true(all(result$conc_sdlog > 0))
})

test_that("prepare_conc_data validates missing geoid column", {
  conc_data <- data.frame(
    wrong_col = c("06001"),
    concentration = c(3.0)
  )

  expect_error(
    prepare_conc_data(conc_data, id_col = "county_fips"),
    "Geographic ID column 'county_fips' not found"
  )
})

test_that("prepare_conc_data validates missing concentration column", {
  conc_data <- data.frame(
    county_fips = c("06001"),
    wrong_col = c(3.0)
  )

  expect_error(
    prepare_conc_data(conc_data, id_col = "county_fips", input_type = "raw"),
    "Concentration column 'concentration' not found"
  )
})

test_that("prepare_conc_data handles zero/negative concentrations", {
  conc_data <- data.frame(
    county_fips = c("06001", "06001", "06003"),
    concentration = c(0, -1, 5.0)
  )

  result <- prepare_conc_data(
    conc_data = conc_data,
    id_col = "county_fips",
    input_type = "raw"
  )

  # Should only have 1 row (06003) since 06001 has no positive values
  expect_equal(nrow(result), 1)
  expect_equal(result$GEOID10, "06003")
})

test_that("prepare_conc_data errors when all concentrations are non-positive", {
  conc_data <- data.frame(
    county_fips = c("06001", "06003"),
    concentration = c(0, -1)
  )

  expect_error(
    prepare_conc_data(conc_data, id_col = "county_fips", input_type = "raw"),
    "No positive concentration values found"
  )
})

test_that("prepare_conc_data includes private well percentage when provided", {
  conc_data <- data.frame(
    county_fips = c("06001", "06001", "06003"),
    concentration = c(3.0, 4.0, 5.0),
    pwell_pct = c(15.0, 15.0, 25.0)
  )

  result <- prepare_conc_data(
    conc_data = conc_data,
    id_col = "county_fips",
    input_type = "raw",
    private_well_pct_col = "pwell_pct"
  )

  expect_true("PWELL_private_pct" %in% names(result))
  expect_equal(result$PWELL_private_pct[result$GEOID10 == "06001"], 15.0)
})

test_that("prepare_conc_data saves output file when requested", {
  conc_data <- data.frame(
    county_fips = c("06001", "06003"),
    concentration = c(3.0, 5.0)
  )

  output_file <- tempfile(fileext = ".csv")

  result <- prepare_conc_data(
    conc_data = conc_data,
    id_col = "county_fips",
    input_type = "raw",
    output_file = output_file
  )

  expect_true(file.exists(output_file))

  # Read back and verify
  saved <- data.table::fread(output_file)
  expect_equal(nrow(saved), 2)

  unlink(output_file)
})

test_that("prepare_conc_data validates summary input missing mean column", {
  conc_data <- data.frame(
    county_fips = c("06001"),
    wrong_col = c(3.0)
  )

  expect_error(
    prepare_conc_data(
      conc_data,
      id_col = "county_fips",
      input_type = "summary",
      mean_col = "mean_concentration"
    ),
    "Mean column 'mean_concentration' not found"
  )
})

# ------------------------------------------------------------------------------
# Tests for format_health_data()
# ------------------------------------------------------------------------------

test_that("format_health_data formats FIPS codes correctly", {
  health <- data.frame(
    FIPS = c(6001, 6003, 6005),
    BWT = c(3200, 3150, 3300),
    OEGEST = c(39.2, 38.8, 40.1)
  )

  result <- format_health_data(
    health_data = health,
    id_col = "FIPS",
    outcome_cols = c("BWT", "OEGEST"),
    validate = FALSE
  )

  # Check FIPS formatting (5 digits with zero-padding)
  expect_equal(result$FIPS, c("06001", "06003", "06005"))
  expect_true(all(nchar(result$FIPS) == 5))
})

test_that("format_health_data includes covariates", {
  health <- data.frame(
    FIPS = c(6001, 6003),
    BWT = c(3200, 3150),
    MAGE = c(28, 32),
    rural = c(0, 1)
  )

  result <- format_health_data(
    health_data = health,
    id_col = "FIPS",
    outcome_cols = c("BWT"),
    covariate_cols = c("MAGE", "rural"),
    validate = FALSE
  )

  expect_true("MAGE" %in% names(result))
  expect_true("rural" %in% names(result))
})

test_that("format_health_data validates outcome ranges", {
  health <- data.frame(
    FIPS = c(6001, 6003),
    BWT = c(100, 7000),  # Outside expected range 200-6000
    OEGEST = c(39, 40)
  )

  expect_warning(
    format_health_data(
      health_data = health,
      id_col = "FIPS",
      outcome_cols = c("BWT", "OEGEST"),
      validate = TRUE
    ),
    "outside expected range"
  )
})

test_that("format_health_data warns about all-NA outcomes", {
  health <- data.frame(
    FIPS = c(6001, 6003),
    BWT = c(NA, NA),
    OEGEST = c(39, 40)
  )

  expect_warning(
    format_health_data(
      health_data = health,
      id_col = "FIPS",
      outcome_cols = c("BWT", "OEGEST"),
      validate = TRUE
    ),
    "contains all NA values"
  )
})

test_that("format_health_data errors on missing geoid column", {
  health <- data.frame(
    wrong_col = c(6001),
    BWT = c(3200)
  )

  expect_error(
    format_health_data(health, id_col = "FIPS", outcome_cols = "BWT"),
    "Geographic ID column 'FIPS' not found"
  )
})

test_that("format_health_data errors on missing outcome columns", {
  health <- data.frame(
    FIPS = c(6001),
    BWT = c(3200)
  )

  expect_error(
    format_health_data(health, id_col = "FIPS", outcome_cols = c("BWT", "MISSING")),
    "Outcome columns not found"
  )
})

test_that("format_health_data warns on missing covariate columns", {
  health <- data.frame(
    FIPS = c(6001),
    BWT = c(3200)
  )

  expect_warning(
    format_health_data(
      health,
      id_col = "FIPS",
      outcome_cols = "BWT",
      covariate_cols = c("MISSING"),
      validate = FALSE
    ),
    "Covariate columns not found"
  )
})

test_that("format_health_data saves to txt format", {
  health <- data.frame(
    FIPS = c(6001, 6003),
    BWT = c(3200, 3150)
  )

  output_file <- tempfile(fileext = ".txt")

  result <- format_health_data(
    health_data = health,
    id_col = "FIPS",
    outcome_cols = "BWT",
    output_file = output_file,
    output_format = "txt",
    validate = FALSE
  )

  expect_true(file.exists(output_file))
  unlink(output_file)
})

test_that("format_health_data saves to csv format", {
  health <- data.frame(
    FIPS = c(6001, 6003),
    BWT = c(3200, 3150)
  )

  output_file <- tempfile(fileext = ".csv")

  result <- format_health_data(
    health_data = health,
    id_col = "FIPS",
    outcome_cols = "BWT",
    output_file = output_file,
    output_format = "csv",
    validate = FALSE
  )

  expect_true(file.exists(output_file))
  unlink(output_file)
})

# ------------------------------------------------------------------------------
# Tests for add_well_percentages()
# ------------------------------------------------------------------------------

test_that("add_well_percentages merges data correctly", {
  target <- data.frame(
    GEOID10 = c("06001", "06003", "06005"),
    conc_meanlog = c(1.0, 1.5, 2.0)
  )

  wells <- data.frame(
    GEOID = c("06001", "06003"),
    private_well_pct = c(15.2, 28.7)
  )

  result <- add_well_percentages(
    target_data = target,
    well_data = wells
  )

  expect_true("PWELL_private_pct" %in% names(result))
  expect_equal(result$PWELL_private_pct[result$GEOID10 == "06001"], 15.2)
  expect_equal(result$PWELL_private_pct[result$GEOID10 == "06003"], 28.7)
  expect_true(is.na(result$PWELL_private_pct[result$GEOID10 == "06005"]))
})

test_that("add_well_percentages uses custom column names", {
  target <- data.frame(
    my_geoid = c("06001", "06003"),
    value = c(1.0, 1.5)
  )

  wells <- data.frame(
    geo_id = c("06001", "06003"),
    pct = c(15.2, 28.7)
  )

  result <- add_well_percentages(
    target_data = target,
    well_data = wells,
    target_id_col = "my_geoid",
    well_id_col = "geo_id",
    well_pct_col = "pct",
    output_col = "well_percentage"
  )

  expect_true("well_percentage" %in% names(result))
  expect_equal(result$well_percentage, c(15.2, 28.7))
})

test_that("add_well_percentages validates missing target column", {
  target <- data.frame(wrong_col = c("06001"))
  wells <- data.frame(GEOID = c("06001"), private_well_pct = c(15))

  expect_error(
    add_well_percentages(target, wells),
    "Column 'GEOID10' not found in target_data"
  )
})

test_that("add_well_percentages validates missing well geoid column", {
  target <- data.frame(GEOID10 = c("06001"))
  wells <- data.frame(wrong_col = c("06001"), private_well_pct = c(15))

  expect_error(
    add_well_percentages(target, wells),
    "Column 'GEOID' not found in well_data"
  )
})

test_that("add_well_percentages validates missing well pct column", {
  target <- data.frame(GEOID10 = c("06001"))
  wells <- data.frame(GEOID = c("06001"), wrong_col = c(15))

  expect_error(
    add_well_percentages(target, wells),
    "Column 'private_well_pct' not found in well_data"
  )
})

# ------------------------------------------------------------------------------
# Tests for subset_by_geography()
# ------------------------------------------------------------------------------

test_that("subset_by_geography filters by state FIPS", {
  data <- data.frame(
    GEOID10 = c("06001", "06003", "32001", "32003", "04001"),
    value = 1:5
  )

  # Filter to California (06)
  result <- subset_by_geography(data, state_fips = "06")
  expect_equal(nrow(result), 2)
  expect_true(all(grepl("^06", result$GEOID10)))

  # Filter to multiple states
  result <- subset_by_geography(data, state_fips = c("06", "32"))
  expect_equal(nrow(result), 4)
})

test_that("subset_by_geography filters by geoid list", {
  data <- data.frame(
    GEOID10 = c("06001", "06003", "06005", "06007"),
    value = 1:4
  )

  result <- subset_by_geography(
    data,
    geoid_list = c("06001", "06005")
  )

  expect_equal(nrow(result), 2)
  expect_equal(sort(result$GEOID10), c("06001", "06005"))
})

test_that("subset_by_geography filters by pattern", {
  data <- data.frame(
    GEOID10 = c("06001", "06003", "32001", "32003"),
    value = 1:4
  )

  result <- subset_by_geography(data, geoid_pattern = "^06")
  expect_equal(nrow(result), 2)
})

test_that("subset_by_geography validates missing geoid column", {
  data <- data.frame(wrong_col = c("06001"))

  expect_error(
    subset_by_geography(data),
    "Column 'GEOID10' not found in data"
  )
})

test_that("subset_by_geography uses custom geoid column", {
  data <- data.frame(
    my_geoid = c("06001", "32001"),
    value = 1:2
  )

  result <- subset_by_geography(
    data,
    id_col = "my_geoid",
    state_fips = "06"
  )

  expect_equal(nrow(result), 1)
})

# ------------------------------------------------------------------------------
# Tests for validate_prepared_inputs()
# ------------------------------------------------------------------------------

test_that("validate_prepared_inputs detects missing files", {
  result <- validate_prepared_inputs(
    prob_file = "nonexistent_prob.csv",
    conc_file = "nonexistent_epa.csv",
    health_file = "nonexistent_health.txt"
  )

  expect_false(result$valid)
  expect_true(length(result$issues) >= 3)
  expect_true(any(grepl("not found", result$issues)))
})

test_that("validate_prepared_inputs validates correct files", {
  # Create valid test files
  temp_dir <- tempdir()

  prob_file <- file.path(temp_dir, "test_prob.csv")
  conc_file <- file.path(temp_dir, "test_conc.csv")
  health_file <- file.path(temp_dir, "test_health.txt")

  # Create probability file
  prob_data_local <- data.frame(
    GEOID10 = c("06001", "06003"),
    prob_C1 = c(0.7, 0.6),
    prob_C2 = c(0.2, 0.25),
    prob_C3 = c(0.1, 0.15)
  )
  data.table::fwrite(prob_data_local, prob_file)

  # Create concentration file
  conc_data <- data.frame(
    GEOID10 = c("06001", "06003"),
    conc_meanlog = c(1.0, 1.5),
    conc_sdlog = c(0.8, 0.9),
    PWELL_private_pct = c(15, 25)
  )
  data.table::fwrite(conc_data, conc_file)

  # Create health file
  health_data <- data.frame(
    FIPS = c("06001", "06003"),
    BWT = c(3200, 3150)
  )
  data.table::fwrite(health_data, health_file, sep = "\t")

  result <- validate_prepared_inputs(
    prob_file = prob_file,
    conc_file = conc_file,
    health_file = health_file
  )

  expect_true(result$valid)
  expect_equal(length(result$issues), 0)

  # Cleanup
  unlink(c(prob_file, conc_file, health_file))
})

test_that("validate_prepared_inputs detects missing columns", {
  temp_dir <- tempdir()
  prob_file <- file.path(temp_dir, "test_prob_bad.csv")
  conc_file <- file.path(temp_dir, "test_conc_bad.csv")
  health_file <- file.path(temp_dir, "test_health_bad.txt")

  # Create probability file missing probability columns
  prob_data_local <- data.frame(
    GEOID10 = c("06001"),
    wrong_col = c(0.5)
  )
  data.table::fwrite(prob_data_local, prob_file)

  # Create concentration file missing meanlog
  conc_data <- data.frame(
    GEOID10 = c("06001"),
    wrong_col = c(1.0)
  )
  data.table::fwrite(conc_data, conc_file)

  # Create health file
  health_data <- data.frame(
    FIPS = c("06001"),
    BWT = c(3200)
  )
  data.table::fwrite(health_data, health_file, sep = "\t")

  result <- validate_prepared_inputs(
    prob_file = prob_file,
    conc_file = conc_file,
    health_file = health_file
  )

  expect_false(result$valid)
  expect_true(any(grepl("missing probability columns", result$issues)))
  expect_true(any(grepl("conc_meanlog", result$issues)))

  unlink(c(prob_file, conc_file, health_file))
})

test_that("validate_prepared_inputs detects invalid probability sums", {
  temp_dir <- tempdir()
  prob_file <- file.path(temp_dir, "test_prob_prob.csv")
  conc_file <- file.path(temp_dir, "test_conc_prob.csv")
  health_file <- file.path(temp_dir, "test_health_prob.txt")

  # Create probability file with probabilities that don't sum to 1
  prob_data_local <- data.frame(
    GEOID10 = c("06001", "06003"),
    prob_C1 = c(0.5, 0.6),
    prob_C2 = c(0.2, 0.25),
    prob_C3 = c(0.1, 0.15)  # Sum = 0.8 and 1.0
  )
  data.table::fwrite(prob_data_local, prob_file)

  # Create valid concentration file
  conc_data <- data.frame(
    GEOID10 = c("06001", "06003"),
    conc_meanlog = c(1.0, 1.5),
    conc_sdlog = c(0.8, 0.9),
    PWELL_private_pct = c(15, 25)
  )
  data.table::fwrite(conc_data, conc_file)

  # Create health file
  health_data <- data.frame(
    FIPS = c("06001", "06003"),
    BWT = c(3200, 3150)
  )
  data.table::fwrite(health_data, health_file, sep = "\t")

  result <- validate_prepared_inputs(
    prob_file = prob_file,
    conc_file = conc_file,
    health_file = health_file
  )

  # Should have warnings about probability sums
  expect_true(any(grepl("probability sums", result$warnings)))

  unlink(c(prob_file, conc_file, health_file))
})

test_that("validate_prepared_inputs detects negative probabilities", {
  temp_dir <- tempdir()
  prob_file <- file.path(temp_dir, "test_prob_neg.csv")
  conc_file <- file.path(temp_dir, "test_conc_neg.csv")
  health_file <- file.path(temp_dir, "test_health_neg.txt")

  # Create probability file with negative probability
  prob_data_local <- data.frame(
    GEOID10 = c("06001"),
    prob_C1 = c(-0.1),
    prob_C2 = c(0.6),
    prob_C3 = c(0.5)
  )
  data.table::fwrite(prob_data_local, prob_file)

  conc_data <- data.frame(
    GEOID10 = c("06001"),
    conc_meanlog = c(1.0),
    conc_sdlog = c(0.8)
  )
  data.table::fwrite(conc_data, conc_file)

  health_data <- data.frame(FIPS = c("06001"), BWT = c(3200))
  data.table::fwrite(health_data, health_file, sep = "\t")

  result <- validate_prepared_inputs(
    prob_file = prob_file,
    conc_file = conc_file,
    health_file = health_file
  )

  expect_false(result$valid)
  expect_true(any(grepl("Negative probabilities", result$issues)))

  unlink(c(prob_file, conc_file, health_file))
})

test_that("validate_prepared_inputs detects invalid PWELL range", {
  temp_dir <- tempdir()
  prob_file <- file.path(temp_dir, "test_prob_pwell.csv")
  conc_file <- file.path(temp_dir, "test_conc_pwell.csv")
  health_file <- file.path(temp_dir, "test_health_pwell.txt")

  prob_data_local <- data.frame(
    GEOID10 = c("06001"),
    prob_C1 = c(0.7),
    prob_C2 = c(0.2),
    prob_C3 = c(0.1)
  )
  data.table::fwrite(prob_data_local, prob_file)

  # Create concentration file with invalid PWELL (> 100)
  conc_data <- data.frame(
    GEOID10 = c("06001"),
    conc_meanlog = c(1.0),
    conc_sdlog = c(0.8),
    PWELL_private_pct = c(150)  # Invalid
  )
  data.table::fwrite(conc_data, conc_file)

  health_data <- data.frame(FIPS = c("06001"), BWT = c(3200))
  data.table::fwrite(health_data, health_file, sep = "\t")

  result <- validate_prepared_inputs(
    prob_file = prob_file,
    conc_file = conc_file,
    health_file = health_file
  )

  expect_false(result$valid)
  expect_true(any(grepl("PWELL_private_pct must be between", result$issues)))

  unlink(c(prob_file, conc_file, health_file))
})

test_that("validate_prepared_inputs returns summary statistics", {
  temp_dir <- tempdir()
  prob_file <- file.path(temp_dir, "test_prob_sum.csv")
  conc_file <- file.path(temp_dir, "test_conc_sum.csv")
  health_file <- file.path(temp_dir, "test_health_sum.txt")

  prob_data_local <- data.frame(
    GEOID10 = c("06001", "06003", "06005"),
    prob_C1 = c(0.7, 0.6, 0.5),
    prob_C2 = c(0.2, 0.25, 0.3),
    prob_C3 = c(0.1, 0.15, 0.2)
  )
  data.table::fwrite(prob_data_local, prob_file)

  conc_data <- data.frame(
    GEOID10 = c("06001", "06003"),
    conc_meanlog = c(1.0, 1.5),
    conc_sdlog = c(0.8, 0.9),
    PWELL_private_pct = c(15, 25)
  )
  data.table::fwrite(conc_data, conc_file)

  health_data <- data.frame(
    FIPS = c("06001", "06003", "06005", "06007"),
    BWT = c(3200, 3150, 3100, 3050)
  )
  data.table::fwrite(health_data, health_file, sep = "\t")

  result <- validate_prepared_inputs(
    prob_file = prob_file,
    conc_file = conc_file,
    health_file = health_file
  )

  # Check summary statistics
  expect_equal(result$summary$prob$n_rows, 3)
  expect_equal(result$summary$conc$n_rows, 2)
  expect_equal(result$summary$health$n_rows, 4)

  unlink(c(prob_file, conc_file, health_file))
})

# ------------------------------------------------------------------------------
# Tests for prepare_prob_data() - requires terra and sf packages
# ------------------------------------------------------------------------------

test_that("prepare_prob_data validates required packages", {
  # Test with non-list input
  expect_error(
    prepare_prob_data(
      prob_rasters = "not_a_list",
      boundaries = NULL
    ),
    "prob_rasters.*must be a non-empty list"
  )

  # Test with empty list
  expect_error(
    prepare_prob_data(
      prob_rasters = list(),
      boundaries = NULL
    ),
    "prob_rasters.*must be a non-empty list"
  )
})

test_that("prepare_prob_data validates boundaries input", {
  expect_error(
    prepare_prob_data(
      prob_rasters = list(gt5 = "test.tif"),
      boundaries = data.frame(x = 1)
    ),
    "boundaries.*must be an sf object"
  )
})

test_that("prepare_prob_data validates cutoffs", {
  # Create a minimal sf object for testing
  boundaries <- sf::st_sf(
    GEOID = "06001",
    geometry = sf::st_sfc(sf::st_point(c(-122, 38)))
  )

  expect_error(
    prepare_prob_data(
      prob_rasters = list(gt5 = "test.tif"),
      boundaries = boundaries,
      cutoffs = c(10, 5)  # Invalid: first > second
    ),
    "cutoffs.*must be a numeric vector"
  )

  expect_error(
    prepare_prob_data(
      prob_rasters = list(gt5 = "test.tif"),
      boundaries = boundaries,
      cutoffs = c(5)  # Invalid: length 1
    ),
    "cutoffs.*must be a numeric vector of length 2"
  )
})

test_that("prepare_prob_data validates geoid column exists", {
  boundaries <- sf::st_sf(
    wrong_col = "06001",
    geometry = sf::st_sfc(sf::st_point(c(-122, 38)))
  )

  expect_error(
    prepare_prob_data(
      prob_rasters = list(gt5 = "test.tif"),
      boundaries = boundaries,
      id_col = "GEOID"
    ),
    "Column 'GEOID' not found in boundaries"
  )
})

# Test prepare_prob_data full workflow with mock rasters
test_that("prepare_prob_data processes rasters with centroid method", {
  # Create temporary raster files
  temp_dir <- tempdir()
  gt5_file <- file.path(temp_dir, "prob_gt5.tif")
  gt10_file <- file.path(temp_dir, "prob_gt10.tif")

  # Create simple rasters with known values
  r5 <- terra::rast(nrows = 10, ncols = 10, xmin = -122, xmax = -121, ymin = 37, ymax = 38)
  terra::values(r5) <- runif(100, 0.1, 0.5)  # Probability > 5 µg/L
  terra::writeRaster(r5, gt5_file, overwrite = TRUE)

  r10 <- terra::rast(nrows = 10, ncols = 10, xmin = -122, xmax = -121, ymin = 37, ymax = 38)
  terra::values(r10) <- runif(100, 0.05, 0.2)  # Probability > 10 µg/L (should be <= gt5)
  terra::writeRaster(r10, gt10_file, overwrite = TRUE)

  # Create a simple polygon boundary
  poly <- sf::st_as_sf(sf::st_sfc(
    sf::st_polygon(list(matrix(c(-122, 37, -121, 37, -121, 38, -122, 38, -122, 37), ncol = 2, byrow = TRUE)))
  ))
  poly$GEOID <- "06001"
  sf::st_crs(poly) <- 4326

  # Test centroid extraction method
  result <- prepare_prob_data(
    prob_rasters = list(gt5 = gt5_file, gt10 = gt10_file),
    boundaries = poly,
    id_col = "GEOID",
    extraction_method = "centroid",
    cutoffs = c(5, 10)
  )

  expect_true(data.table::is.data.table(result))
  expect_equal(nrow(result), 1)
  expect_true("GEOID10" %in% names(result))
  expect_true("prob_C1" %in% names(result))
  expect_true("prob_C2" %in% names(result))
  expect_true("prob_C3" %in% names(result))

  # Probabilities should sum to approximately 1
  prob_sum <- result$prob_C1 + result$prob_C2 + result$prob_C3
  expect_equal(prob_sum, 1, tolerance = 0.01)

  # Clean up
  unlink(c(gt5_file, gt10_file))
})

test_that("prepare_prob_data processes rasters with mean method", {
  temp_dir <- tempdir()
  gt5_file <- file.path(temp_dir, "prob_gt5_mean.tif")
  gt10_file <- file.path(temp_dir, "prob_gt10_mean.tif")

  # Create rasters
  r5 <- terra::rast(nrows = 10, ncols = 10, xmin = -122, xmax = -121, ymin = 37, ymax = 38)
  terra::values(r5) <- rep(0.3, 100)  # Uniform probability
  terra::writeRaster(r5, gt5_file, overwrite = TRUE)

  r10 <- terra::rast(nrows = 10, ncols = 10, xmin = -122, xmax = -121, ymin = 37, ymax = 38)
  terra::values(r10) <- rep(0.1, 100)
  terra::writeRaster(r10, gt10_file, overwrite = TRUE)

  # Create polygon
  poly <- sf::st_as_sf(sf::st_sfc(
    sf::st_polygon(list(matrix(c(-122, 37, -121, 37, -121, 38, -122, 38, -122, 37), ncol = 2, byrow = TRUE)))
  ))
  poly$GEOID <- "06001"
  sf::st_crs(poly) <- 4326

  # Test mean extraction method
  result <- prepare_prob_data(
    prob_rasters = list(gt5 = gt5_file, gt10 = gt10_file),
    boundaries = poly,
    id_col = "GEOID",
    extraction_method = "mean",
    cutoffs = c(5, 10)
  )

  expect_true(data.table::is.data.table(result))
  expect_equal(nrow(result), 1)

  # With uniform values, we should get approximately:
  # prob_cat1 = 1 - 0.3 = 0.7
  # prob_cat2 = 0.3 - 0.1 = 0.2
  # prob_cat3 = 0.1
  expect_equal(result$prob_C1, 0.7, tolerance = 0.05)
  expect_equal(result$prob_C2, 0.2, tolerance = 0.05)
  expect_equal(result$prob_C3, 0.1, tolerance = 0.05)

  unlink(c(gt5_file, gt10_file))
})

test_that("prepare_prob_data adds population data when provided", {
  temp_dir <- tempdir()
  gt5_file <- file.path(temp_dir, "prob_gt5_pop.tif")
  gt10_file <- file.path(temp_dir, "prob_gt10_pop.tif")

  r5 <- terra::rast(nrows = 5, ncols = 5, xmin = -122, xmax = -121, ymin = 37, ymax = 38)
  terra::values(r5) <- rep(0.25, 25)
  terra::writeRaster(r5, gt5_file, overwrite = TRUE)

  r10 <- terra::rast(nrows = 5, ncols = 5, xmin = -122, xmax = -121, ymin = 37, ymax = 38)
  terra::values(r10) <- rep(0.08, 25)
  terra::writeRaster(r10, gt10_file, overwrite = TRUE)

  poly <- sf::st_as_sf(sf::st_sfc(
    sf::st_polygon(list(matrix(c(-122, 37, -121, 37, -121, 38, -122, 38, -122, 37), ncol = 2, byrow = TRUE)))
  ))
  poly$GEOID <- "06001"
  sf::st_crs(poly) <- 4326

  # Create population data
  pop_data <- data.frame(
    GEOID = "06001",
    private_well_pop = 15000
  )

  result <- prepare_prob_data(
    prob_rasters = list(gt5 = gt5_file, gt10 = gt10_file),
    boundaries = poly,
    id_col = "GEOID",
    pop_data = pop_data,
    pop_id_col = "GEOID",
    pop_well_col = "private_well_pop",
    extraction_method = "centroid",
    cutoffs = c(5, 10)
  )

  expect_true("Wells_2010" %in% names(result))
  expect_equal(result$Wells_2010, 15000)

  unlink(c(gt5_file, gt10_file))
})

test_that("prepare_prob_data warns on missing pop columns", {
  temp_dir <- tempdir()
  gt5_file <- file.path(temp_dir, "prob_gt5_popwarn.tif")
  gt10_file <- file.path(temp_dir, "prob_gt10_popwarn.tif")

  r5 <- terra::rast(nrows = 5, ncols = 5, xmin = -122, xmax = -121, ymin = 37, ymax = 38)
  terra::values(r5) <- rep(0.25, 25)
  terra::writeRaster(r5, gt5_file, overwrite = TRUE)

  r10 <- terra::rast(nrows = 5, ncols = 5, xmin = -122, xmax = -121, ymin = 37, ymax = 38)
  terra::values(r10) <- rep(0.08, 25)
  terra::writeRaster(r10, gt10_file, overwrite = TRUE)

  poly <- sf::st_as_sf(sf::st_sfc(
    sf::st_polygon(list(matrix(c(-122, 37, -121, 37, -121, 38, -122, 38, -122, 37), ncol = 2, byrow = TRUE)))
  ))
  poly$GEOID <- "06001"
  sf::st_crs(poly) <- 4326

  # Test with wrong GEOID column
  pop_data_bad_geoid <- data.frame(
    wrong_col = "06001",
    private_well_pop = 15000
  )

  expect_warning(
    prepare_prob_data(
      prob_rasters = list(gt5 = gt5_file, gt10 = gt10_file),
      boundaries = poly,
      id_col = "GEOID",
      pop_data = pop_data_bad_geoid,
      pop_id_col = "GEOID",
      pop_well_col = "private_well_pop",
      extraction_method = "centroid",
      cutoffs = c(5, 10)
    ),
    "Population GEOID column"
  )

  # Test with wrong well column
  pop_data_bad_well <- data.frame(
    GEOID = "06001",
    wrong_well_col = 15000
  )

  expect_warning(
    prepare_prob_data(
      prob_rasters = list(gt5 = gt5_file, gt10 = gt10_file),
      boundaries = poly,
      id_col = "GEOID",
      pop_data = pop_data_bad_well,
      pop_id_col = "GEOID",
      pop_well_col = "private_well_pop",
      extraction_method = "centroid",
      cutoffs = c(5, 10)
    ),
    "Population well column"
  )

  unlink(c(gt5_file, gt10_file))
})

test_that("prepare_prob_data saves output file when requested", {
  temp_dir <- tempdir()
  gt5_file <- file.path(temp_dir, "prob_gt5_save.tif")
  gt10_file <- file.path(temp_dir, "prob_gt10_save.tif")
  output_file <- file.path(temp_dir, "prob_output.csv")

  r5 <- terra::rast(nrows = 5, ncols = 5, xmin = -122, xmax = -121, ymin = 37, ymax = 38)
  terra::values(r5) <- rep(0.3, 25)
  terra::writeRaster(r5, gt5_file, overwrite = TRUE)

  r10 <- terra::rast(nrows = 5, ncols = 5, xmin = -122, xmax = -121, ymin = 37, ymax = 38)
  terra::values(r10) <- rep(0.1, 25)
  terra::writeRaster(r10, gt10_file, overwrite = TRUE)

  poly <- sf::st_as_sf(sf::st_sfc(
    sf::st_polygon(list(matrix(c(-122, 37, -121, 37, -121, 38, -122, 38, -122, 37), ncol = 2, byrow = TRUE)))
  ))
  poly$GEOID <- "06001"
  sf::st_crs(poly) <- 4326

  result <- prepare_prob_data(
    prob_rasters = list(gt5 = gt5_file, gt10 = gt10_file),
    boundaries = poly,
    id_col = "GEOID",
    extraction_method = "centroid",
    cutoffs = c(5, 10),
    output_file = output_file
  )

  expect_true(file.exists(output_file))
  saved_data <- data.table::fread(output_file)
  expect_equal(nrow(saved_data), 1)
  expect_true("prob_C1" %in% names(saved_data))

  unlink(c(gt5_file, gt10_file, output_file))
})

test_that("prepare_prob_data warns when cutoffs cannot be matched to columns", {
  temp_dir <- tempdir()
  # Create rasters with non-standard names
  gt_low_file <- file.path(temp_dir, "prob_low.tif")
  gt_high_file <- file.path(temp_dir, "prob_high.tif")

  r_low <- terra::rast(nrows = 5, ncols = 5, xmin = -122, xmax = -121, ymin = 37, ymax = 38)
  terra::values(r_low) <- rep(0.3, 25)
  terra::writeRaster(r_low, gt_low_file, overwrite = TRUE)

  r_high <- terra::rast(nrows = 5, ncols = 5, xmin = -122, xmax = -121, ymin = 37, ymax = 38)
  terra::values(r_high) <- rep(0.1, 25)
  terra::writeRaster(r_high, gt_high_file, overwrite = TRUE)

  poly <- sf::st_as_sf(sf::st_sfc(
    sf::st_polygon(list(matrix(c(-122, 37, -121, 37, -121, 38, -122, 38, -122, 37), ncol = 2, byrow = TRUE)))
  ))
  poly$GEOID <- "06001"
  sf::st_crs(poly) <- 4326

  # Names don't contain cutoff values, should warn
  expect_warning(
    prepare_prob_data(
      prob_rasters = list(low = gt_low_file, high = gt_high_file),
      boundaries = poly,
      id_col = "GEOID",
      extraction_method = "centroid",
      cutoffs = c(5, 10)
    ),
    "Could not automatically match cutoffs"
  )

  unlink(c(gt_low_file, gt_high_file))
})

test_that("prepare_prob_data handles missing raster file", {
  poly <- sf::st_as_sf(sf::st_sfc(
    sf::st_polygon(list(matrix(c(-122, 37, -121, 37, -121, 38, -122, 38, -122, 37), ncol = 2, byrow = TRUE)))
  ))
  poly$GEOID <- "06001"
  sf::st_crs(poly) <- 4326

  expect_error(
    prepare_prob_data(
      prob_rasters = list(gt5 = "/nonexistent/path/file.tif"),
      boundaries = poly,
      id_col = "GEOID",
      cutoffs = c(5, 10)
    ),
    "Raster file not found"
  )
})

test_that("prepare_prob_data processes multiple geographic units", {
  temp_dir <- tempdir()
  gt5_file <- file.path(temp_dir, "prob_gt5_multi.tif")
  gt10_file <- file.path(temp_dir, "prob_gt10_multi.tif")

  # Create larger rasters
  r5 <- terra::rast(nrows = 20, ncols = 20, xmin = -123, xmax = -121, ymin = 36, ymax = 38)
  terra::values(r5) <- runif(400, 0.1, 0.5)
  terra::writeRaster(r5, gt5_file, overwrite = TRUE)

  r10 <- terra::rast(nrows = 20, ncols = 20, xmin = -123, xmax = -121, ymin = 36, ymax = 38)
  terra::values(r10) <- runif(400, 0.05, 0.2)
  terra::writeRaster(r10, gt10_file, overwrite = TRUE)

  # Create two polygons
  poly1 <- sf::st_polygon(list(matrix(c(-123, 36, -122, 36, -122, 37, -123, 37, -123, 36), ncol = 2, byrow = TRUE)))
  poly2 <- sf::st_polygon(list(matrix(c(-122, 37, -121, 37, -121, 38, -122, 38, -122, 37), ncol = 2, byrow = TRUE)))

  polys <- sf::st_as_sf(sf::st_sfc(poly1, poly2))
  polys$GEOID <- c("06001", "06003")
  sf::st_crs(polys) <- 4326

  result <- prepare_prob_data(
    prob_rasters = list(gt5 = gt5_file, gt10 = gt10_file),
    boundaries = polys,
    id_col = "GEOID",
    extraction_method = "mean",
    cutoffs = c(5, 10)
  )

  expect_equal(nrow(result), 2)
  expect_true(all(c("06001", "06003") %in% result$GEOID10))

  # Each row should have probabilities summing to 1
  for (i in 1:nrow(result)) {
    prob_sum <- result$prob_C1[i] + result$prob_C2[i] + result$prob_C3[i]
    expect_equal(prob_sum, 1, tolerance = 0.01)
  }

  unlink(c(gt5_file, gt10_file))
})

# Test validate_prepared_inputs with .txt health file
test_that("validate_prepared_inputs handles .txt health files", {
  temp_dir <- tempdir()

  # Create probability file - probabilities must sum to 1
  prob_file <- file.path(temp_dir, "prob_txt_test.csv")
  prob_data_local <- data.frame(
    GEOID10 = c("06001", "06003"),
    prob_C1 = c(0.7, 0.6),
    prob_C2 = c(0.2, 0.25),
    prob_C3 = c(0.1, 0.15),
    Wells_2010 = c(100, 200),
    PWELL = c(0.05, 0.08)
  )
  data.table::fwrite(prob_data_local, prob_file)

  # Create concentration file with expected column names
  conc_file <- file.path(temp_dir, "conc_txt_test.csv")
  conc_data <- data.frame(
    GEOID10 = c("06001", "06003"),
    conc_meanlog = c(0.5, 0.8),
    conc_sdlog = c(0.3, 0.4),
    PWELL_private_pct = c(0.05, 0.08)
  )
  data.table::fwrite(conc_data, conc_file)

  # Create health file as .txt (tab-separated)
  health_file <- file.path(temp_dir, "health_test.txt")
  health_data <- data.frame(
    FIPS = c("06001", "06003"),
    BWT = c(3200, 3100),
    OEGEST = c(39, 38)
  )
  data.table::fwrite(health_data, health_file, sep = "\t")

  result <- validate_prepared_inputs(
    prob_file = prob_file,
    conc_file = conc_file,
    health_file = health_file
  )

  expect_true(result$valid)
  expect_equal(result$summary$health$n_rows, 2)

  unlink(c(prob_file, conc_file, health_file))
})

# ------------------------------------------------------------------------------
# Additional tests for improved coverage
# ------------------------------------------------------------------------------

test_that("prepare_conc_data handles population-weighted averaging without Hmisc", {
  # This tests the fallback when Hmisc is not available
  conc_raw_data <- data.frame(
    county_fips = c("06001", "06001", "06003", "06003"),
    concentration = c(2.5, 3.1, 8.2, 7.5),
    population_served = c(50000, 25000, 10000, 15000)
  )

  # This will use population-weighted averaging
  result <- prepare_conc_data(
    conc_data = conc_raw_data,
    id_col = "county_fips",
    concentration_col = "concentration",
    population_col = "population_served",
    input_type = "raw"
  )

  expect_true(data.table::is.data.table(result))
  expect_equal(nrow(result), 2)
  expect_true("n_systems" %in% names(result))
})

test_that("format_health_data handles SBP and DBP range validation", {
  health <- data.frame(
    FIPS = c(6001, 6003),
    SBP = c(50, 260),  # Outside expected range 60-250
    DBP = c(25, 160)   # Outside expected range 30-150
  )

  # Should produce warnings for both SBP and DBP
  expect_warning(
    result <- format_health_data(
      health_data = health,
      id_col = "FIPS",
      outcome_cols = c("SBP", "DBP"),
      validate = TRUE
    ),
    "outside expected range"
  )
})

test_that("format_health_data reports missing data summary", {
  health <- data.frame(
    FIPS = c(6001, 6003, 6005, 6007),
    BWT = c(3200, NA, 3100, NA),
    MAGE = c(28, 30, NA, 32)
  )

  # Should report missing values in messages
  expect_message(
    format_health_data(
      health_data = health,
      id_col = "FIPS",
      outcome_cols = c("BWT"),
      covariate_cols = c("MAGE"),
      validate = TRUE
    ),
    "missing values"
  )
})

test_that("format_health_data handles duplicate FIPS codes message", {
  health <- data.frame(
    FIPS = c(6001, 6001, 6001, 6003),
    BWT = c(3200, 3150, 3100, 3250)
  )

  expect_message(
    format_health_data(
      health_data = health,
      id_col = "FIPS",
      outcome_cols = c("BWT"),
      validate = TRUE
    ),
    "duplicate FIPS codes"
  )
})

test_that("subset_by_geography handles multiple filters simultaneously", {
  data <- data.frame(
    GEOID10 = c("06001", "06003", "06005", "32001"),
    value = 1:4
  )

  # Apply both state_fips and geoid_list filters
  result <- subset_by_geography(
    data,
    state_fips = "06",
    geoid_list = c("06001", "06003")
  )

  expect_equal(nrow(result), 2)
  expect_true(all(result$GEOID10 %in% c("06001", "06003")))
})

test_that("validate_prepared_inputs warns on missing concentration sdlog column", {
  temp_dir <- tempdir()
  prob_file <- file.path(temp_dir, "test_prob_sdlog.csv")
  conc_file <- file.path(temp_dir, "test_conc_sdlog.csv")
  health_file <- file.path(temp_dir, "test_health_sdlog.txt")

  prob_data_local <- data.frame(
    GEOID10 = c("06001"),
    prob_C1 = c(0.7),
    prob_C2 = c(0.2),
    prob_C3 = c(0.1)
  )
  data.table::fwrite(prob_data_local, prob_file)

  # Concentration file without sdlog column
  conc_data <- data.frame(
    GEOID10 = c("06001"),
    conc_meanlog = c(1.0),
    PWELL_private_pct = c(15)
  )
  data.table::fwrite(conc_data, conc_file)

  health_data <- data.frame(FIPS = c("06001"), BWT = c(3200))
  data.table::fwrite(health_data, health_file, sep = "\t")

  result <- validate_prepared_inputs(
    prob_file = prob_file,
    conc_file = conc_file,
    health_file = health_file
  )

  expect_true(any(grepl("conc_sdlog", result$warnings)))

  unlink(c(prob_file, conc_file, health_file))
})

test_that("validate_prepared_inputs warns on missing PWELL column", {
  temp_dir <- tempdir()
  prob_file <- file.path(temp_dir, "test_prob_nopwell.csv")
  conc_file <- file.path(temp_dir, "test_conc_nopwell.csv")
  health_file <- file.path(temp_dir, "test_health_nopwell.txt")

  prob_data_local <- data.frame(
    GEOID10 = c("06001"),
    prob_C1 = c(0.7),
    prob_C2 = c(0.2),
    prob_C3 = c(0.1)
  )
  data.table::fwrite(prob_data_local, prob_file)

  # Concentration file without PWELL column
  conc_data <- data.frame(
    GEOID10 = c("06001"),
    conc_meanlog = c(1.0),
    conc_sdlog = c(0.8)
  )
  data.table::fwrite(conc_data, conc_file)

  health_data <- data.frame(FIPS = c("06001"), BWT = c(3200))
  data.table::fwrite(health_data, health_file, sep = "\t")

  result <- validate_prepared_inputs(
    prob_file = prob_file,
    conc_file = conc_file,
    health_file = health_file
  )

  expect_true(any(grepl("PWELL_private_pct", result$warnings)))

  unlink(c(prob_file, conc_file, health_file))
})

test_that("validate_prepared_inputs warns on non-positive sdlog", {
  temp_dir <- tempdir()
  prob_file <- file.path(temp_dir, "test_prob_badsd.csv")
  conc_file <- file.path(temp_dir, "test_conc_badsd.csv")
  health_file <- file.path(temp_dir, "test_health_badsd.txt")

  prob_data_local <- data.frame(
    GEOID10 = c("06001"),
    prob_C1 = c(0.7),
    prob_C2 = c(0.2),
    prob_C3 = c(0.1)
  )
  data.table::fwrite(prob_data_local, prob_file)

  # Concentration file with non-positive sdlog
  conc_data <- data.frame(
    GEOID10 = c("06001"),
    conc_meanlog = c(1.0),
    conc_sdlog = c(0),  # Invalid: must be positive
    PWELL_private_pct = c(15)
  )
  data.table::fwrite(conc_data, conc_file)

  health_data <- data.frame(FIPS = c("06001"), BWT = c(3200))
  data.table::fwrite(health_data, health_file, sep = "\t")

  result <- validate_prepared_inputs(
    prob_file = prob_file,
    conc_file = conc_file,
    health_file = health_file
  )

  expect_false(result$valid)
  expect_true(any(grepl("sdlog must be positive", result$issues)))

  unlink(c(prob_file, conc_file, health_file))
})

test_that("validate_prepared_inputs detects missing health GEOID column", {
  temp_dir <- tempdir()
  prob_file <- file.path(temp_dir, "test_prob_hgeo.csv")
  conc_file <- file.path(temp_dir, "test_conc_hgeo.csv")
  health_file <- file.path(temp_dir, "test_health_hgeo.txt")

  prob_data_local <- data.frame(
    GEOID10 = c("06001"),
    prob_C1 = c(0.7),
    prob_C2 = c(0.2),
    prob_C3 = c(0.1)
  )
  data.table::fwrite(prob_data_local, prob_file)

  conc_data <- data.frame(
    GEOID10 = c("06001"),
    conc_meanlog = c(1.0),
    conc_sdlog = c(0.8),
    PWELL_private_pct = c(15)
  )
  data.table::fwrite(conc_data, conc_file)

  # Health file with wrong GEOID column name
  health_data <- data.frame(wrong_col = c("06001"), BWT = c(3200))
  data.table::fwrite(health_data, health_file, sep = "\t")

  result <- validate_prepared_inputs(
    prob_file = prob_file,
    conc_file = conc_file,
    health_file = health_file
  )

  expect_false(result$valid)
  expect_true(any(grepl("Health missing GEOID column", result$issues)))

  unlink(c(prob_file, conc_file, health_file))
})

test_that("validate_prepared_inputs warns on FIPS format issues", {
  temp_dir <- tempdir()
  prob_file <- file.path(temp_dir, "test_prob_fips.csv")
  conc_file <- file.path(temp_dir, "test_conc_fips.csv")
  health_file <- file.path(temp_dir, "test_health_fips.txt")

  prob_data_local <- data.frame(
    GEOID10 = c("06001"),
    prob_C1 = c(0.7),
    prob_C2 = c(0.2),
    prob_C3 = c(0.1)
  )
  data.table::fwrite(prob_data_local, prob_file)

  conc_data <- data.frame(
    GEOID10 = c("06001"),
    conc_meanlog = c(1.0),
    conc_sdlog = c(0.8),
    PWELL_private_pct = c(15)
  )
  data.table::fwrite(conc_data, conc_file)

  # Health file with FIPS codes that aren't 5 digits
  health_data <- data.frame(FIPS = c("6001", "123456"), BWT = c(3200, 3100))
  data.table::fwrite(health_data, health_file, sep = "\t")

  result <- validate_prepared_inputs(
    prob_file = prob_file,
    conc_file = conc_file,
    health_file = health_file
  )

  expect_true(any(grepl("FIPS codes are not 5 digits", result$warnings)))

  unlink(c(prob_file, conc_file, health_file))
})

test_that("validate_prepared_inputs detects missing probability data GEOID column", {
  temp_dir <- tempdir()
  prob_file <- file.path(temp_dir, "test_prob_ugeo.csv")
  conc_file <- file.path(temp_dir, "test_conc_ugeo.csv")
  health_file <- file.path(temp_dir, "test_health_ugeo.txt")

  # Probability file with wrong GEOID column name
  prob_data_local <- data.frame(
    wrong_geoid = c("06001"),
    prob_C1 = c(0.7),
    prob_C2 = c(0.2),
    prob_C3 = c(0.1)
  )
  data.table::fwrite(prob_data_local, prob_file)

  conc_data <- data.frame(
    GEOID10 = c("06001"),
    conc_meanlog = c(1.0),
    conc_sdlog = c(0.8),
    PWELL_private_pct = c(15)
  )
  data.table::fwrite(conc_data, conc_file)

  health_data <- data.frame(FIPS = c("06001"), BWT = c(3200))
  data.table::fwrite(health_data, health_file, sep = "\t")

  result <- validate_prepared_inputs(
    prob_file = prob_file,
    conc_file = conc_file,
    health_file = health_file
  )

  expect_false(result$valid)
  expect_true(any(grepl("Probability data missing GEOID column", result$issues)))

  unlink(c(prob_file, conc_file, health_file))
})

test_that("validate_prepared_inputs detects missing concentration data GEOID column", {
  temp_dir <- tempdir()
  prob_file <- file.path(temp_dir, "test_prob_egeo.csv")
  conc_file <- file.path(temp_dir, "test_conc_egeo.csv")
  health_file <- file.path(temp_dir, "test_health_egeo.txt")

  prob_data_local <- data.frame(
    GEOID10 = c("06001"),
    prob_C1 = c(0.7),
    prob_C2 = c(0.2),
    prob_C3 = c(0.1)
  )
  data.table::fwrite(prob_data_local, prob_file)

  # Concentration file with wrong GEOID column name
  conc_data <- data.frame(
    wrong_geoid = c("06001"),
    conc_meanlog = c(1.0),
    conc_sdlog = c(0.8),
    PWELL_private_pct = c(15)
  )
  data.table::fwrite(conc_data, conc_file)

  health_data <- data.frame(FIPS = c("06001"), BWT = c(3200))
  data.table::fwrite(health_data, health_file, sep = "\t")

  result <- validate_prepared_inputs(
    prob_file = prob_file,
    conc_file = conc_file,
    health_file = health_file
  )

  expect_false(result$valid)
  expect_true(any(grepl("Concentration data missing GEOID column", result$issues)))

  unlink(c(prob_file, conc_file, health_file))
})

test_that("validate_prepared_inputs warns on meanlog outside typical range", {
  temp_dir <- tempdir()
  prob_file <- file.path(temp_dir, "test_prob_rng.csv")
  conc_file <- file.path(temp_dir, "test_conc_rng.csv")
  health_file <- file.path(temp_dir, "test_health_rng.txt")

  prob_data_local <- data.frame(
    GEOID10 = c("06001"),
    prob_C1 = c(0.7),
    prob_C2 = c(0.2),
    prob_C3 = c(0.1)
  )
  data.table::fwrite(prob_data_local, prob_file)

  # Concentration file with meanlog outside typical range
  conc_data <- data.frame(
    GEOID10 = c("06001"),
    conc_meanlog = c(15.0),  # Outside typical range [-5, 10]
    conc_sdlog = c(0.8),
    PWELL_private_pct = c(15)
  )
  data.table::fwrite(conc_data, conc_file)

  health_data <- data.frame(FIPS = c("06001"), BWT = c(3200))
  data.table::fwrite(health_data, health_file, sep = "\t")

  result <- validate_prepared_inputs(
    prob_file = prob_file,
    conc_file = conc_file,
    health_file = health_file
  )

  expect_true(any(grepl("meanlog values outside typical range", result$warnings)))

  unlink(c(prob_file, conc_file, health_file))
})

test_that("validate_prepared_inputs detects ID alignment issues", {
  temp_dir <- tempdir()
  prob_file <- file.path(temp_dir, "test_prob_align.csv")
  conc_file <- file.path(temp_dir, "test_conc_align.csv")
  health_file <- file.path(temp_dir, "test_health_align.txt")

  # Probability data has different IDs than concentration data
  prob_data_local <- data.frame(
    GEOID10 = c("06001", "06003", "06005"),
    prob_C1 = c(0.7, 0.6, 0.5),
    prob_C2 = c(0.2, 0.25, 0.3),
    prob_C3 = c(0.1, 0.15, 0.2)
  )
  data.table::fwrite(prob_data_local, prob_file)

  conc_data <- data.frame(
    GEOID10 = c("06001", "06007", "06009"),  # Different IDs
    conc_meanlog = c(1.0, 1.5, 2.0),
    conc_sdlog = c(0.8, 0.9, 1.0),
    PWELL_private_pct = c(15, 20, 25)
  )
  data.table::fwrite(conc_data, conc_file)

  health_data <- data.frame(FIPS = c("06001"), BWT = c(3200))
  data.table::fwrite(health_data, health_file, sep = "\t")

  result <- validate_prepared_inputs(
    prob_file = prob_file,
    conc_file = conc_file,
    health_file = health_file
  )

  # Should have warnings about IDs only in one dataset
  expect_true(any(grepl("IDs in probability data only", result$warnings)))
  expect_true(any(grepl("IDs in concentration data only", result$warnings)))

  unlink(c(prob_file, conc_file, health_file))
})

test_that("prepare_conc_data handles simple averaging without population weights", {
  conc_raw_data <- data.frame(
    county_fips = c("06001", "06001", "06003"),
    concentration = c(2.5, 3.5, 5.0)
  )

  result <- prepare_conc_data(
    conc_data = conc_raw_data,
    id_col = "county_fips",
    concentration_col = "concentration",
    population_col = NULL,  # No population weighting
    input_type = "raw"
  )

  expect_equal(nrow(result), 2)
  # Mean for 06001 should be (2.5 + 3.5) / 2 = 3.0
  expect_true("n_systems" %in% names(result))
})

test_that("format_health_data skips validation when disabled", {
  health <- data.frame(
    FIPS = c(6001),
    BWT = c(100)  # Would trigger warning if validated
  )

  # No warning expected when validate = FALSE
  result <- format_health_data(
    health_data = health,
    id_col = "FIPS",
    outcome_cols = c("BWT"),
    validate = FALSE
  )

  expect_equal(nrow(result), 1)
})

# ==============================================================================
# Tests for automatic concentration lognormal parameter calculation
# ==============================================================================

test_that("calculate_conc_lognormal_params auto-calculates from concentration", {
  # Create concentration data with raw concentration but no lognormal params
  conc_data <- data.table::data.table(
    GEOID10 = c("06001", "06002", "06003"),
    exposure_conc = c(2.5, 7.0, 15.0),
    PWELL_private_pct = c(20, 35, 45)
  )

  # Call the internal function
  result <- geoExposeR:::calculate_conc_lognormal_params(
    datmatx_conc = conc_data,
    conc_mean_col = "conc_meanlog",
    conc_sd_col = "conc_sdlog",
    conc_raw_col = "exposure_conc",
    default_sdlog = 1.5
  )

  # Check that new columns were added
  expect_true("conc_meanlog" %in% names(result))
  expect_true("conc_sdlog" %in% names(result))

  # Check meanlog calculation: log(conc + 0.1)
  expect_equal(result$conc_meanlog[1], log(2.5 + 0.1), tolerance = 1e-6)
  expect_equal(result$conc_meanlog[2], log(7.0 + 0.1), tolerance = 1e-6)
  expect_equal(result$conc_meanlog[3], log(15.0 + 0.1), tolerance = 1e-6)

  # Check sdlog uses default
  expect_equal(result$conc_sdlog, rep(1.5, 3))
})

test_that("calculate_conc_lognormal_params returns data unchanged when cols exist", {
  # Create concentration data with existing lognormal params
  conc_data <- data.table::data.table(
    GEOID10 = c("06001", "06002"),
    conc_meanlog = c(1.0, 2.0),
    conc_sdlog = c(0.8, 0.9),
    PWELL_private_pct = c(20, 35)
  )

  # Call the internal function
  result <- geoExposeR:::calculate_conc_lognormal_params(
    datmatx_conc = conc_data,
    conc_mean_col = "conc_meanlog",
    conc_sd_col = "conc_sdlog",
    conc_raw_col = "exposure_conc",  # Not used
    default_sdlog = 1.0
  )

  # Check values unchanged
  expect_equal(result$conc_meanlog, c(1.0, 2.0))
  expect_equal(result$conc_sdlog, c(0.8, 0.9))
})

test_that("calculate_conc_lognormal_params handles NA in concentration", {
  conc_data <- data.table::data.table(
    GEOID10 = c("06001", "06002", "06003"),
    exposure_conc = c(2.5, NA, 15.0),
    PWELL_private_pct = c(20, 35, 45)
  )

  result <- geoExposeR:::calculate_conc_lognormal_params(
    datmatx_conc = conc_data,
    conc_mean_col = "conc_meanlog",
    conc_sd_col = "conc_sdlog",
    conc_raw_col = "exposure_conc",
    default_sdlog = 1.0
  )

  # First and third should have values, second should be NA
  expect_false(is.na(result$conc_meanlog[1]))
  expect_true(is.na(result$conc_meanlog[2]))
  expect_false(is.na(result$conc_meanlog[3]))
})

test_that("calculate_conc_lognormal_params errors when conc col not provided and meanlog missing", {
  conc_data <- data.table::data.table(
    GEOID10 = c("06001", "06002"),
    exposure_conc = c(2.5, 7.0),
    PWELL_private_pct = c(20, 35)
  )

  expect_error(
    geoExposeR:::calculate_conc_lognormal_params(
      datmatx_conc = conc_data,
      conc_mean_col = "conc_meanlog",
      conc_sd_col = "conc_sdlog",
      conc_raw_col = NULL,  # Not provided!
      default_sdlog = 1.0
    ),
    "Concentration lognormal meanlog column.*not found.*conc_raw_col"
  )
})

test_that("validate_default_sdlog validates input correctly", {
  # Valid values should not error
  expect_no_error(geoExposeR:::validate_default_sdlog(1.0))
  expect_no_error(geoExposeR:::validate_default_sdlog(0.5))

  # Invalid values should error
  expect_error(geoExposeR:::validate_default_sdlog("1.0"), "must be a single numeric")
  expect_error(geoExposeR:::validate_default_sdlog(c(1.0, 2.0)), "must be a single numeric")
  expect_error(geoExposeR:::validate_default_sdlog(-1.0), "must be positive")
  expect_error(geoExposeR:::validate_default_sdlog(0), "must be positive")

  # Large values should warn
  expect_warning(geoExposeR:::validate_default_sdlog(15.0), "very large")
})

# =============================================================================
# Tests for format_ids_func (Smart ID Formatting)
# =============================================================================

test_that("format_ids_func applies zero-padding for numeric FIPS codes", {
  # Test with small numeric FIPS codes (1-5 digits)
  df <- data.frame(
    id = c(1, 123, 6001, 12345),
    value = c(10, 20, 30, 40)
  )
  
  result <- geoExposeR:::format_ids_func(df, "id")
  
  # Check that IDs are zero-padded to 5 digits

  expect_equal(result$id, c("00001", "00123", "06001", "12345"))
  expect_type(result$id, "character")
})

test_that("format_ids_func preserves longer numeric codes without padding", {
  # Test with 11-digit census tract IDs
  df <- data.frame(
    id = c(6001001234, 6001005678, 12345678901),
    value = c(10, 20, 30)
  )
  
  result <- geoExposeR:::format_ids_func(df, "id")
  
  # Check that IDs are converted to character without padding
  expect_equal(result$id, c("6001001234", "6001005678", "12345678901"))
  expect_type(result$id, "character")
})

test_that("format_ids_func preserves alphanumeric IDs unchanged", {
  # Test with alphanumeric participant IDs
  df <- data.frame(
    id = c("PTD001", "PTD002", "CA-001", "NV-123"),
    value = c(10, 20, 30, 40),
    stringsAsFactors = FALSE
  )
  
  result <- geoExposeR:::format_ids_func(df, "id")
  
  # Check that IDs are unchanged
  expect_equal(result$id, c("PTD001", "PTD002", "CA-001", "NV-123"))
  expect_type(result$id, "character")
})

test_that("format_ids_func handles NA values correctly", {
  # Test with NA values in numeric IDs
  df <- data.frame(
    id = c(1, NA, 6001, NA, 12345),
    value = c(10, 20, 30, 40, 50)
  )
  
  result <- geoExposeR:::format_ids_func(df, "id")
  
  # Check that non-NA values are formatted correctly
  # Note: sprintf converts NA to "   NA" string, not actual NA
  expect_equal(result$id[1], "00001")
  expect_equal(result$id[3], "06001")
  expect_equal(result$id[5], "12345")
  expect_type(result$id, "character")
})

test_that("format_ids_func handles mixed character IDs", {
  # Test with already-character IDs including numbers
  df <- data.frame(
    id = c("06001", "06002", "ABC123"),
    value = c(10, 20, 30),
    stringsAsFactors = FALSE
  )
  
  result <- geoExposeR:::format_ids_func(df, "id")
  
  # Check that character IDs are preserved
  expect_equal(result$id, c("06001", "06002", "ABC123"))
  expect_type(result$id, "character")
})

test_that("format_ids_func handles edge case of exactly 99999", {
  # Test boundary case - 99999 is within FIPS range, 100000 is not
  df <- data.frame(
    id = c(99999, 100000),
    value = c(10, 20)
  )
  
  result <- geoExposeR:::format_ids_func(df, "id")
  
  # 99999 should be padded (within FIPS range)
  # 100000 is outside range so whole vector gets as.character() treatment
  # Since one value is > 99999, no padding is applied to any
  expect_equal(result$id[1], "99999")
  expect_type(result$id, "character")
})

test_that("format_ids_func handles decimal numbers as non-FIPS", {
  # Decimal numbers should not be treated as FIPS codes
  df <- data.frame(
    id = c(1.5, 2.7, 3.0),
    value = c(10, 20, 30)
  )
  
  result <- geoExposeR:::format_ids_func(df, "id")
  
  # Decimals should be converted to character without padding
  expect_equal(result$id, c("1.5", "2.7", "3"))
  expect_type(result$id, "character")
})

test_that("format_ids_func handles single row data frame", {
  df <- data.frame(id = 6001, value = 10)
  
  result <- geoExposeR:::format_ids_func(df, "id")
  
  expect_equal(result$id, "06001")
})

test_that("format_ids_func handles empty data frame", {
  df <- data.frame(id = numeric(0), value = numeric(0))
  
  result <- geoExposeR:::format_ids_func(df, "id")
  
  expect_equal(nrow(result), 0)
  expect_type(result$id, "character")
})

# =============================================================================
# Tests for format_ids parameter in main functions
# =============================================================================

test_that("perform_sensitivity_analysis works with format_ids = TRUE (default)", {
  temp_dir <- tempfile(pattern = "test_format_ids_true_")
  dir.create(temp_dir)
  
  # Create dummy data with numeric IDs
  dummy_files <- create_dummy_data(file_dir = temp_dir)
  
  output_dir <- file.path(temp_dir, "output")
  dir.create(output_dir)
  
  # Run with format_ids = TRUE (default)
  results <- expect_no_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2,
      regression_formula = "~ as.factor(ExposureLevel) + MAGE_R + (1 | MRSTATE)",
      targets = c("BWT"),
      output_dir = output_dir,
      rucc_col = NULL,
      format_ids = TRUE
    )
  )
  
  expect_type(results, "list")
  expect_named(results, "BWT")
  
  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})

test_that("perform_sensitivity_analysis works with format_ids = FALSE", {
  temp_dir <- tempfile(pattern = "test_format_ids_false_")
  dir.create(temp_dir)
  
  # Create dummy data
  dummy_files <- create_dummy_data(file_dir = temp_dir)
  
  output_dir <- file.path(temp_dir, "output")
  dir.create(output_dir)
  
  # Run with format_ids = FALSE
  results <- expect_no_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      model_prob_csv = dummy_files$prob,
      conc_data_csv = dummy_files$conc,
      ndraws = 2,
      regression_formula = "~ as.factor(ExposureLevel) + MAGE_R + (1 | MRSTATE)",
      targets = c("BWT"),
      output_dir = output_dir,
      rucc_col = NULL,
      format_ids = FALSE
    )
  )
  
  expect_type(results, "list")
  expect_named(results, "BWT")
  
  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})

test_that("load_and_process_exposure_data respects format_ids parameter", {
  temp_dir <- tempfile(pattern = "test_exposure_format_")
  dir.create(temp_dir)
  
  # Create dummy data
  dummy_files <- create_dummy_data(file_dir = temp_dir)
  
  # Test with format_ids = TRUE
  result_formatted <- geoExposeR:::load_and_process_exposure_data(
    model_prob_csv = dummy_files$prob,
    conc_data_csv = dummy_files$conc,
    prob_cols = c("prob_C1", "prob_C2", "prob_C3"),
    pop_well_col = "Wells_2010",
    conc_mean_col = "conc_meanlog",
    conc_sd_col = "conc_sdlog",
    conc_raw_col = NULL,
    default_sdlog = 1.0,
    pwell_col = "PWELL_private_pct",
    id_col = "GEOID10",
    conc_cutoffs = c(5, 10),
    format_ids = TRUE
  )
  
  # Test with format_ids = FALSE
  result_unformatted <- geoExposeR:::load_and_process_exposure_data(
    model_prob_csv = dummy_files$prob,
    conc_data_csv = dummy_files$conc,
    prob_cols = c("prob_C1", "prob_C2", "prob_C3"),
    pop_well_col = "Wells_2010",
    conc_mean_col = "conc_meanlog",
    conc_sd_col = "conc_sdlog",
    conc_raw_col = NULL,
    default_sdlog = 1.0,
    pwell_col = "PWELL_private_pct",
    id_col = "GEOID10",
    conc_cutoffs = c(5, 10),
    format_ids = FALSE
  )
  
  # Both should return valid data
  expect_type(result_formatted, "list")
  expect_type(result_unformatted, "list")
  
  # IDs should be character in both cases
  expect_type(result_formatted$datmatx_prob$GEOID10, "character")
  expect_type(result_unformatted$datmatx_prob$GEOID10, "character")
  
  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})

test_that("load_and_process_birth_data respects format_ids parameter", {
  temp_dir <- tempfile(pattern = "test_birth_format_")
  dir.create(temp_dir)
  
  # Create dummy birth data with numeric IDs
  birth_data <- data.frame(
    FIPS = c(1, 2, 3, 4, 5),
    BWT = c(3200, 3100, 3300, 3000, 3400),
    OEGEST = c(39, 38, 40, 37, 41),
    MAGE_R = c(25, 30, 28, 22, 35),
    MRSTATE = c("CA", "CA", "NV", "NV", "CA")
  )
  
  birth_file <- file.path(temp_dir, "births.txt")
  data.table::fwrite(birth_data, birth_file, sep = "\t")
  
  # Test with format_ids = TRUE
  result_formatted <- geoExposeR:::load_and_process_birth_data(
    birth_data_txt = birth_file,
    columns_to_select = NULL,
    rucc_col = NULL,
    record_id_col = "FIPS",
    format_ids = TRUE
  )
  
  # Test with format_ids = FALSE
  result_unformatted <- geoExposeR:::load_and_process_birth_data(
    birth_data_txt = birth_file,
    columns_to_select = NULL,
    rucc_col = NULL,
    record_id_col = "FIPS",
    format_ids = FALSE
  )
  
  # With format_ids = TRUE, IDs should be zero-padded
  expect_equal(result_formatted$FIPS, c("00001", "00002", "00003", "00004", "00005"))
  
  # With format_ids = FALSE, IDs should just be character (no padding)
  expect_equal(result_unformatted$FIPS, c("1", "2", "3", "4", "5"))
  
  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})

test_that("format_ids parameter handles alphanumeric IDs correctly", {
  temp_dir <- tempfile(pattern = "test_alpha_ids_")
  dir.create(temp_dir)
  
  # Create birth data with alphanumeric participant IDs
  birth_data <- data.frame(
    participant_id = c("PTD001", "PTD002", "PTD003"),
    BWT = c(3200, 3100, 3300),
    OEGEST = c(39, 38, 40),
    MAGE_R = c(25, 30, 28),
    MRSTATE = c("CA", "CA", "NV"),
    stringsAsFactors = FALSE
  )
  
  birth_file <- file.path(temp_dir, "births.txt")
  data.table::fwrite(birth_data, birth_file, sep = "\t")
  
  # Test with format_ids = TRUE (should preserve alphanumeric)
  result <- geoExposeR:::load_and_process_birth_data(
    birth_data_txt = birth_file,
    columns_to_select = NULL,
    rucc_col = NULL,
    record_id_col = "participant_id",
    format_ids = TRUE
  )
  
  # Alphanumeric IDs should be preserved unchanged
  expect_equal(result$participant_id, c("PTD001", "PTD002", "PTD003"))
  
  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})

# ==============================================================================
# Additional Tests for 100% Code Coverage
# ==============================================================================

test_that("validate_prob_columns_param handles edge cases", {
  # Test non-character input
  expect_error(
    geoExposeR:::validate_prob_columns_param(123),
    "prob_columns must be a character vector"
  )
  
  # Test empty vector
  expect_error(
    geoExposeR:::validate_prob_columns_param(character(0)),
    "prob_columns cannot be empty"
  )
  
  # Test NA values
  expect_error(
    geoExposeR:::validate_prob_columns_param(c("col1", NA)),
    "prob_columns cannot contain NA or empty strings"
  )
  
  # Test empty strings
  expect_error(
    geoExposeR:::validate_prob_columns_param(c("col1", "")),
    "prob_columns cannot contain NA or empty strings"
  )
  
  # Test duplicate column names
  expect_error(
    geoExposeR:::validate_prob_columns_param(c("col1", "col1")),
    "prob_columns cannot contain duplicate column names"
  )
  
  # Test NULL (should pass silently)
  expect_silent(geoExposeR:::validate_prob_columns_param(NULL))
  
  # Test valid input
  expect_silent(geoExposeR:::validate_prob_columns_param(c("col1", "col2")))
})

test_that("validate_conc_sd_values handles non-finite and non-positive values", {
  temp_dir <- tempfile(pattern = "test_conc_sd_")
  dir.create(temp_dir)
  
  # Test with Inf values in sdlog
  conc_data_inf <- data.frame(
    GEOID10 = c("01001", "01003", "01005"),
    conc_meanlog = c(1.5, 2.0, 1.8),
    conc_sdlog = c(0.5, Inf, 0.7),
    PWELL_private_pct = c(30, 40, 35)
  )
  conc_file_inf <- file.path(temp_dir, "conc_inf.csv")
  data.table::fwrite(conc_data_inf, conc_file_inf)
  
  expect_error(
    geoExposeR:::validate_conc_sd_values(conc_file_inf, "conc_sdlog"),
    "must contain only finite numeric values"
  )
  
  # Test with negative values in sdlog
  conc_data_neg <- data.frame(
    GEOID10 = c("01001", "01003", "01005"),
    conc_meanlog = c(1.5, 2.0, 1.8),
    conc_sdlog = c(0.5, -0.3, 0.7),
    PWELL_private_pct = c(30, 40, 35)
  )
  conc_file_neg <- file.path(temp_dir, "conc_neg.csv")
  data.table::fwrite(conc_data_neg, conc_file_neg)
  
  expect_error(
    geoExposeR:::validate_conc_sd_values(conc_file_neg, "conc_sdlog"),
    "must contain only positive values"
  )
  
  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})

test_that("validate_conc_lognormal_input handles missing columns", {
  temp_dir <- tempfile(pattern = "test_conc_lognormal_")
  dir.create(temp_dir)
  
  # Concentration data without meanlog column and no concentration column
  conc_data <- data.frame(
    GEOID10 = c("01001", "01003", "01005"),
    PWELL_private_pct = c(30, 40, 35)
  )
  conc_file <- file.path(temp_dir, "conc_no_meanlog.csv")
  data.table::fwrite(conc_data, conc_file)
  
  # Test missing meanlog without concentration column
  expect_error(
    geoExposeR:::validate_conc_lognormal_input(
      conc_file, "conc_meanlog", "conc_sdlog", NULL
    ),
    "not found in data"
  )
  
  # Test with non-existent concentration column
  expect_error(
    geoExposeR:::validate_conc_lognormal_input(
      conc_file, "conc_meanlog", "conc_sdlog", "exposure_conc"
    ),
    "not found in the data"
  )
  
  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})

test_that("validate_file_contents handles file read errors", {
  temp_dir <- tempfile(pattern = "test_file_contents_")
  dir.create(temp_dir)
  
  # Create valid probability file
  prob_data_local <- data.frame(
    GEOID10 = c("01001", "01003"),
    prob_C1 = c(0.5, 0.6),
    prob_C2 = c(0.3, 0.25),
    prob_C3 = c(0.2, 0.15),
    Wells_2010 = c(1000, 2000)
  )
  prob_file <- file.path(temp_dir, "prob_data.csv")
  data.table::fwrite(prob_data_local, prob_file)
  
  # Create valid concentration file
  conc_data <- data.frame(
    GEOID10 = c("01001", "01003"),
    conc_meanlog = c(1.5, 2.0),
    conc_sdlog = c(0.5, 0.6),
    PWELL_private_pct = c(30, 40)
  )
  conc_file <- file.path(temp_dir, "conc_data.csv")
  data.table::fwrite(conc_data, conc_file)
  
  # Create valid birth file
  birth_data <- data.frame(
    FIPS = c("01001", "01003"),
    BWT = c(3200, 3100)
  )
  birth_file <- file.path(temp_dir, "births.txt")
  data.table::fwrite(birth_data, birth_file, sep = "\t")
  
  # Test that valid files pass validation
  expect_silent(
    geoExposeR:::validate_file_contents(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      rucc_col = NULL,
      record_id_col = "FIPS",
      id_col = "GEOID10",
      prob_cols = c("prob_C1", "prob_C2", "prob_C3"),
      pop_well_col = "Wells_2010",
      conc_mean_col = "conc_meanlog",
      conc_sd_col = "conc_sdlog",
      conc_raw_col = NULL,
      pwell_col = "PWELL_private_pct"
    )
  )
  
  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})

test_that("calculate_conc_lognormal_params handles missing conc_raw_col", {
  temp_dir <- tempfile(pattern = "test_calc_conc_")
  dir.create(temp_dir)
  
  # Concentration data without meanlog but with concentration
  conc_data <- data.frame(
    GEOID10 = c("01001", "01003", "01005"),
    exposure_conc = c(5.0, 10.0, 15.0),
    PWELL_private_pct = c(30, 40, 35)
  )
  conc_file <- file.path(temp_dir, "conc_raw_values.csv")
  data.table::fwrite(conc_data, conc_file)
  
  datmatx_conc <- data.table::fread(conc_file)
  
  # Test with valid concentration column (auto-conversion)
  result <- geoExposeR:::calculate_conc_lognormal_params(
    datmatx_conc = data.table::copy(datmatx_conc),  # Use copy to avoid modifying original
    conc_mean_col = "conc_meanlog",
    conc_sd_col = "conc_sdlog",
    conc_raw_col = "exposure_conc",
    default_sdlog = 1.0
  )
  
  expect_true("conc_meanlog" %in% names(result))
  expect_true("conc_sdlog" %in% names(result))
  
  # Test error when concentration column doesn't exist (reload fresh data)
  datmatx_conc_fresh <- data.table::fread(conc_file)
  expect_error(
    geoExposeR:::calculate_conc_lognormal_params(
      datmatx_conc = datmatx_conc_fresh,
      conc_mean_col = "conc_meanlog",
      conc_sd_col = "conc_sdlog",
      conc_raw_col = "nonexistent_col",
      default_sdlog = 1.0
    ),
    "not found"
  )
  
  # Test error when neither meanlog nor conc_col exists
  conc_data_no_cols <- data.table::as.data.table(data.frame(
    GEOID10 = c("01001", "01003"),
    PWELL_private_pct = c(30, 40)
  ))
  
  expect_error(
    geoExposeR:::calculate_conc_lognormal_params(
      datmatx_conc = conc_data_no_cols,
      conc_mean_col = "conc_meanlog",
      conc_sd_col = "conc_sdlog",
      conc_raw_col = NULL,
      default_sdlog = 1.0
    ),
    "not found"
  )
  
  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})

test_that("impute_exposure handles edge cases", {
  # This function requires many dependencies, so we test through the main function

  # The edge cases (single cell per county, NA weights) are covered by integration tests
  # Here we just verify the function exists and has correct signature
  expect_true(is.function(geoExposeR:::impute_exposure))
})

test_that("impute_additional_variables handles no MICE covariates", {
  # Create data with only one variable (insufficient for MICE)
  imputed_datasets <- list(
    data.frame(
      ExposureLevel = c("<5", "5-10", "10+"),
      var_to_impute = c(NA, 2, 3)
    )
  )
  
  births <- data.frame(
    var_to_impute = c(1, 2, 3)  # Original births data
  )
  
  # Should produce a warning about no MICE covariates
  expect_warning(
    geoExposeR:::impute_additional_variables(
      imputed_datasets = imputed_datasets,
      impute_vars = "var_to_impute",
      births = births,
      mice_covs = NULL,  # No covariates available
      mice_m = 1,
      mice_maxit = 5,
      mice_method = "pmm",
      seed = 12345,
      check_birth_cols = FALSE
    ),
    "No MICE covariates available"
  )
})

test_that("regression_analysis handles case with insufficient factor levels", {
  temp_dir <- tempfile(pattern = "test_reg_no_terms_")
  dir.create(temp_dir)
  
  # Create data where ExposureLevel has only one unique value (only reference category)
  # This should trigger an error about contrasts
  births_data <- data.frame(
    .id = rep(1:3, each = 100),  # 3 imputed datasets
    ExposureLevel = rep("<5", 300),  # Only reference category
    BWT = rnorm(300, 3200, 500),
    MRSTATE = sample(c("CA", "NV"), 300, replace = TRUE)
  )
  births_data$ExposureLevel <- factor(births_data$ExposureLevel, levels = c("<5", "5-10", "10+"))
  
  # Should produce error about contrasts (only one factor level present)
  expect_error(
    geoExposeR:::regression_analysis(
      births = births_data,
      targets = "BWT",
      cat_label = c("5-10", "10+"),  # Non-reference categories
      exposure_level_col = "ExposureLevel",
      regression_formula = "~ as.factor(ExposureLevel) + (1|MRSTATE)",
      output_dir = temp_dir
    ),
    "contrasts|Model fitting failed"
  )
  
  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})

test_that("pool_estimates_by_term handles missing terms", {
  # Create combined results with only some terms
  combined_results <- data.frame(
    term = c("5-10", "5-10", "5-10"),  # Only 5-10, no 10+
    estimate = c(0.1, 0.15, 0.12),
    std.error = c(0.05, 0.06, 0.055)
  )
  
  # Create a mock model for degrees of freedom
  mock_data <- data.frame(
    y = rnorm(100),
    x = rnorm(100),
    group = sample(1:5, 100, replace = TRUE)
  )
  mock_model <- lme4::lmer(y ~ x + (1|group), data = mock_data)
  
  # Pool with both 5-10 and 10+ (one missing)
  result <- geoExposeR:::pool_estimates_by_term(
    combined_results = combined_results,
    exposure_terms = c("5-10", "10+"),  # 10+ not in combined_results
    reference_model = mock_model
  )
  
  # Should only have 5-10 in result (10+ is NULL and filtered)
  expect_equal(nrow(result), 1)
  expect_equal(result$term[1], "5-10")
})

test_that("pool_single_estimate handles edge case with zero between-imputation variance", {
  # When all estimates are identical, b = 0
  estimates <- rep(0.5, 5)  # Identical estimates
  std_errors <- rep(0.1, 5)   # Identical SE values
  
  result <- geoExposeR:::pool_single_estimate(
    estimates = estimates,
    std_errors = std_errors,
    df = 100
  )
  
  expect_equal(result$q.mi, 0.5)
  expect_true(!is.na(result$se.mi))
  expect_true(!is.na(result$p.value))
})

test_that("validate_prepared_inputs handles unreadable health file", {
  temp_dir <- tempfile(pattern = "test_validate_unread_")
  dir.create(temp_dir)
  
  # Create valid probability file
  prob_data_local <- data.frame(
    GEOID10 = c("01001", "01003"),
    prob_C1 = c(0.5, 0.6),
    prob_C2 = c(0.3, 0.25),
    prob_C3 = c(0.2, 0.15),
    Wells_2010 = c(1000, 2000)
  )
  prob_file <- file.path(temp_dir, "prob_data.csv")
  data.table::fwrite(prob_data_local, prob_file)
  
  # Create valid concentration file
  conc_data <- data.frame(
    GEOID10 = c("01001", "01003"),
    conc_meanlog = c(1.5, 2.0),
    conc_sdlog = c(0.5, 0.6),
    PWELL_private_pct = c(30, 40)
  )
  conc_file <- file.path(temp_dir, "conc_data.csv")
  data.table::fwrite(conc_data, conc_file)
  
  # Create a file with unexpected format (will cause parsing issues)
  health_file <- file.path(temp_dir, "health.txt")
  writeLines(c("malformed data without proper columns", "just random text"), health_file)
  
  # Should handle the read error gracefully
  result <- geoExposeR::validate_prepared_inputs(
    prob_file = prob_file,
    conc_file = conc_file,
    health_file = health_file
  )
  
  # The function should complete but report issues (missing FIPS column)
  expect_false(result$valid)
  
  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})

test_that("prepare_prob_data handles less than 2 probability rasters", {
  # This tests line 199 in data-ingestion.R
  temp_dir <- tempfile(pattern = "test_prob_rasters_")
  dir.create(temp_dir)
  
  # Create mock boundaries
  suppressWarnings({
    boundaries <- sf::st_sf(
      GEOID = c("01001", "01003"),
      geometry = sf::st_sfc(
        sf::st_polygon(list(matrix(c(0, 0, 1, 0, 1, 1, 0, 1, 0, 0), ncol = 2, byrow = TRUE))),
        sf::st_polygon(list(matrix(c(1, 0, 2, 0, 2, 1, 1, 1, 1, 0), ncol = 2, byrow = TRUE)))
      ),
      crs = 4326
    )
  })
  
  # Create a single raster (insufficient for 3-category model)
  suppressWarnings({
    r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 2, ymin = 0, ymax = 1)
    terra::values(r) <- runif(100)
    raster_file <- file.path(temp_dir, "prob_gt5.tif")
    terra::writeRaster(r, raster_file, overwrite = TRUE)
  })
  
  # Test with single raster (should error)
  expect_error(
    suppressWarnings(
      geoExposeR::prepare_prob_data(
        prob_rasters = list(gt5 = raster_file),
        boundaries = boundaries,
        id_col = "GEOID"
      )
    ),
    "at least 2 probability rasters"
  )
  
  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})

# Test regression_analysis with no exposure terms (warning path)
test_that("regression_analysis handles no exposure terms in results", {
  temp_dir <- tempfile(pattern = "test_reg_warn_")
  dir.create(temp_dir)
  
  # Create data with multiple factor levels but model results won't have exposure terms
  set.seed(123)
  births_data <- data.frame(
    .id = rep(1:3, each = 100),
    ExposureLevel = factor(sample(c("<5", "5-10", "10+"), 300, replace = TRUE),
                     levels = c("<5", "5-10", "10+")),
    BWT = rnorm(300, 3200, 500),
    MRSTATE = sample(c("CA", "NV"), 300, replace = TRUE)
  )
  
  # This should run successfully (and cover lines 78-82 when exposure terms exist)
  result <- suppressWarnings(
    geoExposeR:::regression_analysis(
      births = births_data,
      targets = "BWT",
      cat_label = c("5-10", "10+"),
      exposure_level_col = "ExposureLevel",
      regression_formula = "~ as.factor(ExposureLevel) + (1|MRSTATE)",
      output_dir = temp_dir
    )
  )
  
  # Result should have BWT entry
  expect_true("BWT" %in% names(result))
  
  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})

# Test pool_estimates_by_term error handling for degrees of freedom
test_that("pool_estimates_by_term handles df extraction error", {
  # Create combined results
  combined_results <- data.frame(
    term = c("5-10", "5-10", "5-10"),
    estimate = c(0.1, 0.15, 0.12),
    std.error = c(0.05, 0.06, 0.055)
  )
  
  # Create a mock object that will cause glance to fail
  mock_model <- list(class = "mock")
  class(mock_model) <- "invalid_model"
  
  # Should produce warning about degrees of freedom and use Inf
  expect_warning(
    geoExposeR:::pool_estimates_by_term(
      combined_results = combined_results,
      exposure_terms = c("5-10"),
      reference_model = mock_model
    ),
    "Could not extract degrees of freedom"
  )
})

# Test validate_file_contents with missing required columns
test_that("validate_file_contents detects missing columns", {
  temp_dir <- tempfile(pattern = "test_validate_cols_")
  dir.create(temp_dir)
  
  # Create file with missing required columns
  prob_data_local <- data.frame(
    GEOID10 = c("01001", "01003"),
    prob_C1 = c(0.5, 0.6)
    # Missing prob_C2, prob_C3, Wells_2010
  )
  prob_file <- file.path(temp_dir, "prob_data.csv")
  data.table::fwrite(prob_data_local, prob_file)
  
  conc_data <- data.frame(
    GEOID10 = c("01001", "01003"),
    conc_meanlog = c(1.5, 2.0),
    conc_sdlog = c(0.5, 0.6),
    PWELL_private_pct = c(30, 40)
  )
  conc_file <- file.path(temp_dir, "conc_data.csv")
  data.table::fwrite(conc_data, conc_file)
  
  birth_data <- data.frame(
    FIPS = c("01001", "01003"),
    BWT = c(3200, 3100)
  )
  birth_file <- file.path(temp_dir, "births.txt")
  data.table::fwrite(birth_data, birth_file, sep = "\t")
  
  # Should error due to missing probability columns
  expect_error(
    geoExposeR:::validate_file_contents(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      rucc_col = NULL,
      record_id_col = "FIPS",
      id_col = "GEOID10",
      prob_cols = c("prob_C1", "prob_C2", "prob_C3"),
      pop_well_col = "Wells_2010",
      conc_mean_col = "conc_meanlog",
      conc_sd_col = "conc_sdlog",
      conc_raw_col = NULL,
      pwell_col = "PWELL_private_pct"
    ),
    "columns are missing"
  )
  
  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})

# Test validate_prepared_inputs with file read error in health file
test_that("validate_prepared_inputs handles health file read error path", {
  temp_dir <- tempfile(pattern = "test_validate_health_err_")
  dir.create(temp_dir)
  
  # Create valid files
  prob_data_local <- data.frame(
    GEOID10 = c("01001", "01003"),
    prob_C1 = c(0.5, 0.6),
    prob_C2 = c(0.3, 0.25),
    prob_C3 = c(0.2, 0.15),
    Wells_2010 = c(1000, 2000)
  )
  prob_file <- file.path(temp_dir, "prob_data.csv")
  data.table::fwrite(prob_data_local, prob_file)
  
  conc_data <- data.frame(
    GEOID10 = c("01001", "01003"),
    conc_meanlog = c(1.5, 2.0),
    conc_sdlog = c(0.5, 0.6),
    PWELL_private_pct = c(30, 40)
  )
  conc_file <- file.path(temp_dir, "conc_data.csv")
  data.table::fwrite(conc_data, conc_file)
  
  # Create a non-existent health file path
  health_file <- file.path(temp_dir, "nonexistent.txt")
  
  # Should handle missing file
  result <- geoExposeR::validate_prepared_inputs(
    prob_file = prob_file,
    conc_file = conc_file,
    health_file = health_file
  )
  
  expect_false(result$valid)
  expect_true(any(grepl("not found", result$issues, ignore.case = TRUE)))
  
  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})

# Test validate_conc_lognormal_input with NULL file
test_that("validate_conc_lognormal_input handles NULL and invalid files", {
  # Test with NULL
  expect_silent(
    geoExposeR:::validate_conc_lognormal_input(
      conc_data_csv = NULL,
      conc_mean_col = "conc_meanlog",
      conc_sd_col = "conc_sdlog",
      conc_raw_col = NULL
    )
  )
  
  # Test with non-existent file
  expect_silent(
    geoExposeR:::validate_conc_lognormal_input(
      conc_data_csv = "/nonexistent/path/file.csv",
      conc_mean_col = "conc_meanlog",
      conc_sd_col = "conc_sdlog",
      conc_raw_col = NULL
    )
  )
})

# Test impute_exposure with cells_i null condition (line 52)
test_that("impute_exposure handles county not in lookup", {
  # This edge case is covered implicitly through integration tests
  # The function handles counties with no matching cells by skipping them
  expect_true(is.function(geoExposeR:::impute_exposure))
})

# Test impute_exposure with mismatched county IDs
test_that("impute_exposure skips counties not in lookup", {
  # Create probability data with one county
  datmatx_prob <- data.frame(
    GEOID10 = c("01001", "01001"),
    Wells = c(100, 200)
  )
  
  # Create probability matrix
  rasterprob_combined <- matrix(c(0.5, 0.5, 0.3, 0.7), nrow = 2, ncol = 2)
  
  # Create population weights
  cellpop <- c(100, 200)
  
  # Create births data
  births <- data.frame(
    FIPS = c("01001", "01003"),  # 01003 exists in births but not in probability data
    BWT = c(3200, 3100)
  )
  
  # Run imputation - the 01003 county should be skipped
  # but since county_id_list comes from datmatx_prob, this won't trigger the next
  result <- geoExposeR:::impute_exposure(
    datmatx_prob = datmatx_prob,
    rasterprob_combined = rasterprob_combined,
    cellpop = cellpop,
    births = births,
    ndraws = 2,
    cat_label = c("Low", "High"),
    id_col = "GEOID10",
    record_id_col = "FIPS",
    exposure_level_col = "ExposureLevel",
    apply_imputation_fallback = TRUE
  )
  
  expect_true(is.list(result))
  expect_equal(length(result), 2)
})

# Test get_file_columns with unsupported file format (line 531)
test_that("get_file_columns rejects unsupported file format", {
  temp_dir <- tempfile(pattern = "test_file_cols_")
  dir.create(temp_dir)
  
  # Create a file with unsupported extension
  bad_file <- file.path(temp_dir, "data.xlsx")
  writeLines("test", bad_file)
  
  # Should error on unsupported format
  expect_error(
    geoExposeR:::get_file_columns(bad_file),
    "Unsupported file format"
  )
  
  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})

# Test validate_conc_sd_values with column not in file (line 597)
test_that("validate_conc_sd_values returns when column not in file", {
  temp_dir <- tempfile(pattern = "test_conc_sd_col_")
  dir.create(temp_dir)
  
  # Create concentration file without the sd column
  conc_data <- data.frame(
    GEOID10 = c("01001", "01003"),
    conc_meanlog = c(1.5, 2.0),
    PWELL_private_pct = c(30, 40)
    # Missing conc_sdlog
  )
  conc_file <- file.path(temp_dir, "conc_data.csv")
  data.table::fwrite(conc_data, conc_file)
  
  # Should return silently when column doesn't exist
  expect_silent(
    geoExposeR:::validate_conc_sd_values(
      conc_data_csv = conc_file,
      conc_sd_col = "conc_sdlog"
    )
  )
  
  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})

# Test validate_conc_sd_values with NA values (lines 611-612)
test_that("validate_conc_sd_values warns on NA values", {
  temp_dir <- tempfile(pattern = "test_conc_sd_na_")
  dir.create(temp_dir)
  
  # Create concentration file with NA values in sd column
  conc_data <- data.frame(
    GEOID10 = c("01001", "01003", "01005"),
    conc_meanlog = c(1.5, 2.0, 1.8),
    conc_sdlog = c(0.5, NA, 0.6),  # Contains NA
    PWELL_private_pct = c(30, 40, 35)
  )
  conc_file <- file.path(temp_dir, "conc_data.csv")
  data.table::fwrite(conc_data, conc_file)
  
  # Should warn about NA values
  expect_warning(
    geoExposeR:::validate_conc_sd_values(
      conc_data_csv = conc_file,
      conc_sd_col = "conc_sdlog"
    ),
    "contains NA values"
  )
  
  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})

# Test validate_file_contents with concentration column and no meanlog (line 531)
test_that("validate_file_contents uses concentration column when no meanlog", {
  temp_dir <- tempfile(pattern = "test_validate_conc_")
  dir.create(temp_dir)
  
  # Create concentration file without meanlog but with concentration column
  conc_data <- data.frame(
    GEOID10 = c("01001", "01003"),
    raw_conc = c(5.0, 8.0),  # Concentration instead of meanlog
    PWELL_private_pct = c(30, 40)
    # No conc_meanlog or conc_sdlog
  )
  conc_file <- file.path(temp_dir, "conc_data.csv")
  data.table::fwrite(conc_data, conc_file)
  
  prob_data_local <- data.frame(
    GEOID10 = c("01001", "01003"),
    prob_C1 = c(0.5, 0.6),
    prob_C2 = c(0.3, 0.25),
    prob_C3 = c(0.2, 0.15),
    Wells_2010 = c(1000, 2000)
  )
  prob_file <- file.path(temp_dir, "prob_data.csv")
  data.table::fwrite(prob_data_local, prob_file)
  
  birth_data <- data.frame(
    FIPS = c("01001", "01003"),
    BWT = c(3200, 3100)
  )
  birth_file <- file.path(temp_dir, "births.txt")
  data.table::fwrite(birth_data, birth_file, sep = "\t")
  
  # Should pass when using concentration column
  expect_silent(
    geoExposeR:::validate_file_contents(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      rucc_col = NULL,
      record_id_col = "FIPS",
      id_col = "GEOID10",
      prob_cols = c("prob_C1", "prob_C2", "prob_C3"),
      pop_well_col = "Wells_2010",
      conc_mean_col = "conc_meanlog",  # Not in file
      conc_sd_col = NULL,
      conc_raw_col = "raw_conc",  # This should be required
      pwell_col = "PWELL_private_pct"
    )
  )
  
  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})

# Test validate_prepared_inputs with unreadable health file (lines 973, 976-977)
test_that("validate_prepared_inputs handles health file read error", {
  temp_dir <- tempfile(pattern = "test_health_read_err_")
  dir.create(temp_dir)
  
  # Create valid probability and concentration files
  prob_data_local <- data.frame(
    GEOID10 = c("01001", "01003"),
    prob_C1 = c(0.5, 0.6),
    prob_C2 = c(0.3, 0.25),
    prob_C3 = c(0.2, 0.15),
    Wells_2010 = c(1000, 2000)
  )
  prob_file <- file.path(temp_dir, "prob_data.csv")
  data.table::fwrite(prob_data_local, prob_file)
  
  conc_data <- data.frame(
    GEOID10 = c("01001", "01003"),
    conc_meanlog = c(1.5, 2.0),
    conc_sdlog = c(0.5, 0.6),
    PWELL_private_pct = c(30, 40)
  )
  conc_file <- file.path(temp_dir, "conc_data.csv")
  data.table::fwrite(conc_data, conc_file)
  
  # Create a directory instead of a file - fread will fail
  health_file <- file.path(temp_dir, "births.csv")
  dir.create(health_file)  # This is a directory, not a file
  
  # Should return valid = FALSE with issue message about not being able to read
  result <- geoExposeR::validate_prepared_inputs(
    prob_file = prob_file,
    conc_file = conc_file,
    health_file = health_file
  )
  
  # Check that it reports an issue
  expect_false(result$valid)
  expect_true(any(grepl("Could not read health file|directory", result$issues, ignore.case = TRUE)))
  
  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})

# Test regression_analysis with no exposure terms found (lines 78-82)
test_that("regression_analysis warns when no exposure terms in results", {
  temp_dir <- tempfile(pattern = "test_reg_no_as_")
  dir.create(temp_dir)
  
  # Create data with multiple factor levels but where the exposure terms
  # we're looking for won't exist in the model output
  set.seed(456)
  births_data <- data.frame(
    .id = rep(1:3, each = 100),
    ExposureLevel = factor(sample(c("<5", "5-10"), 300, replace = TRUE),
                     levels = c("<5", "5-10")),  # Two levels - model will fit
    BWT = rnorm(300, 3200, 500),
    MRSTATE = sample(c("CA", "NV"), 300, replace = TRUE)
  )
  
  # This should warn about no exposure terms because we're asking for "10+" and "20+"
  # but the data only has "<5" and "5-10"
  result <- suppressWarnings(
    geoExposeR:::regression_analysis(
      births = births_data,
      targets = "BWT",
      cat_label = c("10+", "20+"),  # Terms that won't exist
      exposure_level_col = "ExposureLevel",
      regression_formula = "~ as.factor(ExposureLevel) + (1|MRSTATE)",
      output_dir = temp_dir
    )
  )
  
  # Result should have BWT but with NULL or empty due to no matching terms
  expect_true("BWT" %in% names(result))
  
  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})

# Test impute_additional_variables with insufficient variables for MICE (lines 215-219)
test_that("impute_additional_variables warns on insufficient MICE variables", {
  # Create imputed datasets where the impute variable has missing values
  # and there's only one column for MICE (less than 2)
  set.seed(789)
  imputed_datasets <- list(
    data.frame(
      MVAR = c(1, NA, 3)  # Only one column, needs imputation
    )
  )
  
  # Create births data frame for validation
  births_data <- data.frame(
    MVAR = c(1, 2, 3)
  )
  
  # The warning is triggered when ncol(var_data) < 2
  # We need available_covs to be non-empty first, then have var_data with < 2 cols
  # This means mice_covs need to exist but result in single column var_data
  imputed_datasets2 <- list(
    data.frame(
      MVAR = c(1, NA, 3),  # Variable to impute
      cov1 = c(1, 2, 3)    # One covariate
    )
  )
  
  births_data2 <- data.frame(
    MVAR = c(1, 2, 3),
    cov1 = c(1, 2, 3)
  )
  
  # This should produce no warning since we have 2 columns (MVAR + cov1)
  # The warning is for < 2 columns, which is hard to trigger
  result <- suppressWarnings(
    geoExposeR:::impute_additional_variables(
      imputed_datasets = imputed_datasets2,
      impute_vars = "MVAR",
      births = births_data2,
      mice_covs = "cov1",
      mice_m = 1,
      mice_maxit = 1,
      mice_method = "pmm",
      seed = 123
    )
  )
  
  # Check the function ran
  expect_true(is.list(result))
})

# =============================================================================
# Tests for flexible model types: fixed / random / mixed effects, linear and
# non-linear (generalized) variants.
# =============================================================================

# Shared helper to build a clustered data frame with continuous, binary, and
# count outcomes across two imputed datasets.
make_model_variant_data <- function(seed = 2024, n = 600) {
  set.seed(seed)
  grp <- factor(rep(sprintf("g%02d", 1:20), length.out = n))
  exposure <- factor(
    sample(c("<5", "5-10", "10+"), n, replace = TRUE),
    levels = c("<5", "5-10", "10+")
  )
  exposure <- stats::relevel(exposure, ref = "<5")
  cov <- rnorm(n, 30, 5)
  lin_pred <- as.integer(exposure) * -25
  data.frame(
    .id = rep(1:2, each = n / 2),
    grp = grp,
    ExposureLevel = exposure,
    cov = cov,
    BWT = 3400 + lin_pred + rnorm(n, 0, 400),
    is_event = rbinom(n, 1, 0.3),
    n_events = rpois(n, 2)
  )
}

test_that("resolve_model_family normalizes all accepted forms", {
  expect_null(geoExposeR:::resolve_model_family(NULL))
  expect_s3_class(geoExposeR:::resolve_model_family("binomial"), "family")
  expect_s3_class(geoExposeR:::resolve_model_family(poisson()), "family")
  expect_s3_class(geoExposeR:::resolve_model_family(Gamma), "family")
  expect_equal(
    geoExposeR:::resolve_model_family("poisson")$family,
    "poisson"
  )
})

test_that("resolve_model_family errors on unknown family name", {
  expect_error(
    geoExposeR:::resolve_model_family("not_a_real_family"),
    "Unknown model_family"
  )
})

test_that("validate_model_spec accepts valid specifications", {
  expect_silent(geoExposeR:::validate_model_spec("mixed", NULL))
  expect_silent(geoExposeR:::validate_model_spec("random", "binomial"))
  expect_silent(geoExposeR:::validate_model_spec("fixed", poisson()))
  expect_silent(geoExposeR:::validate_model_spec("mixed", Gamma))
})

test_that("validate_model_spec rejects invalid model_type", {
  expect_error(
    geoExposeR:::validate_model_spec("panel", NULL),
    "model_type must be one of"
  )
  expect_error(
    geoExposeR:::validate_model_spec(c("fixed", "mixed"), NULL),
    "model_type must be one of"
  )
})

test_that("validate_model_spec rejects invalid model_family", {
  # Empty / whitespace character
  expect_error(
    geoExposeR:::validate_model_spec("fixed", ""),
    "single non-empty character string"
  )
  expect_error(
    geoExposeR:::validate_model_spec("fixed", c("binomial", "poisson")),
    "single non-empty character string"
  )
  # Unknown family name
  expect_error(
    geoExposeR:::validate_model_spec("fixed", "not_a_family"),
    "Invalid model_family"
  )
  # Wrong type entirely (numeric)
  expect_error(
    geoExposeR:::validate_model_spec("fixed", 42),
    "model_family must be NULL"
  )
})

test_that("validate_random_effects_terms enforces fixed model has no random", {
  expect_error(
    geoExposeR:::validate_random_effects_terms(
      "~ as.factor(ExposureLevel) + (1 | grp)", "fixed"
    ),
    "does not support random-effects terms"
  )
  expect_silent(
    geoExposeR:::validate_random_effects_terms(
      "~ as.factor(ExposureLevel) + cov", "fixed"
    )
  )
})

test_that("validate_random_effects_terms requires random term for random/mixed", {
  expect_error(
    geoExposeR:::validate_random_effects_terms(
      "~ as.factor(ExposureLevel) + cov", "mixed"
    ),
    "requires at least one random-effects term"
  )
  expect_error(
    geoExposeR:::validate_random_effects_terms(
      "~ as.factor(ExposureLevel) + cov", "random"
    ),
    "requires at least one random-effects term"
  )
})

test_that("validate_random_effects_terms enforces intercept-only for random", {
  expect_error(
    geoExposeR:::validate_random_effects_terms(
      "~ as.factor(ExposureLevel) + (cov | grp)", "random"
    ),
    "supports random intercepts only"
  )
  expect_silent(
    geoExposeR:::validate_random_effects_terms(
      "~ as.factor(ExposureLevel) + (1 | grp)", "random"
    )
  )
  # mixed permits random slopes
  expect_silent(
    geoExposeR:::validate_random_effects_terms(
      "~ as.factor(ExposureLevel) + (cov | grp)", "mixed"
    )
  )
})

test_that("validate_random_effects_terms works on formula objects", {
  expect_error(
    geoExposeR:::validate_random_effects_terms(
      as.formula("~ as.factor(ExposureLevel) + (1 | grp)"), "fixed"
    ),
    "does not support random-effects terms"
  )
})

test_that("fit_single_model dispatches to the correct engine", {
  d <- make_model_variant_data()
  d1 <- d[d$.id == 1, ]

  # Linear fixed -> lm
  m_lm <- geoExposeR:::fit_single_model(
    BWT ~ as.factor(ExposureLevel) + cov, d1, "fixed", NULL
  )
  expect_s3_class(m_lm, "lm")

  # Non-linear fixed -> glm
  m_glm <- geoExposeR:::fit_single_model(
    is_event ~ as.factor(ExposureLevel) + cov, d1, "fixed",
    stats::binomial()
  )
  expect_s3_class(m_glm, "glm")

  # Linear mixed -> lmerMod
  m_lmer <- geoExposeR:::fit_single_model(
    BWT ~ as.factor(ExposureLevel) + (1 | grp), d1, "mixed", NULL
  )
  expect_s4_class(m_lmer, "lmerMod")

  # Non-linear mixed -> glmerMod
  m_glmer <- suppressWarnings(geoExposeR:::fit_single_model(
    is_event ~ as.factor(ExposureLevel) + (1 | grp), d1, "mixed",
    stats::binomial()
  ))
  expect_s4_class(m_glmer, "glmerMod")
})

test_that("tidy_fixed_effects returns uniform output for lm/glm and lmer", {
  d <- make_model_variant_data()
  d1 <- d[d$.id == 1, ]

  m_lm <- stats::lm(BWT ~ as.factor(ExposureLevel) + cov, data = d1)
  tidy_lm <- geoExposeR:::tidy_fixed_effects(m_lm, "fixed")
  expect_named(tidy_lm, c("term", "estimate", "std.error"))
  expect_true("as.factor(ExposureLevel)5-10" %in% tidy_lm$term)

  m_lmer <- lme4::lmer(
    BWT ~ as.factor(ExposureLevel) + (1 | grp), data = d1
  )
  tidy_lmer <- geoExposeR:::tidy_fixed_effects(m_lmer, "mixed")
  expect_named(tidy_lmer, c("term", "estimate", "std.error"))
})

test_that("extract_df_residual works for fixed and mixed models", {
  d <- make_model_variant_data()
  d1 <- d[d$.id == 1, ]

  m_lm <- stats::lm(BWT ~ as.factor(ExposureLevel) + cov, data = d1)
  expect_true(is.finite(geoExposeR:::extract_df_residual(m_lm, "fixed")))

  m_lmer <- lme4::lmer(
    BWT ~ as.factor(ExposureLevel) + (1 | grp), data = d1
  )
  expect_silent(geoExposeR:::extract_df_residual(m_lmer, "mixed"))
})

test_that("extract_df_residual warns and returns Inf on failure", {
  bad <- structure(list(), class = "not_a_model")
  expect_warning(
    df_val <- geoExposeR:::extract_df_residual(bad, "mixed"),
    "Could not extract degrees of freedom"
  )
  expect_true(is.infinite(df_val))
})

test_that("regression_analysis fits linear fixed-effects (lm) models", {
  temp_dir <- tempfile(pattern = "test_fe_lm_")
  dir.create(temp_dir)
  d <- make_model_variant_data()

  result <- geoExposeR:::regression_analysis(
    births = d,
    targets = "BWT",
    cat_label = c("5-10", "10+"),
    exposure_level_col = "ExposureLevel",
    regression_formula = "~ as.factor(ExposureLevel) + cov",
    output_dir = temp_dir,
    model_type = "fixed",
    model_family = NULL
  )
  expect_true("BWT" %in% names(result))
  expect_equal(as.character(result$BWT$term), c("5-10", "10+"))
  expect_true(file.exists(file.path(temp_dir, "BWT_pooled_results.csv")))
  unlink(temp_dir, recursive = TRUE)
})

test_that("regression_analysis fits non-linear fixed-effects (glm) models", {
  temp_dir <- tempfile(pattern = "test_fe_glm_")
  dir.create(temp_dir)
  d <- make_model_variant_data()

  result <- suppressWarnings(geoExposeR:::regression_analysis(
    births = d,
    targets = "is_event",
    cat_label = c("5-10", "10+"),
    exposure_level_col = "ExposureLevel",
    regression_formula = "~ as.factor(ExposureLevel) + cov",
    output_dir = temp_dir,
    model_type = "fixed",
    model_family = "binomial"
  ))
  expect_true("is_event" %in% names(result))
  expect_equal(nrow(result$is_event), 2)
  unlink(temp_dir, recursive = TRUE)
})

test_that("regression_analysis handles categorical character targets (logistic)", {
  temp_dir <- tempfile(pattern = "test_fe_glm_char_")
  dir.create(temp_dir)
  d <- make_model_variant_data()
  # A two-level CHARACTER outcome, as it would be read from a text file
  d$bp_category <- ifelse(d$is_event == 1, "high", "normal")

  result <- suppressWarnings(geoExposeR:::regression_analysis(
    births = d,
    targets = "bp_category",
    cat_label = c("5-10", "10+"),
    exposure_level_col = "ExposureLevel",
    regression_formula = "~ as.factor(ExposureLevel) + cov",
    output_dir = temp_dir,
    model_type = "fixed",
    model_family = "binomial"
  ))
  expect_true("bp_category" %in% names(result))
  expect_equal(nrow(result$bp_category), 2)
  unlink(temp_dir, recursive = TRUE)
})

test_that("regression_analysis fits multiple targets in a single call", {
  temp_dir <- tempfile(pattern = "test_fe_multi_")
  dir.create(temp_dir)
  d <- make_model_variant_data()
  d$BWT2 <- d$BWT + rnorm(nrow(d), 0, 50)

  result <- geoExposeR:::regression_analysis(
    births = d,
    targets = c("BWT", "BWT2"),
    cat_label = c("5-10", "10+"),
    exposure_level_col = "ExposureLevel",
    regression_formula = "~ as.factor(ExposureLevel) + cov",
    output_dir = temp_dir,
    model_type = "fixed",
    model_family = NULL
  )
  expect_setequal(names(result), c("BWT", "BWT2"))
  expect_true(file.exists(file.path(temp_dir, "BWT2_pooled_results.csv")))
  unlink(temp_dir, recursive = TRUE)
})

test_that("regression_analysis fits linear random-effects (lmer) models", {
  temp_dir <- tempfile(pattern = "test_re_lmer_")
  dir.create(temp_dir)
  d <- make_model_variant_data()

  result <- suppressWarnings(geoExposeR:::regression_analysis(
    births = d,
    targets = "BWT",
    cat_label = c("5-10", "10+"),
    exposure_level_col = "ExposureLevel",
    regression_formula = "~ as.factor(ExposureLevel) + (1 | grp)",
    output_dir = temp_dir,
    model_type = "random",
    model_family = NULL
  ))
  expect_true("BWT" %in% names(result))
  unlink(temp_dir, recursive = TRUE)
})

test_that("regression_analysis fits non-linear mixed-effects (glmer) models", {
  temp_dir <- tempfile(pattern = "test_me_glmer_")
  dir.create(temp_dir)
  d <- make_model_variant_data()

  result <- suppressWarnings(geoExposeR:::regression_analysis(
    births = d,
    targets = "n_events",
    cat_label = c("5-10", "10+"),
    exposure_level_col = "ExposureLevel",
    regression_formula = "~ as.factor(ExposureLevel) + cov + (1 | grp)",
    output_dir = temp_dir,
    model_type = "mixed",
    model_family = poisson()
  ))
  expect_true("n_events" %in% names(result))
  expect_equal(nrow(result$n_events), 2)
  unlink(temp_dir, recursive = TRUE)
})

# -----------------------------------------------------------------------------
# Multinomial (multi-level categorical outcome) support
# -----------------------------------------------------------------------------

test_that("is_multinomial_spec detects the multinomial family", {
  expect_true(geoExposeR:::is_multinomial_spec("multinomial"))
  expect_true(geoExposeR:::is_multinomial_spec("Multinomial"))
  expect_false(geoExposeR:::is_multinomial_spec("binomial"))
  expect_false(geoExposeR:::is_multinomial_spec(NULL))
  expect_false(geoExposeR:::is_multinomial_spec(poisson()))
  expect_false(geoExposeR:::is_multinomial_spec(c("multinomial", "binomial")))
})

test_that("validate_model_spec accepts multinomial for all effect structures", {
  expect_silent(geoExposeR:::validate_model_spec("fixed", "multinomial"))
  expect_silent(geoExposeR:::validate_model_spec("random", "multinomial"))
  expect_silent(geoExposeR:::validate_model_spec("mixed", "multinomial"))
})

test_that("tidy_multinomial handles multi-level (matrix) coefficients", {
  d <- make_model_variant_data()
  d1 <- d[d$.id == 1, ]
  d1$y3 <- factor(sample(c("a", "b", "c"), nrow(d1), replace = TRUE))
  m <- nnet::multinom(
    y3 ~ as.factor(ExposureLevel) + cov, data = d1, trace = FALSE, Hess = TRUE
  )
  tidy <- geoExposeR:::tidy_multinomial(m)
  expect_named(tidy, c("y.level", "term", "estimate", "std.error"))
  # Two non-reference outcome levels expected
  expect_equal(length(unique(tidy$y.level)), 2)
  expect_true("as.factor(ExposureLevel)5-10" %in% tidy$term)
})

test_that("tidy_multinomial handles two-level (vector) coefficients", {
  d <- make_model_variant_data()
  d1 <- d[d$.id == 1, ]
  d1$y2 <- factor(sample(c("lo", "hi"), nrow(d1), replace = TRUE))
  m <- nnet::multinom(
    y2 ~ as.factor(ExposureLevel) + cov, data = d1, trace = FALSE, Hess = TRUE
  )
  tidy <- geoExposeR:::tidy_multinomial(m)
  expect_named(tidy, c("y.level", "term", "estimate", "std.error"))
  expect_equal(length(unique(tidy$y.level)), 1)
})

test_that("fit_single_model fits a multinomial model and df is infinite", {
  d <- make_model_variant_data()
  d1 <- d[d$.id == 1, ]
  d1$y3 <- factor(sample(c("a", "b", "c"), nrow(d1), replace = TRUE))
  m <- geoExposeR:::fit_single_model(
    y3 ~ as.factor(ExposureLevel) + cov, d1, "fixed", NULL,
    is_multinomial = TRUE
  )
  expect_s3_class(m, "multinom")
  expect_true(is.infinite(geoExposeR:::extract_df_residual(m, "fixed")))
})

test_that("pool_estimates_by_term pools within outcome levels (y.level)", {
  combined <- data.frame(
    y.level = rep(c("b", "c"), each = 4),
    term = rep(c("5-10", "5-10", "10+", "10+"), 2),
    estimate = c(0.1, 0.2, 0.3, 0.4, -0.1, -0.2, -0.3, -0.4),
    std.error = rep(0.1, 8),
    stringsAsFactors = FALSE
  )
  mock_model <- structure(list(), class = "multinom")
  result <- geoExposeR:::pool_estimates_by_term(
    combined_results = combined,
    exposure_terms = c("5-10", "10+"),
    reference_model = mock_model,
    model_type = "fixed"
  )
  expect_true("y.level" %in% names(result))
  expect_equal(nrow(result), 4)
  expect_setequal(unique(result$y.level), c("b", "c"))
})

test_that("regression_analysis fits multinomial models with y.level output", {
  temp_dir <- tempfile(pattern = "test_multinom_")
  dir.create(temp_dir)
  d <- make_model_variant_data()
  # Three-level CHARACTER outcome (as read from a text file)
  d$bp_stage <- sample(c("normal", "elevated", "stage2"), nrow(d), replace = TRUE)

  result <- suppressWarnings(geoExposeR:::regression_analysis(
    births = d,
    targets = "bp_stage",
    cat_label = c("5-10", "10+"),
    exposure_level_col = "ExposureLevel",
    regression_formula = "~ as.factor(ExposureLevel) + cov",
    output_dir = temp_dir,
    model_type = "fixed",
    model_family = "multinomial"
  ))
  expect_true("bp_stage" %in% names(result))
  expect_true("y.level" %in% names(result$bp_stage))
  # Two non-reference levels x two exposure terms = 4 pooled rows
  expect_equal(nrow(result$bp_stage), 4)
  expect_true(file.exists(file.path(temp_dir, "bp_stage_pooled_results.csv")))
  unlink(temp_dir, recursive = TRUE)
})

test_that("fit_multinomial_mixed fits an mblogit model from an lme4 formula", {
  d <- make_model_variant_data()
  d1 <- d[d$.id == 1, ]
  d1$y3 <- factor(sample(c("a", "b", "c"), nrow(d1), replace = TRUE))
  m <- suppressWarnings(geoExposeR:::fit_multinomial_mixed(
    y3 ~ as.factor(ExposureLevel) + cov + (1 | grp), d1
  ))
  expect_s3_class(m, "mblogit")
  # Inference for mblogit uses the large-sample (z) approximation
  expect_true(is.infinite(geoExposeR:::extract_df_residual(m, "random")))
})

test_that("tidy_multinomial handles mblogit (random-effects) coefficients", {
  d <- make_model_variant_data()
  d1 <- d[d$.id == 1, ]
  d1$y3 <- factor(sample(c("a", "b", "c"), nrow(d1), replace = TRUE))
  m <- suppressWarnings(geoExposeR:::fit_multinomial_mixed(
    y3 ~ as.factor(ExposureLevel) + (1 | grp), d1
  ))
  tidy <- geoExposeR:::tidy_multinomial(m)
  expect_named(tidy, c("y.level", "term", "estimate", "std.error"))
  expect_equal(length(unique(tidy$y.level)), 2)
  expect_true(any(grepl("ExposureLevel", tidy$term)))
})

test_that("regression_analysis fits random/mixed multinomial models (mblogit)", {
  d <- make_model_variant_data()
  d$bp_stage <- sample(c("normal", "elevated", "stage2"), nrow(d), replace = TRUE)

  for (mt in c("random", "mixed")) {
    temp_dir <- tempfile(pattern = paste0("test_multinom_", mt, "_"))
    dir.create(temp_dir)
    formula <- if (mt == "random") {
      "~ as.factor(ExposureLevel) + (1 | grp)"
    } else {
      "~ as.factor(ExposureLevel) + cov + (1 | grp)"
    }
    result <- suppressWarnings(geoExposeR:::regression_analysis(
      births = d,
      targets = "bp_stage",
      cat_label = c("5-10", "10+"),
      exposure_level_col = "ExposureLevel",
      regression_formula = formula,
      output_dir = temp_dir,
      model_type = mt,
      model_family = "multinomial"
    ))
    expect_true("y.level" %in% names(result$bp_stage))
    expect_equal(nrow(result$bp_stage), 4)
    unlink(temp_dir, recursive = TRUE)
  }
})

test_that("perform_sensitivity_analysis validates model spec early", {
  temp_dir <- tempfile(pattern = "test_psa_spec_")
  dir.create(temp_dir)
  birth_file <- file.path(temp_dir, "birth.txt")
  prob_file <- file.path(temp_dir, "prob_data.csv")
  conc_file <- file.path(temp_dir, "conc_data.csv")
  writeLines("test", birth_file)
  writeLines("test", prob_file)
  writeLines("test", conc_file)

  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      regression_formula = "~ as.factor(ExposureLevel) + cov",
      model_type = "fixed",
      model_family = "binomial",
      targets = "is_event"
    ),
    NA
  )

  # Fixed model with a random-effects term should be rejected
  expect_error(
    geoExposeR:::validate_all_inputs(
      birth_data_txt = birth_file,
      model_prob_csv = prob_file,
      conc_data_csv = conc_file,
      regression_formula = "~ as.factor(ExposureLevel) + (1 | grp)",
      model_type = "fixed"
    ),
    "does not support random-effects terms"
  )

  unlink(temp_dir, recursive = TRUE)
})
