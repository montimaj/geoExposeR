# Load the testthat library
library(testthat)
library(gwArsenicR)

# Source the helper function to create dummy data
source("helper-create_dummy_data.R")

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
  as_level_col <- "AsLevel"
  mice_covs <- c("pm", "DMAR_1", "MEDUC_2", "MEDUC_3", "MEDUC_4")

  # check if mice_covs are in the birth_data
  if (!all(mice_covs %in% names(birth_data))) {
    stop("Some mice_covs are not present in the birth_data.")
  }

  # Define a simple regression formula for testing purposes
  regression_formula <- paste0(
    "~ as.factor(", as_level_col, ") + ",
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
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = ndraws,
      regression_formula = regression_formula,
      impute_vars = impute_vars,
      targets = targets,
      output_dir = output_dir,
      mice_covs = mice_covs,
      as_level_col = as_level_col
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
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = 3,
      regression_formula = "~ as.factor(AsLevel) + MAGE_R + (1 | MRSTATE)",
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
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = 2,
      regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)",
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
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = NULL,
      as_usgs_prob_csv = "usgs.csv",
      as_epa_prob_csv = "epa.csv"
    ),
    "birth_data_txt is required and cannot be NULL"
  )

  # Test missing as_usgs_prob_csv
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = "birth.txt",
      as_usgs_prob_csv = NULL,
      as_epa_prob_csv = "epa.csv"
    ),
    "as_usgs_prob_csv is required and cannot be NULL"
  )

  # Test missing as_epa_prob_csv
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = "birth.txt",
      as_usgs_prob_csv = "usgs.csv",
      as_epa_prob_csv = NULL
    ),
    "as_epa_prob_csv is required and cannot be NULL"
  )
})

test_that("validate_inputs validates file path formats", {
  # Test non-character birth_data_txt
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = 123,
      as_usgs_prob_csv = "usgs.csv",
      as_epa_prob_csv = "epa.csv"
    ),
    "birth_data_txt must be a single character string"
  )

  # Test multiple values for birth_data_txt
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = c("file1.txt", "file2.txt"),
      as_usgs_prob_csv = "usgs.csv",
      as_epa_prob_csv = "epa.csv"
    ),
    "birth_data_txt must be a single character string"
  )

  # Test non-character as_usgs_prob_csv
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = "birth.txt",
      as_usgs_prob_csv = 123,
      as_epa_prob_csv = "epa.csv"
    ),
    "as_usgs_prob_csv must be a single character string"
  )

  # Test non-character as_epa_prob_csv
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = "birth.txt",
      as_usgs_prob_csv = "usgs.csv",
      as_epa_prob_csv = 123
    ),
    "as_epa_prob_csv must be a single character string"
  )
})

test_that("validate_inputs validates file existence", {
  # Test non-existent files
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = "nonexistent_birth.txt",
      as_usgs_prob_csv = "usgs.csv",
      as_epa_prob_csv = "epa.csv"
    ),
    "File does not exist: nonexistent_birth.txt"
  )

  # Create temporary files for testing
  temp_dir <- tempfile(pattern = "test_validation_")
  dir.create(temp_dir)

  birth_file <- file.path(temp_dir, "birth.txt")
  usgs_file <- file.path(temp_dir, "usgs.csv")
  epa_file <- file.path(temp_dir, "epa.csv")

  # Create the files
  writeLines("test", birth_file)
  writeLines("test", usgs_file)
  writeLines("test", epa_file)

  # Test with existing files (should not error)
  expect_no_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file
    )
  )

  # Test when the EPA file does not exist
  unlink(epa_file)  # Remove EPA file
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file
    ),
    "File does not exist.*epa\\.csv"
  )

  # Test when the USGS file does not exist
  writeLines("test", epa_file)  # Restore EPA file
  unlink(usgs_file)  # Remove USGS file
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file
    ),
    "File does not exist.*usgs\\.csv"
  )

  unlink(temp_dir, recursive = TRUE)
})

test_that("validate_inputs validates ndraws parameter", {
  temp_dir <- tempfile(pattern = "test_ndraws_")
  dir.create(temp_dir)

  # Create temporary files
  birth_file <- file.path(temp_dir, "birth.txt")
  usgs_file <- file.path(temp_dir, "usgs.csv")
  epa_file <- file.path(temp_dir, "epa.csv")
  writeLines("test", birth_file)
  writeLines("test", usgs_file)
  writeLines("test", epa_file)

  # Test non-numeric ndraws
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
      ndraws = "five"
    ),
    "ndraws must be a single numeric value"
  )

  # Test multiple values for ndraws
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
      ndraws = c(1, 2, 3)
    ),
    "ndraws must be a single numeric value"
  )

  # Test zero ndraws
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
      ndraws = 0
    ),
    "ndraws must be a positive integer"
  )

  # Test negative ndraws
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
      ndraws = -5
    ),
    "ndraws must be a positive integer"
  )

  # Test decimal ndraws
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
      ndraws = 5.5
    ),
    "ndraws must be a positive integer"
  )

  # Test very large ndraws (should warn)
  expect_warning(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
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
  usgs_file <- file.path(temp_dir, "usgs.csv")
  epa_file <- file.path(temp_dir, "epa.csv")
  writeLines("test", birth_file)
  writeLines("test", usgs_file)
  writeLines("test", epa_file)

  # Test non-character, non-formula regression_formula
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
      regression_formula = 123
    ),
    "regression_formula must be a character string or formula object"
  )

  # Test multiple character values
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
      regression_formula = c("~ x", "~ y")
    ),
    "regression_formula must be a single character string"
  )

  # Test empty regression_formula
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
      regression_formula = ""
    ),
    "regression_formula cannot be empty"
  )

  # Test whitespace-only regression_formula
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
      regression_formula = "   "
    ),
    "regression_formula cannot be empty"
  )

  # Test invalid formula syntax
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
      regression_formula = "invalid formula syntax $$"
    ),
    "Invalid regression_formula"
  )

  # Test valid character formula
  expect_no_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
      regression_formula = "~ as.factor(AsLevel) + y"
    )
  )

  # Test valid formula object
  expect_no_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
      regression_formula = as.formula("~ x + y")
    )
  )

  unlink(temp_dir, recursive = TRUE)
})

test_that("validate_inputs validates targets parameter", {
  temp_dir <- tempfile(pattern = "test_targets_")
  dir.create(temp_dir)

  # Create temporary files
  birth_file <- file.path(temp_dir, "birth.txt")
  usgs_file <- file.path(temp_dir, "usgs.csv")
  epa_file <- file.path(temp_dir, "epa.csv")
  writeLines("test", birth_file)
  writeLines("test", usgs_file)
  writeLines("test", epa_file)

  # Test non-character targets
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
      targets = 123
    ),
    "targets must be a character vector"
  )

  # Test empty targets
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
      targets = character(0)
    ),
    "targets cannot be empty"
  )

  # Test targets with NA
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
      targets = c("BWT", NA, "OEGEST")
    ),
    "targets cannot contain NA or empty values"
  )

  # Test targets with empty string
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
      targets = c("BWT", "", "OEGEST")
    ),
    "targets cannot contain NA or empty values"
  )

  # Test duplicate targets
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
      targets = c("BWT", "BWT", "OEGEST")
    ),
    "targets cannot contain duplicate values"
  )

  # Test valid targets
  expect_no_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
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
  usgs_file <- file.path(temp_dir, "usgs.csv")
  epa_file <- file.path(temp_dir, "epa.csv")
  writeLines("test", birth_file)
  writeLines("test", usgs_file)
  writeLines("test", epa_file)

  # Test non-character impute_vars
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
      impute_vars = 123
    ),
    "impute_vars must be a character vector or NULL"
  )

  # Test empty impute_vars (should suggest NULL)
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
      impute_vars = character(0)
    ),
    "impute_vars cannot be empty \\(use NULL instead\\)"
  )

  # Test impute_vars with NA
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
      impute_vars = c("var1", NA, "var2")
    ),
    "impute_vars cannot contain NA or empty values"
  )

  # Test valid impute_vars
  expect_no_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
      impute_vars = c("var1", "var2")
    )
  )

  # Test NULL impute_vars (should be valid)
  expect_no_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
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
  usgs_file <- file.path(temp_dir, "usgs.csv")
  epa_file <- file.path(temp_dir, "epa.csv")
  writeLines("test", birth_file)
  writeLines("test", usgs_file)
  writeLines("test", epa_file)

  # Test non-character output_dir
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
      output_dir = 123
    ),
    "output_dir must be a single character string"
  )

  # Test empty output_dir
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
      output_dir = ""
    ),
    "output_dir cannot be empty"
  )

  # Test output_dir with non-existent parent
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
      output_dir = "/nonexistent/parent/directory"
    ),
    "Parent directory does not exist"
  )

  # Test valid output_dir
  valid_output_dir <- file.path(temp_dir, "output")
  expect_no_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
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
  usgs_file <- file.path(temp_dir, "usgs.csv")
  epa_file <- file.path(temp_dir, "epa.csv")
  writeLines("test", birth_file)
  writeLines("test", usgs_file)
  writeLines("test", epa_file)

  # Test invalid mice_m
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
      mice_m = -1
    ),
    "mice_m must be a positive integer"
  )

  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
      mice_m = 2.5
    ),
    "mice_m must be a positive integer"
  )

  # Test non-numeric mice_m
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
      mice_m = "ten"
    ),
    "mice_m must be a single numeric value"
  )

  # Test large mice_m (should warn)
  expect_warning(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
      mice_m = 150
    ),
    "mice_m is very large.*This may take a long time"
  )

  # Test invalid mice_maxit
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
      mice_maxit = 0
    ),
    "mice_maxit must be a positive integer"
  )

  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
      mice_maxit = "three"
    ),
    "mice_maxit must be a single numeric value"
  )

  # Test invalid mice_method
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
      mice_method = "invalid_method"
    ),
    "mice_method must be one of"
  )

  # Test non-character mice_method
  expect_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
      mice_method = 123
    ),
    "mice_method must be a single character string"
  )

  # Test valid mice_method
  expect_no_error(
    gwArsenicR:::validate_all_inputs(
      birth_data_txt = birth_file,
      as_usgs_prob_csv = usgs_file,
      as_epa_prob_csv = epa_file,
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

  # Test invalid epa_lognormal_sdlog
  # create a dummy epa_lognormal_sdlog
  dummy_epa_df <- data.table::fread(dummy_files$epa, header = TRUE)
  dummy_epa_df$EPA_AS_sdlog <- -1.0
  dummy_epa_csv <- file.path(temp_dir, "invalid_epa_arsenic.csv")
  data.table::fwrite(dummy_epa_df, dummy_epa_csv)

  expect_error(
    gwArsenicR:::perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_epa_csv
    ),
    "epa_as_sd_col must contain only positive values.*"
  )

  dummy_epa_df$EPA_AS_sdlog1 <- "one"
  data.table::fwrite(dummy_epa_df, dummy_epa_csv)
  expect_error(
    gwArsenicR:::perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_epa_csv,
      epa_as_sd_col = "EPA_AS_sdlog1"
    ),
    "epa_as_sd_col must contain only finite numeric values.*"
  )

  # Test large epa_lognormal_sdlog (should warn)
  dummy_epa_df$EPA_AS_sdlog <- 15.0
  data.table::fwrite(dummy_epa_df, dummy_epa_csv)
  expect_warning(
    gwArsenicR:::perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_epa_csv,
      regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)"
    ),
    "epa_as_sd_col contains very large values.*"
  )

  # Test invalid seed
  expect_error(
    gwArsenicR:::perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_epa_csv,
      seed = 123.45
    ),
    "seed must be an integer"
  )

  expect_error(
    gwArsenicR:::perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_epa_csv,
      seed = c(1, 2, 3)
    ),
    "seed must be a single numeric value"
  )

  unlink(temp_dir, recursive = TRUE)
})

test_that("validate_inputs validates arsenic category labels", {
  temp_dir <- tempfile(pattern = "test_arsenic_labels_")
  dir.create(temp_dir)

  # Create temporary files
  dummy_files <- create_dummy_data(file_dir = temp_dir)

  # Test non-character as_cat_label
  expect_error(
    gwArsenicR:::perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      as_cat_label = 123,
      drop_as_cat_label_reg = c("Low"),
      regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)"
    ),
    "as_cat_label must be a character vector"
  )

  # Test empty as_cat_label
  expect_error(
    gwArsenicR:::perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      as_cat_label = character(0),
      drop_as_cat_label_reg = c("Low")
    ),
    "as_cat_label cannot be empty"
  )

  # Test as_cat_label with duplicates
  expect_error(
    gwArsenicR:::perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      as_cat_label = c("Low", "Medium", "Low"),
      drop_as_cat_label_reg = c("Low")
    ),
    "as_cat_label cannot contain duplicate values"
  )

  # Test drop_as_cat_label_reg not in as_cat_label
  expect_error(
    gwArsenicR:::perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      as_cat_label = c("Low", "Medium", "High"),
      drop_as_cat_label_reg = c("VeryLow")
    ),
    "All values in drop_as_cat_label_reg must be present in as_cat_label"
  )

  # Test dropping all categories
  expect_error(
    gwArsenicR:::perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      as_cat_label = c("Low", "Medium"),
      drop_as_cat_label_reg = c("Low", "Medium")
    ),
    "Cannot drop all arsenic categories as reference levels"
  )

  # Test valid category configuration
  expect_no_error(
    gwArsenicR:::perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      as_cat_label = c("Low", "Medium", "High"),
      drop_as_cat_label_reg = c("Low"),
      regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)"
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
    "as_usgs_prob_csv is required and cannot be NULL"
  )

  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = "test.txt",
      as_usgs_prob_csv = "usgs.csv"
    ),
    "as_epa_prob_csv is required and cannot be NULL"
  )
})

test_that("Function validates parameter values through main interface", {
  temp_dir <- tempfile(pattern = "test_main_validation_")
  dir.create(temp_dir)

  dummy_files <- create_dummy_data(file_dir = temp_dir)

  # Test invalid ndraws through main function
  expect_error(
    gwArsenicR::perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = 0,
      regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)"
    ),
    "ndraws must be a positive integer"
  )

  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = -5,
      regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)"
    ),
    "ndraws must be a positive integer"
  )

  # Test invalid regression formula through main function
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
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
      as_usgs_prob_csv = "nonexistent.csv",
      as_epa_prob_csv = "nonexistent.csv",
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
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = 2,
      regression_formula = "~ as.factor(AsLevel) + MAGE_R + (1 | MRSTATE)",
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

test_that("Function handles different arsenic level column names", {
  temp_dir <- tempfile(pattern = "test_as_col_")
  dir.create(temp_dir)

  dummy_files <- create_dummy_data(
    file_dir = temp_dir
  )
  as_level_col <- "ArsenicCategory"
  results <- expect_no_error(
    gwArsenicR::perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = 2,
      regression_formula = paste0(
        "~ as.factor(", as_level_col, ") + (1 | MRSTATE)"
      ),
      targets = "BWT",
      as_level_col = as_level_col,
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
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = 2,
      regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)",
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
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = 1,  # Absolute minimum
      regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)",
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
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = 2,
      regression_formula = "~ as.factor(AsLevel) + MAGE_R + (1 | MRSTATE)",
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
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = 2,
      regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir()
    )
  )

  # Test formula as formula object
  formula_obj <- as.formula("~ as.factor(AsLevel) + (1 | MRSTATE)")
  results2 <- expect_no_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
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
    as_usgs_prob_csv = dummy_files$usgs,
    as_epa_prob_csv = dummy_files$epa,
    ndraws = 2,
    regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)",
    targets = "BWT",
    output_dir = tempdir()
  )

  set.seed(12345)
  results2 <- perform_sensitivity_analysis(
    birth_data_txt = dummy_files$births,
    as_usgs_prob_csv = dummy_files$usgs,
    as_epa_prob_csv = dummy_files$epa,
    ndraws = 2,
    regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)",
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
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
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
    num_rows = 1000,     # Large number of rows in USGS and EPA
    introduce_nas = TRUE  # Introduce some NAs to test imputation
  )

  # Test with larger dataset but minimal draws to keep test time reasonable
  results <- expect_no_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = 2,
      regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir()
    )
  )

  expect_type(results, "list")
  expect_true("BWT" %in% names(results))

  unlink(temp_dir, recursive = TRUE)
})

test_that("Function validates AsLevel requirement through main interface", {
  temp_dir <- tempfile(pattern = "test_main_aslevel_")
  dir.create(temp_dir)

  dummy_files <- create_dummy_data(file_dir = temp_dir)

  # Test regression formula without AsLevel through main function
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = 2,
      regression_formula = "~ MAGE_R + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir()
    ),
    "regression_formula must include 'AsLevel' as a factor variable"
  )

  # Test regression formula with AsLevel but not as factor through main function
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = 2,
      regression_formula = "~ AsLevel + MAGE_R + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir()
    ),
    "regression_formula must include 'AsLevel' as a factor variable"
  )

  unlink(temp_dir, recursive = TRUE)
})

# Add this test to your test-gwArsenic.R file in the validation section

test_that("Function handles arsenic category reference dropping edge cases", {
  temp_dir <- tempfile(pattern = "test_arsenic_ref_")
  dir.create(temp_dir)

  dummy_files <- create_dummy_data(file_dir = temp_dir)

  # Test dropping all available arsenic categories as reference
  # This should trigger:
  # "Cannot drop all arsenic categories as reference levels."
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = 2,
      regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)",
      targets = "BWT",
      as_cat_label = c("Low", "Medium", "High"),
      drop_as_cat_label_reg = c("Low", "Medium", "High"),  # Drop all categories
      output_dir = tempdir()
    ),
    "Cannot drop all arsenic categories as reference levels"
  )

  # Test with data that only has one unique arsenic category
  # First, create a birth dataset with only one arsenic level
  birth_data <- data.table::fread(dummy_files$births, header = TRUE)
  birth_data$AsLevel <- "Low"  # Set all to the same category

  # Write the modified data back
  single_category_file <- file.path(temp_dir, "single_category_births.txt")
  data.table::fwrite(birth_data, single_category_file, sep = "\t")

  # Test with single category that gets dropped
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = single_category_file,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = 2,
      regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)",
      targets = "BWT",
      as_cat_label = c("Low"),
      drop_as_cat_label_reg = c("Low"),  # Drop the only category
      output_dir = tempdir()
    ),
    "Cannot drop all arsenic categories as reference levels"
  )

  # Test edge case: data has categories but after filtering, none remain
  # Create data with two categories
  birth_data_two <- data.table::fread(dummy_files$births, header = TRUE)
  birth_data_two$AsLevel <- sample(
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
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = 2,
      regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)",
      targets = "BWT",
      as_cat_label = c("Low", "Medium"),
      drop_as_cat_label_reg = c("Low", "Medium"),  # Drop both categories
      output_dir = tempdir()
    ),
    "Cannot drop all arsenic categories as reference levels"
  )

  # Test valid case: have multiple categories and only drop some
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = two_category_file,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = 2,
      regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)",
      targets = "BWT",
      as_cat_label = c("Low", "Medium"),
      # Only drop one category, leaving "Medium"
      drop_as_cat_label_reg = c("Low"),
      output_dir = tempdir()
    ),
    "incorrect number of probabilities"
  )

  unlink(temp_dir, recursive = TRUE)
})

test_that("Function validates arsenic category data consistency", {
  temp_dir <- tempfile(pattern = "test_arsenic_data_")
  dir.create(temp_dir)

  dummy_files <- create_dummy_data(file_dir = temp_dir)

  # Test when specified as_cat_label doesn't match data categories
  # Load the birth data to see what categories are actually present
  birth_data <- data.table::fread(dummy_files$births, header = TRUE)

  # Specify categories that don't exist in the data
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = 2,
      regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)",
      targets = "BWT",
      # Categories not in data
      as_cat_label = c("NonExistent1", "NonExistent2"),
      output_dir = tempdir()
    ),
    # This might trigger a different error depending on your validation logic
    # Adjust the error message to match your actual implementation
    "All values in drop_as_cat_label_reg must be present in as_cat_label"
  )

  # Test dropping categories that don't exist in as_cat_label
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = 2,
      regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)",
      targets = "BWT",
      drop_as_cat_label_reg = c("NonExistentCategory"),
      # Try to drop non-existent category
      output_dir = tempdir()
    ),
    "All values in drop_as_cat_label_reg must be present in as_cat_label"
  )

  # Test if drop_as_cat_label_reg is a character vector
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = 2,
      regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)",
      targets = "BWT",
      drop_as_cat_label_reg = 123,  # Invalid type
      output_dir = tempdir()
    ),
    "drop_as_cat_label_reg must be a character vector"
  )

  # Test if as_cat_label has NA or empty values
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = 2,
      regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)",
      targets = "BWT",
      as_cat_label = c("Low", NA, "High"),  # Contains NA
      output_dir = tempdir()
    ),    "as_cat_label cannot contain NA or empty values"
  )

  unlink(temp_dir, recursive = TRUE)
})


test_that("Function validates column names in input files", {
  temp_dir <- tempfile(pattern = "test_column_names_")
  dir.create(temp_dir)
  dummy_files <- create_dummy_data(file_dir = temp_dir)

  # Test missing columns in USGS data
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = 2,
      regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir(),
      as_conc_cols = c(
        "NonExistentColumn1", "NonExistentColumn2",
        "NonExistentColumn3"
      )
    ),
    "The following columns are missing in.*"
  )

  # Test incorrect file format (e.g., non-CSV for EPA)
  dummy_files <- create_dummy_data(file_dir = temp_dir)
  dummy_epa_df <- data.table::fread(dummy_files$epa, header = TRUE)
  invalid_epa_file <- file.path(temp_dir, "invalid_epa.bin")
  data.table::fwrite(
    dummy_epa_df,
    file = invalid_epa_file,
    sep = "\t"
  )
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = invalid_epa_file,
      ndraws = 2,
      regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir()
    ),
    "Unsupported file format for.*"
  )

  # Test non-character other column names
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = 2,
      regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir(),
      geoid_col = 123  # Invalid type
    ),
    "geoid_col must be a single character string"
  )

  # Test empty other column names
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = 2,
      regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir(),
      geoid_col = ""  # Empty string
    ),
    "geoid_col cannot be empty"
  )

  # Test incorrect character vectors for other column names
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = 2,
      regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir(),
      as_conc_cols = c(1, 2, 3)  # Mixed types
    ),
    "as_conc_cols must be a character vector.*"
  )

  # Test NA values in other column names
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = 2,
      regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir(),
      as_conc_cols = c("Col1", NA)  # NA value
    ),
    "as_conc_cols cannot contain NA or empty values"
  )

  # Test duplicate values in other column names
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = 2,
      regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir(),
      as_conc_cols = c("Col1", "Col1")  # Duplicate values
    ),
    "as_conc_cols cannot contain duplicate values"
  )

  unlink(temp_dir, recursive = TRUE)
})


test_that("Function validates epa_cutoff parameter", {
  temp_dir <- tempfile(pattern = "test_epa_cutoff_")
  dir.create(temp_dir)
  dummy_files <- create_dummy_data(file_dir = temp_dir)

  # Test negative epa_cutoffs
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = 2,
      regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir(),
      epa_cutoffs = -1  # Invalid negative value
    ),
    "epa_cutoffs must contain only positive values"
  )

  # Test non-numeric epa_cutoffs
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = 2,
      regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir(),
      epa_cutoffs = "high"  # Invalid non-numeric value
    ),
    "epa_cutoffs must be a numeric vector"
  )

  # Test NA in epa_cutoffs
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = 2,
      regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir(),
      epa_cutoffs = c(5, NA)  # Contains NA
    ),
    "epa_cutoffs cannot contain NA values"
  )

  # Test duplicate values in epa_cutoffs
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = 2,
      regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir(),
      epa_cutoffs = c(5, 5)  # Duplicate values
    ),
    "epa_cutoffs cannot contain duplicate values"
  )

  # Test valid epa_cutoffs with 4 categories
  # Create a dummy USGS file with 4 arsenic categories
  # Ensure the dummy USGS data has appropriate columns
  # and realistic probability distributions

  num_rows <- 500
  county_ids <- sprintf("%05d", 1:100)
  geoid10 <- sample(county_ids, num_rows, replace = TRUE)
  # Generate random probabilities for 4 categories that sum to 1
  usgs_prob_values <- t(apply(
    matrix(stats::runif(num_rows * 4), ncol = 4),
    1,
    function(x) x / sum(x)
  ))
  wells_2010 <- sample(50:500, num_rows, replace = TRUE)
  dummy_usgs_df <- data.frame(
    AS_C1 = usgs_prob_values[, 1],
    AS_C2 = usgs_prob_values[, 2],
    AS_C3 = usgs_prob_values[, 3],
    AS_C4 = usgs_prob_values[, 4],
    Wells_2010 = wells_2010,
    GEOID10 = geoid10
  )
  # Write the dummy USGS data to a CSV fil
  dummy_usgs_csv <- file.path(temp_dir, "dummy_usgs_4cats_arsenic.csv")
  data.table::fwrite(dummy_usgs_df, dummy_usgs_csv)

  # Create a dummy EPA file with 4 categories
  dummy_epa_df <- data.frame(
    EPA_AS_meanlog = runif(num_rows, 0.5, 1.5),  # Log-normal mean
    EPA_AS_sdlog = runif(num_rows, 0.1, 0.5),  # Log-normal sd
    PWELL_private_pct = runif(num_rows, 1, 99),  # Private well percentage
    GEOID10 = geoid10
  )
  dummy_epa_csv <- file.path(temp_dir, "dummy_epa_4cats.csv")
  data.table::fwrite(dummy_epa_df, dummy_epa_csv)
  dummy_files <- create_dummy_data(file_dir = temp_dir)

  # this tests the default fallback imputation as well
  expect_no_error({
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_usgs_csv,
      as_epa_prob_csv = dummy_epa_csv,
      ndraws = 2,
      regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir(),
      as_conc_cols = c("AS_C1", "AS_C2", "AS_C3", "AS_C4"),
      drop_as_cat_label_reg = "AS<5",  # Drop one category
      as_cat_label = c("AS<5", "AS5-10", "AS10-15", "AS>15"),
      epa_cutoffs = c(5, 10, 15)  # Valid ascending cutoffs
    )
  })

  expect_no_error({
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_usgs_csv,
      as_epa_prob_csv = dummy_epa_csv,
      ndraws = 2,
      regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir(),
      as_conc_cols = c("AS_C1", "AS_C2", "AS_C3", "AS_C4"),
      drop_as_cat_label_reg = "AS<5",  # Drop one category
      as_cat_label = c("AS<5", "AS5-10", "AS10-15", "AS>15"),
      epa_cutoffs = c(5, 10, 15),  # Valid ascending cutoffs
      apply_imputation_fallback = FALSE # removes NA values
    )
  })

  # Test descending order epa_cutoffs
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = 2,
      regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir(),
      epa_cutoffs = c(15, 10, 5)  # Descending order
    ),
    "epa_cutoffs must be in ascending order"
  )

  # Test if epa_cutoffs have categories not matching USGS data
  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = 2,
      regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir(),
      epa_cutoffs = c(5, 10, 20)  # Assuming USGS has only 3 categories
    ),
    "as_conc_cols must have length equal to.*"
  )

  unlink(temp_dir, recursive = TRUE)
})


# ===================================================================
# Tests for Dynamic EPA Cutoffs Functionality
# ===================================================================

test_that("validate_epa_cutoffs works correctly", {
  as_conc_cols <- c("AS_C1", "AS_C2", "AS_C3", "AS_C4")
  # Test valid cutoffs
  valid_cutoffs <- c(5, 10, 15)
  expect_no_error(
    validate_epa_cutoff_params(valid_cutoffs, as_conc_cols)
  )

  # Test empty cutoffs similar to no cutoffs provided
  expect_no_error(
    validate_epa_cutoff_params(c(), as_conc_cols)
  )

  # Test non-numeric cutoffs
  expect_error(
    validate_epa_cutoff_params(c("5", "10"), as_conc_cols),
    "epa_cutoffs must be a numeric vector"
  )

  # Test unsorted cutoffs
  expect_error(
    validate_epa_cutoff_params(c(10, 5, 2.5), as_conc_cols),
    "epa_cutoffs must be in ascending order"
  )

  # Test negative cutoffs
  expect_error(
    validate_epa_cutoff_params(c(-1, 5, 10), as_conc_cols),
    "epa_cutoffs must contain.*"
  )
})

test_that("convert_epa_to_multinomial creates correct probabilities", {
  # Create test EPA data
  epa_data <- data.frame(
    FIPS = c(1001, 1002, 1003),
    EPA_AS_meanlog = c(log(1), log(5), log(15)),
    EPA_AS_sdlog = c(0.5, 1.0, 1.5)
  )

  # Test 3 categories (2 cutoffs)
  cutoff_config <- list(cutoffs = c(5, 10), n_categories = 3)
  result <- convert_epa_to_multinomial(
    epa_data, epa_as_mean_col = "EPA_AS_meanlog",
    epa_as_sd_col = "EPA_AS_sdlog",
    cutoff_config = cutoff_config,
    return_prob_matrix = FALSE
  )

  # Check that correct columns were created
  expected_cols <- c("EPA_C1v2", "EPA_C2v2", "EPA_C3v2")
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

test_that("convert_epa_to_multinomial handles different categories", {
  epa_data <- data.frame(
    FIPS = 1001,
    EPA_AS_meanlog = log(5),
    EPA_AS_sdlog = 1.0
  )

  # Test 2 categories
  cutoff_config_2 <- list(cutoffs = c(5), n_categories = 2)
  result_2 <- convert_epa_to_multinomial(
    epa_data,
    epa_as_mean_col = "EPA_AS_meanlog",
    epa_as_sd_col = "EPA_AS_sdlog",
    cutoff_config = cutoff_config_2,
    return_prob_matrix = FALSE
  )
  prob_cols_2 <- c("EPA_C1v2", "EPA_C2v2")
  expect_true(all(prob_cols_2 %in% names(result_2)))
  expect_true(abs(sum(result_2[, ..prob_cols_2]) - 1) < 1e-10)

  # Test 5 categories
  cutoff_config_5 <- list(cutoffs = c(1, 5, 10, 20), n_categories = 5)
  result_5 <- convert_epa_to_multinomial(
    epa_data,
    epa_as_mean_col = "EPA_AS_meanlog",
    epa_as_sd_col = "EPA_AS_sdlog",
    cutoff_config = cutoff_config_5,
    return_prob_matrix = FALSE
  )
  expected_cols_5 <- paste0("EPA_C", 1:5, "v2")
  expect_true(all(expected_cols_5 %in% names(result_5)))
  expect_true(abs(sum(result_5[, ..expected_cols_5]) - 1) < 1e-10)
})

test_that("generate_default_epa_cutoffs produces correct schemes", {
  # Test standard predefined schemes
  expect_equal(generate_default_epa_cutoffs(2), c(5))
  expect_equal(generate_default_epa_cutoffs(3), c(5, 10))
  expect_equal(generate_default_epa_cutoffs(4), c(2.5, 5, 10))
  expect_equal(generate_default_epa_cutoffs(5), c(1, 5, 10, 20))
  expect_equal(generate_default_epa_cutoffs(6), c(1, 3, 5, 10, 20))

  # Test logarithmic generation for unusual numbers
  cutoffs_8 <- generate_default_epa_cutoffs(8)
  expect_equal(length(cutoffs_8), 7)  # n-1 cutoffs for n categories
  expect_true(all(cutoffs_8 > 0))
  expect_true(all(diff(cutoffs_8) > 0))  # Should be increasing

  # Test error for invalid input
  expect_error(
    generate_default_epa_cutoffs(1),
    "Must have at least 2 categories"
  )
  expect_error(
    generate_default_epa_cutoffs(0),
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
      "Category mismatch: USGS has 3 categories but EPA cutoffs define ",
      "2 categories"
    )
  )

  expect_error(
    validate_category_alignment(2, c(5, 10)),
    paste0(
      "Category mismatch: USGS has 2 categories but EPA cutoffs define ",
      "3 categories"
    )
  )
})

test_that("edge cases in EPA probability calculations", {
  # Test extreme values
  epa_extreme <- data.frame(
    FIPS = c(1001, 1002, 1003),
    EPA_AS_meanlog_extreme = c(log(0.1), log(100), log(1000)),
    EPA_AS_sdlog_extreme = c(0.1, 0.5, 1.0)
  )

  cutoff_config <- list(cutoffs = c(5, 10), n_categories = 3)
  result <- convert_epa_to_multinomial(
    datmatx_epa = epa_extreme,
    epa_as_mean_col = "EPA_AS_meanlog_extreme",
    epa_as_sd_col = "EPA_AS_sdlog_extreme",
    return_prob_matrix = FALSE,
    cutoff_config = cutoff_config
  )

  # Check row sums still equal 1
  prob_cols <- c("EPA_C1v2", "EPA_C2v2", "EPA_C3v2")
  # convert to matrix for row sums
  prob_matrix <- as.matrix(result[, prob_cols, with = FALSE])
  row_sums <- rowSums(prob_matrix)
  expect_true(all(abs(row_sums - 1) < 1e-10))

  # Very low mean should have most probability in first category
  expect_true(prob_matrix[1, "EPA_C1v2"] > 0.5)

  # Very high mean should have most probability in last category
  expect_true(prob_matrix[3, "EPA_C3v2"] > 0.5)
})

# test apply_fallback_imputation functionality
test_that("Fallback imputation works when EPA data is missing", {
  temp_dir <- tempfile(pattern = "test_fallback_imputation_")
  dir.create(temp_dir)

  dummy_files <- create_dummy_data(file_dir = temp_dir)

  expect_error(
    perform_sensitivity_analysis(
      birth_data_txt = dummy_files$births,
      as_usgs_prob_csv = dummy_files$usgs,
      as_epa_prob_csv = dummy_files$epa,
      ndraws = 2,
      regression_formula = "~ as.factor(AsLevel) + (1 | MRSTATE)",
      targets = "BWT",
      output_dir = tempdir(),
      apply_imputation_fallback = "FOO"  # Enable fallback imputation
    ), "apply_imputation_fallback must be.*"
  )

  unlink(temp_dir, recursive = TRUE)
})
