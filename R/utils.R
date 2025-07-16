#' Format Geographic Identifiers
#'
#' @param data_frame Data frame to format
#' @param geoid_col Column name for geographic identifier
#' @return Data frame with formatted geographic IDs
#'
#' @keywords internal
format_geographic_ids <- function(data_frame, geoid_col) {
  geoid_numeric <- as.numeric(data_frame[[geoid_col]])
  data_frame[[geoid_col]] <- sprintf("%05d", geoid_numeric)
  return(data_frame)
}

#' Validate Input Parameters for Sensitivity Analysis
#'
# @param birth_data_txt Path to birth data text file
#' @param as_usgs_prob_csv Path to USGS arsenic probability CSV file
#' @param as_epa_prob_csv Path to EPA arsenic probability CSV file
#' @param ndraws Number of draws for sensitivity analysis
#' @param regression_formula Regression formula for modeling
#' @param targets Target variables for analysis
#' @param impute_vars Variables to impute
#' @param output_dir Output directory for results
#' @param mice_m Number of multiple imputations
#' @param mice_maxit Maximum iterations for MICE
#' @param mice_method MICE method for imputation
#' @param epa_as_mean_col Column name for EPA arsenic meanlog
#' @param epa_as_sd_col Column name for EPA arsenic sdlog
#' @param seed Random seed for reproducibility
#' @param as_cat_label Arsenic category label
#' @param as_level_col Column name for arsenic level
#' @param drop_as_cat_label_reg Regular expression to drop arsenic categories
#' @param epa_cutoffs EPA cutoffs for arsenic categories
#' @param usgs_columns USGS columns to include in analysis
#' @param validate_cutoffs Whether to validate cutoffs
#' @param columns_to_select Columns to select from the data
#' @param rucc_col Column name for RUCC classification
#' @param geoid_col Column name for geographic identifier
#' @param as_conc_cols Arsenic concentration columns
#' @param pop_well_col Column name for population well
#' @param birth_county_col Column name for birth county
#' @param epa_pwell_col Column name for EPA population well
#' @param mice_covs Covariates for MICE imputation
#' @param apply_imputation_fallback Whether to apply imputation fallback
#' @return NULL (throws error if validation fails)
#' @keywords internal
validate_all_inputs <- function(
  birth_data_txt = NULL,
  as_usgs_prob_csv = NULL,
  as_epa_prob_csv = NULL,
  ndraws = NULL,
  regression_formula = NULL,
  targets = NULL,
  impute_vars = NULL,
  output_dir = NULL,
  mice_m = NULL,
  mice_maxit = NULL,
  mice_method = NULL,
  epa_as_mean_col = NULL,
  epa_as_sd_col = NULL,
  seed = NULL,
  as_cat_label = NULL,
  as_level_col = "AsLevel",
  drop_as_cat_label_reg = NULL,
  epa_cutoffs = NULL,
  usgs_columns = NULL,
  validate_cutoffs = TRUE,
  columns_to_select = NULL,
  rucc_col = NULL,
  geoid_col = NULL,
  as_conc_cols = NULL,
  pop_well_col = NULL,
  birth_county_col = NULL,
  epa_pwell_col = NULL,
  mice_covs = NULL,
  apply_imputation_fallback = TRUE
) {

  # Validate required file paths
  validate_file_paths(birth_data_txt, as_usgs_prob_csv, as_epa_prob_csv)

  # Validate numeric parameters
  validate_numeric_params(ndraws, mice_m, mice_maxit, seed)

  # Validate imputation fallback parameter
  validate_imputation_fallback(apply_imputation_fallback)

  # Validate character vectors
  validate_character_vectors(
    targets = targets,
    impute_vars = impute_vars,
    as_cat_label = as_cat_label,
    mice_covs = mice_covs,
    columns_to_select = columns_to_select,
    as_conc_cols = as_conc_cols
  )

  # Validate character scalars (column names)
  validate_column_names(
    as_level_col, rucc_col, geoid_col, pop_well_col,
    birth_county_col, epa_as_mean_col, epa_as_sd_col,
    epa_pwell_col
  )

  # Validate MICE parameters
  validate_mice_params(mice_method)

  # Validate output directory
  validate_output_directory(output_dir)

  # Validate arsenic category parameters
  validate_arsenic_categories(as_cat_label, drop_as_cat_label_reg)

  # Validate EPA cutoffs
  validate_epa_cutoff_params(epa_cutoffs, as_conc_cols)

  # Validate EPA standard deviation column
  validate_epa_sd_values(as_epa_prob_csv, epa_as_sd_col)

  # Validate USGS columns
  validate_usgs_columns_param(usgs_columns)

  # Validate file contents and column existence
  validate_file_contents(
    birth_data_txt, as_usgs_prob_csv, as_epa_prob_csv,
    rucc_col, birth_county_col, geoid_col, as_conc_cols,
    pop_well_col, epa_as_mean_col, epa_as_sd_col,
    epa_pwell_col
  )

  # Validate regression formula
  validate_regression_formula(regression_formula, as_level_col)

  invisible(NULL)
}

#' Validate File Paths
#' @keywords internal
validate_file_paths <- function(
    birth_data_txt, as_usgs_prob_csv,
    as_epa_prob_csv) {
  # Check required parameters
  if (is.null(birth_data_txt)) {
    stop("birth_data_txt is required and cannot be NULL")
  }
  if (is.null(as_usgs_prob_csv)) {
    stop("as_usgs_prob_csv is required and cannot be NULL")
  }
  if (is.null(as_epa_prob_csv)) {
    stop("as_epa_prob_csv is required and cannot be NULL")
  }

  # Check parameter types
  file_params <- list(
    birth_data_txt = birth_data_txt,
    as_usgs_prob_csv = as_usgs_prob_csv,
    as_epa_prob_csv = as_epa_prob_csv
  )

  for (param_name in names(file_params)) {
    param_value <- file_params[[param_name]]
    if (!is.character(param_value) || length(param_value) != 1) {
      stop(param_name, " must be a single character string")
    }
  }

  # Check file existence
  for (param_name in names(file_params)) {
    file_path <- file_params[[param_name]]
    if (!file.exists(file_path)) {
      stop("File does not exist: ", file_path)
    }
  }
}

#' Validate Numeric Parameters
#' @keywords internal
validate_numeric_params <- function(ndraws, mice_m, mice_maxit, seed) {
  # Validate ndraws
  if (!is.null(ndraws)) {
    if (!is.numeric(ndraws) || length(ndraws) != 1) {
      stop("ndraws must be a single numeric value")
    }
    if (ndraws <= 0 || ndraws != as.integer(ndraws)) {
      stop("ndraws must be a positive integer")
    }
    if (ndraws > 10000) {
      warning("ndraws is very large (", ndraws, 
              "). This may take a long time to compute.")
    }
  }

  # Validate mice_m
  if (!is.null(mice_m)) {
    if (!is.numeric(mice_m) || length(mice_m) != 1) {
      stop("mice_m must be a single numeric value")
    }
    if (mice_m <= 0 || mice_m != as.integer(mice_m)) {
      stop("mice_m must be a positive integer")
    }
    if (mice_m > 100) {
      warning("mice_m is very large (", mice_m, 
              "). This may take a long time to compute.")
    }
  }

  # Validate mice_maxit
  if (!is.null(mice_maxit)) {
    if (!is.numeric(mice_maxit) || length(mice_maxit) != 1) {
      stop("mice_maxit must be a single numeric value")
    }
    if (mice_maxit <= 0 || mice_maxit != as.integer(mice_maxit)) {
      stop("mice_maxit must be a positive integer")
    }
  }

  # Validate seed
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1) {
      stop("seed must be a single numeric value")
    }
    if (seed != as.integer(seed)) {
      stop("seed must be an integer")
    }
  }
}

#' Validate Regression Formula
#' @keywords internal
validate_regression_formula <- function(regression_formula, as_level_col) {
  if (!is.null(regression_formula)) {
    if (!is.character(regression_formula) &&
          !inherits(regression_formula, "formula")) {
      stop("regression_formula must be a character string or formula object")
    }

    if (is.character(regression_formula)) {
      if (length(regression_formula) != 1) {
        stop("regression_formula must be a single character string")
      }

      if (nchar(trimws(regression_formula)) == 0) {
        stop("regression_formula cannot be empty")
      }

      # Try to convert to formula to check validity
      tryCatch(
        as.formula(regression_formula),
        error = function(e) {
          stop("Invalid regression_formula: ", e$message)
        }
      )

      pattern <- paste0("as\\.factor\\(", as_level_col,
                        "\\)|factor\\(", as_level_col, "\\)")
      if (!grepl(pattern, regression_formula)) {
        stop("regression_formula must include '", as_level_col,
             "' as a factor variable")
      }
    }
  }
}

#' Validate Character Vector Parameters
#' @keywords internal
validate_character_vectors <- function(
    targets, impute_vars, as_cat_label,
    mice_covs, columns_to_select,
    as_conc_cols) {

  char_vector_params <- list(
    targets = targets,
    impute_vars = impute_vars,
    as_cat_label = as_cat_label,
    mice_covs = mice_covs,
    columns_to_select = columns_to_select,
    as_conc_cols = as_conc_cols
  )

  for (param_name in names(char_vector_params)) {
    param_value <- char_vector_params[[param_name]]

    if (!is.null(param_value)) {
      if (!is.character(param_value)) {
        stop(param_name, " must be a character vector or NULL")
      }

      if (length(param_value) == 0) {
        if (param_name == "impute_vars") {
          stop("impute_vars cannot be empty (use NULL instead)")
        } else {
          stop(param_name, " cannot be empty")
        }
      }

      if (any(is.na(param_value)) ||
            any(nchar(trimws(param_value)) == 0)) {
        stop(param_name, " cannot contain NA or empty values")
      }

      if (length(unique(param_value)) != length(param_value)) {
        stop(param_name, " cannot contain duplicate values")
      }
    }
  }
}

#' Validate Column Name Parameters
#' @keywords internal
validate_column_names <- function(
    as_level_col, rucc_col, geoid_col,
    pop_well_col, birth_county_col,
    epa_as_mean_col, epa_as_sd_col,
    epa_pwell_col) {

  char_scalar_params <- list(
    as_level_col = as_level_col,
    rucc_col = rucc_col,
    geoid_col = geoid_col,
    pop_well_col = pop_well_col,
    birth_county_col = birth_county_col,
    epa_as_mean_col = epa_as_mean_col,
    epa_as_sd_col = epa_as_sd_col,
    epa_pwell_col = epa_pwell_col
  )

  for (param_name in names(char_scalar_params)) {
    param_value <- char_scalar_params[[param_name]]
    if (!is.null(param_value)) {
      if (!is.character(param_value) || length(param_value) != 1) {
        stop(param_name, " must be a single character string")
      }
      if (nchar(trimws(param_value)) == 0) {
        stop(param_name, " cannot be empty")
      }
    }
  }
}

#' Validate MICE Parameters
#' @keywords internal
validate_mice_params <- function(mice_method) {
  if (!is.null(mice_method)) {
    if (!is.character(mice_method) || length(mice_method) != 1) {
      stop("mice_method must be a single character string")
    }

    valid_methods <- c(
      "pmm", "logreg", "polyreg", "polr", "mean", "norm", "norm.nob",
      "norm.boot", "norm.predict"
    )
    if (!mice_method %in% valid_methods) {
      stop("mice_method must be one of: ",
           paste(valid_methods, collapse = ", "))
    }
  }
}

#' Validate Output Directory
#' @keywords internal
validate_output_directory <- function(output_dir) {
  if (!is.null(output_dir)) {
    if (!is.character(output_dir) || length(output_dir) != 1) {
      stop("output_dir must be a single character string")
    }

    if (nchar(trimws(output_dir)) == 0) {
      stop("output_dir cannot be empty")
    }

    # Check if parent directory exists (if output_dir doesn't exist yet)
    if (!dir.exists(output_dir)) {
      parent_dir <- dirname(output_dir)
      if (!dir.exists(parent_dir)) {
        stop("Parent directory does not exist: ", parent_dir)
      }
    }
  }
}

#' Validate Arsenic Category Parameters
#' @keywords internal
validate_arsenic_categories <- function(as_cat_label, drop_as_cat_label_reg) {
  if (!is.null(drop_as_cat_label_reg)) {
    if (!is.character(drop_as_cat_label_reg)) {
      stop("drop_as_cat_label_reg must be a character vector")
    }

    # Check that drop_as_cat_label_reg values are in as_cat_label
    if (!is.null(as_cat_label)) {
      if (!all(drop_as_cat_label_reg %in% as_cat_label)) {
        stop("All values in drop_as_cat_label_reg must be present in ",
             "as_cat_label")
      }

      # Check that we don't drop all categories
      remaining_cats <- setdiff(as_cat_label, drop_as_cat_label_reg)
      if (length(remaining_cats) == 0) {
        stop("Cannot drop all arsenic categories as reference levels")
      }
    }
  }
}

#' Validate EPA Cutoffs Parameter
#' @keywords internal
validate_epa_cutoff_params <- function(epa_cutoffs, as_conc_cols) {
  if (!is.null(epa_cutoffs)) {
    if (!is.numeric(epa_cutoffs)) {
      stop("epa_cutoffs must be a numeric vector")
    }
    if (any(is.na(epa_cutoffs))) {
      stop("epa_cutoffs cannot contain NA values")
    }
    if (any(epa_cutoffs <= 0)) {
      stop("epa_cutoffs must contain only positive values")
    }
    if (is.unsorted(epa_cutoffs)) {
      stop("epa_cutoffs must be in ascending order")
    }
    if (length(unique(epa_cutoffs)) != length(epa_cutoffs)) {
      stop("epa_cutoffs cannot contain duplicate values")
    }

    # Check alignment with as_conc_cols if both are provided
    if (!is.null(as_conc_cols)) {
      epa_categories <- length(epa_cutoffs) + 1
      if (length(as_conc_cols) != epa_categories) {
        stop("as_conc_cols must have length equal to length(epa_cutoffs) + 1")
      }
    }
  }
}

#' Validate USGS Columns Parameter
#' @keywords internal
validate_usgs_columns_param <- function(usgs_columns) {
  if (!is.null(usgs_columns)) {
    if (!is.character(usgs_columns)) {
      stop("usgs_columns must be a character vector")
    }
    if (length(usgs_columns) == 0) {
      stop("usgs_columns cannot be empty")
    }
    if (any(is.na(usgs_columns)) || any(nchar(usgs_columns) == 0)) {
      stop("usgs_columns cannot contain NA or empty strings")
    }
    if (length(unique(usgs_columns)) != length(usgs_columns)) {
      stop("usgs_columns cannot contain duplicate column names")
    }
  }
}

#' Validate File Contents and Column Existence
#' @keywords internal
validate_file_contents <- function(
    birth_data_txt, as_usgs_prob_csv,
    as_epa_prob_csv, rucc_col, birth_county_col,
    geoid_col, as_conc_cols, pop_well_col,
    epa_as_mean_col, epa_as_sd_col,
    epa_pwell_col) {

  # Define required columns for each file
  birth_cols <- c(rucc_col, birth_county_col)
  if (is.null(rucc_col)) {
    birth_cols <- c(birth_county_col)
  }
  usgs_cols <- c(geoid_col, as_conc_cols, pop_well_col)
  epa_cols <- c(geoid_col, epa_as_mean_col, epa_as_sd_col, epa_pwell_col)

  data_list <- list(
    birth_data_txt = birth_data_txt,
    as_usgs_prob_csv = as_usgs_prob_csv,
    as_epa_prob_csv = as_epa_prob_csv
  )

  data_col_checks <- list(
    birth_data_txt = birth_cols,
    as_usgs_prob_csv = usgs_cols,
    as_epa_prob_csv = epa_cols
  )

  # Check column existence in each file
  for (data_param in names(data_list)) {
    file_path <- data_list[[data_param]]
    required_cols <- data_col_checks[[data_param]]

    # Read only the header to get column names
    if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
      cols_in_file <- names(
        data.table::fread(
          file_path, nrows = 0, header = TRUE,
          showProgress = FALSE, nThread = parallel::detectCores()
        )
      )
    } else if (grepl("\\.txt$", file_path, ignore.case = TRUE)) {
      cols_in_file <- names(
        data.table::fread(
          file_path, nrows = 0, header = TRUE, sep = "\t",
          showProgress = FALSE, nThread = parallel::detectCores()
        )
      )
    } else {
      stop("Unsupported file format for: ", file_path)
    }

    missing_cols <- setdiff(required_cols, cols_in_file)
    if (length(missing_cols) > 0) {
      stop("The following columns are missing in ", data_param, ": ",
           paste(missing_cols, collapse = ", "))
    }
  }

  # Validate EPA standard deviation column values
  validate_epa_sd_values(as_epa_prob_csv, epa_as_sd_col)
}

#' Validate EPA Standard Deviation Values
#' @keywords internal
validate_epa_sd_values <- function(as_epa_prob_csv, epa_as_sd_col) {
  if (!is.null(epa_as_sd_col)) {
    # Read the EPA data to check the column
    epa_data <- data.table::fread(
      as_epa_prob_csv, header = TRUE,
      showProgress = FALSE,
      nThread = parallel::detectCores()
    )

    epa_sd_log_values <- epa_data[[epa_as_sd_col]]

    # Check non-numeric values
    if (!all(is.finite(epa_sd_log_values))) {
      stop("epa_as_sd_col must contain only finite numeric values. ",
           "Check for non-numeric values or extreme outliers.")
    }

    if (any(epa_sd_log_values <= 0)) {
      stop("epa_as_sd_col must contain only positive values ")
    }

    # Check large values
    if (any(epa_sd_log_values > 10)) {
      warning("epa_as_sd_col contains very large values (> 10). This may ",
              "lead to extreme probabilities in the lognormal distribution.")
    }
  }
}

#' Validate Imputation Fallback Parameter
#' @keywords internal
validate_imputation_fallback <- function(apply_imputation_fallback) {
  if (!is.logical(apply_imputation_fallback) ||
        length(apply_imputation_fallback) != 1) {
    stop("apply_imputation_fallback must be a single logical value ",
         "(TRUE or FALSE)")
  }
}