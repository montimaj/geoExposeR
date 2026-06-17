#' @importFrom data.table fread fwrite data.table as.data.table rbindlist
#'   set setnames setDT copy `:=` .SD .N
#' @importFrom dplyr `%>%` select mutate group_by do all_of
#' @importFrom rlang .data
#' @importFrom stats as.formula na.omit plnorm pt qt update
#' @importFrom stats lm glm df.residual coef
#' @importFrom stats gaussian binomial poisson Gamma inverse.gaussian
#'   quasibinomial quasipoisson quasi
#' @importFrom terra rast vect extract
#' @importFrom sf st_centroid st_crs st_transform st_read st_as_sf
#' @importFrom exactextractr exact_extract
#' @importFrom Hmisc wtd.mean wtd.var
#' @import lme4
#' @importFrom nnet multinom
#' @importFrom mclogit mblogit
#' @import rms
#' @import broom.mixed
#' @import Amelia
#' @import mice
#' @keywords internal
"_PACKAGE"

#' Format Identifiers (Smart Detection)
#'
#' Intelligently formats identifiers based on their detected type. For numeric
#' IDs that appear to be FIPS codes (1-5 digits), applies zero-padding to
#' 5 digits. For other formats (alphanumeric, longer codes, etc.), converts
#' to character strings without modification.
#'
#' @param data_frame Data frame to format
#' @param id_col Column name for the identifier column
#' @return Data frame with formatted identifiers as character strings
#'
#' @details
#' The function uses the following logic:
#' \itemize{
#'   \item If all non-NA values are numeric and <= 99999, applies 5-digit
#'         zero-padding (typical for county FIPS codes)
#'   \item If values are already character strings or contain non-numeric
#'         values, converts to character without modification
#'   \item NA values are preserved as NA
#' }
#'
#' @keywords internal
format_ids_func <- function(data_frame, id_col) {
  ids <- data_frame[[id_col]]
  
  # Check if IDs appear to be numeric FIPS-like codes
  # (all numeric values between 1 and 99999)
  if (is.numeric(ids)) {
    non_na_ids <- ids[!is.na(ids)]
    if (length(non_na_ids) > 0 && 
        all(non_na_ids == floor(non_na_ids)) &&  # All integers
        all(non_na_ids >= 1 & non_na_ids <= 99999)) {  # FIPS-like range
      # Apply 5-digit zero-padding for FIPS codes
      data_frame[[id_col]] <- sprintf("%05d", as.integer(ids))
    } else {
      # Longer numeric codes or decimals - just convert to character
      data_frame[[id_col]] <- as.character(ids)
    }
  } else {
    # Already character or mixed - ensure character type
    data_frame[[id_col]] <- as.character(ids)
  }
  
  return(data_frame)
}

#' Validate Input Parameters for Sensitivity Analysis
#'
# @param birth_data_txt Path to birth data text file
#' @param model_prob_csv Path to modeled probability data CSV file
#' @param conc_data_csv Path to measured concentration data CSV file
#' @param ndraws Number of draws for sensitivity analysis
#' @param regression_formula Regression formula for modeling
#' @param targets Target variables for analysis
#' @param impute_vars Variables to impute
#' @param output_dir Output directory for results
#' @param mice_m Number of multiple imputations
#' @param mice_maxit Maximum iterations for MICE
#' @param mice_method MICE method for imputation
#' @param conc_mean_col Column name for concentration lognormal meanlog
#' @param conc_sd_col Column name for concentration lognormal sdlog
#' @param seed Random seed for reproducibility
#' @param cat_label Exposure category label
#' @param exposure_level_col Column name for exposure level
#' @param drop_cat_label_ref Regular expression to drop exposure categories
#' @param conc_cutoffs Concentration cutoffs for exposure categories
#' @param prob_columns Probability columns to include in analysis
#' @param validate_cutoffs Whether to validate cutoffs
#' @param columns_to_select Columns to select from the data
#' @param rucc_col Column name for RUCC classification
#' @param id_col Column name for geographic identifier
#' @param prob_cols Exposure probability columns
#' @param pop_well_col Column name for population well
#' @param record_id_col Column name for birth county
#' @param pwell_col Column name for private well percentage
#' @param mice_covs Covariates for MICE imputation
#' @param apply_imputation_fallback Whether to apply imputation fallback
#' @param model_type Effects structure of the regression model: one of
#'   `"mixed"`, `"random"`, or `"fixed"`
#' @param model_family Error distribution / link for non-linear (generalized)
#'   models, or `NULL` for a linear model
#' @return NULL (throws error if validation fails)
#' @keywords internal
validate_all_inputs <- function(
  birth_data_txt = NULL,
  model_prob_csv = NULL,
  conc_data_csv = NULL,
  ndraws = NULL,
  regression_formula = NULL,
  model_type = "mixed",
  model_family = NULL,
  targets = NULL,
  impute_vars = NULL,
  output_dir = NULL,
  mice_m = NULL,
  mice_maxit = NULL,
  mice_method = NULL,
  conc_mean_col = NULL,
  conc_sd_col = NULL,
  conc_raw_col = NULL,
  default_sdlog = 1.0,
  seed = NULL,
  cat_label = NULL,
  exposure_level_col = "ExposureLevel",
  drop_cat_label_ref = NULL,
  conc_cutoffs = NULL,
  prob_columns = NULL,
  validate_cutoffs = TRUE,
  columns_to_select = NULL,
  rucc_col = NULL,
  id_col = NULL,
  prob_cols = NULL,
  pop_well_col = NULL,
  record_id_col = NULL,
  pwell_col = NULL,
  mice_covs = NULL,
  apply_imputation_fallback = TRUE
) {

  # Validate required file paths
  validate_file_paths(birth_data_txt, model_prob_csv, conc_data_csv)

  # Validate numeric parameters
  validate_numeric_params(ndraws, mice_m, mice_maxit, seed)

  # Validate imputation fallback parameter
  validate_imputation_fallback(apply_imputation_fallback)

  # Validate default_sdlog
  validate_default_sdlog(default_sdlog)

  # Validate character vectors
  validate_character_vectors(
    targets = targets,
    impute_vars = impute_vars,
    cat_label = cat_label,
    mice_covs = mice_covs,
    columns_to_select = columns_to_select,
    prob_cols = prob_cols
  )

  # Validate character scalars (column names)
  validate_column_names(
    exposure_level_col, rucc_col, id_col, pop_well_col,
    record_id_col, conc_mean_col, conc_sd_col,
    pwell_col
  )

  # Validate MICE parameters
  validate_mice_params(mice_method)

  # Validate output directory
  validate_output_directory(output_dir)

  # Validate exposure category parameters
  validate_exposure_categories(cat_label, drop_cat_label_ref)

  # Validate concentration cutoffs
  validate_conc_cutoff_params(conc_cutoffs, prob_cols)

  # Validate concentration lognormal columns (either lognormal cols or concentration col)
  validate_conc_lognormal_input(
    conc_data_csv, conc_mean_col, conc_sd_col, conc_raw_col
  )

  # Validate probability columns
  validate_prob_columns_param(prob_columns)

  # Validate file contents and column existence
  validate_file_contents(
    birth_data_txt, model_prob_csv, conc_data_csv,
    rucc_col, record_id_col, id_col, prob_cols,
    pop_well_col, conc_mean_col, conc_sd_col,
    conc_raw_col, pwell_col
  )

  # Validate model specification (effects structure and family)
  validate_model_spec(model_type, model_family)

  # Validate regression formula
  validate_regression_formula(
    regression_formula, exposure_level_col, model_type
  )

  invisible(NULL)
}

#' Validate Model Specification (effects structure and family)
#'
#' Checks that `model_type` is one of the supported effects structures and that
#' `model_family` is a valid family specification (or `NULL` for linear models).
#'
#' @param model_type Effects structure: one of `"mixed"`, `"random"`, `"fixed"`
#' @param model_family `NULL`, a family name, the special value `"multinomial"`,
#'   a family object, or a family-generating function
#' @return NULL (throws error if validation fails)
#' @keywords internal
validate_model_spec <- function(model_type, model_family) {
  valid_types <- c("mixed", "random", "fixed")
  if (!is.character(model_type) || length(model_type) != 1 ||
        !(model_type %in% valid_types)) {
    stop(
      "model_type must be one of: ",
      paste(valid_types, collapse = ", ")
    )
  }
  if (!is.null(model_family)) {
    if (is.character(model_family)) {
      if (length(model_family) != 1 || nchar(trimws(model_family)) == 0) {
        stop("model_family must be a single non-empty character string")
      }
      if (!is_multinomial_spec(model_family)) {
        # "multinomial" is a recognized special value (fitted with
        # nnet::multinom for fixed effects or mclogit::mblogit for
        # random-/mixed-effects); any other character value must resolve to a
        # stats family object.
        tryCatch(
          resolve_model_family(model_family),
          error = function(e) {
            stop("Invalid model_family: ", conditionMessage(e))
          }
        )
      }
    } else if (!inherits(model_family, "family") &&
                 !is.function(model_family)) {
      stop(
        "model_family must be NULL, a character string, a family object, ",
        "or a family-generating function"
      )
    }
  }
  invisible(NULL)
}

#' Validate File Paths
#' @keywords internal
validate_file_paths <- function(
    birth_data_txt, model_prob_csv,
    conc_data_csv) {
  # Check required parameters
  if (is.null(birth_data_txt)) {
    stop("birth_data_txt is required and cannot be NULL")
  }
  if (is.null(model_prob_csv)) {
    stop("model_prob_csv is required and cannot be NULL")
  }
  if (is.null(conc_data_csv)) {
    stop("conc_data_csv is required and cannot be NULL")
  }

  # Check parameter types
  file_params <- list(
    birth_data_txt = birth_data_txt,
    model_prob_csv = model_prob_csv,
    conc_data_csv = conc_data_csv
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
#'
#' @param regression_formula Regression formula (character string or formula
#'   object)
#' @param exposure_level_col Column name for the exposure level factor
#' @param model_type Effects structure of the model: one of `"mixed"`,
#'   `"random"`, or `"fixed"`. Controls whether random-effects terms are
#'   required or forbidden. Defaults to `"mixed"`.
#' @return NULL (throws error if validation fails)
#' @keywords internal
validate_regression_formula <- function(
  regression_formula,
  exposure_level_col,
  model_type = "mixed"
) {
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

      pattern <- paste0("as\\.factor\\(", exposure_level_col,
                        "\\)|factor\\(", exposure_level_col, "\\)")
      if (!grepl(pattern, regression_formula)) {
        stop("regression_formula must include '", exposure_level_col,
             "' as a factor variable")
      }
    }

    # Check random-effects terms against the requested model_type
    validate_random_effects_terms(regression_formula, model_type)
  }
}

#' Validate Random-Effects Terms Against Model Type
#'
#' Ensures the formula's random-effects structure matches `model_type`:
#' fixed-effects models must not contain random terms, random- and mixed-effects
#' models must contain at least one, and random-effects models permit
#' intercept-only random terms (e.g. `(1 | group)`).
#'
#' @param regression_formula Regression formula (character or formula object)
#' @param model_type One of `"mixed"`, `"random"`, or `"fixed"`
#' @return NULL (throws error if validation fails)
#' @keywords internal
validate_random_effects_terms <- function(regression_formula, model_type) {
  formula_string <- if (is.character(regression_formula)) {
    regression_formula
  } else {
    paste(deparse(regression_formula), collapse = " ")
  }
  # A random-effects term looks like "( ... | ... )"
  random_term_pattern <- "\\([^|()]*\\|[^()]*\\)"
  has_random <- grepl(random_term_pattern, formula_string)

  if (model_type == "fixed") {
    if (has_random) {
      stop(
        "model_type = 'fixed' does not support random-effects terms such as ",
        "'(1 | group)'. Use model_type = 'random' or 'mixed' instead."
      )
    }
    return(invisible(NULL))
  }

  # model_type is "random" or "mixed": require at least one random-effects term
  if (!has_random) {
    stop(
      "model_type = '", model_type, "' requires at least one random-effects ",
      "term such as '(1 | group)'. Use model_type = 'fixed' for models ",
      "without random effects."
    )
  }

  if (model_type == "random") {
    # Random-intercept-only: the term to the left of '|' must be just "1"
    random_terms <- regmatches(
      formula_string,
      gregexpr(random_term_pattern, formula_string)
    )[[1]]
    slopes <- vapply(random_terms, function(term) {
      lhs <- trimws(sub("\\(([^|]*)\\|.*", "\\1", term))
      lhs != "1"
    }, logical(1))
    if (any(slopes)) {
      stop(
        "model_type = 'random' supports random intercepts only ",
        "(e.g. '(1 | group)'). Use model_type = 'mixed' for random slopes."
      )
    }
  }
  invisible(NULL)
}

#' Validate Character Vector Parameters
#' @keywords internal
validate_character_vectors <- function(
    targets, impute_vars, cat_label,
    mice_covs, columns_to_select,
    prob_cols) {

  char_vector_params <- list(
    targets = targets,
    impute_vars = impute_vars,
    cat_label = cat_label,
    mice_covs = mice_covs,
    columns_to_select = columns_to_select,
    prob_cols = prob_cols
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
    exposure_level_col, rucc_col, id_col,
    pop_well_col, record_id_col,
    conc_mean_col, conc_sd_col,
    pwell_col) {

  char_scalar_params <- list(
    exposure_level_col = exposure_level_col,
    rucc_col = rucc_col,
    id_col = id_col,
    pop_well_col = pop_well_col,
    record_id_col = record_id_col,
    conc_mean_col = conc_mean_col,
    conc_sd_col = conc_sd_col,
    pwell_col = pwell_col
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

#' Validate Exposure Category Parameters
#' @keywords internal
validate_exposure_categories <- function(cat_label, drop_cat_label_ref) {
  if (!is.null(drop_cat_label_ref)) {
    if (!is.character(drop_cat_label_ref)) {
      stop("drop_cat_label_ref must be a character vector")
    }

    # Check that drop_cat_label_ref values are in cat_label
    if (!is.null(cat_label)) {
      if (!all(drop_cat_label_ref %in% cat_label)) {
        stop("All values in drop_cat_label_ref must be present in ",
             "cat_label")
      }

      # Check that we don't drop all categories
      remaining_cats <- setdiff(cat_label, drop_cat_label_ref)
      if (length(remaining_cats) == 0) {
        stop("Cannot drop all exposure categories as reference levels")
      }
    }
  }
}

#' Validate Concentration Cutoffs Parameter
#' @keywords internal
validate_conc_cutoff_params <- function(conc_cutoffs, prob_cols) {
  if (!is.null(conc_cutoffs)) {
    if (!is.numeric(conc_cutoffs)) {
      stop("conc_cutoffs must be a numeric vector")
    }
    if (any(is.na(conc_cutoffs))) {
      stop("conc_cutoffs cannot contain NA values")
    }
    if (any(conc_cutoffs <= 0)) {
      stop("conc_cutoffs must contain only positive values")
    }
    if (is.unsorted(conc_cutoffs)) {
      stop("conc_cutoffs must be in ascending order")
    }
    if (length(unique(conc_cutoffs)) != length(conc_cutoffs)) {
      stop("conc_cutoffs cannot contain duplicate values")
    }

    # Check alignment with prob_cols if both are provided
    if (!is.null(prob_cols)) {
      conc_categories <- length(conc_cutoffs) + 1
      if (length(prob_cols) != conc_categories) {
        stop("prob_cols must have length equal to length(conc_cutoffs) + 1")
      }
    }
  }
}

#' Validate Probability Columns Parameter
#' @keywords internal
validate_prob_columns_param <- function(prob_columns) {
  if (!is.null(prob_columns)) {
    if (!is.character(prob_columns)) {
      stop("prob_columns must be a character vector")
    }
    if (length(prob_columns) == 0) {
      stop("prob_columns cannot be empty")
    }
    if (any(is.na(prob_columns)) || any(nchar(prob_columns) == 0)) {
      stop("prob_columns cannot contain NA or empty strings")
    }
    if (length(unique(prob_columns)) != length(prob_columns)) {
      stop("prob_columns cannot contain duplicate column names")
    }
  }
}

#' Validate File Contents and Column Existence
#' @keywords internal
validate_file_contents <- function(
    birth_data_txt, model_prob_csv,
    conc_data_csv, rucc_col, record_id_col,
    id_col, prob_cols, pop_well_col,
    conc_mean_col, conc_sd_col,
    conc_raw_col, pwell_col) {

  # Define required columns for each file
  birth_cols <- c(rucc_col, record_id_col)
  if (is.null(rucc_col)) {
    birth_cols <- c(record_id_col)
  }
  prob_required_cols <- c(id_col, prob_cols, pop_well_col)

  # For concentration columns, check what's available in the file first
  conc_cols_in_file <- tryCatch(
    get_file_columns(conc_data_csv),
    error = function(e) character(0)
  )

  # Concentration required columns depend on whether lognormal cols exist or conc col
  # Handle NULL values properly
  has_meanlog <- !is.null(conc_mean_col) && conc_mean_col %in% conc_cols_in_file
  has_sdlog <- !is.null(conc_sd_col) && conc_sd_col %in% conc_cols_in_file
  has_conc <- !is.null(conc_raw_col) && conc_raw_col %in% conc_cols_in_file

  # Build concentration required columns based on available data
  conc_cols <- c(id_col, pwell_col)
  if (has_meanlog) {
    conc_cols <- c(conc_cols, conc_mean_col)
  }
  if (has_sdlog) {
    conc_cols <- c(conc_cols, conc_sd_col)
  }
  if (!has_meanlog && has_conc) {
    conc_cols <- c(conc_cols, conc_raw_col)
  }

  data_list <- list(
    birth_data_txt = birth_data_txt,
    model_prob_csv = model_prob_csv,
    conc_data_csv = conc_data_csv
  )

  data_col_checks <- list(
    birth_data_txt = birth_cols,
    model_prob_csv = prob_required_cols,
    conc_data_csv = conc_cols
  )

  # Check column existence in each file
  for (data_param in names(data_list)) {
    file_path <- data_list[[data_param]]
    required_cols <- data_col_checks[[data_param]]

    cols_in_file <- get_file_columns(file_path)

    missing_cols <- setdiff(required_cols, cols_in_file)
    if (length(missing_cols) > 0) {
      stop("The following columns are missing in ", data_param, ": ",
           paste(missing_cols, collapse = ", "))
    }
  }

  # Validate concentration standard deviation column values only if column exists
  if (has_sdlog) {
    validate_conc_sd_values(conc_data_csv, conc_sd_col)
  }
}

#' Get Column Names from File
#' @keywords internal
get_file_columns <- function(file_path) {
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
  return(cols_in_file)
}

#' Validate Concentration Standard Deviation Values
#' @keywords internal
validate_conc_sd_values <- function(conc_data_csv, conc_sd_col) {
  if (!is.null(conc_sd_col)) {
    # Check if the column exists first
    cols_in_file <- get_file_columns(conc_data_csv)
    if (!conc_sd_col %in% cols_in_file) {
      # Column doesn't exist, will be auto-generated - skip validation
      return(invisible(NULL))
    }

    # Read the concentration data to check the column
    conc_data <- data.table::fread(
      conc_data_csv, header = TRUE,
      showProgress = FALSE,
      nThread = parallel::detectCores()
    )

    conc_sd_log_values <- conc_data[[conc_sd_col]]

    # Check for NA values first
    if (any(is.na(conc_sd_log_values))) {
      warning("conc_sd_col contains NA values. These rows will use ",
              "fallback probability calculations.")
    }

    # Check non-numeric values (excluding NAs)
    non_na_values <- conc_sd_log_values[!is.na(conc_sd_log_values)]
    if (length(non_na_values) > 0 && !all(is.finite(non_na_values))) {
      stop("conc_sd_col must contain only finite numeric values. ",
           "Check for non-numeric values or extreme outliers.")
    }

    if (length(non_na_values) > 0 && any(non_na_values <= 0)) {
      stop("conc_sd_col must contain only positive values ")
    }

    # Check large values
    if (length(non_na_values) > 0 && any(non_na_values > 10)) {
      warning("conc_sd_col contains very large values (> 10). This may ",
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

#' Validate Default Sdlog Parameter
#' @keywords internal
validate_default_sdlog <- function(default_sdlog) {
  if (!is.numeric(default_sdlog) || length(default_sdlog) != 1) {
    stop("default_sdlog must be a single numeric value")
  }
  if (default_sdlog <= 0) {
    stop("default_sdlog must be positive")
  }
  if (default_sdlog > 10) {
    warning("default_sdlog is very large (> 10). This may lead to ",
            "extreme probabilities in the lognormal distribution.")
  }
}

#' Validate Concentration Lognormal Input Parameters
#'
#' Validates that either lognormal columns exist or concentration column
#' is provided for automatic calculation.
#' @keywords internal
validate_conc_lognormal_input <- function(
    conc_data_csv,
    conc_mean_col,
    conc_sd_col,
    conc_raw_col
) {
  # Skip validation if concentration file path is NULL or doesn't exist
  if (is.null(conc_data_csv) || !file.exists(conc_data_csv)) {
    return(invisible(NULL))
  }

  # Try to get columns from file, skip if parsing fails
  cols_in_file <- tryCatch(
    get_file_columns(conc_data_csv),
    error = function(e) NULL
  )

  # If we couldn't parse the file, or file is not a proper CSV (e.g., test
  # files with single "test" content), skip this validation - other
  # validation functions will catch file format issues
  if (is.null(cols_in_file) || length(cols_in_file) < 2) {
    return(invisible(NULL))
  }

  has_meanlog <- conc_mean_col %in% cols_in_file
  has_sdlog <- conc_sd_col %in% cols_in_file

  # If meanlog doesn't exist and no concentration column provided, error
  if (!has_meanlog && is.null(conc_raw_col)) {
    stop(
      "Concentration lognormal meanlog column '", conc_mean_col, "' not found in data. ",
      "Either provide data with lognormal parameters, or specify ",
      "'conc_raw_col' to enable automatic conversion from raw ",
      "concentration values."
    )
  }

  # If concentration column specified but doesn't exist in file
  if (!is.null(conc_raw_col) && !conc_raw_col %in% cols_in_file) {
    stop(
      "Concentration column '", conc_raw_col,
      "' not found in the data."
    )
  }

  # Validate sdlog values if column exists
  if (has_sdlog) {
    validate_conc_sd_values(conc_data_csv, conc_sd_col)
  }
}
