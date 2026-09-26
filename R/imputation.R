#' Impute Exposure Levels
#'
#' Internal function to create multiple imputed datasets with 
#' exposure levels based on county-level probability distributions.
#'
#' @param datmatx_prob Data frame containing probability data with geographic
#'   identifiers
#' @param rasterprob_combined Matrix of combined probability distributions
#' @param cellpop Vector of cell population weights
#' @param births Data frame containing birth records
#' @param ndraws Number of imputed datasets to create
#' @param cat_label Character vector of exposure category labels
#' @param id_col Column name for geographic/participant identifier in probability data
#' @param record_id_col Column name for matching identifier in health data
#' @param exposure_level_col Column name for exposure level variable
#' @param apply_imputation_fallback Logical indicating whether to apply fallback
#' imputation for missing exposure levels. The fallback fills missing
#' values with the most common category in the respective dataset. If
#' enabled, this fallback will be applied to any missing exposure levels
#' after the initial imputation step. Otherwise, it will remove any rows
#' with missing exposure levels.
#' @return List of imputed datasets
#'
#' @keywords internal
impute_exposure <- function(
  datmatx_prob,
  rasterprob_combined,
  cellpop,
  births,
  ndraws,
  cat_label,
  id_col,
  record_id_col,
  exposure_level_col,
  apply_imputation_fallback
) {
  # Get unique counties and setup
  county_id_list <- na.omit(unique(datmatx_prob[[id_col]]))
  n_counties <- length(county_id_list)
  ncats <- length(cat_label)
  # Pre-allocate matrix with appropriate type
  exposure_temp <- matrix(NA_integer_, nrow = n_counties, ncol = ndraws)
  # Create efficient county lookup using split()
  county_lookup <- split(
    seq_len(nrow(datmatx_prob)),
    datmatx_prob[[id_col]]
  )
  # Vectorized probability calculation and sampling
  for (i in seq_len(n_counties)) {
    county_id <- county_id_list[i]
    cells_i <- county_lookup[[county_id]]
    # Optimized probability calculation
    if (length(cells_i) == 1) {
      # Single cell - no aggregation needed
      prob_values <- rasterprob_combined[cells_i, ]
    } else {
      # Multiple cells - use efficient weighted mean
      weights <- cellpop[cells_i]
      # Handle NA values in weights - use na.rm=TRUE for sum
      weight_sum <- sum(weights, na.rm = TRUE)
      if (!is.na(weight_sum) && weight_sum > 0) {
        # Replace NA weights with 0 for weighted calculation
        weights[is.na(weights)] <- 0
        prob_values <- colSums(
          rasterprob_combined[cells_i, , drop = FALSE] * weights
        ) / weight_sum
      } else {
        prob_values <- colMeans(rasterprob_combined[cells_i, , drop = FALSE],
                                na.rm = TRUE)
      }
    }
    # Normalize probabilities (handle numerical precision)
    prob_sum <- sum(prob_values, na.rm = TRUE)
    if (!is.na(prob_sum) && prob_sum > 0) {
      prob_values <- prob_values / prob_sum
      # Sample all draws at once for this county
      exposure_temp[i, ] <- sample.int(
        n = ncats,
        size = ndraws,
        replace = TRUE,
        prob = prob_values
      )
    }
  }
  # Report health records that cannot receive an imputed exposure: those whose
  # identifier has no exposure unit, and those in units without a valid
  # probability vector
  record_ids <- as.character(births[[record_id_col]])
  valid_units <- as.character(
    county_id_list[rowSums(!is.na(exposure_temp)) > 0]
  )
  no_unit <- !record_ids %in% as.character(county_id_list)
  no_prob <- !no_unit & !record_ids %in% valid_units
  n_missing <- sum(no_unit | no_prob)
  if (n_missing > 0) {
    warning(
      n_missing, " of ", length(record_ids), " health records (",
      sprintf("%.1f%%", 100 * n_missing / length(record_ids)),
      ") have no imputed exposure: ", sum(no_unit), " with a '",
      record_id_col, "' value absent from the exposure data and ",
      sum(no_prob), " in units without a valid probability vector. ",
      if (apply_imputation_fallback) {
        paste0(
          "They are assigned the most common category of each imputed ",
          "dataset (apply_imputation_fallback = TRUE)."
        )
      } else {
        "They are dropped (apply_imputation_fallback = FALSE)."
      },
      call. = FALSE
    )
  }
  # Pre-allocate results list
  imputed_datasets <- vector("list", ndraws)
  # Create imputed datasets efficiently
  for (j in seq_len(ndraws)) {
    # Create exposure assignments for this draw
    exposure_draw_df <- data.frame(
      countyID = county_id_list,
      temp_as_level = factor(
        exposure_temp[, j],
        levels = seq_len(ncats),
        labels = cat_label
      ),
      stringsAsFactors = FALSE
    )
    # Set proper column name
    names(exposure_draw_df)[2] <- exposure_level_col
    # Merge with births data
    imputed_datasets[[j]] <- merge(
      births,
      exposure_draw_df,
      by.x = record_id_col,
      by.y = "countyID",
      all.x = TRUE,
      sort = FALSE
    )
  }
  # Check if fallback imputation is needed
  if (apply_imputation_fallback) {
    if (
      any(
        sapply(
          imputed_datasets,
          function(x) any(is.na(x[[exposure_level_col]]))
        )
      )
    ) {
      message(
        "Some exposure levels are NA after initial imputation, ",
        "filling with most common category"
      )
      # fill missing exposure levels with most common category
      for (j in seq_along(imputed_datasets)) {
        # Find the most common category in the current dataset
        most_common_cat <- names(
          sort(
            table(imputed_datasets[[j]][[exposure_level_col]]),
            decreasing = TRUE
          )
        )[1]
        # Fill NAs with the most common category
        na_idx <- is.na(imputed_datasets[[j]][[exposure_level_col]])
        imputed_datasets[[j]][[exposure_level_col]][na_idx] <- most_common_cat
      }
      # Validate that exposure levels are now filled
    }
  } else {
    # Remove rows with missing exposure levels
    for (j in seq_along(imputed_datasets)) {
      imputed_datasets[[j]] <- imputed_datasets[[j]][
        !is.na(imputed_datasets[[j]][[exposure_level_col]]),
        ,
        drop = FALSE
      ]
    }
  }
  return(imputed_datasets)
}

#' Impute Additional Variables Using MICE
#'
#' Internal function to impute missing values in additional variables using
#' MICE.
#'
#' @param imputed_datasets List of datasets with exposure already
#'   imputed
#' @param impute_vars Character vector of variables to impute
#' @param births Original birth data frame
#' @param mice_covs Character vector of covariates for MICE
#' @param mice_m Number of MICE imputations
#' @param mice_maxit Maximum MICE iterations
#' @param mice_method MICE imputation method
#' @param seed Random seed
#' @param check_birth_cols Logical indicating whether to check if
#'   impute_vars exist in the birth data
#' @return List of datasets with additional variables imputed
#'
#' @keywords internal
impute_additional_variables <- function(
  imputed_datasets,
  impute_vars,
  births,
  mice_covs,
  mice_m,
  mice_maxit,
  mice_method,
  seed,
  check_birth_cols = TRUE
) {
  if (check_birth_cols) {
    # Validate input variables
    missing_vars <- setdiff(impute_vars, names(births))
    if (length(missing_vars) > 0) {
      stop(
        "Variables not found in birth data: ",
        paste(missing_vars, collapse = ", ")
      )
    }
  }

  # Process each variable to impute
  for (var in impute_vars) {
    for (j in seq_along(imputed_datasets)) {
      # Check if covariates exist in current dataset
      available_covs <- intersect(mice_covs, names(imputed_datasets[[j]]))
      if (length(available_covs) == 0) {
        warning(
          "No MICE covariates available for variable ", var,
          " in draw ", j
        )
        next
      }

      # Create subset for imputation - fix the data.table issue
      var_data <- as.data.frame(imputed_datasets[[j]])[c(var, available_covs)]
      # Only impute if there are missing values
      if (any(is.na(var_data[[var]]))) {
        # Run MICE imputation with error handling
        mice_result <- tryCatch({
          suppressWarnings(
            mice::mice(
              as.data.frame(var_data),
              m = mice_m,
              maxit = mice_maxit,
              method = mice_method,
              seed = seed,
              printFlag = FALSE
            )
          )
        }, error = function(e) {
          warning(
            "MICE failed for variable ", var,
            " in draw ", j, ": ", e$message
          )
          return(NULL)
        })
        # Extract imputed values if successful
        if (!is.null(mice_result)) {
          imputed_datasets[[j]][[var]] <- mice::complete(mice_result)[[var]]
        }
      }
    }
  }
  return(imputed_datasets)
}

#' Validate Imputed Datasets
#'
#' Internal function to validate that exposure imputation was
#' successful.
#'
#' @param imputed_datasets List of imputed datasets
#' @param exposure_level_col Column name for exposure level variable
#' @param cat_label Character vector of valid exposure category labels
#'
#' @keywords internal
validate_imputed_datasets <- function(
  imputed_datasets,
  exposure_level_col,
  cat_label
) {
  for (j in seq_along(imputed_datasets)) {
    # Check for missing exposure levels
    if (any(is.na(imputed_datasets[[j]][[exposure_level_col]]))) {
      stop("Missing values found in exposure level column in draw ", j)
    }
    # Check for invalid exposure levels
    dataset_levels <- as.character(imputed_datasets[[j]][[exposure_level_col]])
    if (!all(dataset_levels %in% cat_label)) {
      invalid_levels <- setdiff(dataset_levels, cat_label)
      stop(
        "Invalid exposure levels found in draw ", j, ": ",
        paste(invalid_levels, collapse = ", ")
      )
    }
  }
  message(
    "All ", length(imputed_datasets),
    " imputed datasets validated successfully"
  )
}
