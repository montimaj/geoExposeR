#' Load and Process Exposure Data
#'
#' Internal function to load probability and concentration data, convert
#' concentration lognormal to multinomial, and create weighted combined
#' probability matrix.
#'
#' This function is used within the sensitivity analysis workflow.
#' @param model_prob_csv Path to modeled probability data CSV
#' @param conc_data_csv Path to measured concentration data CSV
#' @param prob_cols Column names for exposure probability categories
#' @param pop_well_col Column name for population weights
#' @param conc_mean_col Column name for concentration lognormal meanlog
#' @param conc_sd_col Column name for concentration lognormal sdlog
#' @param conc_raw_col Column name for raw concentration values (optional)
#' @param default_sdlog Default sdlog value when not provided (default: 1.0)
#' @param pwell_col Column name for private well percentage
#' @param id_col Column name for geographic identifier
#' @param conc_cutoffs A numeric vector of cutoffs for categorizing exposure
#' levels based on concentration data
#' @param format_ids Logical. Whether to apply smart identifier formatting.
#'   If TRUE (default), numeric FIPS-like codes (1-5 digits) are zero-padded,
#'   while alphanumeric identifiers are preserved unchanged.
#' @return List containing processed exposure data
#'
#' @keywords internal
load_and_process_arsenic_data <- function(
  model_prob_csv,
  conc_data_csv,
  prob_cols,
  pop_well_col,
  conc_mean_col,
  conc_sd_col,
  conc_raw_col = NULL,
  default_sdlog = 1.0,
  pwell_col,
  id_col,
  conc_cutoffs,
  format_ids = TRUE
) {
  # Load probability data (returns a list)
  prob_input <- load_prob_data(model_prob_csv, prob_cols, pop_well_col)

  # Load concentration data
  datmatx_conc <- data.table::fread(
    conc_data_csv,
    header = TRUE,
    showProgress = FALSE,
    nThread = parallel::detectCores()
  )

  # Auto-calculate lognormal parameters if needed
  datmatx_conc <- calculate_conc_lognormal_params(
    datmatx_conc = datmatx_conc,
    conc_mean_col = conc_mean_col,
    conc_sd_col = conc_sd_col,
    conc_raw_col = conc_raw_col,
    default_sdlog = default_sdlog
  )

  # Validate data alignment - use the actual data frames
  validate_data_alignment(prob_input$datmatx_prob, datmatx_conc)
  n_prob_categories <- length(prob_cols)

  # Auto-generate concentration cutoffs based on probability column count
  if (is.null(conc_cutoffs)) {
    conc_cutoffs <- generate_default_conc_cutoffs(n_prob_categories)

    # Validate that concentration categories match probability categories
    validate_category_alignment(n_prob_categories, conc_cutoffs)

    # Process cutoffs
    cutoff_config <- validate_conc_cutoffs(conc_cutoffs)

    message("Detected ", n_prob_categories, " exposure categories\n")
    message("Using auto-generated cutoffs:",
            paste(cutoff_config$cutoffs), "\n")
  } else {
    # User provided cutoffs are already validated and aligned
    cutoff_config <- validate_conc_cutoffs(
      conc_cutoffs, validate_cutoffs = FALSE
    )
    message("Detected ", n_prob_categories, " exposure categories\n")
    message("Using user-provided cutoffs:",
            paste(cutoff_config$cutoffs), "\n")
  }

  # Convert concentration lognormal to multinomial
  rasterprob_conc <- convert_conc_to_multinomial(
    datmatx_conc = datmatx_conc,
    conc_mean_col = conc_mean_col,
    conc_sd_col = conc_sd_col,
    cutoff_config = cutoff_config,
    return_prob_matrix = TRUE
  )

  # Create weighted combined probability matrix
  rasterprob_combined <- create_weighted_prob_matrix(
    rasterprob_model = prob_input$rasterprob_model,  # Use the extracted matrix
    rasterprob_conc = rasterprob_conc,
    datmatx_conc = datmatx_conc,
    pwell_col = pwell_col
  )

  # replace negative values with zero
  rasterprob_combined[rasterprob_combined < 0] <- 0

  # Format identifiers if requested
  if (format_ids) {
    formatted_prob <- format_ids_func(
      prob_input$datmatx_prob,
      id_col
    )
  } else {
    # Just ensure character type without formatting
    formatted_prob <- prob_input$datmatx_prob
    formatted_prob[[id_col]] <- as.character(formatted_prob[[id_col]])
  }

  return(list(
    datmatx_prob = formatted_prob,
    rasterprob_combined = rasterprob_combined,
    cellpop = prob_input$cellpop
  ))
}

#' Load Probability Data
#'
#' @param model_prob_csv Path to modeled probability data CSV file
#' @param prob_cols Column names for exposure probability categories
#' @param pop_well_col Column name for population weights
#' @return List with probability data, probability matrix, and cell populations
#'
#' @keywords internal
load_prob_data <- function(model_prob_csv, prob_cols, pop_well_col) {
  datmatx_prob <- data.table::fread(
    model_prob_csv,
    header = TRUE,
    showProgress = FALSE,
    nThread = parallel::detectCores()
  )

  # Validate required columns
  missing_cols <- setdiff(prob_cols, names(datmatx_prob))
  if (length(missing_cols) > 0) {
    stop(
      "The following columns are missing in probability data: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  # Extract probability matrix and population weights
  rasterprob_model <- as.matrix(datmatx_prob[, prob_cols, with = FALSE])
  cellpop <- datmatx_prob[[pop_well_col]]

  return(list(
    datmatx_prob = datmatx_prob,
    rasterprob_model = rasterprob_model,
    cellpop = cellpop
  ))
}

#' Load and Process Birth Data
#'
#' Internal function to load birth data, select columns, create rural indicator,
#' and format county identifiers.
#'
#' @param birth_data_txt Path to birth data file
#' @param columns_to_select Column names to select (NULL for all)
#' @param rucc_col Column name for Rural-Urban Continuum Code
#' @param record_id_col Column name for birth county identifier
#' @param format_ids Logical. Whether to apply smart identifier formatting.
#'   If TRUE (default), numeric FIPS-like codes (1-5 digits) are zero-padded,
#'   while alphanumeric identifiers are preserved unchanged.
#' @return Processed birth data frame
#'
#' @keywords internal
load_and_process_birth_data <- function(
  birth_data_txt,
  columns_to_select,
  rucc_col,
  record_id_col,
  format_ids = TRUE
) {
  # Load birth data
  births <- data.table::fread(
    birth_data_txt,
    header = TRUE,
    nThread = parallel::detectCores()
  )

  # Select specified columns if provided
  if (!is.null(columns_to_select)) {
    births <- births %>% dplyr::select(dplyr::all_of(columns_to_select))
  }

  # Create rural indicator if RUCC column exists
  if (!is.null(rucc_col) && rucc_col %in% names(births)) {
    births <- births %>%
      dplyr::mutate(rural = .data[[rucc_col]] > 3)
  }

  # Format identifiers if requested
  if (format_ids) {
    births <- format_ids_func(births, record_id_col)
  } else {
    # Just ensure character type without formatting
    births[[record_id_col]] <- as.character(births[[record_id_col]])
  }

  return(births)
}

#' Calculate Concentration Lognormal Parameters
#'
#' Automatically calculates lognormal parameters (meanlog, sdlog) from raw
#' concentration values if the lognormal columns don't exist.
#'
#' @param datmatx_conc Concentration data frame
#' @param conc_mean_col Column name for concentration lognormal meanlog
#' @param conc_sd_col Column name for concentration lognormal sdlog
#' @param conc_raw_col Column name for raw concentration values (optional)
#' @param default_sdlog Default sdlog value (default: 1.0)
#' @return Modified concentration data frame with lognormal parameters
#'
#' @keywords internal
calculate_conc_lognormal_params <- function(
  datmatx_conc,
  conc_mean_col,
  conc_sd_col,
  conc_raw_col = NULL,
  default_sdlog = 1.0
) {
  # Check if lognormal columns already exist
 has_meanlog <- conc_mean_col %in% names(datmatx_conc)
 has_sdlog <- conc_sd_col %in% names(datmatx_conc)

  # If both columns exist, return as-is
  if (has_meanlog && has_sdlog) {
    return(datmatx_conc)
  }

  # If meanlog doesn't exist, we need to calculate it
  if (!has_meanlog) {
    if (is.null(conc_raw_col)) {
      stop(
        "Concentration lognormal meanlog column '", conc_mean_col,
        "' not found. Please provide 'conc_raw_col' to specify the raw ",
        "concentration column for automatic conversion."
      )
    }

    if (!conc_raw_col %in% names(datmatx_conc)) {
      stop(
        "Concentration column '", conc_raw_col,
        "' not found in the data."
      )
    }

    # Calculate meanlog from concentration: log(conc + 0.1) to avoid log(0)
    conc_values <- datmatx_conc[[conc_raw_col]]

    # Use data.table::set() to avoid shallow copy warnings
    data.table::set(
      datmatx_conc,
      j = conc_mean_col,
      value = ifelse(
        !is.na(conc_values) & conc_values > 0,
        log(conc_values + 0.1),
        NA_real_
      )
    )

    message(
      "Calculated lognormal meanlog from '", conc_raw_col,
      "' column using: log(concentration + 0.1)"
    )
  }

  # If sdlog doesn't exist, use default value
  if (!has_sdlog) {
    # Use data.table::set() to avoid shallow copy warnings
    data.table::set(
      datmatx_conc,
      j = conc_sd_col,
      value = default_sdlog
    )

    message(
      "Using default lognormal sdlog = ", default_sdlog
    )
  }

  return(datmatx_conc)
}

#' Validate Data Alignment
#'
#' @param datmatx_prob Probability data frame
#' @param datmatx_conc Concentration data frame
#'
#' @keywords internal
validate_data_alignment <- function(datmatx_prob, datmatx_conc) {
  if (nrow(datmatx_prob) != nrow(datmatx_conc)) {
    stop("Probability and concentration data files have different numbers of rows.")
  }
}

#' Convert Concentration Lognormal to Multinomial
#'
#' @param datmatx_conc Concentration data frame
#' @param conc_mean_col Column name for concentration lognormal meanlog
#' @param conc_sd_col Column name for concentration lognormal sdlog
#' @param cutoff_config List with cutoff values and number of categories
#' @param return_prob_matrix Whether to return the probability matrix. Default
#'   TRUE. If FALSE, returns the modified data frame.
#' @return Probability matrix or modified data frame
#'
#' @keywords internal
convert_conc_to_multinomial <- function(
  datmatx_conc,
  conc_mean_col,
  conc_sd_col,
  cutoff_config,
  return_prob_matrix = TRUE
) {
  cutoffs <- cutoff_config$cutoffs
  n_categories <- cutoff_config$n_categories

  # Create column names dynamically
  prob_cols <- paste0("conc_C", seq_len(n_categories), "v2")

  # Initialize probability columns
  for (col in prob_cols) {
    datmatx_conc[[col]] <- 0
  }

  # If data.frame, convert to data.table
  if (!data.table::is.data.table(datmatx_conc)) {
    datmatx_conc <- data.table::as.data.table(datmatx_conc, header = TRUE)
  }

  # Calculate probabilities for each category
  for (i in seq_len(n_categories)) {
    if (i == 1) {
      # First category: 0 to first cutoff
      lower_bound <- 0
      upper_bound <- cutoffs[1]
    } else if (i == n_categories) {
      # Last category: last cutoff to infinity
      lower_bound <- cutoffs[length(cutoffs)]
      upper_bound <- Inf
    } else {
      # Middle categories: between cutoffs
      lower_bound <- cutoffs[i - 1]
      upper_bound <- cutoffs[i]
    }

    # Calculate cumulative probabilities
    if (upper_bound == Inf) {
      prob_upper <- 1
    } else {
      prob_upper <- plnorm(
        upper_bound,
        meanlog = datmatx_conc[[conc_mean_col]],
        sdlog = datmatx_conc[[conc_sd_col]]
      )
    }

    prob_lower <- plnorm(
      lower_bound,
      meanlog = datmatx_conc[[conc_mean_col]],
      sdlog = datmatx_conc[[conc_sd_col]]
    )

    # Assign probability for this category
    datmatx_conc[[prob_cols[i]]] <- prob_upper - prob_lower
  }

  # After calculating probabilities in convert_conc_to_multinomial
  prob_cols <- paste0("conc_C", seq_len(n_categories), "v2")
  row_sums <- rowSums(datmatx_conc[, prob_cols, with = FALSE])
  for (col in prob_cols) {
    # Normalize probabilities
    datmatx_conc[[col]] <- datmatx_conc[[col]] / row_sums
  }

  # Check if normalization is needed
  # tolerance <- 1e-12
  # needs_normalization <- any(abs(row_sums - 1) > tolerance)

  # if (needs_normalization) {
  #   warning(
  #     "Row sums deviate from 1.0, applying normalization. Max deviation: ",
  #     max(abs(row_sums - 1))
  #   )

  #   # Normalize only if needed
  #   for (col in prob_cols) {
  #   # Normalize probabilities
  #   datmatx_conc[[col]] <- datmatx_conc[[col]] / row_sums
  # }
  # } else {
  #   message("Row sums are within tolerance, no normalization needed")
  # }

  if (!return_prob_matrix) {
    return(datmatx_conc)
  }
  # Extract and return the probability matrix
  rasterprob_conc <- as.matrix(datmatx_conc[, prob_cols, with = FALSE])
  return(rasterprob_conc)
}

#' Create Weighted Probability Matrix
#'
#' @param rasterprob_model Modeled probability matrix
#' @param rasterprob_conc Concentration-derived probability matrix
#' @param datmatx_conc Concentration data frame
#' @param pwell_col Column name for private well percentage
#' @return Combined weighted probability matrix
#'
#' @keywords internal
create_weighted_prob_matrix <- function(
  rasterprob_model,
  rasterprob_conc,
  datmatx_conc,
  pwell_col
) {
  w_private <- datmatx_conc[[pwell_col]] / 100
  w_public <- 1 - w_private

  # Element-wise multiplication and addition for final weighted probability
  result <- (rasterprob_model * w_private) + (rasterprob_conc * w_public)

  # Handle NA values in concentration probabilities by falling back to model probabilities
  # This happens when concentration meanlog/sdlog are NA
  na_rows <- which(apply(result, 1, function(x) any(is.na(x))))
  if (length(na_rows) > 0) {
    result[na_rows, ] <- rasterprob_model[na_rows, , drop = FALSE]
  }

  return(result)
}


#' Generate Default Concentration Cutoffs Based on Probability Categories
#'
#' @param n_categories Number of probability categories
#' @return Numeric vector of concentration cutoffs
#' 
#' @keywords internal
generate_default_conc_cutoffs <- function(n_categories) {

  # Define cutoff schemes for different numbers of categories
  cutoff_schemes <- list(
    "2" = c(5),                    # 2 categories: <5, >=5
    "3" = c(5, 10),                # 3 categories: <5, 5-10, >10 (default)
    "4" = c(2.5, 5, 10),           # 4 categories: <2.5, 2.5-5, 5-10, >10
    "5" = c(1, 5, 10, 20),         # 5 categories: <1, 1-5, 5-10, 10-20, >20
    "6" = c(1, 3, 5, 10, 20),      # 6 categories: <1, 1-3, 3-5, 5-10, 10-20, >20
    "7" = c(1, 2.5, 5, 7.5, 10, 20) # 7 categories: fine-grained
  )

  if (as.character(n_categories) %in% names(cutoff_schemes)) {
    cutoffs <- cutoff_schemes[[as.character(n_categories)]]
  } else {
    # For unusual numbers of categories, generate logarithmic spacing
    cutoffs <- generate_logarithmic_cutoffs(n_categories)
  }

  return(cutoffs)
}

#' Generate Logarithmic Cutoffs for Unusual Category Counts
#'
#' @param n_categories Number of categories
#' @return Numeric vector of logarithmically spaced cutoffs
#' 
#' @keywords internal
generate_logarithmic_cutoffs <- function(n_categories) {

  if (n_categories < 2) {
    stop("Must have at least 2 categories")
  }

  # Generate logarithmically spaced cutoffs from 0.5 to 50 ug/L
  min_cutoff <- 0.5
  max_cutoff <- 50

  # Generate n_categories-1 cutoffs (for n_categories bins)
  log_cutoffs <- seq(
    log(min_cutoff),
    log(max_cutoff),
    length.out = n_categories - 1
  )
  cutoffs <- exp(log_cutoffs)

  # Round to reasonable precision
  cutoffs <- round(cutoffs, 1)

  return(cutoffs)
}

#' Validate and Process Concentration Cutoffs
#' @param conc_cutoffs Numeric vector of cutoff values
#' @param validate_cutoffs Whether to perform validation
#' @return List with processed cutoffs
#'
#' @keywords internal
validate_conc_cutoffs <- function(conc_cutoffs, validate_cutoffs = TRUE) {

  if (validate_cutoffs) {
    # Check cutoffs are numeric and sorted
    if (!is.numeric(conc_cutoffs)) {
      stop("conc_cutoffs must be numeric")
    }

    if (!identical(conc_cutoffs, sort(conc_cutoffs))) {
      stop("conc_cutoffs must be in ascending order")
    }

    # Check for negative values
    if (any(conc_cutoffs < 0)) {
      stop("conc_cutoffs cannot contain negative values")
    }
  }

  return(list(
    cutoffs = conc_cutoffs,
    n_categories = length(conc_cutoffs) + 1
  ))
}

#' Validate Category Alignment Between Probability and Concentration Data
#'
#' @param n_prob_categories Number of probability categories
#' @param conc_cutoffs Concentration cutoff values
#' @return TRUE if aligned, otherwise stops with error
#'
#' @keywords internal
validate_category_alignment <- function(n_prob_categories, conc_cutoffs) {

  n_conc_categories <- length(conc_cutoffs) + 1

  if (n_prob_categories != n_conc_categories) {
    stop(sprintf(
      paste0(
        "Category mismatch: probability data has %d categories but ",
        "concentration cutoffs define %d categories. ",
        "Please adjust conc_cutoffs or prob_cols to have ",
        "matching numbers of categories."
      ),
      n_prob_categories, n_conc_categories
    ))
  }

  return(TRUE)
}
