#' Load and Process Arsenic Data
#'
#' Internal function to load USGS and EPA arsenic data, convert EPA lognormal
#' to multinomial, and create weighted combined probability matrix.
#'
#' This function is used within the sensitivity analysis workflow.
#' @param as_usgs_prob_csv Path to USGS probability CSV
#' @param as_epa_prob_csv Path to EPA probability CSV
#' @param as_conc_cols Column names for arsenic concentration probabilities
#' @param pop_well_col Column name for population weights
#' @param epa_as_mean_col Column name for EPA arsenic meanlog
#' @param epa_as_sd_col Column name for EPA arsenic sdlog
#' @param epa_as_conc_col Column name for raw EPA arsenic concentration (optional)
#' @param epa_default_sdlog Default sdlog value when not provided (default: 1.0)
#' @param epa_pwell_col Column name for private well percentage
#' @param id_col Column name for geographic identifier
#' @param epa_cutoffs A numeric vector of cutoffs for categorizing arsenic
#' levels based on EPA data
#' @param format_ids Logical. Whether to apply smart identifier formatting.
#'   If TRUE (default), numeric FIPS-like codes (1-5 digits) are zero-padded,
#'   while alphanumeric identifiers are preserved unchanged.
#' @return List containing processed arsenic data
#'
#' @keywords internal
load_and_process_arsenic_data <- function(
  as_usgs_prob_csv,
  as_epa_prob_csv,
  as_conc_cols,
  pop_well_col,
  epa_as_mean_col,
  epa_as_sd_col,
  epa_as_conc_col = NULL,
  epa_default_sdlog = 1.0,
  epa_pwell_col,
  id_col,
  epa_cutoffs,
  format_ids = TRUE
) {
  # Load USGS data (returns a list)
  usgs_data <- load_usgs_data(as_usgs_prob_csv, as_conc_cols, pop_well_col)

  # Load EPA data
  datmatx_epa <- data.table::fread(
    as_epa_prob_csv,
    header = TRUE,
    showProgress = FALSE,
    nThread = parallel::detectCores()
  )

  # Auto-calculate lognormal parameters if needed
  datmatx_epa <- calculate_epa_lognormal_params(
    datmatx_epa = datmatx_epa,
    epa_as_mean_col = epa_as_mean_col,
    epa_as_sd_col = epa_as_sd_col,
    epa_as_conc_col = epa_as_conc_col,
    epa_default_sdlog = epa_default_sdlog
  )

  # Validate data alignment - use the actual data frames
  validate_data_alignment(usgs_data$datmatx_usgs, datmatx_epa)
  n_usgs_categories <- length(as_conc_cols)

  # Auto-generate EPA cutoffs based on USGS column count
  if (is.null(epa_cutoffs)) {
    epa_cutoffs <- generate_default_epa_cutoffs(n_usgs_categories)

    # Validate that EPA categories match USGS categories
    validate_category_alignment(n_usgs_categories, epa_cutoffs)

    # Process cutoffs
    cutoff_config <- validate_epa_cutoffs(epa_cutoffs)

    message("Detected ", n_usgs_categories, " exposure categories\n")
    message("Using auto-generated cutoffs:",
            paste(cutoff_config$cutoffs), "\n")
  } else {
    # User provided cutoffs are already validated and aligned with USGS
    cutoff_config <- validate_epa_cutoffs(
      epa_cutoffs, validate_cutoffs = FALSE
    )
    message("Detected ", n_usgs_categories, " exposure categories\n")
    message("Using user-provided cutoffs:",
            paste(cutoff_config$cutoffs), "\n")
  }

  # Convert EPA lognormal to multinomial
  rasterprob_epa <- convert_epa_to_multinomial(
    datmatx_epa = datmatx_epa,
    epa_as_mean_col = epa_as_mean_col,
    epa_as_sd_col = epa_as_sd_col,
    cutoff_config = cutoff_config,
    return_prob_matrix = TRUE
  )

  # Create weighted combined probability matrix
  rasterprob_combined <- create_weighted_prob_matrix(
    rasterprob_usgs = usgs_data$rasterprob_usgs,  # Use the extracted matrix
    rasterprob_epa = rasterprob_epa,
    datmatx_epa = datmatx_epa,
    epa_pwell_col = epa_pwell_col
  )

  # replace negative values with zero
  rasterprob_combined[rasterprob_combined < 0] <- 0

  # Format identifiers if requested
  if (format_ids) {
    formatted_usgs <- format_ids_func(
      usgs_data$datmatx_usgs,
      id_col
    )
  } else {
    # Just ensure character type without formatting
    formatted_usgs <- usgs_data$datmatx_usgs
    formatted_usgs[[id_col]] <- as.character(formatted_usgs[[id_col]])
  }

  return(list(
    datmatx_usgs = formatted_usgs,
    rasterprob_combined = rasterprob_combined,
    cellpop = usgs_data$cellpop
  ))
}

#' Load USGS Data
#'
#' @param as_usgs_prob_csv Path to USGS CSV file
#' @param as_conc_cols Column names for concentration probabilities
#' @param pop_well_col Column name for population weights
#' @return List with USGS data, probability matrix, and cell populations
#'
#' @keywords internal
load_usgs_data <- function(as_usgs_prob_csv, as_conc_cols, pop_well_col) {
  datmatx_usgs <- data.table::fread(
    as_usgs_prob_csv,
    header = TRUE,
    showProgress = FALSE,
    nThread = parallel::detectCores()
  )

  # Validate required columns
  missing_cols <- setdiff(as_conc_cols, names(datmatx_usgs))
  if (length(missing_cols) > 0) {
    stop(
      "The following columns are missing in USGS data: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  # Extract probability matrix and population weights
  rasterprob_usgs <- as.matrix(datmatx_usgs[, as_conc_cols, with = FALSE])
  cellpop <- datmatx_usgs[[pop_well_col]]

  return(list(
    datmatx_usgs = datmatx_usgs,
    rasterprob_usgs = rasterprob_usgs,
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

#' Calculate EPA Lognormal Parameters
#'
#' Automatically calculates lognormal parameters (meanlog, sdlog) from raw
#' arsenic concentration values if the lognormal columns don't exist.
#'
#' @param datmatx_epa EPA data frame
#' @param epa_as_mean_col Column name for EPA arsenic meanlog
#' @param epa_as_sd_col Column name for EPA arsenic sdlog
#' @param epa_as_conc_col Column name for raw arsenic concentration (optional)
#' @param epa_default_sdlog Default sdlog value (default: 1.0)
#' @return Modified EPA data frame with lognormal parameters
#'
#' @keywords internal
calculate_epa_lognormal_params <- function(
  datmatx_epa,
  epa_as_mean_col,
  epa_as_sd_col,
  epa_as_conc_col = NULL,
  epa_default_sdlog = 1.0
) {
  # Check if lognormal columns already exist
 has_meanlog <- epa_as_mean_col %in% names(datmatx_epa)
 has_sdlog <- epa_as_sd_col %in% names(datmatx_epa)

  # If both columns exist, return as-is
  if (has_meanlog && has_sdlog) {
    return(datmatx_epa)
  }

  # If meanlog doesn't exist, we need to calculate it
  if (!has_meanlog) {
    if (is.null(epa_as_conc_col)) {
      stop(
        "EPA lognormal meanlog column '", epa_as_mean_col, "' not found. ",
        "Please provide 'epa_as_conc_col' to specify the raw arsenic ",
        "concentration column for automatic conversion."
      )
    }

    if (!epa_as_conc_col %in% names(datmatx_epa)) {
      stop(
        "EPA arsenic concentration column '", epa_as_conc_col,
        "' not found in the data."
      )
    }

    # Calculate meanlog from concentration: log(conc + 0.1) to avoid log(0)
    conc_values <- datmatx_epa[[epa_as_conc_col]]

    # Use data.table::set() to avoid shallow copy warnings
    data.table::set(
      datmatx_epa,
      j = epa_as_mean_col,
      value = ifelse(
        !is.na(conc_values) & conc_values > 0,
        log(conc_values + 0.1),
        NA_real_
      )
    )

    message(
      "Calculated EPA lognormal meanlog from '", epa_as_conc_col,
      "' column using: log(concentration + 0.1)"
    )
  }

  # If sdlog doesn't exist, use default value
  if (!has_sdlog) {
    # Use data.table::set() to avoid shallow copy warnings
    data.table::set(
      datmatx_epa,
      j = epa_as_sd_col,
      value = epa_default_sdlog
    )

    message(
      "Using default EPA lognormal sdlog = ", epa_default_sdlog
    )
  }

  return(datmatx_epa)
}

#' Validate Data Alignment
#'
#' @param datmatx_usgs USGS data frame
#' @param datmatx_epa EPA data frame
#'
#' @keywords internal
validate_data_alignment <- function(datmatx_usgs, datmatx_epa) {
  if (nrow(datmatx_usgs) != nrow(datmatx_epa)) {
    stop("USGS and EPA data files have different numbers of rows.")
  }
}

#' Convert EPA Lognormal to Multinomial
#'
#' @param datmatx_epa EPA data frame
#' @param epa_as_mean_col Column name for EPA arsenic meanlog
#' @param epa_as_sd_col Column name for EPA arsenic sdlog
#' @param cutoff_config List with cutoff values and number of categories
#' @param return_prob_matrix Whether to return the probability matrix. Default
#'   TRUE. If FALSE, returns the modified data frame.
#' @return Probability matrix or modified data frame
#'
#' @keywords internal
convert_epa_to_multinomial <- function(
  datmatx_epa,
  epa_as_mean_col,
  epa_as_sd_col,
  cutoff_config,
  return_prob_matrix = TRUE
) {
  cutoffs <- cutoff_config$cutoffs
  n_categories <- cutoff_config$n_categories

  # Create column names dynamically
  prob_cols <- paste0("EPA_C", seq_len(n_categories), "v2")

  # Initialize probability columns
  for (col in prob_cols) {
    datmatx_epa[[col]] <- 0
  }

  # If datamatx_epa is a data.frame, convert to data.table
  if (!data.table::is.data.table(datmatx_epa)) {
    datmatx_epa <- data.table::as.data.table(datmatx_epa, header = TRUE)
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
        meanlog = datmatx_epa[[epa_as_mean_col]],
        sdlog = datmatx_epa[[epa_as_sd_col]]
      )
    }

    prob_lower <- plnorm(
      lower_bound,
      meanlog = datmatx_epa[[epa_as_mean_col]],
      sdlog = datmatx_epa[[epa_as_sd_col]]
    )

    # Assign probability for this category
    datmatx_epa[[prob_cols[i]]] <- prob_upper - prob_lower
  }

  # After calculating probabilities in convert_epa_to_multinomial
  prob_cols <- paste0("EPA_C", seq_len(n_categories), "v2")
  row_sums <- rowSums(datmatx_epa[, prob_cols, with = FALSE])
  for (col in prob_cols) {
    # Normalize probabilities
    datmatx_epa[[col]] <- datmatx_epa[[col]] / row_sums
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
  #   datmatx_epa[[col]] <- datmatx_epa[[col]] / row_sums
  # }
  # } else {
  #   message("Row sums are within tolerance, no normalization needed")
  # }

  if (!return_prob_matrix) {
    return(datmatx_epa)
  }
  # Extract and return the probability matrix
  rasterprob_epa <- as.matrix(datmatx_epa[, prob_cols, with = FALSE])
  return(rasterprob_epa)
}

#' Create Weighted Probability Matrix
#'
#' @param rasterprob_usgs USGS probability matrix
#' @param rasterprob_epa EPA probability matrix
#' @param datmatx_epa EPA data frame
#' @param epa_pwell_col Column name for private well percentage
#' @return Combined weighted probability matrix
#'
#' @keywords internal
create_weighted_prob_matrix <- function(
  rasterprob_usgs,
  rasterprob_epa,
  datmatx_epa,
  epa_pwell_col
) {
  w_private <- datmatx_epa[[epa_pwell_col]] / 100
  w_public <- 1 - w_private

  # Element-wise multiplication and addition for final weighted probability
  result <- (rasterprob_usgs * w_private) + (rasterprob_epa * w_public)

  # Handle NA values in EPA probabilities by falling back to USGS probabilities
  # This happens when EPA meanlog/sdlog are NA
  na_rows <- which(apply(result, 1, function(x) any(is.na(x))))
  if (length(na_rows) > 0) {
    result[na_rows, ] <- rasterprob_usgs[na_rows, , drop = FALSE]
  }

  return(result)
}


#' Generate Default EPA Cutoffs Based on USGS Categories
#'
#' @param n_categories Number of USGS categories
#' @return Numeric vector of EPA cutoffs
#' 
#' @keywords internal
generate_default_epa_cutoffs <- function(n_categories) {

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

#' Generate Logarithmic EPA Cutoffs for Unusual Category Counts
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

#' Validate and Process EPA Cutoffs
#' @param epa_cutoffs Numeric vector of cutoff values
#' @param validate_cutoffs Whether to perform validation
#' @return List with processed cutoffs
#'
#' @keywords internal
validate_epa_cutoffs <- function(epa_cutoffs, validate_cutoffs = TRUE) {

  if (validate_cutoffs) {
    # Check cutoffs are numeric and sorted
    if (!is.numeric(epa_cutoffs)) {
      stop("epa_cutoffs must be numeric")
    }

    if (!identical(epa_cutoffs, sort(epa_cutoffs))) {
      stop("epa_cutoffs must be in ascending order")
    }

    # Check for negative values
    if (any(epa_cutoffs < 0)) {
      stop("epa_cutoffs cannot contain negative values")
    }
  }

  return(list(
    cutoffs = epa_cutoffs,
    n_categories = length(epa_cutoffs) + 1
  ))
}

#' Validate Category Alignment Between USGS and EPA
#'
#' @param n_usgs_categories Number of USGS categories
#' @param epa_cutoffs EPA cutoff values
#' @return TRUE if aligned, otherwise stops with error
#'
#' @keywords internal
validate_category_alignment <- function(n_usgs_categories, epa_cutoffs) {

  n_epa_categories <- length(epa_cutoffs) + 1

  if (n_usgs_categories != n_epa_categories) {
    stop(sprintf(
      paste0(
        "Category mismatch: USGS has %d categories but EPA cutoffs define ",
        "%d categories. Please adjust epa_cutoffs or usgs_columns to have ",
        "matching numbers of categories."
      ),
      n_usgs_categories, n_epa_categories
    ))
  }

  return(TRUE)
}
