#' Data Ingestion Functions for geoExposeR
#'
#' Functions for preparing input data from raw sources including:
#' GeoTIFF probability rasters, measured concentration data,
#' health outcome data formatting, and Census/population data integration.
#'
#' @name data-ingestion
#' @keywords internal
NULL

# Declare global variables for data.table/dplyr NSE
utils::globalVariables(c(
  ".", ".N", "FIPS", "mean_as", "sd_as", "conc_sdlog", "conc_meanlog",
  "prob_cat1", "prob_cat2", "prob_cat3", "total", "..cols_to_keep",
  "..prob_data_cols"
))

#' Prepare Modeled Probability Data from GeoTIFFs
#'
#' Extracts exposure probability estimates from GeoTIFF rasters
#' and aggregates them by geographic unit (e.g., county) to create the
#' CSV format required by geoExposeR.
#'
#' @param prob_rasters A named list of file paths to GeoTIFF rasters containing
#'   exceedance probabilities. Names should indicate thresholds (e.g.,
#'   `list(gt1 = "prob_gt1.tif", gt5 = "prob_gt5.tif", gt10 = "prob_gt10.tif")`).
#' @param boundaries An sf object containing polygon boundaries (e.g., counties)
#'   for spatial aggregation.
#' @param id_col Character. Column name in `boundaries` containing the
#'   geographic identifier. Supports any identifier format (FIPS codes, census
#'   tract IDs, ZIP codes, custom IDs). Default is "GEOID".
#' @param pop_data Optional data frame with population data. Should have columns
#'   for geographic ID and private well population.
#' @param pop_id_col Character. Column name in `pop_data` for geographic ID.
#'   Default is "GEOID".
#' @param pop_well_col Character. Column name in `pop_data` for private well
#'   population. Default is "private_well_pop".
#' @param extraction_method Character. Method for extracting raster values:
#'   "centroid" (faster) or "mean" (more accurate). Default is "mean".
#' @param cutoffs Numeric vector of length 2 specifying arsenic concentration
#'   cutoffs in ug/L for creating categories. Default is `c(5, 10)` for
#'   categories: <5, 5-10, >=10.
#' @param output_file Optional file path to save the resulting CSV.
#'
#' @return A data.table with columns for geographic ID and multinomial
#'   probabilities for each arsenic category.
#'
#' @details
#' The function converts exceedance probabilities (P(As > threshold)) to
#' multinomial probabilities that sum to 1:
#' \itemize{
#'   \item Category 1: P(As < lower_cutoff) = 1 - P(As > lower_cutoff)
#'   \item Category 2: P(lower <= As < upper) = P(As > lower) - P(As > upper)
#'   \item Category 3: P(As >= upper_cutoff) = P(As > upper_cutoff)
#' }
#'
#' @section Required Packages:
#' This function requires the `terra` and `sf` packages. For the "mean"
#' extraction method, `exactextractr` is also recommended for better accuracy.
#'
#' @examples
#' \dontrun{
#' # Load required packages
#' library(terra)
#' library(sf)
#'
#' # Define raster paths
#' rasters <- list(
#'   gt5 = "path/to/prob_gt5.tif",
#'   gt10 = "path/to/prob_gt10.tif"
#' )
#'
#' # Load county boundaries
#' counties <- st_read("path/to/counties.shp")
#'
#' # Extract and format probability data
#' prob_data <- prepare_prob_data(
#'   prob_rasters = rasters,
#'   boundaries = counties,
#'   id_col = "GEOID",
#'   cutoffs = c(5, 10),
#'   output_file = "prob_model_output.csv"
#' )
#' }
#'
#' @export
prepare_prob_data <- function(
  prob_rasters,
  boundaries,
  id_col = "GEOID",
  pop_data = NULL,
  pop_id_col = "GEOID",
  pop_well_col = "private_well_pop",
  extraction_method = c("mean", "centroid"),
  cutoffs = c(5, 10),
  output_file = NULL
) {
  # Validate inputs
 extraction_method <- match.arg(extraction_method)

  # Validate rasters input
  if (!is.list(prob_rasters) || length(prob_rasters) < 1) {
    stop("'prob_rasters' must be a non-empty list of file paths")
  }

  # Validate boundaries
  if (!inherits(boundaries, "sf")) {
    stop("'boundaries' must be an sf object")
  }

  if (!id_col %in% names(boundaries)) {
    stop(
      "Column '", id_col, "' not found in boundaries. ",
      "Available columns: ", paste(names(boundaries), collapse = ", ")
    )
  }

  # Validate cutoffs
  if (length(cutoffs) != 2 || !is.numeric(cutoffs) || cutoffs[1] >= cutoffs[2]) {
    stop("'cutoffs' must be a numeric vector of length 2 with cutoffs[1] < cutoffs[2]")
  }

  # Load rasters
  message("Loading GeoTIFF rasters...")
  raster_list <- lapply(names(prob_rasters), function(name) {
    path <- prob_rasters[[name]]
    if (!file.exists(path)) {
      stop("Raster file not found: ", path)
    }
    r <- terra::rast(path)
    message("  Loaded: ", name, " (", path, ")")
    return(r)
  })
  names(raster_list) <- names(prob_rasters)

  # Extract values based on method
  message("Extracting probabilities using '", extraction_method, "' method...")

  if (extraction_method == "centroid") {
    # Use centroids (faster) - suppress sf warning about constant attributes
    centroids <- suppressWarnings(sf::st_centroid(boundaries))
    centroids_vect <- terra::vect(centroids)

    extracted <- lapply(names(raster_list), function(name) {
      vals <- terra::extract(raster_list[[name]], centroids_vect)
      return(vals[, 2])  # Second column contains the values
    })
    names(extracted) <- names(raster_list)

  } else {
    # Use polygon means (more accurate)
    extracted <- lapply(names(raster_list), function(name) {
      exactextractr::exact_extract(raster_list[[name]], boundaries, fun = "mean")
    })
    names(extracted) <- names(raster_list)
  }

  # Create output data.table - use set operations to avoid shallow copy warnings
  message("Converting to multinomial probabilities...")
  result <- data.table::data.table(
    GEOID = boundaries[[id_col]]
  )

  # Add extracted probabilities using set() to avoid shallow copy warning
  for (name in names(extracted)) {
    data.table::set(result, j = paste0("prob_", name), value = extracted[[name]])
  }

  # Determine which columns to use for cutoffs
  # Expecting names like "gt5", "gt10" or similar
  lower_cutoff <- cutoffs[1]
  upper_cutoff <- cutoffs[2]

  # Find probability columns corresponding to cutoffs
  prob_cols <- names(result)[grepl("^prob_", names(result))]

  # Try to match cutoffs to probability columns
  lower_col <- NULL
  upper_col <- NULL

  for (col in prob_cols) {
    if (grepl(as.character(lower_cutoff), col)) {
      lower_col <- col
    }
    if (grepl(as.character(upper_cutoff), col)) {
      upper_col <- col
    }
  }

  if (is.null(lower_col) || is.null(upper_col)) {
    warning(
      "Could not automatically match cutoffs to probability columns. ",
      "Using first two probability columns in order."
    )
    if (length(prob_cols) >= 2) {
      lower_col <- prob_cols[1]
      upper_col <- prob_cols[2]
    } else {
      stop("Need at least 2 probability rasters for 3-category model")
    }
  }

  # Calculate multinomial probabilities
  # Category 1: As < lower_cutoff
  # Category 2: lower_cutoff <= As < upper_cutoff
  # Category 3: As >= upper_cutoff
  result[, prob_cat1 := pmax(0, 1 - get(lower_col))]
  result[, prob_cat2 := pmax(0, get(lower_col) - get(upper_col))]
  result[, prob_cat3 := pmax(0, get(upper_col))]

  # Normalize to ensure sum = 1
  result[, total := prob_cat1 + prob_cat2 + prob_cat3]
  result[total > 0, `:=`(
    prob_cat1 = prob_cat1 / total,
    prob_cat2 = prob_cat2 / total,
    prob_cat3 = prob_cat3 / total
  )]
  result[, total := NULL]

  # Rename to standard column names
  data.table::setnames(
    result,
    c("prob_cat1", "prob_cat2", "prob_cat3"),
    c("prob_C1", "prob_C2", "prob_C3")
  )

  # Remove intermediate probability columns
  cols_to_remove <- prob_cols
  result[, (cols_to_remove) := NULL]

  # Add population data if provided
  if (!is.null(pop_data)) {
    message("Adding population data...")
    pop_dt <- data.table::as.data.table(pop_data)

    if (!pop_id_col %in% names(pop_dt)) {
      warning(
        "Population GEOID column '", pop_id_col, "' not found. ",
        "Skipping population data."
      )
    } else if (!pop_well_col %in% names(pop_dt)) {
      warning(
        "Population well column '", pop_well_col, "' not found. ",
        "Skipping population data."
      )
    } else {
      # Merge population data
      pop_subset <- pop_dt[, c(pop_id_col, pop_well_col), with = FALSE]
      data.table::setnames(pop_subset, pop_id_col, "GEOID")
      data.table::setnames(pop_subset, pop_well_col, "Wells_2010")
      result <- merge(result, pop_subset, by = "GEOID", all.x = TRUE)
    }
  }

  # Rename GEOID to GEOID10 for compatibility
  data.table::setnames(result, "GEOID", "GEOID10")

  # Save to file if requested
  if (!is.null(output_file)) {
    message("Saving to: ", output_file)
    data.table::fwrite(result, output_file)
  }

  message("Done! Created probability data with ", nrow(result), " geographic units.")
  return(result)
}


#' Prepare Measured Concentration Data with Lognormal Parameters
#'
#' Processes measured concentration data and converts to lognormal
#' distribution parameters required by geoExposeR.
#'
#' @param conc_data A data frame containing concentration measurements. Can be
#'   raw measurements or summary statistics.
#' @param id_col Character. Column name containing geographic identifiers.
#'   Default is "county_fips".
#' @param concentration_col Character. Column name containing arsenic
#'   concentrations in ug/L. Default is "arsenic_concentration".
#' @param population_col Optional character. Column name for population served
#'   by each water system, used for weighted averaging.
#' @param private_well_pct_col Optional character. Column name containing
#'   private well percentage. If NULL, must be provided separately or merged
#'   from another source.
#' @param input_type Character. Type of input data: "raw" for individual
#'   measurements, "summary" for pre-aggregated statistics. Default is "raw".
#' @param mean_col Character. For summary input, column name for mean
#'   concentration. Default is "mean_concentration".
#' @param sd_col Character. For summary input, column name for standard
#'   deviation. Default is "sd_concentration".
#' @param default_sdlog Numeric. Default sdlog value to use when SD cannot be
#'   calculated. Default is 1.0.
#' @param min_concentration Numeric. Minimum concentration value to avoid log(0).
#'   Default is 0.1.
#' @param output_file Optional file path to save the resulting CSV.
#'
#' @return A data.table with columns for geographic ID, conc_meanlog
#'   (lognormal mean), conc_sdlog (lognormal SD), and optionally
#'   PWELL_private_pct.
#'
#' @details
#' The function converts arithmetic mean and standard deviation to lognormal
#' parameters using the following formulas:
#' \deqn{meanlog = \log(\mu^2 / \sqrt{\sigma^2 + \mu^2})}
#' \deqn{sdlog = \sqrt{\log(1 + \sigma^2/\mu^2)}}
#'
#' where \eqn{\mu} is the arithmetic mean and \eqn{\sigma} is the standard
#' deviation of the concentrations.
#'
#' @examples
#' \dontrun{
#' # Example with raw concentration data
#' conc_raw <- data.frame(
#'   county_fips = c("06001", "06001", "06003", "06003"),
#'   arsenic_concentration = c(2.5, 3.1, 8.2, 7.5),
#'   population_served = c(50000, 25000, 10000, 15000)
#' )
#'
#' conc_processed <- prepare_conc_data(
#'   conc_data = conc_raw,
#'   id_col = "county_fips",
#'   concentration_col = "arsenic_concentration",
#'   population_col = "population_served",
#'   input_type = "raw",
#'   output_file = "conc_lognormal.csv"
#' )
#' }
#'
#' @export
prepare_conc_data <- function(
  conc_data,
  id_col = "county_fips",
  concentration_col = "arsenic_concentration",
  population_col = NULL,
  private_well_pct_col = NULL,
  input_type = c("raw", "summary"),
  mean_col = "mean_concentration",
  sd_col = "sd_concentration",
  default_sdlog = 1.0,
  min_concentration = 0.1,
  output_file = NULL
) {
  input_type <- match.arg(input_type)

  # Convert to data.table
  dt <- data.table::as.data.table(conc_data)

  # Validate geoid column
  if (!id_col %in% names(dt)) {
    stop(
      "Geographic ID column '", id_col, "' not found. ",
      "Available columns: ", paste(names(dt), collapse = ", ")
    )
  }

  if (input_type == "raw") {
    # Validate concentration column
    if (!concentration_col %in% names(dt)) {
      stop(
        "Concentration column '", concentration_col, "' not found. ",
        "Available columns: ", paste(names(dt), collapse = ", ")
      )
    }

    message("Processing raw concentration measurements...")

    # Filter positive concentrations
    dt <- dt[get(concentration_col) > 0]

    if (nrow(dt) == 0) {
      stop("No positive concentration values found")
    }

    # Aggregate by geographic unit
    if (!is.null(population_col) && population_col %in% names(dt)) {
      # Population-weighted aggregation
      message("Using population-weighted averaging...")
      result <- dt[, .(
        mean_as = stats::weighted.mean(
          get(concentration_col),
          get(population_col),
          na.rm = TRUE
        ),
        sd_as = sqrt(Hmisc::wtd.var(
          get(concentration_col),
          get(population_col),
          na.rm = TRUE
        )),
        n_systems = .N,
        total_pop = sum(get(population_col), na.rm = TRUE)
      ), by = id_col]
    } else {
      # Simple averaging
      result <- dt[, .(
        mean_as = mean(get(concentration_col), na.rm = TRUE),
        sd_as = stats::sd(get(concentration_col), na.rm = TRUE),
        n_systems = .N
      ), by = id_col]
    }

  } else {
    # Summary input type
    message("Processing summary concentration statistics...")

    # Validate required columns
    if (!mean_col %in% names(dt)) {
      stop(
        "Mean column '", mean_col, "' not found. ",
        "Available columns: ", paste(names(dt), collapse = ", ")
      )
    }

    result <- data.table::copy(dt)
    result[, mean_as := get(mean_col)]

    if (sd_col %in% names(dt)) {
      result[, sd_as := get(sd_col)]
    } else {
      message("SD column not found. Using default sdlog = ", default_sdlog)
      result[, sd_as := NA_real_]
    }
  }

  # Convert to lognormal parameters
  message("Converting to lognormal parameters...")

  result[, `:=`(
    # Ensure minimum concentration
    mean_as = pmax(mean_as, min_concentration),
    # Handle missing/zero SD
    sd_as = ifelse(is.na(sd_as) | sd_as <= 0, mean_as * 0.5, sd_as)
  )]

  # Calculate lognormal parameters
  # meanlog = log(mu^2 / sqrt(sigma^2 + mu^2))
  # sdlog = sqrt(log(1 + sigma^2/mu^2))
  result[, `:=`(
    conc_meanlog = log(mean_as^2 / sqrt(sd_as^2 + mean_as^2)),
    conc_sdlog = sqrt(log(1 + (sd_as^2 / mean_as^2)))
  )]

  # Handle edge cases
  result[is.na(conc_sdlog) | !is.finite(conc_sdlog), conc_sdlog := default_sdlog]
  result[is.na(conc_meanlog) | !is.finite(conc_meanlog),
         conc_meanlog := log(min_concentration)]

  # Add private well percentage if provided
  if (!is.null(private_well_pct_col) && private_well_pct_col %in% names(dt)) {
    # Get first value per geographic unit (should be same for all rows)
    pwell_data <- dt[, .(PWELL_private_pct = get(private_well_pct_col)[1]),
                     by = id_col]
    result <- merge(result, pwell_data, by = id_col, all.x = TRUE)
  }

  # Select and rename columns
  cols_to_keep <- c(id_col, "conc_meanlog", "conc_sdlog")
  if ("PWELL_private_pct" %in% names(result)) {
    cols_to_keep <- c(cols_to_keep, "PWELL_private_pct")
  }
  if ("n_systems" %in% names(result)) {
    cols_to_keep <- c(cols_to_keep, "n_systems")
  }

  result <- result[, ..cols_to_keep]

  # Rename geoid column to standard name
  data.table::setnames(result, id_col, "GEOID10")

  # Save to file if requested
  if (!is.null(output_file)) {
    message("Saving to: ", output_file)
    data.table::fwrite(result, output_file)
  }

  message("Done! Created concentration data with ", nrow(result), " geographic units.")
  return(result)
}


#' Format Health Outcome Data for geoExposeR
#'
#' Formats and validates health outcome data to ensure compatibility with
#' geoExposeR's requirements.
#'
#' @param health_data A data frame containing health outcomes and covariates.
#' @param id_col Character. Column name containing geographic identifiers
#'   (e.g., FIPS codes). Default is "FIPS".
#' @param outcome_cols Character vector. Column names of health outcome
#'   variables to include.
#' @param covariate_cols Optional character vector. Column names of covariate
#'   variables to include.
#' @param geoid_digits Integer. Number of digits for FIPS code formatting.
#'   Default is 5 (county level).
#' @param validate Logical. Whether to perform validation checks. Default TRUE.
#' @param output_file Optional file path to save the resulting file.
#' @param output_format Character. Output format: "txt" (tab-delimited) or
#'   "csv". Default is "txt".
#'
#' @return A data.table with properly formatted health outcome data.
#'
#' @details
#' The function performs the following operations:
#' \itemize{
#'   \item Formats FIPS codes to specified number of digits with zero-padding
#'   \item Validates that outcome columns exist and contain valid data
#'   \item Standardizes missing value representation to NA
#'   \item Optionally validates data ranges for common health outcomes
#' }
#'
#' @examples
#' \dontrun{
#' health <- data.frame(
#'   FIPS = c(6001, 6003, 6005),
#'   BWT = c(3200, 3150, NA),
#'   OEGEST = c(39.2, 38.8, 40.1),
#'   MAGE = c(28, 32, 25)
#' )
#'
#' formatted <- format_health_data(
#'   health_data = health,
#'   id_col = "FIPS",
#'   outcome_cols = c("BWT", "OEGEST"),
#'   covariate_cols = c("MAGE"),
#'   output_file = "health_outcomes.txt"
#' )
#' }
#'
#' @export
format_health_data <- function(
  health_data,
  id_col = "FIPS",
  outcome_cols,
  covariate_cols = NULL,
  geoid_digits = 5,
  validate = TRUE,
  output_file = NULL,
  output_format = c("txt", "csv")
) {
  output_format <- match.arg(output_format)

  # Convert to data.table
  dt <- data.table::as.data.table(health_data)

  # Validate geoid column
  if (!id_col %in% names(dt)) {
    stop(
      "Geographic ID column '", id_col, "' not found. ",
      "Available columns: ", paste(names(dt), collapse = ", ")
    )
  }

  # Validate outcome columns
  missing_outcomes <- outcome_cols[!outcome_cols %in% names(dt)]
  if (length(missing_outcomes) > 0) {
    stop(
      "Outcome columns not found: ", paste(missing_outcomes, collapse = ", "),
      ". Available columns: ", paste(names(dt), collapse = ", ")
    )
  }

  # Validate covariate columns
  if (!is.null(covariate_cols)) {
    missing_covs <- covariate_cols[!covariate_cols %in% names(dt)]
    if (length(missing_covs) > 0) {
      warning(
        "Covariate columns not found (will be skipped): ",
        paste(missing_covs, collapse = ", ")
      )
      covariate_cols <- covariate_cols[covariate_cols %in% names(dt)]
    }
  }

  message("Formatting health outcome data...")

  # Format FIPS codes
  dt[, FIPS := sprintf(paste0("%0", geoid_digits, "d"), as.integer(get(id_col)))]

  # Select columns to keep
  cols_to_keep <- c("FIPS", outcome_cols)
  if (!is.null(covariate_cols) && length(covariate_cols) > 0) {
    cols_to_keep <- c(cols_to_keep, covariate_cols)
  }

  result <- dt[, ..cols_to_keep]

  # Validation checks
  if (validate) {
    message("Performing validation checks...")

    # Check for all-NA outcomes
    for (col in outcome_cols) {
      if (all(is.na(result[[col]]))) {
        warning("Outcome column '", col, "' contains all NA values")
      }
    }

    # Check common outcome ranges
    common_outcomes <- list(
      BWT = c(200, 6000),      # Birth weight in grams
      OEGEST = c(20, 45),      # Gestational age in weeks
      SBP = c(60, 250),        # Systolic blood pressure
      DBP = c(30, 150)         # Diastolic blood pressure
    )

    for (col in outcome_cols) {
      if (col %in% names(common_outcomes)) {
        range_check <- common_outcomes[[col]]
        vals <- result[[col]][!is.na(result[[col]])]
        if (length(vals) > 0) {
          if (min(vals) < range_check[1] || max(vals) > range_check[2]) {
            warning(
              "Outcome '", col, "' has values outside expected range [",
              range_check[1], ", ", range_check[2], "]. ",
              "Actual range: [", min(vals), ", ", max(vals), "]"
            )
          }
        }
      }
    }

    # Check for duplicate FIPS codes (warn if present)
    n_unique <- length(unique(result$FIPS))
    if (n_unique < nrow(result)) {
      message(
        "Note: ", nrow(result) - n_unique, " duplicate FIPS codes found. ",
        "This is expected for individual-level data."
      )
    }

    # Report missing data summary
    for (col in c(outcome_cols, covariate_cols)) {
      n_missing <- sum(is.na(result[[col]]))
      if (n_missing > 0) {
        pct_missing <- round(100 * n_missing / nrow(result), 1)
        message("  Column '", col, "': ", n_missing, " missing values (", pct_missing, "%)")
      }
    }
  }

  # Save to file if requested
  if (!is.null(output_file)) {
    message("Saving to: ", output_file)
    if (output_format == "txt") {
      data.table::fwrite(result, output_file, sep = "\t")
    } else {
      data.table::fwrite(result, output_file)
    }
  }

  message(
    "Done! Formatted health data with ", nrow(result), " records and ",
    length(outcome_cols), " outcomes."
  )
  return(result)
}


#' Add Private Well Population Percentages
#'
#' Merges private well percentage data into probability or concentration datasets.
#'
#' @param target_data A data.table or data frame to add well percentages to.
#' @param well_data A data frame containing private well percentages by
#'   geographic unit.
#' @param target_id_col Character. Column name for geographic ID in
#'   target_data. Default is "GEOID10".
#' @param well_id_col Character. Column name for geographic ID in
#'   well_data. Default is "GEOID".
#' @param well_pct_col Character. Column name for private well percentage.
#'   Default is "private_well_pct".
#' @param output_col Character. Name for the output column. Default is
#'   "PWELL_private_pct".
#'
#' @return The target_data with an added column for private well percentage.
#'
#' @examples
#' \dontrun{
#' # Add well percentages to concentration data
#' census_wells <- data.frame(
#'   GEOID = c("06001", "06003"),
#'   private_well_pct = c(15.2, 28.7)
#' )
#'
#' conc_with_wells <- add_well_percentages(
#'   target_data = conc_data,
#'   well_data = census_wells,
#'   well_pct_col = "private_well_pct"
#' )
#' }
#'
#' @export
add_well_percentages <- function(
  target_data,
  well_data,
  target_id_col = "GEOID10",
  well_id_col = "GEOID",
  well_pct_col = "private_well_pct",
  output_col = "PWELL_private_pct"
) {
  # Convert to data.tables
  target_dt <- data.table::as.data.table(target_data)
  well_dt <- data.table::as.data.table(well_data)

  # Validate columns
  if (!target_id_col %in% names(target_dt)) {
    stop("Column '", target_id_col, "' not found in target_data")
  }
  if (!well_id_col %in% names(well_dt)) {
    stop("Column '", well_id_col, "' not found in well_data")
  }
  if (!well_pct_col %in% names(well_dt)) {
    stop("Column '", well_pct_col, "' not found in well_data")
  }

  # Prepare well data for merge
  well_subset <- well_dt[, c(well_id_col, well_pct_col), with = FALSE]
  data.table::setnames(well_subset, well_id_col, target_id_col)
  data.table::setnames(well_subset, well_pct_col, output_col)

  # Merge
  result <- merge(target_dt, well_subset, by = target_id_col, all.x = TRUE)

  # Report merge results
  n_matched <- sum(!is.na(result[[output_col]]))
  n_total <- nrow(result)
  message(
    "Merged well percentages: ", n_matched, "/", n_total,
    " (", round(100 * n_matched / n_total, 1), "%) records matched"
  )

  return(result)
}


#' Subset Data by Geographic Area
#'
#' Filters prepared datasets to a specific geographic subset (state, region,
#' or custom list of geographic IDs).
#'
#' @param data A data.table or data frame to subset.
#' @param id_col Character. Column name for geographic IDs. Default is "GEOID10".
#' @param state_fips Optional character vector. Two-digit state FIPS codes to
#'   include. Data will be filtered to counties within these states.
#' @param geoid_list Optional character vector. Specific geographic IDs to
#'   include.
#' @param geoid_pattern Optional character. Regular expression pattern to match
#'   geographic IDs.
#'
#' @return A subset of the input data matching the specified criteria.
#'
#' @examples
#' \dontrun{
#' # Subset to California counties (state FIPS = "06")
#' ca_data <- subset_by_geography(
#'   data = prob_data,
#'   state_fips = "06"
#' )
#'
#' # Subset to specific counties
#' selected <- subset_by_geography(
#'   data = prob_data,
#'   geoid_list = c("06001", "06003", "06005")
#' )
#' }
#'
#' @export
subset_by_geography <- function(
  data,
  id_col = "GEOID10",
  state_fips = NULL,
  geoid_list = NULL,
  geoid_pattern = NULL
) {
  dt <- data.table::as.data.table(data)

  if (!id_col %in% names(dt)) {
    stop("Column '", id_col, "' not found in data")
  }

  original_rows <- nrow(dt)

  if (!is.null(state_fips)) {
    # Filter by state FIPS (first 2 digits of county FIPS)
    state_pattern <- paste0("^(", paste(state_fips, collapse = "|"), ")")
    dt <- dt[grepl(state_pattern, get(id_col))]
    message("Filtered to states: ", paste(state_fips, collapse = ", "))
  }

  if (!is.null(geoid_list)) {
    dt <- dt[get(id_col) %in% geoid_list]
    message("Filtered to ", length(geoid_list), " specified geographic IDs")
  }

  if (!is.null(geoid_pattern)) {
    dt <- dt[grepl(geoid_pattern, get(id_col))]
    message("Filtered by pattern: ", geoid_pattern)
  }

  message(
    "Subset complete: ", nrow(dt), "/", original_rows,
    " records retained (", round(100 * nrow(dt) / original_rows, 1), "%)"
  )

  return(dt)
}


#' Validate Prepared Input Files
#'
#' Performs comprehensive validation of prepared input files to ensure they
#' are ready for use with geoExposeR.
#'
#' @param prob_file Path to the prepared probability data CSV file.
#' @param conc_file Path to the prepared concentration lognormal parameters CSV file.
#' @param health_file Path to the prepared health outcome file.
#' @param prob_id_col Character. GEOID column in probability data. Default "GEOID10".
#' @param conc_id_col Character. GEOID column in concentration data. Default "GEOID10".
#' @param health_id_col Character. GEOID column in health data. Default "FIPS".
#' @param prob_data_cols Character vector. Probability column names in probability data.
#' @param fix_issues Logical. Whether to attempt automatic fixes for common
#'   issues. Default is FALSE.
#'
#' @return A list with validation results including:
#'   \itemize{
#'     \item valid: Logical indicating if all checks passed
#'     \item issues: Character vector of identified issues
#'     \item warnings: Character vector of warnings (non-critical)
#'     \item summary: Summary statistics for each dataset
#'   }
#'
#' @examples
#' \dontrun{
#' validation <- validate_prepared_inputs(
#'   prob_file = "prob_model.csv",
#'   conc_file = "conc_params.csv",
#'   health_file = "health_outcomes.txt"
#' )
#'
#' if (validation$valid) {
#'   message("All validations passed!")
#' } else {
#'   cat("Issues found:\n")
#'   cat(paste(validation$issues, collapse = "\n"))
#' }
#' }
#'
#' @export
validate_prepared_inputs <- function(
  prob_file,
  conc_file,
  health_file,
  prob_id_col = "GEOID10",
  conc_id_col = "GEOID10",
  health_id_col = "FIPS",
  prob_data_cols = c("prob_C1", "prob_C2", "prob_C3"),
  fix_issues = FALSE
) {
  issues <- character(0)
  warnings <- character(0)
  summary_list <- list()

  message("Validating prepared input files...\n")

  # ===== Probability Data Validation =====
  message("Checking probability data...")
  if (!file.exists(prob_file)) {
    issues <- c(issues, paste("Probability file not found:", prob_file))
  } else {
    prob_dt <- data.table::fread(prob_file)
    summary_list$prob <- list(
      n_rows = nrow(prob_dt),
      n_cols = ncol(prob_dt),
      columns = names(prob_dt)
    )

    # Check GEOID column
    if (!prob_id_col %in% names(prob_dt)) {
      issues <- c(issues, paste("Probability data missing GEOID column:", prob_id_col))
    }

    # Check probability columns
    missing_prob <- prob_data_cols[!prob_data_cols %in% names(prob_dt)]
    if (length(missing_prob) > 0) {
      issues <- c(issues, paste(
        "Probability data missing probability columns:",
        paste(missing_prob, collapse = ", ")
      ))
    } else {
      # Check probability sums
      prob_sums <- rowSums(prob_dt[, ..prob_data_cols], na.rm = TRUE)
      bad_sums <- sum(abs(prob_sums - 1) > 0.01, na.rm = TRUE)
      if (bad_sums > 0) {
        warnings <- c(warnings, paste(
          "Probability data:", bad_sums, "rows have probability sums != 1"
        ))
      }

      # Check for negative probabilities
      has_negative <- any(prob_dt[, ..prob_data_cols] < 0, na.rm = TRUE)
      if (has_negative) {
        issues <- c(issues, "Probability data: Negative probabilities found")
      }
    }

    message("  ", nrow(prob_dt), " geographic units")
  }

  # ===== Concentration Data Validation =====
  message("Checking concentration data...")
  if (!file.exists(conc_file)) {
    issues <- c(issues, paste("Concentration file not found:", conc_file))
  } else {
    conc_dt <- data.table::fread(conc_file)
    summary_list$conc <- list(
      n_rows = nrow(conc_dt),
      n_cols = ncol(conc_dt),
      columns = names(conc_dt)
    )

    # Check GEOID column
    if (!conc_id_col %in% names(conc_dt)) {
      issues <- c(issues, paste("Concentration data missing GEOID column:", conc_id_col))
    }

    # Check lognormal columns
    if (!"conc_meanlog" %in% names(conc_dt)) {
      issues <- c(issues, "Concentration data missing 'conc_meanlog' column")
    } else {
      # Check meanlog range
      meanlog_range <- range(conc_dt$conc_meanlog, na.rm = TRUE)
      if (meanlog_range[1] < -5 || meanlog_range[2] > 10) {
        warnings <- c(warnings, paste(
          "Concentration data: meanlog values outside typical range [-5, 10]:",
          "actual [", meanlog_range[1], ",", meanlog_range[2], "]"
        ))
      }
    }

    if (!"conc_sdlog" %in% names(conc_dt)) {
      warnings <- c(warnings, "Concentration data missing 'conc_sdlog' column (will use default)")
    } else {
      # Check sdlog is positive
      if (any(conc_dt$conc_sdlog <= 0, na.rm = TRUE)) {
        issues <- c(issues, "Concentration data: sdlog must be positive")
      }
    }

    # Check private well percentage
    if (!"PWELL_private_pct" %in% names(conc_dt)) {
      warnings <- c(warnings, "Concentration data missing 'PWELL_private_pct' column")
    } else {
      pwell_range <- range(conc_dt$PWELL_private_pct, na.rm = TRUE)
      if (pwell_range[1] < 0 || pwell_range[2] > 100) {
        issues <- c(issues, paste(
          "Concentration data: PWELL_private_pct must be between 0 and 100:",
          "actual [", pwell_range[1], ",", pwell_range[2], "]"
        ))
      }
    }

    message("  ", nrow(conc_dt), " geographic units")
  }

  # ===== Health Data Validation =====
  message("Checking health outcome data...")
  if (!file.exists(health_file)) {
    issues <- c(issues, paste("Health file not found:", health_file))
  } else {
    # Try to detect format
    health <- tryCatch({
      if (grepl("\\.txt$", health_file, ignore.case = TRUE)) {
        data.table::fread(health_file, sep = "\t")
      } else {
        data.table::fread(health_file)
      }
    }, error = function(e) {
      issues <<- c(issues, paste("Could not read health file:", e$message))
      return(NULL)
    })

    if (!is.null(health)) {
      summary_list$health <- list(
        n_rows = nrow(health),
        n_cols = ncol(health),
        columns = names(health)
      )

      # Check GEOID column
      if (!health_id_col %in% names(health)) {
        issues <- c(issues, paste("Health missing GEOID column:", health_id_col))
      } else {
        # Check FIPS format
        fips_vals <- health[[health_id_col]]
        if (any(nchar(as.character(fips_vals)) != 5, na.rm = TRUE)) {
          warnings <- c(warnings, "Health: Some FIPS codes are not 5 digits")
        }
      }

      message("  ", nrow(health), " records")
    }
  }

  # ===== Cross-file validation =====
  message("Checking geographic ID alignment...")
  if (exists("prob_dt") && exists("conc_dt") &&
      prob_id_col %in% names(prob_dt) && conc_id_col %in% names(conc_dt)) {
    prob_ids <- unique(prob_dt[[prob_id_col]])
    conc_ids <- unique(conc_dt[[conc_id_col]])

    common_ids <- intersect(prob_ids, conc_ids)
    prob_only <- setdiff(prob_ids, conc_ids)
    conc_only <- setdiff(conc_ids, prob_ids)

    message("  ", length(common_ids), " geographic IDs in both datasets")
    if (length(prob_only) > 0) {
      warnings <- c(warnings, paste(length(prob_only), "IDs in probability data only"))
    }
    if (length(conc_only) > 0) {
      warnings <- c(warnings, paste(length(conc_only), "IDs in concentration data only"))
    }
  }

  # ===== Report results =====
  message("\n========== Validation Results ==========")

  valid <- length(issues) == 0

  if (valid) {
    message("[OK] All validation checks passed!")
  } else {
    message("[X] ", length(issues), " issue(s) found:")
    for (issue in issues) {
      message("  ERROR: ", issue)
    }
  }

  if (length(warnings) > 0) {
    message("\n", length(warnings), " warning(s):")
    for (warn in warnings) {
      message("  WARNING: ", warn)
    }
  }

  return(list(
    valid = valid,
    issues = issues,
    warnings = warnings,
    summary = summary_list
  ))
}
