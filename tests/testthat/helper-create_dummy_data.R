#' Create Dummy Data for Testing gwArsenicR
#'
#' @description
#' This helper function generates a complete set of dummy data files required to
#' run `perform_sensitivity_analysis()`. It creates data that mimics the
#' structure of real-world inputs, including USGS multinomial probabilities,
#' EPA lognormal parameters, and individual-level birth data.
#'
#' @param file_dir The directory where the dummy data files will be saved.
#' @param num_rows The number of rows of raster-level data to generate.
#' @param num_births The number of rows of birth data to generate.
#' @param seed An integer for the random number generator to ensure
#'   reproducibility.
#' @param introduce_nas A logical indicating whether to introduce NAs in the
#'   birth data.
#'
#' @return A list containing the file paths to the three generated data files:
#'   1. USGS arsenic probabilities CSV.
#'   2. EPA arsenic parameters CSV.
#'   3. Birth data TXT file.
#'
#' @keywords internal
create_dummy_data <- function(
  file_dir = "test_data",
  num_rows = 500, # Raster-level data
  num_births = 2000, # Individual-level data
  seed = 12345,
  introduce_nas = TRUE # Whether to introduce NAs in birth data
) {
  # Set seed for reproducibility
  set.seed(seed)

  # Ensure the directory exists
  if (!dir.exists(file_dir)) {
    dir.create(file_dir, recursive = TRUE)
  }

  # --- Generate Common Geographic Identifiers ---
  # Simulate 100 counties
  county_ids <- sprintf("%05d", 1:100)
  geoid10 <- sample(county_ids, num_rows, replace = TRUE)
  birth_fips <- sample(county_ids, num_births, replace = TRUE)

  # --- 1. USGS Arsenic Probability Data (Private Wells) ---
  usgs_prob_values <- t(apply(
    matrix(stats::runif(num_rows * 3), ncol = 3),
    1,
    function(x) x / sum(x)
  ))
  wells_2010 <- sample(50:500, num_rows, replace = TRUE)
  dummy_usgs_df <- data.frame(
    RFC3_C1v2 = usgs_prob_values[, 1],
    RFC3_C2v2 = usgs_prob_values[, 2],
    RFC3_C3v2 = usgs_prob_values[, 3],
    Wells_2010 = wells_2010,
    GEOID10 = geoid10
  )
  dummy_usgs_csv <- file.path(file_dir, "dummy_usgs_arsenic.csv")
  data.table::fwrite(dummy_usgs_df, dummy_usgs_csv)

  # --- 2. EPA Arsenic Data (Public Wells) ---
  epa_lognormal_meanlog <- stats::runif(num_rows, 0.5, 1.5)
  epa_lognormal_sdlog <- stats::runif(num_rows, 0.1, 0.5)
  # You can use a constant value for epa_lognormal_sdlog if desired
  # For example, epa_lognormal_sdlog <- 1.2
  private_well_percent <- stats::runif(num_rows, min = 1, max = 99)
  dummy_epa_df <- data.frame(
    EPA_AS_meanlog = epa_lognormal_meanlog,
    EPA_AS_sdlog = epa_lognormal_sdlog,
    PWELL_private_pct = private_well_percent,
    GEOID10 = geoid10
  )
  dummy_epa_csv <- file.path(file_dir, "dummy_epa_arsenic.csv")
  data.table::fwrite(dummy_epa_df, dummy_epa_csv)

  # --- 3. Birth Data ---
  births_df <- data.frame(
    MAGE_R = sample(18:45, num_births, replace = TRUE),
    MRACEHISP_F = factor(
      sample(
        c("White", "Black", "Hispanic", "Other"),
        num_births,
        replace = TRUE
      )
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
    # Introduce ~10% NAs into MAGE_R and smoke columns
    n_na <- floor(0.1 * num_births)
    births_df$MAGE_R[sample.int(num_births, n_na)] <- NA
    births_df$smoke[sample.int(num_births, n_na)] <- NA
  }
  dummy_births_txt <- file.path(file_dir, "dummy_births.txt")
  data.table::fwrite(births_df, dummy_births_txt, sep = "\t")

  return(list(
    usgs = dummy_usgs_csv,
    epa = dummy_epa_csv,
    births = dummy_births_txt
  ))
}
