#' Perform Sensitivity Analysis for Environmental Exposure
#'
#' @description
#' This function performs a sensitivity analysis to assess the impact of
#' environmental exposure on specified health outcomes. It integrates modeled
#' probability data and measured concentration data to create a probabilistic
#' exposure model, weighted by the proportion of the population using each
#' source. It then uses multiple imputation to generate a series of datasets
#' with imputed exposure levels.
#'
#' @param ndraws An integer specifying the number of imputed datasets to
#'   generate.
#' @param model_prob_csv A file path to the CSV containing modeled multinomial
#'   probabilities of exposure levels (e.g., from predictive models).
#' @param conc_data_csv A file path to the CSV containing measured concentration
#'   lognormal parameters and the percentage of private well users.
#' @param birth_data_txt A file path to the text file containing the birth data.
#' @param regression_formula A string or formula object for the regression
#'   model.
#' @param output_dir A file path to the directory where output files will be
#'   saved.
#' @param targets A character vector of the dependent variables (health
#'   outcomes) for the regression analysis. Defaults to c("OEGEST", "BWT").
#' @param impute_vars A character vector of column names to be imputed using
#'   MICE. If NULL, no additional imputation is performed. Defaults to NULL.
#' @param as_cat_label A character vector of labels for the arsenic
#'   concentration categories. Defaults to c("As<5", "As5-10", "As10+").
#' @param drop_as_cat_label_reg A character vector of arsenic categories to be
#'   used as the reference level in the regression. Defaults to c("As<5").
#' @param columns_to_select A character vector of column names to be selected
#'   from the birth data. If NULL, all columns are used.
#' @param rucc_col A character string specifying the column name for the
#'   Rural-Urban Continuum Code. Defaults to NULL.
#' @param id_col A character string specifying the column name for the
#'   Geographic Identifier in the arsenic data. Supports any identifier format
#'   including county FIPS codes, census tract IDs, ZIP codes, or participant-
#'   specific IDs. The only requirement is that identifiers match between
#'   datasets. Defaults to "GEOID10".
#' @param prob_cols A character vector of column names from the probability data
#'   for the exposure probability categories. Defaults to
#'   c("prob_C1", "prob_C2", "prob_C3").
#' @param conc_cutoffs A numeric vector of cutoffs for categorizing exposure
#' levels based on concentration data. If NULL, will be automatically determined
#' based on number of probability columns. Defaults to c(5, 10)
#' @param pop_well_col A character string specifying the column name from the
#'   probability data for the population of well users. Defaults to "Wells_2010".
#'   This has to be at the same scale as the `id_col` and `record_id_col`.
#' @param record_id_col A character string specifying the column name for
#'   the geographic/participant identifier in the health data. Supports any
#'   identifier format (FIPS codes, census tracts, ZIP codes, participant IDs).
#'   Must have matching values with `id_col` in the arsenic data.
#'   Defaults to "FIPS".
#' @param as_level_col A character string for the name of the imputed arsenic
#'   level column. Defaults to "AsLevel".
#' @param conc_mean_col A character string for the column name of the
#'   concentration lognormal meanlog parameter. Defaults to "conc_meanlog".
#'   If this column doesn't exist and `conc_raw_col` is provided, it will
#'   be calculated automatically.
#' @param conc_sd_col A character string for the column name of the
#'   concentration lognormal sdlog parameter. Defaults to "conc_sdlog".
#'   If this column doesn't exist, `default_sdlog` will be used.
#' @param conc_raw_col A character string for the column name of raw
#'   concentration values (in measurement units). If provided and the
#'   meanlog/sdlog columns don't exist, lognormal parameters will be calculated
#'   automatically as: meanlog = log(concentration + 0.1), sdlog = default_sdlog.
#'   Defaults to NULL.
#' @param default_sdlog A numeric value for the default lognormal sdlog
#'   when the sdlog column is not provided. Defaults to 1.0, which is typical
#'   for environmental concentration data.
#' @param pwell_col A character string for the column name of the
#'   percentage of private well users. Defaults to "PWELL_private_pct".
#' @param seed An integer for setting the random number generator seed for
#'   reproducibility. Defaults to 12345.
#' @param mice_m An integer for the number of multiple imputations for
#'   covariates. Defaults to 1.
#' @param mice_maxit An integer for the maximum number of iterations for mice.
#'   Defaults to 5.
#' @param mice_method A character string specifying the imputation method for
#'   covariates. Defaults to "pmm" (predictive mean matching).
#' @param mice_covs A character vector of column names to be used as
#'   covariates in the mice imputation process. Defaults to c("pm").
#'   Make sure these columns exist in the birth data.
#' @param apply_imputation_fallback Logical indicating whether to apply fallback
#'  imputation for missing arsenic levels. The fallback fills missing
#'  values with the most common category in the respective dataset. If
#'  enabled, this fallback will be applied to any missing arsenic levels
#'  after the initial imputation step. Otherwise, it will remove any rows
#'  with missing arsenic levels.
#' @param format_ids Logical indicating whether to apply smart identifier
#'   formatting. If TRUE (default), numeric FIPS-like codes (1-5 digits) are
#'   zero-padded to 5 digits for consistency, while alphanumeric identifiers
#'   (census tracts, participant IDs, etc.) are preserved unchanged. Set to
#'   FALSE to skip formatting and use identifiers as-is (converted to character).
#'
#' @return A list of data frames, each containing the pooled regression results
#'   for a target health outcome.
#'
#' @examples
#' \dontrun{
#' # Basic usage with minimal parameters
#' results <- perform_sensitivity_analysis(
#'   ndraws = 10,
#'   model_prob_csv = "prob_model_data.csv",
#'   conc_data_csv = "conc_measured_data.csv",
#'   birth_data_txt = "birth_outcomes.txt",
#'   regression_formula = "~ as.factor(AsLevel) + maternal_age + (1|county)",
#'   output_dir = "results/",
#'   targets = c("birth_weight", "gestational_age")
#' )
#'
#' # Advanced usage with custom parameters and MICE imputation
#' results <- perform_sensitivity_analysis(
#'   ndraws = 100,
#'   model_prob_csv = "data/prob_data.csv",
#'   conc_data_csv = "data/conc_data.csv",
#'   birth_data_txt = "data/births.txt",
#'   regression_formula = "~ as.factor(AsLevel) + MAGE_R + rural + (1|FIPS)",
#'   output_dir = "sensitivity_results/",
#'   targets = c("OEGEST", "BWT"),
#'   impute_vars = c("MAGE_R", "education"),
#'   mice_m = 5,
#'   mice_maxit = 10,
#'   seed = 42
#' )
#' }
#'
#' @export
perform_sensitivity_analysis <- function(
  ndraws = 10,
  model_prob_csv = NULL,
  conc_data_csv = NULL,
  birth_data_txt = NULL,
  regression_formula = NULL,
  output_dir = "sensitivity_results/",
  targets = c("OEGEST", "BWT"),
  impute_vars = NULL,
  as_cat_label = c("As<5", "As5-10", "As10+"),
  drop_as_cat_label_reg = c("As<5"),
  columns_to_select = NULL,
  rucc_col = NULL,
  id_col = "GEOID10",
  prob_cols = c("prob_C1", "prob_C2", "prob_C3"),
  conc_cutoffs = c(5, 10),
  pop_well_col = "Wells_2010",
  record_id_col = "FIPS",
  as_level_col = "AsLevel",
  conc_mean_col = "conc_meanlog",
  conc_sd_col = "conc_sdlog",
  conc_raw_col = NULL,
  default_sdlog = 1.0,
  pwell_col = "PWELL_private_pct",
  seed = 12345,
  mice_m = 1,
  mice_maxit = 5,
  mice_method = "pmm",
  mice_covs = c("pm"),
  apply_imputation_fallback = TRUE,
  format_ids = TRUE
) {

  # --- 0. Validate all input parameters ---
  validate_all_inputs(
    birth_data_txt = birth_data_txt,
    model_prob_csv = model_prob_csv,
    conc_data_csv = conc_data_csv,
    ndraws = ndraws,
    regression_formula = regression_formula,
    targets = targets,
    impute_vars = impute_vars,
    output_dir = output_dir,
    mice_m = mice_m,
    mice_maxit = mice_maxit,
    mice_method = mice_method,
    conc_mean_col = conc_mean_col,
    conc_sd_col = conc_sd_col,
    conc_raw_col = conc_raw_col,
    default_sdlog = default_sdlog,
    seed = seed,
    as_cat_label = as_cat_label,
    as_level_col = as_level_col,
    drop_as_cat_label_reg = drop_as_cat_label_reg,
    conc_cutoffs = conc_cutoffs,
    columns_to_select = columns_to_select,
    rucc_col = rucc_col,
    id_col = id_col,
    prob_cols = prob_cols,
    pop_well_col = pop_well_col,
    record_id_col = record_id_col,
    pwell_col = pwell_col,
    mice_covs = mice_covs,
    apply_imputation_fallback = apply_imputation_fallback
  )

  set.seed(seed)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # --- 1. Load and Process Arsenic Data ---
  arsenic_data <- load_and_process_arsenic_data(
    model_prob_csv = model_prob_csv,
    conc_data_csv = conc_data_csv,
    prob_cols = prob_cols,
    pop_well_col = pop_well_col,
    conc_mean_col = conc_mean_col,
    conc_sd_col = conc_sd_col,
    conc_raw_col = conc_raw_col,
    default_sdlog = default_sdlog,
    pwell_col = pwell_col,
    id_col = id_col,
    conc_cutoffs = conc_cutoffs,
    format_ids = format_ids
  )

  # --- 2. Load and Process Birth Data ---
  births <- load_and_process_birth_data(
    birth_data_txt = birth_data_txt,
    columns_to_select = columns_to_select,
    rucc_col = rucc_col,
    record_id_col = record_id_col,
    format_ids = format_ids
  )

  # --- 3. Multiple Imputation ---
  imputed_datasets <- impute_arsenic_exposure(
    datmatx_prob = arsenic_data$datmatx_prob,
    rasterprob_combined = arsenic_data$rasterprob_combined,
    cellpop = arsenic_data$cellpop,
    births = births,
    ndraws = ndraws,
    as_cat_label = as_cat_label,
    id_col = id_col,
    record_id_col = record_id_col,
    as_level_col = as_level_col,
    apply_imputation_fallback = apply_imputation_fallback
  )
  # Impute additional variables if specified
  if (!is.null(impute_vars)) {
    imputed_datasets <- impute_additional_variables(
      imputed_datasets = imputed_datasets,
      impute_vars = impute_vars,
      births = births,
      mice_covs = mice_covs,
      mice_m = mice_m,
      mice_maxit = mice_maxit,
      mice_method = mice_method,
      seed = seed
    )
  }
  # Validate imputed datasets
  validate_imputed_datasets(imputed_datasets, as_level_col, as_cat_label)

  # --- 4. Regression and Pooling ---

  combined_births <- data.table::rbindlist(imputed_datasets, idcol = ".id")
  # Ensure the reference level is set correctly
  filtered_as_cat_label <- base::setdiff(as_cat_label, drop_as_cat_label_reg)
  combined_births[[as_level_col]] <- stats::relevel(
    combined_births[[as_level_col]],
    ref = drop_as_cat_label_reg
  )
  res <- regression_analysis(
    births = combined_births,
    targets = targets,
    as_cat_label = filtered_as_cat_label,
    as_level_col = as_level_col,
    regression_formula = regression_formula,
    output_dir = output_dir
  )
  return(res)
}
