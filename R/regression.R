#' Perform Regression Analysis
#'
#' This function performs regression analysis on multiply imputed data using
#' linear mixed-effects models and pools the results using Rubin's Rules.
#'
#' @param births A data frame containing the birth data with imputation
#'   identifiers.
#' @param targets A character vector of target variables (health outcomes).
#' @param as_cat_label A character vector of arsenic categories (excluding
#'   reference).
#' @param as_level_col A character string for the arsenic level column name.
#' @param regression_formula A formula object or character string for the
#'   regression model.
#' @param output_dir A character string for the output directory path.
#' @return A named list of data frames with pooled regression results for each
#'   target.
#'
#' @keywords internal
regression_analysis <- function(
  births,
  targets,
  as_cat_label,
  as_level_col,
  regression_formula,
  output_dir
) {
  # Validate inputs
  if (!".id" %in% names(births)) {
    stop("births data must contain '.id' column for imputation grouping")
  }
  if (!all(targets %in% names(births))) {
    missing_targets <- setdiff(targets, names(births))
    stop(
      "Target variables not found in data: ",
      paste(missing_targets, collapse = ", ")
    )
  }
  # Pre-compute arsenic factor terms to avoid repeated string operations
  arsenic_terms <- paste0("as.factor(", as_level_col, ")", as_cat_label)
  # Initialize results list
  results <- vector("list", length(targets))
  names(results) <- targets  
  # Process each target variable
  for (var in targets) {
    message("Processing target variable: ", var)
    # Create formula for current target
    full_formula <- if (is.character(regression_formula)) {
      as.formula(paste(var, regression_formula))
    } else {
      update(regression_formula, paste(var, "~ ."))
    }
    # Fit models for each imputed dataset with error handling
    model_fits <- tryCatch({
      births %>%
        dplyr::group_by(.data$.id) %>%
        dplyr::do(
          fitmodel = lme4::lmer(
            full_formula,
            data = .,
            control = lme4::lmerControl(calc.derivs = FALSE)
          )
        )
    }, error = function(e) {
      stop("Model fitting failed for ", var, ": ", e$message)
    })
    # Extract coefficients and standard errors more efficiently
    model_results <- lapply(model_fits$fitmodel, function(model) {
      model_summary <- broom.mixed::tidy(model, effects = "fixed")
      # Filter for arsenic terms only
      model_summary[
        model_summary$term %in% arsenic_terms,
        c("term", "estimate", "std.error")
      ]
    })
    # Combine results across imputations
    combined_results <- do.call(rbind, model_results)
    if (nrow(combined_results) == 0) {
      warning(
        "No arsenic terms found in model for ", var,
        ". Check factor levels."
      )
      next
    }
    # Pool results using Rubin's Rules for each arsenic term
    pooled_results <- pool_estimates_by_term(
      combined_results,
      arsenic_terms,
      model_fits$fitmodel[[1]]
    )
    # Save results
    output_file <- file.path(output_dir, paste0(var, "_pooled_results.csv"))
    data.table::fwrite(pooled_results, output_file)

    # Rename pooled results terms to match arsenic categories
    pooled_results$term <- factor(
      pooled_results$term,
      levels = arsenic_terms,
      labels = as_cat_label
    )
    # Store in results list
    results[[var]] <- pooled_results
  }
  return(results)
}

#' Pool regression estimates using Rubin's Rules
#'
#' @param combined_results Data frame with term, estimate, and std.error columns
#' @param arsenic_terms Character vector of arsenic terms to pool
#' @param reference_model A fitted model to extract degrees of freedom
#' @return Data frame with pooled estimates and confidence intervals
#'
#' @keywords internal
pool_estimates_by_term <- function(
  combined_results,
  arsenic_terms,
  reference_model
) {
  # Get degrees of freedom from reference model
  deg_freedom <- tryCatch({
    as.integer(broom.mixed::glance(reference_model)$df.residual)
  }, error = function(e) {
    warning(
      "Could not extract degrees of freedom, using large sample approximation"
    )
    Inf
  })
  # Pool estimates for each arsenic term
  pooled_list <- lapply(arsenic_terms, function(term) {
    term_data <- combined_results[combined_results$term == term, ]
    if (nrow(term_data) == 0) {
      return(NULL)
    }
    # Apply Rubin's Rules
    pooled_estimate <- pool_single_estimate(
      estimates = term_data$estimate,
      std_errors = term_data$std.error,
      df = deg_freedom
    )
    # Add term identifier
    pooled_estimate$term <- term
    return(pooled_estimate)
  })
  # Combine results and remove NULL entries
  filtered_list <- Filter(Negate(is.null), pooled_list)
  if (length(filtered_list) == 0) {
    return(NULL)
  }
  pooled_results <- do.call(rbind, filtered_list)
  # Reorder columns for clarity
  if (!is.null(pooled_results) && nrow(pooled_results) > 0) {
    pooled_results <- pooled_results[, c(
      "term", "q.mi", "se.mi", "statistic",
      "conf.low", "conf.high", "p.value"
    )]
  }
  return(pooled_results)
}

#' Apply Rubin's Rules to pool a single parameter estimate
#'
#' @param estimates Numeric vector of parameter estimates across imputations
#' @param std_errors Numeric vector of standard errors across imputations
#' @param df Degrees of freedom for t-distribution
#' @return Data frame with pooled estimate and inference
#'
#' @keywords internal
pool_single_estimate <- function(estimates, std_errors, df) {
  m <- length(estimates)
  if (m == 1) {
    # Single imputation case
    return(data.frame(
      q.mi = estimates[1],
      se.mi = std_errors[1],
      statistic = estimates[1] / std_errors[1],
      conf.low = estimates[1] + std_errors[1] * qt(0.025, df),
      conf.high = estimates[1] + std_errors[1] * qt(0.975, df),
      p.value = 2 * pt(
        abs(estimates[1] / std_errors[1]), df, lower.tail = FALSE
      )
    ))
  }
  # Multiple imputation case - apply Rubin's Rules
  # Pooled estimate
  q_bar <- mean(estimates)
  # Within-imputation variance
  u_bar <- mean(std_errors^2)
  # Between-imputation variance
  b <- sum((estimates - q_bar)^2) / (m - 1)
  # Total variance
  t_var <- u_bar + (1 + 1 / m) * b
  se_pooled <- sqrt(t_var)
  # Test statistic
  t_stat <- q_bar / se_pooled
  # Degrees of freedom for t-distribution (if finite)
  if (is.finite(df) && b > 0) {
    r <- (1 + 1 / m) * b / u_bar
    df_adjusted <- (m - 1) * (1 + 1 / r)^2
    df_final <- min(df, df_adjusted)
  } else {
    df_final <- df
  }
  # Confidence intervals and p-value
  t_critical <- qt(0.975, df_final)
  conf_low <- q_bar - t_critical * se_pooled
  conf_high <- q_bar + t_critical * se_pooled
  p_value <- 2 * pt(abs(t_stat), df_final, lower.tail = FALSE)
  return(data.frame(
    q.mi = q_bar,
    se.mi = se_pooled,
    statistic = t_stat,
    conf.low = conf_low,
    conf.high = conf_high,
    p.value = p_value
  ))
}
