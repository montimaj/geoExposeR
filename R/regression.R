#' Perform Regression Analysis
#'
#' This function performs regression analysis on multiply imputed data and pools
#' the results using Rubin's Rules. It supports fixed-effects, random-effects,
#' and mixed-effects models, each available in a linear variant (Gaussian /
#' identity link) or a non-linear variant (generalized models with a `family`,
#' e.g. logistic or Poisson).
#'
#' The fitting engine is selected from the combination of `model_type` and
#' `model_family`:
#'
#' \tabular{lll}{
#'   \strong{model_type} \tab \strong{linear (model_family = NULL)} \tab
#'     \strong{non-linear (model_family set)} \cr
#'   `"fixed"`  \tab [stats::lm()]   \tab [stats::glm()]  \cr
#'   `"random"` \tab [lme4::lmer()]  \tab [lme4::glmer()] \cr
#'   `"mixed"`  \tab [lme4::lmer()]  \tab [lme4::glmer()]
#' }
#'
#' Multi-level categorical outcomes are supported via multinomial logistic
#' regression (`model_family = "multinomial"`): fixed-effects fits use
#' [nnet::multinom()], while random- and mixed-effects fits use
#' [mclogit::mblogit()]. Multinomial results include a `y.level` column
#' identifying the outcome category (relative to the reference level) that each
#' pooled estimate refers to.
#'
#' @param births A data frame containing the birth data with imputation
#'   identifiers.
#' @param targets A character vector of target variables (health outcomes).
#' @param cat_label A character vector of exposure categories (excluding
#'   reference).
#' @param exposure_level_col A character string for the exposure level column
#'   name.
#' @param regression_formula A formula object or character string for the
#'   regression model.
#' @param output_dir A character string for the output directory path.
#' @param model_type A character string selecting the effects structure. One of
#'   `"mixed"` (default, fixed + random effects), `"random"` (random-intercept
#'   only), or `"fixed"` (no random effects).
#' @param model_family The error distribution / link for non-linear
#'   (generalized) models. Accepts `NULL` (default, a linear/Gaussian model), a
#'   family name such as `"binomial"` or `"poisson"`, the special value
#'   `"multinomial"` (multi-level categorical outcome, fitted with
#'   [nnet::multinom()] for fixed effects or [mclogit::mblogit()] for
#'   random-/mixed-effects), a [stats::family()] object, or a family-generating
#'   function. When non-`NULL`, generalized models are fitted.
#' @return A named list of data frames with pooled regression results for each
#'   target.
#'
#' @keywords internal
regression_analysis <- function(
  births,
  targets,
  cat_label,
  exposure_level_col,
  regression_formula,
  output_dir,
  model_type = "mixed",
  model_family = NULL
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
  # Multinomial (multi-level categorical) outcomes use a dedicated engine
  # (nnet::multinom) rather than a stats family object.
  is_multinomial <- is_multinomial_spec(model_family)
  # Resolve the (optional) family once for all targets
  family_obj <- if (is_multinomial) NULL else resolve_model_family(model_family)
  # Pre-compute exposure factor terms to avoid repeated string operations
  exposure_terms <- paste0("as.factor(", exposure_level_col, ")", cat_label)
  # Exposure coefficients are contrasts with the reference (first) category.
  # In an imputed dataset where that category was drawn for no record, the
  # model uses a different baseline, so its coefficients are not comparable.
  exposure_values <- births[[exposure_level_col]]
  if (!is.null(exposure_values)) {
    ref_level <- if (is.factor(exposure_values)) {
      levels(exposure_values)[1]
    } else {
      sort(unique(stats::na.omit(as.character(exposure_values))))[1]
    }
    has_ref <- tapply(
      as.character(exposure_values) == ref_level, births$.id,
      function(x) any(x, na.rm = TRUE)
    )
    if (!all(has_ref)) {
      warning(
        "The reference exposure category '", ref_level, "' was drawn for no ",
        "record in ", sum(!has_ref), " of ", length(has_ref), " imputed ",
        "datasets. Their exposure coefficients are contrasts with a different ",
        "category, so the pooled estimates mix contrasts. Choose a more ",
        "common reference category (drop_cat_label_ref) or merge rare ",
        "categories.",
        call. = FALSE
      )
    }
  }
  # Initialize results list
  results <- vector("list", length(targets))
  names(results) <- targets
  # Process each target variable
  for (var in targets) {
    message("Processing target variable: ", var)
    # For categorical models (binomial logistic or multinomial), coerce a
    # character response (e.g. "high"/"normal" read from a text file) into a
    # factor so that glm()/glmer()/multinom() can model it. The first factor
    # level is treated as the reference category.
    is_binomial <- !is.null(family_obj) && family_obj$family == "binomial"
    if ((is_binomial || is_multinomial) && is.character(births[[var]])) {
      births[[var]] <- factor(births[[var]])
    }
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
          fitmodel = fit_single_model(
            formula = full_formula,
            data = .,
            model_type = model_type,
            family_obj = family_obj,
            is_multinomial = is_multinomial
          )
        )
    }, error = function(e) {
      stop("Model fitting failed for ", var, ": ", e$message)
    })
    # Extract coefficients and standard errors more efficiently. Multinomial
    # models additionally carry a "y.level" column (the outcome category).
    model_results <- lapply(model_fits$fitmodel, function(model) {
      model_summary <- tidy_fixed_effects(model, model_type, is_multinomial)
      keep_cols <- intersect(
        c("term", "estimate", "std.error", "y.level"),
        names(model_summary)
      )
      # Filter for exposure terms only
      model_summary[model_summary$term %in% exposure_terms, keep_cols]
    })
    # Combine results across imputations
    combined_results <- do.call(rbind, model_results)
    if (nrow(combined_results) == 0) {
      warning(
        "No exposure terms found in model for ", var,
        ". Check factor levels."
      )
      next
    }
    # Pool results using Rubin's Rules for each exposure term
    pooled_results <- pool_estimates_by_term(
      combined_results,
      exposure_terms,
      model_fits$fitmodel[[1]],
      model_type,
      n_imputations = length(model_fits$fitmodel),
      target = var
    )
    # Save results
    output_file <- file.path(output_dir, paste0(var, "_pooled_results.csv"))
    data.table::fwrite(pooled_results, output_file)

    # Rename pooled results terms to match exposure categories
    pooled_results$term <- factor(
      pooled_results$term,
      levels = exposure_terms,
      labels = cat_label
    )
    # Store in results list
    results[[var]] <- pooled_results
  }
  return(results)
}

#' Resolve a model family specification to a family object
#'
#' Normalizes the `model_family` argument into either `NULL` (linear model) or a
#' [stats::family()] object usable by [stats::glm()] / [lme4::glmer()].
#'
#' @param model_family `NULL`, a family name (character), a family object, or a
#'   family-generating function.
#' @return `NULL` for linear models, otherwise a family object.
#'
#' @keywords internal
resolve_model_family <- function(model_family) {
  if (is.null(model_family)) {
    return(NULL)
  }
  if (inherits(model_family, "family")) {
    return(model_family)
  }
  if (is.function(model_family)) {
    return(model_family())
  }
  # Character family name (e.g. "binomial", "poisson")
  family_fun <- tryCatch(
    get(model_family, mode = "function"),
    error = function(e) {
      stop("Unknown model_family: '", model_family, "'")
    }
  )
  family_fun()
}

#' Test whether a model_family requests a multinomial model
#'
#' @param model_family A `model_family` specification.
#' @return `TRUE` if `model_family` is the string `"multinomial"` (case
#'   insensitive), otherwise `FALSE`.
#'
#' @keywords internal
is_multinomial_spec <- function(model_family) {
  is.character(model_family) &&
    length(model_family) == 1 &&
    tolower(model_family) == "multinomial"
}

#' Fit a single regression model
#'
#' Dispatches to the appropriate fitting engine based on `model_type` and
#' whether a generalized `family_obj` was supplied.
#'
#' @param formula A model formula.
#' @param data A data frame for a single imputed dataset.
#' @param model_type One of `"fixed"`, `"random"`, or `"mixed"`.
#' @param family_obj `NULL` for a linear model, or a family object for a
#'   generalized (non-linear) model.
#' @param is_multinomial Logical; if `TRUE`, fit a multinomial logistic model
#'   with [nnet::multinom()] (requires `model_type = "fixed"`).
#' @return A fitted model object.
#'
#' @keywords internal
fit_single_model <- function(formula, data, model_type, family_obj,
                             is_multinomial = FALSE) {
  if (is_multinomial) {
    # Multi-level categorical outcome: multinomial logistic regression.
    if (model_type == "fixed") {
      # Fixed effects only: nnet::multinom. Hess = TRUE stores the Hessian so
      # that standard errors can be computed without re-evaluating the
      # (grouped) data frame at summary() time. Unlike the other engines,
      # multinom keeps unused factor levels and reports a zero estimate for an
      # exposure category drawn for no record, so drop them first.
      return(nnet::multinom(
        formula, data = droplevels(data), trace = FALSE, Hess = TRUE
      ))
    }
    # Random- or mixed-effects multinomial: mclogit::mblogit, which takes the
    # random-effects structure through a separate `random` argument.
    return(fit_multinomial_mixed(formula, data))
  }
  if (model_type == "fixed") {
    # No random effects: ordinary linear / generalized linear model
    if (is.null(family_obj)) {
      stats::lm(formula, data = data)
    } else {
      stats::glm(formula, data = data, family = family_obj)
    }
  } else {
    # Random-effects or mixed-effects model via lme4
    if (is.null(family_obj)) {
      lme4::lmer(
        formula,
        data = data,
        control = lme4::lmerControl(calc.derivs = FALSE)
      )
    } else {
      lme4::glmer(
        formula,
        data = data,
        family = family_obj,
        control = lme4::glmerControl(calc.derivs = FALSE)
      )
    }
  }
}

#' Fit a random-/mixed-effects multinomial model
#'
#' Splits an `lme4`-style formula (with inline random-effects terms such as
#' `(1 | group)`) into a fixed-effects formula plus a `random` specification
#' (using [reformulas::nobars()] and [reformulas::findbars()]), and fits a
#' baseline-category multinomial logistic model with [mclogit::mblogit()].
#'
#' @param formula A model formula containing at least one random-effects term.
#' @param data A data frame for a single imputed dataset.
#' @return A fitted [mclogit::mblogit()] model object.
#'
#' @keywords internal
fit_multinomial_mixed <- function(formula, data) {
  # Separate the fixed-effects part from the random-effects terms.
  fixed_formula <- reformulas::nobars(formula)
  random_bars <- reformulas::findbars(formula)
  random_spec <- lapply(random_bars, function(bar) {
    stats::as.formula(paste0("~ ", deparse(bar)))
  })
  # mblogit accepts a single random formula or a list of them.
  if (length(random_spec) == 1) {
    random_spec <- random_spec[[1]]
  }
  mclogit::mblogit(
    fixed_formula,
    random = random_spec,
    data = data,
    trace = FALSE
  )
}

#' Tidy the fixed-effect coefficients of a fitted model
#'
#' Returns a uniform data frame of fixed-effect terms regardless of the fitting
#' engine used.
#'
#' @param model A fitted model object (`lm`, `glm`, `lmerMod`, `glmerMod`, or
#'   `multinom`).
#' @param model_type One of `"fixed"`, `"random"`, or `"mixed"`.
#' @param is_multinomial Logical; if `TRUE`, the model is a multinomial
#'   [nnet::multinom()] fit and the returned data frame additionally contains a
#'   `y.level` column.
#' @return A data frame with `term`, `estimate`, and `std.error` columns (plus
#'   `y.level` for multinomial models).
#'
#' @keywords internal
tidy_fixed_effects <- function(model, model_type, is_multinomial = FALSE) {
  if (is_multinomial) {
    return(tidy_multinomial(model))
  }
  if (model_type == "fixed") {
    # lm / glm: pull the coefficient table from the model summary
    coef_table <- stats::coef(summary(model))
    data.frame(
      term = rownames(coef_table),
      estimate = coef_table[, "Estimate"],
      std.error = coef_table[, "Std. Error"],
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  } else {
    # lmer / glmer: use broom.mixed to extract fixed effects
    model_summary <- broom.mixed::tidy(model, effects = "fixed")
    as.data.frame(
      model_summary[, c("term", "estimate", "std.error")]
    )
  }
}

#' Tidy the coefficients of a multinomial model
#'
#' Converts the coefficient and standard-error tables of a multinomial model
#' into a long data frame with one row per (outcome level, term) combination.
#' Handles both fixed-effects [nnet::multinom()] fits and random-/mixed-effects
#' [mclogit::mblogit()] fits.
#'
#' @param model A fitted [nnet::multinom()] or [mclogit::mblogit()] model object.
#' @return A data frame with `y.level`, `term`, `estimate`, and `std.error`
#'   columns.
#'
#' @keywords internal
tidy_multinomial <- function(model) {
  if (inherits(model, "mblogit")) {
    # mclogit::mblogit summary rownames are "<level>~<term>".
    coef_table <- summary(model)$coefficients
    row_labels <- rownames(coef_table)
    return(data.frame(
      y.level = sub("~.*", "", row_labels),
      term = sub("^[^~]*~", "", row_labels),
      estimate = coef_table[, "Estimate"],
      std.error = coef_table[, "Std. Error"],
      stringsAsFactors = FALSE,
      row.names = NULL
    ))
  }
  coefs <- stats::coef(model)
  std_errors <- summary(model)$standard.errors
  if (is.matrix(coefs)) {
    # Three or more outcome levels: rows are non-reference outcome levels,
    # columns are model terms.
    grid <- expand.grid(
      y.level = rownames(coefs),
      term = colnames(coefs),
      stringsAsFactors = FALSE
    )
    grid$estimate <- as.vector(coefs)
    grid$std.error <- as.vector(std_errors)
    grid[, c("y.level", "term", "estimate", "std.error")]
  } else {
    # Two outcome levels: coefficients are a single named vector, like a binary
    # logistic model. The single non-reference level is the second level.
    data.frame(
      y.level = model$lev[2],
      term = names(coefs),
      estimate = as.numeric(coefs),
      std.error = as.numeric(std_errors),
      stringsAsFactors = FALSE
    )
  }
}

#' Extract residual degrees of freedom from a fitted model
#'
#' @param model A fitted model object.
#' @param model_type One of `"fixed"`, `"random"`, or `"mixed"`.
#' @return An integer number of residual degrees of freedom, or `Inf` if it
#'   cannot be determined.
#'
#' @keywords internal
extract_df_residual <- function(model, model_type) {
  if (inherits(model, c("multinom", "mblogit", "mclogit"))) {
    # Multinomial fits have no residual df: use Inf as the complete-data df.
    # pool_single_estimate() still applies Rubin's df when B > 0.
    return(Inf)
  }
  tryCatch({
    if (model_type == "fixed") {
      as.integer(stats::df.residual(model))
    } else {
      as.integer(broom.mixed::glance(model)$df.residual)
    }
  }, error = function(e) {
    warning(
      "Could not extract degrees of freedom; using Rubin's degrees of ",
      "freedom without a complete-data cap"
    )
    Inf
  })
}

#' Pool regression estimates using Rubin's Rules
#'
#' @param combined_results Data frame with term, estimate, and std.error columns
#' @param exposure_terms Character vector of exposure terms to pool
#' @param reference_model A fitted model to extract degrees of freedom
#' @param model_type One of `"fixed"`, `"random"`, or `"mixed"`. Used to select
#'   how residual degrees of freedom are extracted. Defaults to `"mixed"`.
#' @param n_imputations Number of imputed datasets to which a model was fitted,
#'   or `NULL` to skip the check. When supplied, a warning lists every exposure
#'   term (per outcome level for multinomial models) that could be estimated in
#'   fewer datasets, because such terms are pooled over fewer imputations, or
#'   omitted when no dataset estimates them.
#' @param target Optional name of the outcome, used in that warning.
#' @return Data frame with pooled estimates and confidence intervals. When
#'   `combined_results` contains a `y.level` column (multinomial models),
#'   estimates are pooled per outcome level and the column is retained.
#'
#' @keywords internal
pool_estimates_by_term <- function(
  combined_results,
  exposure_terms,
  reference_model,
  model_type = "mixed",
  n_imputations = NULL,
  target = NULL
) {
  # Get degrees of freedom from reference model
  deg_freedom <- extract_df_residual(reference_model, model_type)
  # Multinomial models carry an outcome-level dimension; pool within each level.
  has_levels <- "y.level" %in% names(combined_results)
  outcome_levels <- if (has_levels) {
    unique(combined_results$y.level)
  } else {
    NA_character_
  }
  pooled_list <- list()
  # Terms estimable in fewer than all imputed datasets (e.g. a category drawn
  # for no record in some datasets)
  shortfalls <- character(0)
  for (lev in outcome_levels) {
    for (term in exposure_terms) {
      term_data <- if (has_levels) {
        combined_results[
          combined_results$y.level == lev & combined_results$term == term,
        ]
      } else {
        combined_results[combined_results$term == term, ]
      }
      if (!is.null(n_imputations) && nrow(term_data) < n_imputations) {
        shortfalls <- c(shortfalls, paste0(
          if (has_levels) paste0(lev, ": ") else "",
          term, " (", nrow(term_data), " of ", n_imputations, ")"
        ))
      }
      if (nrow(term_data) == 0) {
        next
      }
      # Apply Rubin's Rules
      pooled_estimate <- pool_single_estimate(
        estimates = term_data$estimate,
        std_errors = term_data$std.error,
        df = deg_freedom
      )
      # Add identifiers
      pooled_estimate$term <- term
      if (has_levels) {
        pooled_estimate$y.level <- lev
      }
      pooled_list[[length(pooled_list) + 1]] <- pooled_estimate
    }
  }
  if (length(shortfalls) > 0) {
    warning(
      if (!is.null(target)) paste0("For target '", target, "', some") else "Some",
      " exposure terms could not be estimated in every imputed dataset ",
      "(for example, because a category was drawn for no record), so they ",
      "are pooled over fewer imputations: ",
      paste(shortfalls, collapse = "; "),
      ". Terms estimable in no dataset are omitted from the results. ",
      "Check category counts across imputations or merge rare categories.",
      call. = FALSE
    )
  }
  # Combine results
  if (length(pooled_list) == 0) {
    return(NULL)
  }
  pooled_results <- do.call(rbind, pooled_list)
  # Reorder columns for clarity
  if (!is.null(pooled_results) && nrow(pooled_results) > 0) {
    base_cols <- c(
      "term", "q.mi", "se.mi", "statistic",
      "conf.low", "conf.high", "p.value"
    )
    ordered_cols <- if (has_levels) c("y.level", base_cols) else base_cols
    pooled_results <- pooled_results[, ordered_cols]
  }
  return(pooled_results)
}

#' Apply Rubin's Rules to pool a single parameter estimate
#'
#' @param estimates Numeric vector of parameter estimates across imputations
#' @param std_errors Numeric vector of standard errors across imputations
#' @param df Complete-data (residual) degrees of freedom; `Inf` for
#'   large-sample models such as multinomial fits. With between-imputation
#'   variance, inference uses Rubin's (1987) degrees of freedom
#'   \eqn{(m - 1)(1 + 1/r)^2}, capped at `df`.
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
  # Rubin's (1987) degrees of freedom, capped at the complete-data value. This
  # also applies when df is infinite (multinomial models), where it is the
  # only small-m correction.
  if (b > 0) {
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
