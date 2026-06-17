#' =============================================================================
#' Generate simulated demographics/health data for the geoExposeR examples
#' =============================================================================
#'
#' Produces `Demographics_dom.txt`: a tab-delimited file of *individual-level*
#' health records (one row per participant) keyed to the county FIPS codes used
#' in `prob_dom.csv` / `conc_dom.csv`.
#'
#' Unlike an earlier version of this file (which held a single row per county),
#' each county here contains MULTIPLE participants. This individual-within-
#' county structure is what makes a geographic random intercept `(1 | id)`
#' identifiable: with repeated observations per county, `lme4` can estimate the
#' between-county variance separately from the residual (within-county)
#' variance. The `age_decade` grouping (only 5 levels) is retained so the
#' examples can also demonstrate using age as a *contextual* grouping variable.
#'
#' All values are randomly generated dummy data for demonstration only.
#'
#' Reproducible: run `Rscript examples/input_data/generate_demographics.R`
#' from the package root.
#' =============================================================================

set.seed(2026)

prob_path <- file.path("examples", "input_data", "prob_dom.csv")
out_path  <- file.path("examples", "input_data", "Demographics_dom.txt")

prob <- data.table::fread(prob_path)
ids  <- prob$id

# County-level arsenic-risk proxy: modeled probability of high arsenic
# (P(As >= 10 ug/L)). Used only to give simulated SBP a faint, plausible
# county-level signal so the example output is not pure noise.
risk <- prob$RFCprobC3
risk[is.na(risk)] <- mean(risk, na.rm = TRUE)
risk_centered <- risk - mean(risk)

# Latent county-level baseline shift (the quantity a random intercept targets).
county_intercept <- rnorm(length(ids), mean = 0, sd = 5)

age_levels <- c("30-39", "40-49", "50-59", "60-69", "70-79")
age_mid    <- c(35, 45, 55, 65, 75)
age_probs  <- c(0.15, 0.20, 0.25, 0.25, 0.15)

# SBP generative model (mmHg):
#   118 baseline + ~0.4 mmHg per year of age above 50
#   + a small county-level arsenic-risk effect + county random intercept
#   + individual residual noise.
sbp_per_year_age <- 0.4
sbp_risk_effect  <- 12

records <- vector("list", length(ids))
for (i in seq_along(ids)) {
  n_i  <- sample(4:12, size = 1)
  ages <- sample(age_levels, size = n_i, replace = TRUE, prob = age_probs)
  amid <- age_mid[match(ages, age_levels)]
  mu <- 118 +
    sbp_per_year_age * (amid - 50) +
    sbp_risk_effect * risk_centered[i] +
    county_intercept[i]
  sbp <- round(rnorm(n_i, mean = mu, sd = 9))
  sbp <- pmin(pmax(sbp, 90), 200)  # keep within a plausible physiologic range
  records[[i]] <- data.frame(
    id = ids[i],
    SBP = sbp,
    age_decade = ages,
    stringsAsFactors = FALSE
  )
}

demographics <- do.call(rbind, records)

data.table::fwrite(demographics, out_path, sep = "\t")

cat("Wrote", nrow(demographics), "participant records across",
    length(unique(demographics$id)), "counties to", out_path, "\n")
cat("Participants per county: range",
    paste(range(table(demographics$id)), collapse = "-"),
    "| mean", round(mean(table(demographics$id)), 1), "\n")
