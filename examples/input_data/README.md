# Example Input Data Documentation

> **⚠️ All data files in this directory are simulated dummy data** created for demonstration and testing purposes only. The county-level exposure files (`prob_dom.csv`, `conc_dom.csv`) contain 50 synthetic counties; the health file (`Demographics_dom.txt`) contains ~400 synthetic individual-level records nested within those 50 counties. They have realistic column structures but randomly generated values. Actual research data cannot be shared due to data use agreements and privacy restrictions.
>
> `Demographics_dom.txt` is produced by `generate_demographics.R` (run `Rscript examples/input_data/generate_demographics.R` from the package root).

## Core Arsenic Data

### prob_dom.csv — Modeled Arsenic Probability Data

Multinomial probability distributions of arsenic concentration categories derived from modeled predictions (e.g., USGS machine learning predictions) for private domestic wells. This represents the **modeled probability** input.

| Column | Description | Units / Format |
|--------|-------------|----------------|
| `id` | Geographic identifier (county FIPS code) | Integer (e.g., 1001 = FIPS 01001) |
| `pub_water_pop` | Population served by public water supply in the geographic unit | Count (persons); some values may be `NA` |
| `RFCprobC1` | P(arsenic < 5 µg/L) — Probability of low arsenic concentration | Proportion (0–1) |
| `RFCprobC2` | P(5 ≤ arsenic < 10 µg/L) — Probability of moderate arsenic concentration | Proportion (0–1) |
| `RFCprobC3` | P(arsenic ≥ 10 µg/L) — Probability of high arsenic concentration (above EPA MCL) | Proportion (0–1) |
| `RFC_maxprob` | Name of the column with the highest probability for this geographic unit | Character (e.g., "RFCprobC1") |
| `RFC_maxprob_num` | Value of the highest probability | Proportion (0–1) |
| `Dom_Pop_Percent` | Percentage of population using domestic (private) wells | Percentage (0–100) |

**Notes:** Probabilities `RFCprobC1 + RFCprobC2 + RFCprobC3` should sum to approximately 1.0 for each row. The EPA Maximum Contaminant Level (MCL) for arsenic is 10 µg/L.

---

### conc_dom.csv — Measured Arsenic Concentration Data

Arsenic concentrations measured in public water supply systems, aggregated at the county level. This represents the **measured concentration** input.

| Column | Description | Units / Format |
|--------|-------------|----------------|
| `id` | Geographic identifier (county FIPS code), matching probability data | Integer (e.g., 1001 = FIPS 01001) |
| `weighted_as20062008` | Population-weighted mean arsenic concentration for 2006–2008 | µg/L (micrograms per liter); some values may be `NA` |
| `county_pop_served20062008` | Total county population served by public water systems during 2006–2008 | Count (persons) |
| `well` | Well indicator variable | Numeric (0.5 in dummy data) |
| `Dom_Pop_Percent` | Percentage of population using domestic (private) wells | Percentage (0–100) |

**Notes:** Missing `weighted_as20062008` values indicate counties without monitoring data. The package handles these by falling back to probability-only estimates.

---

### Demographics_dom.txt — Health Outcome Data

Tab-delimited file containing **individual-level** health outcome measurements and demographic covariates for the study population. Each row is one participant, and **multiple participants are nested within each county `id`**. This individual-within-county structure is what makes a geographic random intercept `(1 | id)` identifiable (it gives `lme4` repeated observations per county, so between-county variance can be estimated separately from within-county residual variance).

| Column | Description | Units / Format |
|--------|-------------|----------------|
| `id` | Geographic identifier (county FIPS code), matching probability and concentration data; repeats across the participants in a county | Integer (e.g., 1001 = FIPS 01001) |
| `SBP` | Systolic Blood Pressure — the primary health outcome | mmHg (millimeters of mercury) |
| `age_decade` | Age group of the participant | Character (e.g., "30-39", "40-49", "50-59", "60-69", "70-79") |

**Notes on using these columns in a model:** The recommended random effect is a geographic random intercept by county, `(1 | id)` — exposure is assigned at the county level while `SBP` is measured on individuals nested within counties. `age_decade` can instead be used as a *contextual grouping variable* via `(1 | age_decade)` (a separate baseline per age band, with **no random slope**), but it has only 5 levels, so the between-group variance is estimated from very few groups and should be interpreted with caution. To simply adjust for age, use a fixed effect, `+ as.factor(age_decade)`. The `examples/run_geoExposeR_model_variants.R` script demonstrates these options.

---

## Nitrate Data

### prob_nitrate.csv — Modeled Nitrate Probability Data

Multinomial probability distributions of nitrate concentration categories. Same column structure as `prob_dom.csv`.

| Column | Description | Units / Format |
|--------|-------------|----------------|
| `id` | Geographic identifier (county FIPS code) | Integer |
| `pub_water_pop` | Population served by public water supply | Count (persons) |
| `RFCprobC1` | P(nitrate < 5 mg/L NO₃-N) — Probability of low nitrate | Proportion (0–1) |
| `RFCprobC2` | P(5 ≤ nitrate < 10 mg/L NO₃-N) — Probability of moderate nitrate | Proportion (0–1) |
| `RFCprobC3` | P(nitrate ≥ 10 mg/L NO₃-N) — Probability exceeding EPA MCL | Proportion (0–1) |
| `RFC_maxprob` | Column name with highest probability | Character |
| `RFC_maxprob_num` | Value of the highest probability | Proportion (0–1) |
| `Dom_Pop_Percent` | Percentage of population using domestic wells | Percentage (0–100) |

**Notes:** The EPA MCL for nitrate is 10 mg/L as NO₃-N.

### conc_nitrate.csv — Measured Nitrate Concentration Data

| Column | Description | Units / Format |
|--------|-------------|----------------|
| `id` | Geographic identifier (county FIPS code) | Integer |
| `weighted_nitrate` | Population-weighted mean nitrate concentration | mg/L NO₃-N (milligrams per liter as nitrate-nitrogen) |
| `county_pop_served` | County population served by public water | Count (persons) |
| `well` | Well indicator variable | Numeric |
| `Dom_Pop_Percent` | Percentage of population using domestic wells | Percentage (0–100) |

### Health_nitrate.txt — Infant Health Outcome Data

Tab-delimited file with infant health outcomes relevant to nitrate exposure.

| Column | Description | Units / Format |
|--------|-------------|----------------|
| `id` | Geographic identifier (county FIPS code) | Integer |
| `methemoglobin_pct` | Methemoglobin percentage — primary outcome for nitrate exposure | % (methemoglobin >2% is elevated in infants) |
| `infant_age_group` | Age group of the infant, used as a random effect | Character (e.g., "0-3m", "3-6m") |
| `birth_weight_kg` | Birth weight of the infant | kg (kilograms) |

---

## Uranium Data

### prob_uranium.csv — Modeled Uranium Probability Data

Multinomial probability distributions of uranium concentration categories. Same column structure as `prob_dom.csv`.

| Column | Description | Units / Format |
|--------|-------------|----------------|
| `id` | Geographic identifier (county FIPS code) | Integer |
| `pub_water_pop` | Population served by public water supply | Count (persons) |
| `RFCprobC1` | P(uranium < 15 µg/L) — Probability of low uranium | Proportion (0–1) |
| `RFCprobC2` | P(15 ≤ uranium < 30 µg/L) — Probability of moderate uranium | Proportion (0–1) |
| `RFCprobC3` | P(uranium ≥ 30 µg/L) — Probability exceeding EPA MCL | Proportion (0–1) |
| `RFC_maxprob` | Column name with highest probability | Character |
| `RFC_maxprob_num` | Value of the highest probability | Proportion (0–1) |
| `Dom_Pop_Percent` | Percentage of population using domestic wells | Percentage (0–100) |

**Notes:** The EPA MCL for uranium is 30 µg/L.

### conc_uranium.csv — Measured Uranium Concentration Data

| Column | Description | Units / Format |
|--------|-------------|----------------|
| `id` | Geographic identifier (county FIPS code) | Integer |
| `weighted_uranium` | Population-weighted mean uranium concentration | µg/L (micrograms per liter) |
| `county_pop_served` | County population served by public water | Count (persons) |
| `well` | Well indicator variable | Numeric |
| `Dom_Pop_Percent` | Percentage of population using domestic wells | Percentage (0–100) |

### Health_uranium.txt — Kidney Function Data

Tab-delimited file with kidney function outcomes relevant to uranium exposure.

| Column | Description | Units / Format |
|--------|-------------|----------------|
| `id` | Geographic identifier (county FIPS code) | Integer |
| `eGFR` | Estimated Glomerular Filtration Rate — primary outcome for uranium exposure | mL/min/1.73m² (eGFR <60 indicates chronic kidney disease) |
| `age_years` | Age of the participant | Years |
| `sex` | Sex of the participant, used as a random effect | Character ("M" or "F") |

---

## PM2.5 (Air Quality) Data

### prob_pm25.csv — Modeled PM2.5 Probability Data

Multinomial probability distributions of PM2.5 concentration categories from satellite/hybrid models. Same column structure as `prob_dom.csv`. This example demonstrates that the probability data slot can accept any model-based probability data, not just predictions from a specific agency.

| Column | Description | Units / Format |
|--------|-------------|----------------|
| `id` | Geographic identifier (county FIPS code) | Integer |
| `pub_water_pop` | Population in the geographic unit (analogous to well population for weighting) | Count (persons) |
| `RFCprobC1` | P(PM2.5 < 9 µg/m³) — Probability of meeting EPA annual standard | Proportion (0–1) |
| `RFCprobC2` | P(9 ≤ PM2.5 < 35 µg/m³) — Probability of moderate air quality | Proportion (0–1) |
| `RFCprobC3` | P(PM2.5 ≥ 35 µg/m³) — Probability exceeding EPA 24-hour standard | Proportion (0–1) |
| `RFC_maxprob` | Column name with highest probability | Character |
| `RFC_maxprob_num` | Value of the highest probability | Proportion (0–1) |
| `Dom_Pop_Percent` | Percentage used for population weighting | Percentage (0–100) |

**Notes:** EPA annual PM2.5 standard is 9 µg/m³ (updated 2024). EPA 24-hour standard is 35 µg/m³. WHO annual guideline is 5 µg/m³.

### conc_pm25.csv — Ground Monitor PM2.5 Data

| Column | Description | Units / Format |
|--------|-------------|----------------|
| `id` | Geographic identifier (county FIPS code) | Integer |
| `weighted_pm25` | Population-weighted mean PM2.5 concentration from ground monitors | µg/m³ (micrograms per cubic meter) |
| `county_pop` | County population | Count (persons) |
| `monitoring_stations` | Number of air quality monitoring stations in the county | Count |
| `Dom_Pop_Percent` | Percentage used for population weighting | Percentage (0–100) |

### Health_pm25.txt — Respiratory Function Data

Tab-delimited file with respiratory health outcomes relevant to PM2.5 exposure.

| Column | Description | Units / Format |
|--------|-------------|----------------|
| `id` | Geographic identifier (county FIPS code) | Integer |
| `FEV1_pct_predicted` | Forced Expiratory Volume in 1 second as % of predicted value — primary outcome | % predicted (FEV1 <80% often indicates obstructive lung disease) |
| `age_years` | Age of the participant | Years |
| `smoking_status` | Smoking status, used as a random effect | Character (e.g., "never", "former", "current") |
