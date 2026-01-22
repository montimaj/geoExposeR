---
title: 'gwArsenicR: An R package for modeling health effects of drinking water arsenic exposure'
tags:
  - R
  - arsenic
  - groundwater
  - epidemiology
  - environmental health
  - multiple imputation
  - exposure assessment
  - birth outcomes
  - geospatial
  - mixed-effects models
authors:
  - name: Sayantan Majumdar
    orcid: 0000-0002-3539-0147
    affiliation: 1
    corresponding: true
  - name: Scott M. Bartell
    orcid: 0000-0001-7797-2906
    affiliation: 2
  - name: Melissa A. Lombard
    orcid: 0000-0001-5924-6556
    affiliation: 3
  - name: Ryan G. Smith
    orcid: 0000-0002-3747-6868
    affiliation: 4
  - name: Matthew O. Gribble
    orcid: 0000-0002-1614-2981
    affiliation: 5
affiliations:
  - name: Division of Hydrologic Sciences, Desert Research Institute, Reno, NV, USA
    index: 1
  - name: Joe C. Wen School of Population & Public Health, University of California, Irvine, CA, USA
    index: 2
  - name: New England Water Science Center, U.S. Geological Survey, Pembroke, NH, USA
    index: 3
  - name: Department of Civil and Environmental Engineering, Colorado State University, Fort Collins, CO, USA
    index: 4
  - name: Division of Occupational, Environmental, and Climate Medicine, Department of Medicine, University of California, San Francisco, CA, USA
    index: 5
date: 17 January 2026
bibliography: paper.bib
---

# Summary

gwArsenicR is an R package [@R-base] that provides a streamlined workflow for epidemiologists and public health researchers to estimate arsenic exposure from private and public well water and to assess its association with various health outcomes. The package implements a sophisticated statistical approach that combines geospatial arsenic prediction models with a multiple imputation framework to account for exposure uncertainty. The core functionality integrates predicted arsenic concentration probabilities from U.S. Geological Survey (USGS) models for private wells [@Lombard-USGS-2021] with U.S. Environmental Protection Agency (EPA) data for public water systems. It then uses a weighting scheme based on the proportion of the population served by each water source [@Johnson-USGS-2019] to create a unified, county-level exposure probability distribution. By encapsulating these complex methodologies into a single, user-friendly function, gwArsenicR makes this type of analysis more accessible, reproducible, and standardized. While designed for drinking water arsenic exposure assessments, the underlying statistical framework is broadly applicable to other environmental contaminants, such as nitrate, uranium, manganese, per- and polyfluoroalkyl substances (PFAS), or air pollutants like particulate matter (PM~2.5~), where spatially-predicted exposure probabilities and health outcome data are available.

# Statement of Need

Chronic exposure to arsenic in drinking water is a significant public health problem, linked to a range of adverse health outcomes. Assessing this risk is challenging, particularly for the estimated 40 million Americans who rely on unregulated private wells [@DeSimone-2015]. Foundational studies by @Bulka-2022 and @Lombard-2021 have established robust methodologies for this purpose, but their implementation requires specialized expertise in data integration, statistical modeling, and computational methods. This presents a significant barrier for many environmental health researchers who may lack the necessary programming skills to replicate or adapt these complex analyses.

gwArsenicR directly addresses this need by packaging the entire analytical pipeline into an easy-to-use R function. It automates the process of:

1. Converting lognormal EPA arsenic data to a multinomial distribution.
2. Combining USGS and EPA models based on population-use weights.
3. Performing multiple imputation to manage exposure uncertainty.
4. Fitting mixed-effects regression models to health data.
5. Pooling results according to Rubin's Rules [@Rubin-1987].

By providing this accessible tool, gwArsenicR lowers the barrier to entry for conducting high-quality research on the health impacts of groundwater contaminants, promoting reproducibility and wider application of these state-of-the-art methods in environmental epidemiology.

# Design and Implementation

\autoref{fig:workflow} illustrates the complete analytical workflow implemented in gwArsenicR. The pipeline begins with four primary input data sources: USGS arsenic probability rasters [@Lombard-USGS-2021], EPA public water system data, private well population estimates [@Johnson-USGS-2019], and health outcome data. These inputs flow through a data ingestion stage that leverages geospatial packages (terra, sf, exactextractr) [@terra-2024; @sf-2024; @exactextractr-2024] for raster extraction and data.table [@data.table-2025] for efficient data manipulation. The probability integration stage converts EPA lognormal parameters to multinomial distributions compatible with USGS probability categories, then creates population-weighted probability matrices reflecting actual water source usage. Multiple imputation generates replicate datasets with probabilistically assigned arsenic exposure levels, with optional Multivariate Imputation by Chained Equations (MICE) imputation [@mice-2011] for missing covariates. Finally, mixed-effects regression models [@lme4-2015] are fit to each imputed dataset and results are pooled according to Rubin's Rules to produce final estimates with proper uncertainty quantification.

![gwArsenicR analytical workflow showing the data processing pipeline from input sources through statistical modeling to pooled estimates. Packages used at each stage are indicated in brackets. Pooled estimates include regression coefficients ($\beta$), standard errors (SE), 95% confidence intervals (CI), and p-values, all computed using Rubin's Rules to properly account for imputation uncertainty. \label{fig:workflow}](workflow.png)

## Methods

The gwArsenicR package implements the following analytical workflow:

### Step 1: Data Ingestion and Validation

The package loads three input datasets and validates their structure:

1. **USGS probability data** — Multinomial probabilities for arsenic concentration categories (e.g., <5, 5–10, ≥10 µg/L) in private wells, derived from machine learning models [@Lombard-2021; @Lombard-USGS-2021]
2. **EPA lognormal parameters** — Mean and standard deviation (log-scale) of arsenic concentrations in public water systems, plus the percentage of the population using private wells
3. **Health outcome data** — Individual or aggregated health records with geographic identifiers and covariates

A common geographic identifier (e.g., county FIPS code, census tract ID, or participant-specific ID) links all three datasets. This identifier serves as the key for merging arsenic exposure probabilities with health records and defines the geographic unit at which exposure is characterized.

### Step 2: EPA Lognormal-to-Multinomial Conversion

EPA arsenic concentrations are characterized by lognormal distributions. If raw concentration values are provided instead of pre-calculated lognormal parameters, the package automatically computes them using the transformation $\mu = \log(\text{concentration} + 0.1)$ with a configurable default standard deviation ($\sigma = 1.0$ by default). The package then converts these continuous distributions into multinomial probabilities matching the USGS categories by computing:

$$P(\text{As} < c_i) = \Phi\left(\frac{\log(c_i) - \mu}{\sigma}\right)$$

where $\mu$ and $\sigma$ are the lognormal parameters (meanlog, sdlog), $c_i$ are the concentration cutoffs, and $\Phi$ is the standard normal CDF.

### Step 3: Probability Integration via Population Weighting

For each geographic unit, the package creates a unified exposure probability distribution by weighting USGS (private well) and EPA (public water) probabilities according to the population served:

$$P_{\text{combined}}(k) = w \cdot P_{\text{USGS}}(k) + (1-w) \cdot P_{\text{EPA}}(k)$$

where $w$ is the proportion of the population using private wells and $k$ indexes arsenic categories.

### Step 4: Multiple Imputation of Arsenic Exposure

Rather than assigning a single exposure category, the package generates $m$ imputed datasets, where $m$ is controlled by the `ndraws` parameter (default = 10). The arsenic concentration cutoffs defining exposure categories are customizable via the `epa_cutoffs` parameter (default: 5 and 10 µg/L, yielding three categories: <5, 5–10, and ≥10 µg/L).

The arsenic imputation uses a domain-specific probabilistic sampling approach. For each geographic unit, the package retrieves the combined probability distribution $P_{\text{combined}}(k)$ from the weighted USGS+EPA matrix. When multiple raster cells exist within a geographic unit, a population-weighted average of probabilities is computed. The package then samples from this multinomial distribution using R's `sample.int()` function, generating all $m$ draws for each geographic unit simultaneously. Each health record is subsequently merged with its geographic unit's sampled arsenic category based on the matching identifier. Importantly, all individuals within the same geographic unit receive the same arsenic category within a given imputed dataset, but different imputed datasets may assign different categories to the same unit. This approach propagates exposure uncertainty from the underlying geospatial models through all subsequent analyses.

### Step 5: Covariate Imputation (Optional)

If covariates contain missing values, the package applies Multiple Imputation by Chained Equations (MICE) to impute missing data within each arsenic-imputed dataset, using the mice package [@mice-2011] with user-specified parameters.

Unlike the arsenic imputation which draws from pre-computed external probability distributions, MICE imputation operates at the individual level and leverages internal data relationships. The method constructs regression models using available covariates and imputes missing values based on the observed data patterns. MICE iteratively cycles through variables with missing data, drawing from the posterior predictive distribution of each variable conditional on all other variables. Users can configure the number of MICE imputations (`mice_m`), maximum iterations for convergence (`mice_maxit`), and the imputation method (e.g., predictive mean matching via `mice_method = "pmm"`). This two-stage imputation strategy—first arsenic exposure from external probability models, then covariates from internal data relationships—ensures that all sources of uncertainty are properly propagated through the analysis.

### Step 6: Mixed-Effects Regression

For each imputed dataset and each health outcome target, the package fits a linear mixed-effects regression model using lme4 [@lme4-2015]:

$$Y_{ij} = \beta_0 + \beta_1 \text{AsLevel}_{ij} + \mathbf{X}_{ij}\boldsymbol{\gamma} + u_j + \varepsilon_{ij}$$

where $Y_{ij}$ is the health outcome for individual $i$ in geographic unit $j$, $\text{AsLevel}$ is the imputed arsenic category, $\mathbf{X}$ contains covariates, $\boldsymbol{\gamma}$ is the vector of regression coefficients for the covariates, $u_j$ is the random intercept for geographic unit $j$ (typically indexed by the GEOID), and $\varepsilon_{ij}$ is the residual error. The random effects and residuals are assumed to be normally distributed: $u_j \sim N(0, \sigma_u^2)$ and $\varepsilon_{ij} \sim N(0, \sigma^2)$, with $u_j$ and $\varepsilon_{ij}$ independent. The random intercept accounts for within-area correlation among individuals sharing the same geographic identifier, providing proper standard errors for clustered data.

### Step 7: Pooling Results via Rubin's Rules

Regression estimates from all imputed datasets are combined using Rubin's Rules [@Rubin-1987] to obtain final pooled estimates with proper variance that accounts for both within-imputation and between-imputation uncertainty:

$$\bar{Q} = \frac{1}{m}\sum_{i=1}^{m} \hat{Q}_i$$

$$B = \frac{1}{m-1}\sum_{i=1}^{m} (\hat{Q}_i - \bar{Q})^2$$

$$T = \bar{U} + \left(1 + \frac{1}{m}\right)B$$

where $\bar{Q}$ is the pooled estimate, $\bar{U}$ is the mean within-imputation variance, $B$ is the between-imputation variance, $T$ is the total variance, and $m$ is the number of imputations (i.e., `ndraws`). Confidence intervals and p-values are computed using the appropriate degrees of freedom adjustment.

### Step 8: Output Generation

The package returns pooled regression coefficients ($\beta$), standard errors, 95% confidence intervals, and p-values for each exposure category and health outcome, along with optional saved intermediate datasets for diagnostics and reproducibility.

## Software Architecture

gwArsenicR follows a modular design philosophy centered around the main exported function `perform_sensitivity_analysis()`, which orchestrates the complete analytical workflow. The package architecture consists of six main modules:

* **Main Analysis Module** (`gwArsenic.R`): Contains the primary `perform_sensitivity_analysis()` function that orchestrates the complete workflow.
* **Data Ingestion Module** (`data-ingestion.R`): Provides six exported functions for preparing input data from raw sources, including USGS raster extraction, EPA data processing, and health data formatting.
* **Data Loading Module** (`data-loading.R`): Handles import and preprocessing of USGS probability data, EPA lognormal parameters, and health outcome datasets with comprehensive input validation.
* **Imputation Module** (`imputation.R`): Implements multiple imputation for arsenic exposure levels and optional MICE [@mice-2011] imputation for missing covariates.
* **Regression Module** (`regression.R`): Fits linear mixed-effects models using lme4 [@lme4-2015] and pools results according to Rubin's Rules using Amelia [@Amelia-2011].
* **Utilities Module** (`utils.R`): Provides data formatting, validation, and helper functions.

## Key Features and Capabilities

**Flexible Data Integration**: The package accommodates the different probability structures of USGS and EPA data by converting EPA lognormal parameters to multinomial probabilities compatible with USGS discrete probability categories. Population-weighted averaging creates unified exposure probability distributions that reflect the actual water source usage patterns in each geographic area. The package supports any geographic identifier format—including county FIPS codes, census tract IDs, ZIP codes, or even participant-specific identifiers—enabling applications at multiple spatial scales. Smart identifier formatting (enabled by default via `format_ids = TRUE`) automatically detects ID types: numeric FIPS-like codes (1-5 digits) receive zero-padding for consistency, while alphanumeric identifiers (e.g., census tracts, participant IDs) are preserved unchanged. This feature can be disabled with `format_ids = FALSE`.

**Automatic EPA Parameter Calculation**: For users who have raw arsenic concentration data instead of pre-calculated lognormal parameters, gwArsenicR can automatically compute meanlog and sdlog values. By specifying the `epa_as_conc_col` parameter, the package applies the transformation `meanlog = log(concentration + 0.1)` with a configurable `epa_default_sdlog` (default: 1.0), eliminating the need for manual preprocessing.

**Comprehensive Data Ingestion Pipeline**: The package exports seven functions total—one main analysis function and six data preparation functions:

* `perform_sensitivity_analysis()`: Main analysis function that orchestrates the complete workflow including data integration, imputation, regression, and result pooling.
* `prepare_usgs_data()`: Prepares USGS arsenic probability data from raster files, extracting county-level statistics with optional population weighting.
* `prepare_epa_data()`: Processes EPA public water system data, converting raw concentrations to lognormal parameters or handling pre-summarized data.
* `format_health_data()`: Formats health outcome datasets with proper geographic identifiers and validates required columns.
* `add_well_percentages()`: Merges private well usage percentages from census data to target datasets for population weighting.
* `subset_by_geography()`: Filters datasets by state FIPS codes, specific geographic IDs, or pattern matching.
* `validate_prepared_inputs()`: Validates all prepared datasets for consistency, checking column alignment and identifier matching before analysis.

**Robust Multiple Imputation**: gwArsenicR implements a two-stage imputation strategy. First, it performs probabilistic assignment of arsenic exposure categories based on the integrated probability distributions. Second, it optionally implements MICE [@mice-2011] for missing demographic and health covariates, ensuring that all sources of uncertainty are properly propagated through the analysis.

**Hierarchical Modeling**: The package leverages lme4 [@lme4-2015] to fit linear mixed-effects models that account for the hierarchical structure typical of population health data (e.g., individuals nested within counties, counties within states). This approach provides more accurate standard errors and better accounts for geographic clustering.

**Statistical Rigor**: Results from multiple imputed datasets are pooled using Rubin's Rules [@Rubin-1987], which properly combines point estimates and accounts for both within-imputation and between-imputation variance. This ensures that confidence intervals and p-values correctly reflect the uncertainty introduced by the imputation process.

**Performance Optimization**: The package uses data.table [@data.table-2025] for efficient handling of large datasets and implements vectorized operations for probability calculations, making it practical for national-scale epidemiological studies.

**Quality Assurance**: gwArsenicR includes comprehensive testing with 476 tests achieving 100% code coverage, automated continuous integration across multiple R versions and operating systems (macOS, Ubuntu, Windows), and extensive input validation to prevent common user errors.

# Research Enabled by gwArsenicR

By simplifying a complex analytical workflow, gwArsenicR can accelerate and broaden the scope of research in environmental health. The package enables researchers to:

* **Rapidly Test Hypotheses**: Easily investigate the health effects of arsenic exposure on a wide range of outcomes, such as birth weight, gestational age, cancer incidence, or cardiovascular disease, using large, nationally representative datasets.
* **Conduct Subgroup Analyses**: Modify the regression formula to explore whether associations differ across specific subpopulations (e.g., by race/ethnicity, socioeconomic status, or geographic region), which can help identify vulnerable groups.
* **Inform Public Health Policy**: Generate robust evidence that can inform policy decisions, such as establishing state-level guidelines for private well testing or targeting public health interventions in high-risk areas.
* **Serve as a Methodological Template**: The package framework can be adapted to study other environmental contaminants where exposure data is uncertain and requires integration from multiple sources.

Ultimately, gwArsenicR empowers a wider research community to contribute to our understanding of the health risks posed by groundwater contaminants, fostering more timely and impactful science.

# Broader Applicability

While gwArsenicR is designed for arsenic exposure assessment, the underlying framework is scale-agnostic and contaminant-agnostic, making it adaptable to a wide range of environmental epidemiology applications. The package operates at any geographic scale—county, census tract, ZIP code, or even individual participant level—provided that identifiers match across input datasets. This flexibility enables analyses ranging from national-scale studies using county FIPS codes to fine-grained urban analyses using census tracts. Beyond arsenic, the core methodology can be applied to other groundwater contaminants (e.g., uranium, nitrate, manganese, PFAS), air quality exposures (e.g., PM2.5, ozone, wildfire smoke), or any environmental exposure assessment problem that shares similar characteristics: probabilistic exposure estimates from multiple data sources, population-weighted integration of those sources, and the need to propagate uncertainty through statistical inference. Users adapting the framework to new domains would prepare analogous probability and distributional parameter datasets for their contaminant of interest, modify the concentration cutoffs to match relevant exposure thresholds, and adjust the regression formula for appropriate health outcomes. By providing a generalizable template for integrating uncertain exposure data with health outcomes, gwArsenicR can serve as a foundation for rigorous, reproducible environmental health research across multiple domains.

# Availability

gwArsenicR is available on GitLab at https://code.usgs.gov/envirohealth/gwArsenicR under a CC0 1.0 Universal license. Comprehensive documentation, including a detailed vignette with worked examples and full API reference, is available at https://envirohealth.code-pages.usgs.gov/gwArsenicR/. The package requires R >= 4.4.0 and can be installed via:

```r
# Install remotes if not already installed
install.packages("remotes")

# Install gwArsenicR from GitLab (use a specific tag or branch)
remotes::install_gitlab(
  repo = "envirohealth/gwArsenicR@v1.0.0",
  host = "code.usgs.gov"
)
```

# Example Usage

The following example demonstrates a basic analysis workflow:

```r
library(gwArsenicR)

# Run sensitivity analysis with 5 imputed datasets
results <- perform_sensitivity_analysis(
  ndraws = 5,
  as_usgs_prob_csv = "usgs_arsenic_probabilities.csv",
  as_epa_prob_csv = "epa_lognormal_parameters.csv",
  birth_data_txt = "birth_outcomes.txt",
  regression_formula = "BWT ~ As_cat + MAGE + MRACE + (1|state)",
  output_dir = "analysis_results/",
  targets = c("BWT", "OEGEST"),
  impute_vars = c("MAGE", "MRACE")
)

# Access pooled regression results
results$pooled_results
```

The package also provides data preparation functions for users starting with raw GeoTIFF rasters or unprocessed EPA data. Complete examples demonstrating the full data pipeline are available in the package vignette and the `examples/` directory on GitLab.

# Related Software

While R provides excellent packages for multiple imputation [@mice-2011], mixed-effects modeling [@lme4-2015], and geospatial analysis [@terra-2024; @sf-2024], no existing package provides an integrated workflow specifically designed for arsenic exposure assessment that combines USGS and EPA data sources. General-purpose imputation packages like mice and Amelia [@Amelia-2011] handle missing data but do not address the domain-specific challenge of integrating probabilistic exposure estimates from multiple environmental models. Similarly, packages like exactextractr [@exactextractr-2024] provide efficient raster extraction but require substantial custom code to implement the full analytical pipeline described by @Bulka-2022 and @Lombard-2021. gwArsenicR fills this gap by providing a domain-specific solution that handles the unique data integration and uncertainty quantification challenges in groundwater contaminant epidemiology, reducing implementation time from weeks to hours.

# Limitations

The current implementation has the following limitations that users should consider. First, the package is designed primarily for U.S.-based studies using USGS probability rasters and EPA public water system data; international applications would require adaptation of the data ingestion functions to accommodate different data formats and geographic coding systems. Second, the default analysis assumes arsenic exposure can be approximated at the county level using population-weighted averages; finer-scale exposure assessment (e.g., census tract or individual well level) is supported through flexible identifier handling but may require additional data preparation. Third, the mixed-effects regression framework assumes linear relationships; users investigating non-linear exposure-response relationships should consider extensions or post-hoc analyses. Finally, computational requirements scale with the number of imputed datasets and the size of the health outcome data; national-scale analyses with many imputations may require substantial memory and processing time.

# Author Contributions

**Sayantan Majumdar**: Software, Methodology, Validation, Writing - Original Draft, Funding Acquisition. **Scott M. Bartell**: Methodology, Writing - Review & Editing. **Melissa A. Lombard**: Data Curation, Methodology, Resources, Writing - Review & Editing. **Ryan G. Smith**: Methodology, Funding Acquisition, Writing - Review & Editing, Project Administration. **Matthew O. Gribble**: Conceptualization, Funding Acquisition, Project Administration, Methodology, Supervision, Writing - Review & Editing.

# Acknowledgments

This work was supported by the National Heart, Lung, and Blood Institute (R21HL159574) and funding from the United States Geological Survey's John Wesley Powell Center for Analysis and Synthesis. Any use of trade, firm, or product names is for descriptive purposes only and does not imply endorsement by the U.S. Government.

# References
