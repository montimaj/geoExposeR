---
title: 'geoExposeR: An R package for modeling health effects of environmental exposures'
tags:
  - R
  - contaminants
  - pollutants
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
date: March 20 2026
bibliography: paper.bib
---

# SUMMARY

geoExposeR is an R package [@R-base] that provides a streamlined workflow for epidemiologists and public health researchers to assess associations between environmental exposures and health outcomes while properly accounting for exposure uncertainty. The package implements a general statistical framework that combines two complementary types of exposure data, namely modeled probability estimates (e.g., from machine learning or geostatistical models) and measured concentration data (e.g., from monitoring networks), using population-weighted integration. It then performs multiple imputation to propagate exposure uncertainty through the analysis, fits mixed-effects regression models to health data, and pools results using Rubin's Rules [@Rubin-1987]. This framework is applicable to any environmental contaminant for which probabilistic exposure estimates and health outcome data are available, including drinking water contaminants (arsenic, nitrate, uranium, manganese, PFAS), air pollutants (PM2.5, ozone), and other environmental exposures. The package supports flexible geographic scales (e.g., county, census tract, ZIP code, or individual level) and includes worked examples for arsenic, nitrate, uranium, and PM2.5 to demonstrate its generalizability. The motivating application is drinking water arsenic exposure assessment, integrating predicted arsenic concentration probabilities from U.S. Geological Survey (USGS) machine learning models for private wells [@Lombard-USGS-2021] with U.S. Environmental Protection Agency (EPA) monitoring data [@Nigra2025; @Angley2026; @Ravalli2022] for public water systems and a weighting scheme based on private well usage [@Johnson-USGS-2019]. By encapsulating complex data integration, imputation, and regression methodologies into user-friendly functions, geoExposeR makes rigorous environmental health research more accessible, reproducible, and standardized.

# STATEMENT OF NEED

Environmental epidemiologists face a fundamental challenge: exposure to environmental contaminants is rarely measured directly for study participants. Instead, researchers must estimate exposures using predictive models, monitoring data, or biomarkers, each with inherent uncertainty. Failing to account for this uncertainty can lead to biased effect estimates and misleading conclusions about health risks.

This challenge is particularly acute for drinking water contaminants like arsenic, due to widespread lack of monitoring. An estimated 40 million Americans rely on unregulated private wells for their drinking water [@DeSimone-2014]. Foundational studies by @Bulka-2022 and @Lombard-2021 developed sophisticated methods to address exposure uncertainty by combining USGS machine learning predictions with EPA monitoring data through population-weighted integration and multiple imputation. However, implementing these methods requires specialized expertise in geospatial analysis, statistical modeling, and R programming, which presents a significant barrier for many researchers.

geoExposeR directly addresses this need by packaging the complete analytical pipeline into accessible R functions. The package automates:

1. Integration of probabilistic exposure estimates from multiple data sources
2. Population-weighted combination of exposure probabilities
3. Multiple imputation to propagate exposure uncertainty
4. Mixed-effects regression modeling with proper clustering adjustment
5. Result pooling via Rubin's Rules [@Rubin-1987] with appropriate variance estimation

Beyond arsenic, researchers studying other environmental exposures face similar methodological challenges. The package includes worked examples for nitrate (infant methemoglobinemia), uranium (kidney function), and PM2.5 (respiratory outcomes), demonstrating how the framework generalizes to other contaminants and health endpoints. By lowering technical barriers to rigorous uncertainty quantification, geoExposeR promotes reproducibility and enables wider application of state-of-the-art methods in environmental epidemiology.

# STATE OF THE FIELD

Several R packages address components of the environmental epidemiology workflow, but none provides an integrated solution for probabilistic exposure assessment with uncertainty propagation. The mice package [@mice-2011] offers flexible multiple imputation but lacks domain-specific functionality for integrating geospatial exposure models. Amelia [@Amelia-2011] provides efficient imputation for cross-sectional and time-series data but does not handle multinomial probability distributions from environmental models. Geospatial packages like terra [@terra-2024], sf [@sf-2024], and exactextractr [@exactextractr-2024] enable sophisticated raster extraction but require substantial custom code to implement population-weighted probability integration. The lme4 package [@lme4-2015] fits mixed-effects models but does not address the upstream challenge of exposure assignment under uncertainty.

We chose to build geoExposeR rather than contribute to existing packages because the core challenge of integrating probabilistic exposure estimates from multiple complementary data sources and propagating that uncertainty through regression requires a cohesive workflow spanning data ingestion, probability theory, imputation, and hierarchical modeling. This integration represents our unique scholarly contribution: operationalizing the methodological advances of @Bulka-2022 and @Lombard-2021 into a reusable, validated software framework. Existing packages excel at individual components but cannot be easily composed to implement this specific analytical paradigm. geoExposeR bridges this gap by providing a domain-specific solution while leveraging established packages (mice, lme4, data.table) for specialized computations.

# DESIGN AND IMPLEMENTATION

We describe each step of the geoExposeR workflow in general terms, applicable to any environmental contaminant. Each step concludes with a paragraph illustrating the step using the drinking water arsenic case study, which is the motivating application for which geoExposeR was originally developed. The same framework applies to other contaminants (nitrate, uranium, PM2.5, etc.) by substituting the appropriate data sources and concentration thresholds.

## METHODS

\autoref{fig:workflow} illustrates the complete analytical workflow implemented in geoExposeR. The pipeline begins with four primary input data sources: modeled probability estimates for the contaminant of interest, measured concentration data from monitoring networks, population weighting data reflecting the fraction of the population exposed through each source (e.g., private well usage estimates from @Johnson-USGS-2019), and health outcome data. These inputs flow through a data ingestion stage that leverages geospatial packages (terra, sf, exactextractr) [@terra-2024; @sf-2024; @exactextractr-2024] for raster extraction and data.table [@data.table-2025] for efficient data manipulation. The probability integration stage converts measured concentration lognormal parameters to multinomial distributions compatible with modeled probability categories, then creates population-weighted probability matrices reflecting actual exposure source usage. If the concentration data includes raw measurements rather than lognormal parameters, the package automatically converts these to a lognormal distribution using a default standard deviation of 1.0 unless the user specifies a different value. Next, multiple imputation generates replicate datasets with probabilistically assigned exposure levels, with optional Multivariate Imputation by Chained Equations (MICE) [@mice-2011] for missing covariates. Finally, mixed-effects regression models [@lme4-2015] are fit to each imputed dataset and the results are pooled according to Rubin's Rules to produce final estimates with proper uncertainty quantification. The complete step-by-step workflow is as follows:

![geoExposeR analytical workflow showing the data processing pipeline from input sources through statistical modeling to pooled estimates. Packages used at each stage are indicated in brackets. Pooled estimates include regression coefficients ($\beta$), standard errors (SE), 95% confidence intervals (CI), and p-values, all computed using Rubin's Rules to properly account for imputation uncertainty. \label{fig:workflow}](workflow.png)


### STEP 1: DATA INGESTION AND VALIDATION

The package loads four input datasets and validates their structure:

**1. Modeled probability data** (`model_prob_csv`) --- Multinomial probabilities for contaminant concentration categories, typically derived from predictive or machine learning models applied to spatial covariates.

**2. Measured concentration data** (`conc_data_csv`) --- Lognormal distribution parameters (meanlog, sdlog) describing contaminant concentrations from monitoring networks.

**3. Population weighting data** --- The proportion of the population exposed through each source (e.g., private wells vs. public water systems) for each geographic unit. This dataset is merged into the concentration data during data preparation using the `add_well_percentages()` function, which joins on geographic identifiers and produces the weighting column (specified via the `pwell_col` parameter) used in the probability integration step.

**4. Health outcome data** (`birth_data_txt`) --- Individual or aggregated health records with geographic identifiers and covariates.

A common geographic identifier (e.g., county FIPS code, census tract ID, ZIP code, or participant-specific ID) links all four datasets. This identifier serves as the key for merging exposure probabilities with health records and defines the geographic unit at which exposure is characterized.

In the case of the arsenic example, modeled probability data are multinomial probabilities for arsenic concentration categories (<5, 5–10, ≥10 µg/L) in private wells, derived from USGS machine learning models [@Lombard-2021; @Lombard-USGS-2021]. Measured concentration data are log-scale mean and standard deviation of arsenic concentrations in public water systems from EPA monitoring [@Nigra2025; @Angley2026; @Ravalli2022], with the percentage of the population using private wells as the weighting variable. The common identifier is the county FIPS code.

### STEP 2: CONCENTRATION-TO-MULTINOMIAL CONVERSION

Measured concentrations are typically characterized by lognormal distributions. If raw concentration values are provided instead of pre-calculated lognormal parameters, the package automatically computes them using the transformation $\mu = \log(\text{concentration} + 0.1)$ with a configurable standard deviation ($\sigma = 1.0$ by default). The package then converts these continuous distributions into multinomial probabilities matching the modeled probability categories by computing:

$$P(\text{conc} < c_i) = \Phi\left(\frac{\log(c_i) - \mu}{\sigma}\right)$$

where $\mu$ and $\sigma$ are the lognormal parameters (meanlog, sdlog), $c_i$ are the concentration cutoffs defined by the user via the `conc_cutoffs` parameter, and $\Phi$ is the standard normal CDF.

In the case of the arsenic example, concentration data from the EPA [@Nigra2025; @Angley2026; @Ravalli2022] were converted to multinomial probabilities matching the USGS model categories using the concentration cutoffs 5 and 10 µg/L, yielding three categories: <5 (below detection/low exposure), 5–10 (moderate), and ≥10 µg/L (above the EPA Maximum Contaminant Level).

### STEP 3: PROBABILITY INTEGRATION VIA POPULATION WEIGHTING

For each geographic unit, the package creates a unified exposure probability distribution by weighting modeled and measured probabilities according to the population served by each exposure source:

$$P_{\text{combined}}(k) = w \cdot P_{\text{model}}(k) + (1-w) \cdot P_{\text{measured}}(k)$$

where $w$ is the proportion of the population exposed through the source represented by the modeled data (obtained from the population weighting data ingested in Step 1 via `add_well_percentages()`), $P_{\text{model}}(k)$ are the modeled multinomial probabilities, $P_{\text{measured}}(k)$ are the probabilities derived from measured concentration data, and $k$ indexes exposure categories.

In the case of the arsenic example, $w$ represents the proportion of the population using private wells [@Johnson-USGS-2019], $P_{\text{model}}(k)$ are the USGS machine learning predictions for private well arsenic concentrations [@Lombard-2021; @Lombard-USGS-2021], and $P_{\text{measured}}(k)$ are the EPA-derived probabilities for public water systems [@Nigra2025; @Angley2026; @Ravalli2022].

### STEP 4: MULTIPLE IMPUTATION OF EXPOSURE

Rather than assigning a single exposure category, the package generates $m$ imputed datasets, where $m$ is controlled by the `ndraws` parameter (default = 10). The concentration cutoffs defining exposure categories are customizable via the `conc_cutoffs` parameter.

The exposure imputation uses a domain-specific probabilistic sampling approach. For each geographic unit, the package retrieves the combined probability distribution $P_{\text{combined}}(k)$ from the weighted probability matrix. When multiple raster cells exist within a geographic unit, a population-weighted average of probabilities is computed. The package then samples from this multinomial distribution using R's `sample.int()` function, generating all $m$ draws for each geographic unit simultaneously. Each health record is subsequently merged with its geographic unit's sampled exposure category based on the matching identifier. Importantly, all individuals within the same geographic unit receive the same exposure category within a given imputed dataset, but different imputed datasets may assign different categories to the same unit. This approach propagates exposure uncertainty from the underlying geospatial models through all subsequent analyses.

In the case of the arsenic example, the default cutoffs of 5 and 10 µg/L yield three arsenic categories (<5, 5–10, ≥10 µg/L). Each imputed dataset assigns one of these categories to every geographic unit by sampling from the combined USGS/EPA probability distribution, so that uncertainty in arsenic exposure classification is fully propagated.

### STEP 5: COVARIATE IMPUTATION (OPTIONAL)

If covariates contain missing values, the package applies Multiple Imputation by Chained Equations (MICE) to impute missing data within each exposure-imputed dataset, using the mice package [@mice-2011] with user-specified parameters.

Unlike the exposure imputation which draws from pre-computed external probability distributions, MICE imputation operates at the individual level and leverages internal data relationships. The method constructs regression models using available covariates and imputes missing values based on the observed data patterns. MICE iteratively cycles through variables with missing data, drawing from the posterior predictive distribution of each variable conditional on all other variables. Users can configure the number of MICE imputations (`mice_m`), maximum iterations for convergence (`mice_maxit`), and the imputation method (e.g., predictive mean matching via `mice_method = "pmm"`). This two-stage imputation strategy (first exposure levels from external probability models, then covariates from internal data relationships) ensures that all sources of uncertainty are properly propagated through the analysis.

### STEP 6: MIXED-EFFECTS REGRESSION

For each imputed dataset and each health outcome target, the package fits a linear mixed-effects regression model using lme4 [@lme4-2015]:

$$Y_{ij} = \beta_0 + \beta_1 \text{ExpLevel}_{ij} + \mathbf{X}_{ij}\boldsymbol{\gamma} + u_j + \varepsilon_{ij}$$

where $Y_{ij}$ is the health outcome for individual $i$ in geographic unit $j$, $\beta_0$ is the intercept, $\text{ExpLevel}$ is the imputed exposure category, $\beta_1$ is the coefficient for the exposure category, $\mathbf{X}$ contains covariates, $\boldsymbol{\gamma}$ is the vector of regression coefficients for the covariates, $u_j$ is the random intercept for geographic unit $j$ (typically indexed by the GEOID), and $\varepsilon_{ij}$ is the residual error. The random effects and residuals are assumed to be normally distributed: $u_j \sim N(0, \sigma_u^2)$ and $\varepsilon_{ij} \sim N(0, \sigma^2)$, with $u_j$ and $\varepsilon_{ij}$ independent. The random intercept accounts for within-area correlation among individuals sharing the same geographic identifier, providing proper standard errors for clustered data.

### STEP 7: POOLING RESULTS VIA RUBIN'S RULES

Regression estimates from all imputed datasets are combined using Rubin's Rules [@Rubin-1987] to obtain final pooled estimates with proper variance that accounts for both within-imputation and between-imputation uncertainty:

$$\bar{Q} = \frac{1}{m}\sum_{i=1}^{m} \hat{Q}_i$$

$$B = \frac{1}{m-1}\sum_{i=1}^{m} (\hat{Q}_i - \bar{Q})^2$$

$$T = \bar{U} + \left(1 + \frac{1}{m}\right)B$$

where $\hat{Q}_i$ is the point estimate from the $i$-th imputed dataset, $\bar{Q}$ is the pooled estimate (i.e., the average of the point estimates across all imputed datasets), $\bar{U}$ is the mean within-imputation variance, $B$ is the between-imputation variance, $T$ is the total variance, and $m$ is the number of imputations (i.e., `ndraws`). Confidence intervals and p-values are computed using the appropriate degrees of freedom adjustment.

### STEP 8: OUTPUT GENERATION

The package returns pooled regression coefficients ($\beta$), standard errors, 95% confidence intervals, and p-values for each exposure category and health outcome, along with optional saved intermediate datasets for diagnostics and reproducibility.

## SOFTWARE ARCHITECTURE

geoExposeR follows a modular design philosophy centered around the main exported function `perform_sensitivity_analysis()`, which orchestrates the complete analytical workflow. The package architecture consists of six main modules:

**1. Main Analysis Module** (`geoExposeR.R`): Contains the primary `perform_sensitivity_analysis()` function that orchestrates the complete workflow.

**2. Data Ingestion Module** (`data-ingestion.R`): Provides six exported functions for preparing input data from raw sources, including probability raster extraction, concentration data processing, and health data formatting.

**3. Data Loading Module** (`data-loading.R`): Handles import and preprocessing of modeled probability data, measured concentration parameters, and health outcome datasets with comprehensive input validation.

**4. Imputation Module** (`imputation.R`): Implements multiple imputation for exposure levels and optional MICE [@mice-2011] imputation for missing covariates.

**5. Regression Module** (`regression.R`): Fits linear mixed-effects models using lme4 [@lme4-2015] and pools results according to Rubin's Rules using Amelia [@Amelia-2011].

**6. Utilities Module** (`utils.R`): Provides data formatting, validation, and helper functions using Hmisc [@Hmisc-2024].

## SOFTWARE DESIGN

We made several deliberate design trade-offs to balance flexibility, usability, and statistical rigor. First, we chose a "pipeline" architecture with a single entry point (`perform_sensitivity_analysis()`) rather than exposing granular functions, reducing the learning curve for epidemiologists who need results, not programming flexibility. However, we also export six data preparation functions for advanced users requiring custom preprocessing workflows.

Second, we implemented a two-stage imputation strategy (external probability-based imputation for exposures, followed by optional MICE for covariates) rather than a unified joint imputation model. This separation reflects the distinct data-generating processes: exposure probabilities derive from validated external models, while covariate missingness arises from internal data collection. Joint imputation would conflate these mechanisms and complicate interpretation.

Third, we chose data.table over tidyverse for internal data manipulation, prioritizing computational efficiency for national-scale datasets (millions of birth records) over code readability. This enables analyses that would be impractical with memory-intensive alternatives.

Finally, we invested heavily in input validation (52 validation functions) and comprehensive testing (476 tests, 100% coverage) to fail early with informative errors. Environmental health analyses inform policy decisions; silent errors in exposure assignment could have real-world consequences. This defensive design philosophy ensures that users receive correct results or clear diagnostics, never spurious outputs.

## KEY FEATURES AND CAPABILITIES

**1. Flexible Data Integration**: The package accommodates the different probability structures of modeled and measured data by converting measured concentration lognormal parameters to multinomial probabilities compatible with modeled discrete probability categories. Population-weighted averaging creates unified exposure probability distributions that reflect the actual exposure source usage patterns in each geographic area. The package supports any geographic identifier format, including county FIPS codes, census tract IDs, ZIP codes, or even participant-specific identifiers, enabling applications at multiple spatial scales. Smart identifier formatting (enabled by default via `format_ids = TRUE`) automatically detects ID types: numeric FIPS-like codes (1-5 digits) receive zero-padding for consistency, while alphanumeric identifiers (e.g., census tracts, participant IDs) are preserved unchanged. This feature can be disabled with `format_ids = FALSE`.

**2. Automatic Concentration Parameter Calculation**: For users who have raw concentration data instead of pre-calculated lognormal parameters, geoExposeR can automatically compute meanlog and sdlog values. By specifying the `conc_raw_col` parameter, the package applies the transformation `meanlog = log(concentration + 0.1)` with a configurable `default_sdlog` (default: 1.0), eliminating the need for manual preprocessing.

**3. Comprehensive Data Ingestion Pipeline**: The package exports seven functions total: one main analysis function and six data preparation functions.

* `perform_sensitivity_analysis()`: Main analysis function that orchestrates the complete workflow including data integration, imputation, regression, and result pooling.
* `prepare_prob_data()`: Prepares modeled probability data from raster files, extracting geographic-unit-level statistics with optional population weighting.
* `prepare_conc_data()`: Processes measured concentration data, converting raw concentrations to lognormal parameters or handling pre-summarized data.
* `format_health_data()`: Formats health outcome datasets with proper geographic identifiers and validates required columns.
* `add_well_percentages()`: Merges private well usage percentages from census data to target datasets for population weighting.
* `subset_by_geography()`: Filters datasets by state FIPS codes, specific geographic IDs, or pattern matching.
* `validate_prepared_inputs()`: Validates all prepared datasets for consistency, checking column alignment and identifier matching before analysis.

**4. Robust Multiple Imputation**: geoExposeR implements a two-stage imputation strategy. First, it performs probabilistic assignment of exposure categories based on the integrated probability distributions. Second, it optionally implements MICE [@mice-2011] for missing demographic and health covariates, ensuring that all sources of uncertainty are properly propagated through the analysis.

**5. Hierarchical Modeling**: The package leverages lme4 [@lme4-2015] to fit linear mixed-effects models that account for the hierarchical structure typical of population health data (e.g., individuals nested within counties, counties within states). This approach provides more accurate standard errors and better accounts for geographic clustering.

**6. Statistical Rigor**: Results from multiple imputed datasets are pooled using Rubin's Rules [@Rubin-1987], which properly combines point estimates and accounts for both within-imputation and between-imputation variance. This ensures that confidence intervals and p-values correctly reflect the uncertainty introduced by the imputation process.

**7. Performance Optimization**: The package uses data.table [@data.table-2025] for efficient handling of large datasets and implements vectorized operations for probability calculations, making it practical for national-scale epidemiological studies. The package is designed to run locally on a standard desktop or laptop; county-level analyses with tens of thousands of health records and 5–10 imputations typically complete within minutes. Larger analyses, such as national-scale studies with millions of records or a high number of imputations, may benefit from machines with additional memory (16+ GB RAM) but do not require high-performance computing (HPC) resources.

**8. Quality Assurance**: geoExposeR includes comprehensive testing with 476 tests achieving 100% code coverage using testthat [@testthat-2024], automated continuous integration across multiple R versions and operating systems (macOS, Ubuntu, Windows), and extensive input validation to prevent common user errors. The package follows standard R package development conventions [@Wickham-2023].

# RESEARCH ENABLED BY geoExposeR

By simplifying a complex analytical workflow, geoExposeR can accelerate and broaden the scope of research in environmental health. The package enables researchers to:

**1. Rapidly Test Hypotheses**: Investigate the health effects of environmental exposures on a wide range of outcomes (birth weight, gestational age, cancer incidence, cardiovascular disease, kidney function, respiratory health) using large, nationally representative datasets without extensive custom programming.

**2. Conduct Subgroup Analyses**: Modify the regression formula to explore whether associations differ across specific subpopulations (e.g., by race/ethnicity, socioeconomic status, or geographic region), helping identify vulnerable groups who may bear disproportionate exposure burdens.

**3. Perform Multi-Contaminant Studies**: Apply the same analytical framework across multiple environmental exposures (e.g., arsenic, nitrate, uranium, PM2.5) to compare effect sizes, identify co-exposure patterns, or prioritize contaminants for regulatory attention.

**4. Quantify Exposure Uncertainty**: Propagate uncertainty from probabilistic exposure models through regression analyses, providing confidence intervals that honestly reflect the limitations of indirect exposure assessment rather than treating estimated exposures as known values.

**5. Inform Public Health Policy**: Generate robust evidence that can inform policy decisions, such as establishing guidelines for private well testing, setting regulatory standards, or targeting public health interventions in high-risk areas.

Ultimately, geoExposeR empowers a wider research community to contribute to our understanding of the health risks posed by environmental contaminants, fostering more timely, rigorous, and reproducible science.

# RESEARCH IMPACT STATEMENT

geoExposeR operationalizes methods from two foundational publications in environmental epidemiology. @Bulka-2022 introduced the statistical framework for integrating USGS private well predictions with EPA public water system data [@Ravalli2022] using population weighting and multiple imputation, while @Lombard-2021 developed the machine learning models that generate probabilistic arsenic concentration estimates. These works have collectively been cited over 100 times and inform ongoing policy discussions about private well contamination. geoExposeR makes these validated methods accessible to the broader research community.

The package demonstrates community-readiness through several signals: comprehensive documentation including a detailed vignette with worked examples, 476 automated tests achieving 100% code coverage, continuous integration across three operating systems (macOS, Ubuntu, Windows) and multiple R versions, and worked examples spanning four contaminants (arsenic, nitrate, uranium, PM2.5) with complete reproducible workflows. The package is hosted on USGS GitLab with full version control history and follows standard R package conventions for installation and use.

Analyses using geoExposeR methodology are currently being prepared for submission to peer-reviewed epidemiology journals, examining associations between drinking water arsenic and birth outcomes in U.S. populations. The package reduces implementation time from weeks of custom programming to hours, enabling researchers to focus on scientific questions rather than software engineering.

# BROADER APPLICABILITY

While geoExposeR was originally developed for arsenic exposure assessment, the underlying framework is scale-agnostic and contaminant-agnostic, making it directly applicable to a wide range of environmental epidemiology applications. The package operates at any geographic scale (county, census tract, ZIP code, or even individual participant level) provided that identifiers match across input datasets. This flexibility enables analyses ranging from national-scale studies using county FIPS codes to fine-grained urban analyses using census tracts. Beyond arsenic, the core methodology can be applied to other groundwater contaminants (e.g., uranium, nitrate, manganese, PFAS), air quality exposures (e.g., PM2.5, ozone, wildfire smoke), or any environmental exposure assessment problem that shares similar characteristics: probabilistic exposure estimates from multiple data sources, population-weighted integration of those sources, and the need to propagate uncertainty through statistical inference.

To demonstrate this broader applicability, the package includes three additional worked examples in the `examples/` directory with dummy datasets:

**1. Nitrate and infant methemoglobinemia** (`run_geoExposeR_nitrate.R`): Analyzes drinking water nitrate exposure using cutoffs of 5 and 10 mg/L NO3-N (EPA Maximum Contaminant Level [MCL] = 10 mg/L) and its association with infant methemoglobin levels [@EPA-DWR-2024].

**2. Uranium and kidney function** (`run_geoExposeR_uranium.R`): Assesses uranium nephrotoxicity using cutoffs of 15 and 30 µg/L (EPA MCL = 30 µg/L) with estimated glomerular filtration rate (eGFR) as the health outcome [@EPA-DWR-2024].

**3. PM2.5 air quality and respiratory function** (`run_geoExposeR_pm25.R`): Extends the framework beyond groundwater to air quality, using PM2.5 cutoffs of 9 and 35 µg/m³ (EPA annual standard = 9 µg/m³; 24-hour standard = 35 µg/m³) with FEV1 percent predicted (forced expiratory volume in 1 second as a percentage of the expected value for a healthy person of the same age, sex, and height) as the respiratory health outcome. This example demonstrates how the modeled probability and measured concentration data parameters can be reinterpreted as satellite/model predictions and ground monitoring data, respectively [@EPA-NAAQS-PM-2024].

Users adapting the framework to new domains would prepare analogous probability and distributional parameter datasets for their contaminant of interest, modify the concentration cutoffs via `conc_cutoffs` to match relevant exposure thresholds (e.g., regulatory standards), customize category labels via `cat_label`, and adjust the regression formula for appropriate health outcomes. By providing a generalizable template for integrating uncertain exposure data with health outcomes, geoExposeR can serve as a foundation for rigorous, reproducible environmental health research across multiple domains.

# AVAILABILITY

geoExposeR is available on GitLab at https://code.usgs.gov/envirohealth/geoExposeR under a CC0 1.0 Universal license. Comprehensive documentation, including a detailed vignette with worked examples and full API reference, is available at https://envirohealth.code-pages.usgs.gov/geoExposeR/. The package requires R >= 4.4.0 and can be installed via:

```r
# Install remotes if not already installed
install.packages("remotes")

# Install geoExposeR from GitLab (use a specific tag or branch)
remotes::install_gitlab(
  repo = "envirohealth/geoExposeR@v1.0.0",
  host = "code.usgs.gov"
)
```

# EXAMPLE USAGE

The following example demonstrates a basic analysis workflow:

```r
library(geoExposeR)

# Run sensitivity analysis with 5 imputed datasets
results <- perform_sensitivity_analysis(
  ndraws = 5,
  model_prob_csv = "modeled_probabilities.csv",
  conc_data_csv = "concentration_parameters.csv",
  birth_data_txt = "birth_outcomes.txt",
  regression_formula = "BWT ~ As_cat + MAGE + MRACE + (1|state)",
  output_dir = "analysis_results/",
  targets = c("BWT", "OEGEST"),
  impute_vars = c("MAGE", "MRACE")
)

# Access pooled regression results for the first target
results$BWT
```

The package also provides data preparation functions for users starting with raw GeoTIFF rasters or unprocessed concentration data. Complete examples demonstrating the full data pipeline, including additional examples for nitrate, uranium, and PM2.5 exposure assessments, are available in the package vignette and the `examples/` directory on GitLab.

# LIMITATIONS

The current implementation has the following limitations that users should consider. First, while the package is designed to be flexible, the default workflows are optimized for U.S.-based studies; international applications may require adaptation of data ingestion functions to accommodate different data formats and geographic coding systems. Second, the default analysis assumes exposure can be approximated at the county or geographic unit level using population-weighted averages; finer-scale exposure assessment (e.g., census tract or individual level) is supported through flexible identifier handling but may require additional data preparation. Third, the mixed-effects regression framework assumes linear relationships between exposure categories and health outcomes; users investigating non-linear exposure-response relationships should consider extensions or post-hoc analyses. Fourth, the uncertainty quantification approach relies on multiple imputation, which assumes that missing data mechanisms are ignorable (missing at random); violations of this assumption may affect inference. Finally, computational requirements scale with the number of imputed datasets and the size of the health outcome data; large-scale analyses with many imputations may require substantial memory and processing time.

# AUTHOR CONTRIBUTIONS

**Sayantan Majumdar**: Software, Methodology, Validation, Writing - Original Draft, Funding Acquisition. **Scott M. Bartell**: Methodology, Writing - Review & Editing. **Melissa A. Lombard**: Data Curation, Methodology, Resources, Writing - Review & Editing. **Ryan G. Smith**: Methodology, Funding Acquisition, Writing - Review & Editing, Project Administration. **Matthew O. Gribble**: Conceptualization, Funding Acquisition, Project Administration, Methodology, Supervision, Writing - Review & Editing.

# AI USAGE DISCLOSURE

Generative AI tools (GitHub Copilot and Claude) were used during the development of geoExposeR to assist with code generation, documentation writing, and manuscript preparation. Specifically, AI assistance was employed for: (1) generating boilerplate code for input validation functions, (2) drafting initial roxygen2 documentation strings, (3) creating unit test scaffolding, and (4) editing and refining text in this manuscript. All AI-generated content was critically reviewed, manually verified for correctness, and substantially revised by the authors. Statistical methods, software architecture decisions, and scientific interpretations were developed entirely by the human authors. Code correctness was verified through comprehensive unit testing (476 tests, 100% coverage), evaluating the package using actual data, and manual inspection of intermediate outputs. The final software and manuscript reflect the authors' expertise and judgment; AI tools served as productivity aids rather than intellectual contributors.

# ACKNOWLEDGMENTS

This work was supported by the National Heart, Lung, and Blood Institute (R21HL159574) and funding from the United States Geological Survey's John Wesley Powell Center for Analysis and Synthesis. Any use of trade, firm, or product names is for descriptive purposes only and does not imply endorsement by the U.S. Government.

# REFERENCES
