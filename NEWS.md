# geoExposeR 1.1.0

## Exposure Input Alignment

The modeled probability data (`model_prob_csv`) and the measured concentration data (`conc_data_csv`) are combined element-wise by row position, so row *i* of one file must describe the same unit as row *i* of the other. Previously the only check was that the two files had the same number of rows; if they were sorted differently, or covered different units with a coincidentally equal row count, every unit silently received another unit's concentration and private-well weight, and the analysis completed without any error or warning.

### Bug Fixes
* **Row alignment of the two exposure files is now verified** using `id_col`, which is present in both files but was previously read only from the probability data. `validate_data_alignment()` now:
  * returns the concentration data **reordered to match** the probability data when both contain the same identifiers in a different row order, emitting a `message()` describing the correction;
  * raises an **error** when the two files cover different sets of identifiers, naming the identifiers unique to each;
  * raises an **error** when the identifiers match but at least one file repeats them, since the rows cannot then be paired unambiguously;
  * emits a **warning** and falls back to the previous positional behavior when `id_col` is absent from either file, so existing workflows without a shared identifier continue to run.

Files that are already correctly aligned are unaffected: the output is unchanged, and no message or warning is emitted.

## Flexible Regression Model Types

This release adds flexible regression model options to `perform_sensitivity_analysis()`, so users can choose among fixed-effects, random-effects, and mixed-effects models in both linear and non-linear (generalized) variants.

### New Features
* **`model_type` argument**: Select the effects structure of the regression model:
  * `"fixed"` - fixed effects only, no random effects (fitted with `stats::lm()` / `stats::glm()`).
  * `"random"` - random-intercept-only models (fitted with `lme4::lmer()` / `lme4::glmer()`).
  * `"mixed"` (default) - fixed plus random effects (fitted with `lme4::lmer()` / `lme4::glmer()`). This preserves the previous default behavior, so existing code is unaffected.
* **`model_family` argument**: Select the response distribution for non-linear (generalized) models. Accepts `NULL` (default, a linear/Gaussian model), a family name such as `"binomial"` or `"poisson"`, a `stats::family()` object, or a family-generating function. When set, generalized models are fitted (`stats::glm()` for fixed effects, otherwise `lme4::glmer()`).
* **Engine selection** is automatic from the `model_type` × `model_family` combination; exposure-term coefficients are pooled across imputed datasets with Rubin's Rules regardless of the engine used.
* **Categorical outcomes**: For logistic (`model_family = "binomial"`) models, a two-level *character* target column (e.g. `"high"`/`"normal"`, as read from a text file) is automatically coerced to a factor, so categorical outcomes can be modeled alongside 0/1 indicators without manual preprocessing. Multiple targets of the same type can be analysed in a single call.
* **Multi-level categorical outcomes**: New `model_family = "multinomial"` fits a multinomial logistic regression for outcomes with three or more categories. Fixed-effects models use `nnet::multinom()`; random- and mixed-effects models use `mclogit::mblogit()`, where the inline `(1 | group)` random terms in `regression_formula` are translated to the engine's random-effects specification automatically. Pooled results gain a `y.level` column identifying the outcome category. Adds `nnet` and `mclogit` to package dependencies.

### Validation
* `regression_formula` is now validated against `model_type`: fixed-effects models reject random-effects terms, random- and mixed-effects models require at least one random-effects term, and random-effects models require intercept-only random terms.
* `model_type` and `model_family` are validated up front with clear error messages.

### Documentation and Examples
* New example script `examples/run_geoExposeR_model_variants.R` demonstrating all six fixed/random/mixed × linear/non-linear combinations.
* New example script `examples/run_geoExposeR_fixed_multi_target.R` demonstrating fixed-effects models over multiple targets at once, including categorical (character and 0/1) outcomes.
* New "Regression Model Options" section in the README.
* **Corrected random-effects interpretation in the examples.** The example health data (`examples/input_data/Demographics_dom.txt`) is now individual-level with multiple participants nested within each county (~400 records across 50 counties), generated reproducibly by the new `examples/input_data/generate_demographics.R`. The main example and vignette now use a geographic random intercept `(1 | id)` instead of `(1 | age_decade)`; documentation clarifies that grouping on `age_decade` is a *contextual* random intercept (no random slope) over only a handful of levels, and that a fixed effect (`+ as.factor(age_decade)`) is preferred for simply adjusting for age.

### Quality
* Maintains **100% test coverage** (579 tests, up from 476), with new tests covering every model-fitting and validation branch.

# geoExposeR 1.0.0

## First Stable Release

This is the first stable release of geoExposeR, a feature-complete R package for modeling health effects of environmental exposures. Originally designed for drinking water contaminants but broadly applicable to other environmental exposures.

### Highlights
* **100% test coverage** with 476 comprehensive tests
* **12 required dependencies** for robust functionality
* **7 exported functions** for data ingestion, analysis, and validation
* **Complete documentation** including vignettes and pkgdown site
* **CC0 1.0 Universal license** for unrestricted use

### New Features

#### Core Functionality
* **Main Function**: `perform_sensitivity_analysis()` - Complete workflow for exposure analysis (e.g., arsenic in drinking water)
* **Data Integration**: Combines modeled probability data with measured concentration data
* **Multiple Imputation**: Implements probabilistic exposure assignment across multiple datasets
* **Statistical Analysis**: Mixed-effects regression with proper pooling using Rubin's Rules
* **MICE Integration**: Optional multiple imputation for missing covariates

#### Data Loading and Processing
* **Probability Data Support**: Load and process modeled probability data from GeoTIFF rasters or tabular sources
* **Concentration Data Integration**: Convert measured concentration data to lognormal-derived multinomial probabilities
* **Weighted Combination**: Population-weighted integration of probability and concentration models
* **Flexible Identifier Support**: Works with any geographic or participant-level identifier (FIPS codes, census tracts, ZIP codes, or custom IDs)
* **Smart ID Formatting**: Automatic detection of identifier type—applies zero-padding for numeric FIPS codes (1-5 digits) while preserving alphanumeric identifiers unchanged. Controlled via `format_ids` parameter (enabled by default)
* **Birth Data Processing**: Specialized handling of health outcome datasets
* **Input Validation**: `validate_prepared_inputs()` for comprehensive pre-analysis data checks

#### Spatial Data Support
* **Raster Extraction**: Extract well water usage percentages from raster data using `add_well_percentages()`
* **Geographic Subsetting**: `subset_by_geography()` for filtering data by spatial boundaries
* **Coordinate Handling**: CRS transformation and centroid computation via sf and terra

#### Multiple Imputation Framework
* **Exposure Imputation**: Geographic-level probabilistic assignment based on combined models
* **Covariate Imputation**: MICE-based imputation for missing demographic and health variables
* **Validation Tools**: Comprehensive checks for imputation quality and convergence
* **Flexible Configuration**: Customizable imputation parameters and methods

#### Statistical Analysis
* **Mixed-Effects Models**: Support for complex nested data structures
* **Custom Formulas**: Flexible regression formula specification
* **Multiple Outcomes**: Simultaneous analysis of multiple health endpoints
* **Rubin's Rules**: Proper pooling of estimates across imputed datasets
* **Confidence Intervals**: Accurate uncertainty quantification

### Package Infrastructure

#### Testing and Quality Assurance
* **Comprehensive Test Suite**: 100% code coverage with automated testing
* **Synthetic Data Testing**: All tests use generated dummy data for reproducibility
* **GitLab CI/CD**: Automated testing across multiple R versions and operating systems
* **Performance Optimization**: Fast execution with `ndraws = 2` for development testing

#### Documentation and Usability
* **Extensive Documentation**: Complete roxygen2 documentation for all functions
* **Package Website**: Automated pkgdown site generation and deployment
* **Vignettes**: Step-by-step tutorials with working examples
* **README**: Comprehensive installation and usage instructions

#### Development Tools
* **Modular Code Organization**: Separate files for data ingestion, data loading, imputation, regression, and utilities
* **Internal Function Architecture**: Well-organized internal functions with clear responsibilities
* **Contribution Guidelines**: Detailed CONTRIBUTING.md with development setup instructions
* **Code Style**: Follows tidyverse style guidelines with consistent formatting

### Dependencies

#### Required Packages
* **Amelia**: Multiple imputation using bootstrapped EM algorithm
* **broom.mixed**: Tidy statistical output formatting
* **data.table**: High-performance data manipulation
* **dplyr**: Data transformation and summarization
* **exactextractr**: Fast raster value extraction for spatial analysis
* **Hmisc**: Weighted statistical functions
* **lme4**: Mixed-effects modeling
* **mice**: Multiple imputation by chained equations
* **rlang**: Tidy evaluation framework
* **rms**: Regression modeling strategies
* **sf**: Simple features for spatial vector data
* **terra**: Raster and spatial data processing

#### System Requirements
* **R** (>= 4.4.0): Minimum R version for compatibility
* **Pandoc**: Required for vignette building and documentation

### Performance and Scalability

#### Optimization Features
* **Vectorized Operations**: Efficient matrix operations for probability calculations
* **Memory Management**: Optimized data structures for large datasets
* **Parallel Processing**: Multi-core support for data loading operations
* **Caching**: Efficient dependency caching in CI/CD workflows

#### Practical Considerations
* **Geographic-Level Analysis**: Designed for exposure assessment at any geographic scale (counties, census tracts, ZIP codes, or custom regions)
* **Flexible Sample Sizes**: Supports datasets from small studies to large epidemiological cohorts
* **Configurable Parameters**: Adjustable imputation counts and convergence criteria
* **Output Management**: Comprehensive result saving and structured output formats

### Academic Integration

#### Methodological Foundation
* **Published Methods**: Implements approaches from Bulka et al. (2022) and Lombard et al. (2021)
* **Statistical Rigor**: Follows best practices for uncertainty quantification in exposure assessment
* **Reproducible Research**: Seed control and comprehensive output logging
* **Validation**: Extensive testing against known statistical properties

#### Research Applications
* **Birth Outcomes**: Specialized support for pregnancy and birth outcome studies
* **Environmental Epidemiology**: Designed for environmental health research workflows
* **Risk Assessment**: Tools for population-level exposure estimation
* **Policy Research**: Support for regulatory and public health decision-making

### Known Limitations

* **Geographic Scope**: Designed for US-based studies but supports any identifier format for international applications
* **Data Format Requirements**: Specific column naming conventions required for input data
* **Memory Usage**: Large datasets may require substantial memory for multiple imputation
* **Convergence**: MICE convergence warnings expected for complex imputation models

### Future Directions

* **International Support**: Expansion to non-US geographic coding systems
* **Additional Applications**: Expanding worked examples for other contaminants and exposure pathways
* **Advanced Modeling**: Integration of spatial autocorrelation and temporal trends
* **Visualization Tools**: Enhanced plotting and diagnostic visualization capabilities

---

## Development Team

Dr. Sayantan Majumdar (DRI), Dr. Scott M. Bartell (UC Irvine), Dr. Melissa A. Lombard (USGS), Dr. Ryan G. Smith (Colorado State University), Dr. Matthew O. Gribble (UCSF)

## Funding

This work was supported by the National Heart, Lung, and Blood Institute (R21HL159574) and funding from the United States Geological Survey’s John Wesley Powell Center for Analysis and Synthesis.

## Citation

Software release:

Majumdar, S., Bartell, S. M., Lombard, M. A., Smith, R. G., and Gribble, M. O., 2026, geoExposeR: An R package for modeling health effects of environmental exposures: U.S. Geological Survey software release, https://doi.org/10.5066/P1JGUKMD

Accompanying paper (under review at JOSS):

Majumdar, S., Bartell, S. M., Lombard, M. A., Smith, R. G., and Gribble, M. O., 2026, geoExposeR: An R package for modeling health effects of environmental exposures, Journal of Open Source Software* (under review).